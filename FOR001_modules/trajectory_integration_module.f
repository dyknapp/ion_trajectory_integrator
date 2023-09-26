#define MAX_TRAJECTORY_POINTS 1048576
C     ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE OTHER FILES

      module trajectory_integration
      use iso_c_binding, only: c_int, c_double
      implicit none

      contains

C     Function for linear interpolation on a 4x4x4 lattice
      function lininterpolate3D(matrix, xin, yin, zin, d_grid) 
     & bind(c, name = "lininterpolate3D")
      real(c_double), intent(in) :: xin, yin, zin, d_grid
      real(c_double), dimension(64), intent(in) :: matrix
      real(c_double) :: lininterpolate3D
      real(c_double) :: x_rel, y_rel, z_rel
      real(c_double) :: x, y, z
      real(c_double) :: a, b, c, d, e, f, g, h
      real(c_double) :: cg, ad, eh, bf
      real(c_double) :: fx, xp
      integer :: x_grid, y_grid, z_grid

C     Correct to grid spacing
      x = xin / d_grid
      y = yin / d_grid
      z = zin / d_grid

C     Grid coordinates
      x_grid = FLOOR(x)
      y_grid = FLOOR(y)
      z_grid = FLOOR(z)

C     Find the coordinate relative to the calculated index
      x_rel = x - x_grid
      y_rel = y - y_grid
      z_rel = z - z_grid

C     Store 8 of the nearest gridpoints to interpolate between
      a = matrix(x_grid     + 4*(y_grid - 1) + 16*(z_grid - 1))
      b = matrix(x_grid + 1 + 4*(y_grid - 1) + 16*(z_grid - 1))
      c = matrix(x_grid     + 4*(y_grid    ) + 16*(z_grid - 1))
      d = matrix(x_grid     + 4*(y_grid - 1) + 16*(z_grid    ))
      e = matrix(x_grid + 1 + 4*(y_grid    ) + 16*(z_grid - 1))
      f = matrix(x_grid + 1 + 4*(y_grid - 1) + 16*(z_grid    ))
      g = matrix(x_grid     + 4*(y_grid    ) + 16*(z_grid    ))
      h = matrix(x_grid + 1 + 4*(y_grid    ) + 16*(z_grid    ))

C     interpolate matrix along the gridlines
      cg = z_rel * g + (1 - z_rel) * c
      ad = z_rel * d + (1 - z_rel) * a
      eh = z_rel * h + (1 - z_rel) * e
      bf = z_rel * f + (1 - z_rel) * b

C     interplate matrix along faces of the evaluated grid
      fx  = y_rel * cg + (1 - y_rel) * ad;
      xp = y_rel * eh + (1 - y_rel) * bf;

      lininterpolate3D = (x_rel * xp) + ((1 - x_rel) * fx)
      end function

C     Calculate the electric field by computing the gradient of a point 
C           of a given 4x4x4 potential lattice.
      subroutine potEfunc(potential, x, y, z, d, Ex, Ey, Ez)
     & bind(c, name = "potEfunc")
      real(c_double), intent(in) :: x, y, z, d
      real(c_double), dimension(4, 4, 4), intent(in) :: potential
      real(c_double), intent(out) :: Ex, Ey, Ez

      Ex = -(linInterpolate3D(potential, x + d/2.0, y, z, d) 
     &          - linInterpolate3D(potential, x - d/2.0, y, z, d)) / (d)
      Ey = -(linInterpolate3D(potential, x, y + d/2.0, z, d) 
     &          - linInterpolate3D(potential, x, y - d/2.0, z, d)) / (d)
      Ez = -(linInterpolate3D(potential, x, y, z + d/2.0, d)
     &          - linInterpolate3D(potential, x, y, z - d/2.0, d)) / (d)
      end subroutine

C     Subroutine for interpolating electrode voltages
C     expects a size [n_electrodes, time_steps] sized matrix
      subroutine interpolate_voltages(t, voltages, step_times, 
     &      n_electrodes, time_steps, interpolated_voltages)
     & bind(c, name = "interpolate_voltages")
      integer(c_int), intent(in) :: n_electrodes, time_steps
      real(c_double), intent(in) :: t
      real(c_double), dimension(time_steps, n_electrodes)
     &      , intent(in) :: voltages
      real(c_double), dimension(time_steps), intent(in)
     &      :: step_times
      real(c_double), dimension(n_electrodes), intent(out)
     &      :: interpolated_voltages
      real(c_double) :: start_time
      real(c_double) ::   end_time
      integer :: prev_point_idx, next_point_idx, idx

      start_time =  step_times(1)
      end_time   = step_times(time_steps)

      prev_point_idx = 1 + 
     &      FLOOR((time_steps - 1) 
     &      * (t - start_time) / (end_time - start_time))
      next_point_idx = 1 + prev_point_idx

      do idx = 1, n_electrodes
            interpolated_voltages(idx) 
     &            = ((t - step_times(prev_point_idx))
     &            * (voltages(next_point_idx, idx) 
     &                  - voltages(prev_point_idx, idx))
     &            / (step_times(next_point_idx) 
     &                  - step_times(prev_point_idx)))
     &            + voltages(prev_point_idx, idx)
      end do
      end subroutine

C     Subroutine for computing the electric field at a point
      subroutine field_at(t, voltages, step_times, 
     &      n_electrodes, time_steps,
     &      x, y, z, d, potential_maps, dimensions,
     &      ex, ey, ez)
     & bind(c, name = "field_at")
C     Variable declarations:      
      integer(c_int), intent(in) :: n_electrodes, time_steps
      real(c_double), intent(in) :: t, x, y, z, d
      real(c_double), dimension(time_steps, n_electrodes)
     &      , intent(in) :: voltages
      real(c_double), dimension(time_steps), intent(in)
     &      :: step_times
      integer(c_int), dimension(3), intent(in) :: dimensions
      real(c_double), 
     &      dimension(n_electrodes, dimensions(1), 
     &      dimensions(2), dimensions(3)), 
     &      intent(in) :: potential_maps
      real(c_double), intent(out) :: ex, ey, ez
      real(c_double), dimension(n_electrodes) :: interpolated_voltages
      real(c_double), dimension(4, 4, 4) :: potential
      real(c_double) :: tx, ty, tz
      integer :: cx, cy, cz, ix, iy, iz, idx


C     What are the elctrode voltages right now?
      call interpolate_voltages(t, voltages, step_times, 
     &            n_electrodes, time_steps, interpolated_voltages)

C     Initialize potential array
      do ix = 1, 4
            do iy = 1, 4
                  do iz = 1, 4
                        potential(ix, iy, iz) = 0.0
                  end do
            end do
      end do

C     Cut out a 4x4x4 area of the potential maps and compute linear sum
      cx = FLOOR(x / d)
      cy = FLOOR(y / d)
      cz = FLOOR(z / d)             ! Center coordinates
      do idx = 1, n_electrodes      ! Array slicing:
            do ix = 1, 4
                  do iy = 1, 4
                        do iz = 1, 4
                              potential(ix, iy, iz) =
     &                                             potential(ix, iy, iz)
     &                                     + (interpolated_voltages(idx)
     &     * potential_maps(idx, cx + ix - 2, cy + iy - 2, cz + iz - 2))
                        end do
                  end do
            end do
      end do

C     Correct the coordinates' offset and compute the E field
      tx = 2*d + x - cx*d
      ty = 2*d + y - cy*d
      tz = 2*d + z - cz*d
      call potEfunc(potential, tx, ty, tz, d, ex, ey, ez)
      end subroutine

      function is_dead(dimensions, is_electrode, x, y, z, d)
      integer(c_int), dimension(3), intent(in) :: dimensions
      integer(c_int), intent(in),
     &      dimension(dimensions(1), dimensions(2), dimensions(3))
     &      :: is_electrode
      real(c_double), intent(in) :: x, y, z, d
      logical :: is_dead
            if ((x < 2 * d) .or. (y < 2 * d) .or. (z < 2 * d)
     &            .or. (x > (dimensions(1) - 2) * d)
     &            .or. (y > (dimensions(2) - 2) * d)
     &            .or. (z > (dimensions(3) - 2) * d)) then
                  is_dead = .true.
            else if (is_electrode(NINT(x/d), NINT(y/d), NINT(z/d)) == 1)
     &       then
                  is_dead = .true.
            else
                  is_dead = .false.
            end if
      end function


C     Subroutine for Verlet integration of ion trajectory.
      subroutine integrate_trajectory(xx, yy, zz, vxx, vyy, vzz,
     &                        potential_maps, voltages, step_times_in,
     &                        time_steps, dimensions, is_electrode,
     &                        n_electrodes, m, q, din, maxdist, maxt,
     &                        x_traj,y_traj,z_traj,ts,exs,eys,ezs,its)
     & bind(c, name = "integrate_trajectory")
C     Variable declarations:
C     Dummy variables:
      real(c_double), intent(in) :: xx, yy, zz, vxx, vyy, vzz, m, q, din
      real(c_double), intent(in) :: maxdist, maxt
      integer(c_int), intent(in) :: time_steps, n_electrodes
      real(c_double), dimension(time_steps, n_electrodes)
     &      , intent(in) :: voltages
      integer(c_int), dimension(3), intent(in) :: dimensions
      integer(c_int), intent(in),
     &      dimension(dimensions(1), dimensions(2), dimensions(3))
     &      :: is_electrode
      real(c_double), 
     &      dimension(n_electrodes, dimensions(1), 
     &      dimensions(2), dimensions(3)), 
     &      intent(in) :: potential_maps
      real(c_double), dimension(time_steps), intent(in) :: step_times_in
      real(c_double), dimension(MAX_TRAJECTORY_POINTS), intent(out)
     &      :: x_traj, y_traj, z_traj, ts, exs, eys, ezs
      real(c_double), intent(out) :: its

C     Local variables
      real(c_double), dimension(time_steps) :: step_times
      real(c_double) :: x, y, z, vx, vy, vz, t, mdist, mt, cmr,
     &      a, v, tv, ta, tstep, ex, ey, ez, ex_new, ey_new, ez_new, d
      integer :: idx, iter
      logical :: dead


C C     Unit conversions.  For simplicity, we work with SI units within 
C C           this function.
C       x = xx * 1.0e-3
C       y = yy * 1.0e-3
C       z = zz * 1.0e-3                           ! mm -> m
C       vx = vxx * 1.0e+3
C       vy = vyy * 1.0e+3
C       vz = vzz * 1.0e+3                         ! mm/us -> m/s
C       do idx = 1, time_steps
C             step_times(idx) = step_times_in(idx) * 1.0e-6
C       end do
C       t = step_times(1)                         ! us -> s
C       d = din * 1e-3                            ! mm -> m
C       mdist = maxdist * 1.0e-3
C       mt = maxt * 1.0e-6                        ! us -> s

C C     Charge-to-mass ratio
C       cmr = ((1.602176634e-19) * q) / ((1.660539067e-27) * m);

C       call field_at(t, voltages, step_times, 
C      &       n_electrodes, time_steps,
C      &       x, y, z, d, potential_maps, dimensions,
C      &       ex, ey, ez)    

C C     Main loop of integration
C       iter = 0
C C     Check if particle is alive
C       dead = is_dead(dimensions, is_electrode, x, y, z, d)
C       if (dead) then
C             t = mt;
C       end if
C       do while ((t < mt) .and. (iter < MAX_TRAJECTORY_POINTS))
C             iter = iter + 1

C C           Record current state before it gets modified by the 
C C                 integration step
C             x_traj(iter) = x  * 1.0e+3
C             y_traj(iter) = y  * 1.0e+3
C             z_traj(iter) = z  * 1.0e+3
C             exs(iter)    = ex * 1.0e-3
C             eys(iter)    = ey * 1.0e-3
C             ezs(iter)    = ez * 1.0e-3
C             ts(iter)     = t  * 1.0e+6

C C           Integration timestep calculation
C             a = SQRT(ex*ex + ey*ey + ez*ez) * cmr
C             v = SQRT(vx*vx + vy*vy + vz*vz)

C             a = a + 1.0e-15
C             v = v + 1.0e-15
C             tv = mdist / v
C             ta = SQRT(2.0 * mdist / a)
C             tstep = tv * ta / (tv + ta)
C             tstep = MIN(tstep, mt * 1.0e-3)
C             t = t + tstep

C C           Integration step
C             x = x + tstep * vx + tstep * tstep * ex * cmr / 2
C             y = y + tstep * vy + tstep * tstep * ey * cmr / 2
C             z = z + tstep * vz + tstep * tstep * ez * cmr / 2


C             call field_at(t, voltages, step_times, 
C      &            n_electrodes, time_steps,
C      &            x, y, z, d, potential_maps, dimensions,
C      &            ex_new, ey_new, ez_new)

C C           Check if particle is alive
C             dead = is_dead(dimensions, is_electrode, x, y, z, d)
C             if (dead) then
C                   t = mt;
C             end if

C             vx = vx + tstep * (ex_new + ex) * cmr / 2
C             vy = vy + tstep * (ey_new + ey) * cmr / 2
C             vz = vz + tstep * (ez_new + ez) * cmr / 2
C             ex = ex_new
C             ey = ey_new
C             ez = ez_new
C       end do
C       its = dble(iter)
      end subroutine

      subroutine fly_ensemble(particles, xs,ys,zs,vxs,vys,vzs,
     &                        potential_maps, voltages, step_times_in,
     &                        time_steps, dimensions, is_electrode,
     &                        n_electrodes, m, q, din, maxdist, maxt,
     &                        x_trajs, y_trajs, z_trajs, 
     &                        tss, exss, eyss, ezss, itss)
     & bind(c, name = "fly_ensemble")
C     Local variables:
      real(c_double), dimension(MAX_TRAJECTORY_POINTS)
     &      :: x_traj, y_traj, z_traj, ts, exs, eys, ezs
      integer :: i, j
      real(c_double) :: its

C     Variable declarations:
C     Dummy variables:
      integer(c_int), intent(in) :: particles!, interps

      real(c_double), dimension(particles), intent(in) 
     &      :: xs, ys, zs, vxs, vys, vzs

      real(c_double), intent(in) :: m, q, din

      real(c_double), intent(in) :: maxdist, maxt
     
      integer(c_int), intent(in) :: time_steps, n_electrodes

      real(c_double), dimension(time_steps, n_electrodes), intent(in) 
     &      :: voltages

      integer(c_int), dimension(3), intent(in) :: dimensions

      integer(c_int), intent(in),
     &      dimension(dimensions(1), dimensions(2), dimensions(3))
     &      :: is_electrode

      real(c_double), 
     &      dimension(n_electrodes, dimensions(1), 
     &      dimensions(2), dimensions(3)), intent(in) 
     &      :: potential_maps

      real(c_double), dimension(time_steps), intent(in) 
     &      :: step_times_in

      real(c_double), dimension(particles, MAX_TRAJECTORY_POINTS), 
     &      intent(out) :: x_trajs,y_trajs,z_trajs,tss,exss,eyss,ezss

      real(c_double), dimension(particles),
     &      intent(out) :: itss

C     Allocate the output variables
      !interps = MAX_TRAJECTORY_POINTS ! For now, no interpolation
C       allocate(x_trajs(particles, MAX_TRAJECTORY_POINTS))
C       allocate(y_trajs(particles, MAX_TRAJECTORY_POINTS))
C       allocate(z_trajs(particles, MAX_TRAJECTORY_POINTS))
C       allocate(tss(particles, MAX_TRAJECTORY_POINTS))
C       allocate(exss(particles, MAX_TRAJECTORY_POINTS))
C       allocate(eyss(particles, MAX_TRAJECTORY_POINTS))
C       allocate(ezss(particles, MAX_TRAJECTORY_POINTS))
C       allocate(itss(particles))

C     Run the trajectory simulations
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j) SCHEDULE(STATIC,16)
      do i = 1,particles
            call integrate_trajectory(xs(i), ys(i), zs(i), 
     &                          vxs(i), vys(i), vzs(i),
     &                          potential_maps, voltages, step_times_in,
     &                          time_steps, dimensions, is_electrode,
     &                          n_electrodes, m, q, din, maxdist, maxt,
     &                          x_traj,y_traj,z_traj,ts,exs,eys,ezs,its)
            do j = 1,NINT(its)
                  x_trajs(i, j) = x_traj(j)
                  y_trajs(i, j) = y_traj(j)
                  z_trajs(i, j) = z_traj(j)
                  tss(i, j)     = ts(j)
                  exss(i, j)    = exs(j)
                  eyss(i, j)    = eys(j)
                  ezss(i, j)    = ezs(j)
            end do
            itss(i) = its
            continue
      end do
      !$OMP END PARALLEL DO

      end subroutine

      end module
