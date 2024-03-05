#define DEBUG 0
#define PI 3.14159265358979323846
#define ECHARGE 1.602176634e-19
#define VACEPSILON 8.8541878128e-12
#define PROTONMASS 1.67262192369e-27
#define MAX_TRAJECTORY_POINTS 1048576
C     ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE OTHER FILES

      module trajectory_integration
      use iso_c_binding, only: c_int, c_double
      contains

      function minimal_test(in) bind(c, name="minimal_test")
      real(c_double), intent(in) :: in
      real(c_double) :: minimal_test
      minimal_test = 2 * in
      end function

C     Assumes that the independent variable increases monotonically
C     Also assumes that the input is a function
      subroutine interpolate_monotonic_1d(in_length, out_length, 
     &                                    t, f, new_t, new_f)
        integer(c_int), intent(in) :: in_length, out_length
        real(c_double),  dimension(in_length), intent(in)  :: f, t 
        real(c_double), dimension(out_length), intent(in)  :: new_t
        real(c_double), dimension(out_length), intent(out) :: new_f
        integer, dimension(out_length) :: idcs
        integer :: new_t_idx, old_t_idx, idx

        ! Initialize idcs
        idcs = -1
        new_t_idx = 1

        ! If new_t starts too early, we need to set it up for linear
        !   extrapolation at the start
        do while (new_t(new_t_idx) < t(1))
          idcs(new_t_idx) = 1
          new_t_idx = new_t_idx + 1
        end do

        do old_t_idx = 2, in_length - 1
          ! If we passed by an old_t value, record the index
          do while (t(old_t_idx) .GT. new_t(new_t_idx))
            idcs(new_t_idx) = old_t_idx - 1
            new_t_idx = new_t_idx + 1
            ! break if we have all of the idcs we need
            if (new_t_idx .GT. out_length) then
              goto 1
            end if
          end do
        end do
 1      continue

        ! Fill the remaining for linear extrapolation at the end
        do idx = new_t_idx, out_length
          if (idcs(idx) .EQ. -1) then
            idcs(idx) = in_length - 1
          end if
        end do

        ! Now we can do the interpolation in parallel
        new_f = -1
        do idx = 1,out_length
          if (idcs(idx) .GT. 0) then
            f1 = f(idcs(idx))
            f2 = f(idcs(idx) + 1)
            t1 = t(idcs(idx))
            t2 = t(idcs(idx) + 1)
            rel_t = new_t(idx) - t1
            new_f(idx) = f1 + rel_t * (f2 - f1) / (t2 - t1)
          end if
        end do
      end subroutine interpolate_monotonic_1d


      subroutine linspace_fortran(a, b, n, linspace)
        integer(c_int), intent(in) :: n
        real(c_double), intent(in) :: a, b
        real(c_double), dimension(n), intent(out) :: linspace
        if (b .EQ. a) then
          linspace = a
        else
          dx = (b - a) / (n - 1)
          do idx = 1,n
            linspace(idx) = dx * dble(idx - 1) + a
          end do
        end if
      end subroutine linspace_fortran

C     Function for linear interpolation on a 4x4x4 lattice
      pure function lininterpolate3D(matrix, xin, yin, zin, d_grid) 
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
      pure subroutine potEfunc(potential, x, y, z, d, Ex, Ey, Ez)
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
      pure subroutine interpolate_voltages(t, voltages, step_times, 
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
      potential = 0.0

C     Cut out a 4x4x4 area of the potential maps and compute linear sum
      cx = FLOOR(x / d)
      cy = FLOOR(y / d)
      cz = FLOOR(z / d)             ! Center coordinates

      if (DEBUG .eq. 1) then
            if     ((cx .le. 2) 
     &         .or. (cy .le. 2) 
     &         .or. (cz .le. 2)
     &         .or. (cx .ge. (dimensions(1) - 2))
     &         .or. (cy .ge. (dimensions(2) - 2))
     &         .or. (cz .ge. (dimensions(3) - 2)))    then
                  ex = 0.0
                  ey = 0.0
                  ez = 0.0
                  open (10, file='output_file.txt', position="append",
     &                        status='unknown')
                  write(10, *) "field_at"
                  write(10, *) t, x, y, z, d, cx, cy, cz
                  close(10)
                  goto 5
            end if
      end if

      do idx = 1,n_electrodes      ! Array slicing:
            do ix = 1,4
                  do iy = 1,4
                        do iz = 1,4
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
 5    continue
      end subroutine

      logical function is_dead(dimensions, is_electrode, x, y, z, d)
      integer(c_int), dimension(3), intent(in) :: dimensions
      integer(c_int), intent(in),
     &      dimension(dimensions(1), dimensions(2), dimensions(3))
     &      :: is_electrode
      real(c_double), intent(in) :: x, y, z, d
      integer(c_int) :: gx, gy, gz
      gx = NINT(x / d)
      gy = NINT(y / d)
      gz = NINT(z / d)
      is_dead = .false.
      if (    (gx .le. 2) 
     &   .or. (gy .le. 2) 
     &   .or. (gz .le. 2)
     &   .or. (gx .ge. (dimensions(1) - 2))
     &   .or. (gy .ge. (dimensions(2) - 2))
     &   .or. (gz .ge. (dimensions(3) - 2))
     &   .or. ISNAN(x)
     &   .or. ISNAN(y)
     &   .or. ISNAN(z)                        ) then
            is_dead = .true.
      else if (is_electrode(gx, gy, gz) .eq. 1) then
            is_dead = .true.
      end if
      end function

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
      
C     Unit conversions.  For simplicity, we work with SI units within 
C           this function.
      x = xx * 1.0e-3
      y = yy * 1.0e-3
      z = zz * 1.0e-3                           ! mm -> m
      vx = vxx * 1.0e+3
      vy = vyy * 1.0e+3
      vz = vzz * 1.0e+3                         ! mm/us -> m/s
      do idx = 1, time_steps
            step_times(idx) = step_times_in(idx) * 1.0e-6
      end do
      t = step_times(1)                         ! us -> s
      d = din * 1e-3                            ! mm -> m
      mdist = maxdist * 1.0e-3
      mt = maxt * 1.0e-6                        ! us -> s

C     Charge-to-mass ratio
      cmr = ((ECHARGE) * q) / ((PROTONMASS) * m);

C     Main loop of integration
      iter = 0
C     Check if particle is alive
      dead = is_dead(dimensions, is_electrode, x, y, z, d)
      if (dead) then
            goto 3
      end if
      call field_at(t, voltages, step_times, 
     &       n_electrodes, time_steps,
     &       x, y, z, d, potential_maps, dimensions,
     &       ex, ey, ez)
      do while ((t < mt) .and. (iter < MAX_TRAJECTORY_POINTS))
            iter = iter + 1

C           Record current state before it gets modified by the 
C                 integration step
            x_traj(iter) = x  * 1.0e+3
            y_traj(iter) = y  * 1.0e+3
            z_traj(iter) = z  * 1.0e+3
            exs(iter)    = ex * 1.0e-3
            eys(iter)    = ey * 1.0e-3
            ezs(iter)    = ez * 1.0e-3
            ts(iter)     = t  * 1.0e+6

C           Integration timestep calculation
            a = SQRT(ex*ex + ey*ey + ez*ez) * cmr
            v = SQRT(vx*vx + vy*vy + vz*vz)

            a = a + 1.0e-15
            v = v + 1.0e-15
            tv = mdist / v
            ta = SQRT(2.0 * mdist / a)
            tstep = tv * ta / (tv + ta)
            tstep = MIN(tstep, mt * 1.0e-3)
            t = t + tstep

C           Integration step
            x = x + tstep * vx + tstep * tstep * ex * cmr / 2
            y = y + tstep * vy + tstep * tstep * ey * cmr / 2
            z = z + tstep * vz + tstep * tstep * ez * cmr / 2


            call field_at(t, voltages, step_times, 
     &            n_electrodes, time_steps,
     &            x, y, z, d, potential_maps, dimensions,
     &            ex_new, ey_new, ez_new)

C           Check if particle is alive
            dead = is_dead(dimensions, is_electrode, x, y, z, d)
            if (dead) then
                  goto 3
            end if

            vx = vx + tstep * (ex_new + ex) * cmr / 2
            vy = vy + tstep * (ey_new + ey) * cmr / 2
            vz = vz + tstep * (ez_new + ez) * cmr / 2
            ex = ex_new
            ey = ey_new
            ez = ez_new
      end do
 3    its = dble(iter)
      end subroutine


C     Subroutine for Verlet integration of ion trajectory.
C     Data is recorded with time interval t_int
      subroutine integrate_trajectory_lite(pts, xx,yy,zz, vxx,vyy,vzz,
     &                        potential_maps, voltages, step_times_in,
     &                        time_steps, dimensions, is_electrode,
     &                        n_electrodes, m, q, din, maxdist, maxt,
     &                        x_traj, y_traj, z_traj, ts, its)
     & bind(c, name = "integrate_trajectory_lite")
C     Variable declarations:
C     Dummy variables:
      real(c_double), intent(in) :: xx, yy, zz, vxx, vyy, vzz, m, q, din
      real(c_double), intent(in) 
     &      :: maxdist, maxt
      integer(c_int), intent(in) :: time_steps, n_electrodes, pts
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
      integer(c_int) :: out_length
      real(c_double), dimension(pts), 
     &      intent(out)
     &      :: x_traj, y_traj, z_traj, ts
      real(c_double), intent(out) :: its

      real(c_double), dimension(MAX_TRAJECTORY_POINTS) :: 
     &      xs, ys, zs, exs, eys, ezs, t

      call integrate_trajectory(xx, yy, zz, vxx, vyy, vzz,
     &                        potential_maps, voltages, step_times_in,
     &                        time_steps, dimensions, is_electrode,
     &                        n_electrodes, m, q, din, maxdist, maxt,
     &                        xs, ys, zs, t, exs, eys, ezs, its)

      call linspace_fortran(real(0.0, c_double), t(MAX(1, NINT(its)-1)),
     &      pts, ts)

      call interpolate_monotonic_1d(MAX_TRAJECTORY_POINTS, pts, 
     &                                    t, xs, ts, x_traj)
      call interpolate_monotonic_1d(MAX_TRAJECTORY_POINTS, pts, 
     &                                    t, ys, ts, y_traj)
      call interpolate_monotonic_1d(MAX_TRAJECTORY_POINTS, pts, 
     &                                    t, zs, ts, z_traj)
      end subroutine

      subroutine fly_ensemble(interps, particles, xs,ys,zs,vxs,vys,vzs,
     &                        potential_maps, voltages, step_times_in,
     &                        time_steps, dimensions, is_electrode,
     &                        n_electrodes, m, q, din, maxdist, maxt,
     &                        x_trajs, y_trajs, z_trajs, tss, itss)
     & bind(c, name = "fly_ensemble")
C     Variable declarations:
C     Dummy variables:
      integer(c_int), intent(in) :: particles, interps
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
      real(c_double), dimension(particles, interps), 
     &      intent(out) :: x_trajs, y_trajs, z_trajs, tss
      real(c_double), dimension(particles),
     &      intent(out) :: itss

C     Local variables:
      real(c_double), dimension(interps)
     &      :: x, y, z, t
      real(c_double) :: its, increment
      integer :: i, j, k

      if (DEBUG .eq. 1) then
            open (10, file='output_file.txt',
     &            status='unknown')
            write(10, *) "fly_ensemble"
            close(10)
      end if

C C     Run the trajectory simulations
      do i = 1,particles
            if (DEBUG .eq. 1) then
                  open (10, file='output_file.txt', position="append",
     &                  status='unknown')
                  write(10, *) "PARTICLE: ", i
                  write(10, *) "INITIAL CONDITIONS:"
                  write(10, *) xs(i),ys(i),zs(i),vxs(i),vys(i),vzs(i)
                  close(10)
            end if
            call integrate_trajectory_lite(interps, xs(i), ys(i), zs(i),
     &                          vxs(i), vys(i), vzs(i),
     &                          potential_maps, voltages, step_times_in,
     &                          time_steps, dimensions, is_electrode,
     &                          n_electrodes, m, q, din, maxdist, maxt,
     &                          x, y, z, t, its)

            do j = 1,interps
                  x_trajs(i, j) = x(j)
                  y_trajs(i, j) = y(j)
                  z_trajs(i, j) = z(j)
                  tss(i, j)     = t(j)
            end do
            itss(i) = dble(its)
      end do
      end subroutine

      subroutine check_dead_cloud(particles, dead, xs, ys, zs, d,
     &                            dimensions, is_electrode)
      integer(c_int), intent(in)
     &      :: particles
      real(c_double), intent(in) 
     &      :: d
      real(c_double), dimension(particles), intent(in)
     &      :: xs, ys, zs
      integer(c_int), dimension(3), intent(in) 
     &      :: dimensions
      integer(c_int), intent(in),
     &      dimension(dimensions(1), dimensions(2), dimensions(3))
     &      :: is_electrode
      logical, dimension(particles), intent(inout)
     &      :: dead

      do i = 1,particles
            if (.NOT. dead(i)) then
                  dead(i) =
     &                  is_dead(dimensions, is_electrode, 
     &                          xs(i), ys(i), zs(i), d)
            end if
      end do
      end subroutine check_dead_cloud

      subroutine coulomb(xdiff, ydiff, zdiff, q1, q2, m1, m2,
     &      ax1, ay1, az1, ax2, ay2, az2)
      real(c_double), intent(in)
     &      :: xdiff, ydiff, zdiff, q1, q2, m1, m2
      real(c_double), intent(inout)
     &      :: ax1, ay1, az1, ax2, ay2, az2
      distance_sqrd =
     &        (xdiff * xdiff)
     &      + (ydiff * ydiff)
     &      + (zdiff * zdiff) ! m^2
      distance = SQRT(distance_sqrd)
      coulomb_force =
     &        (q1 * q2 * (ECHARGE**2.0) 
     &        / ( VACEPSILON * 4.0 * PI ))
     &            / distance_sqrd

      ax1 = ax1
     &      + (coulomb_force / (PROTONMASS * m1))
     &      * xdiff / distance
      ay1 = ay1
     &      + (coulomb_force / (PROTONMASS * m1))
     &      * ydiff / distance
      az1 = az1
     &      + (coulomb_force / (PROTONMASS * m1))
     &      * zdiff / distance
     
      ax2 = ax2
     &      - (coulomb_force / (PROTONMASS * m2))
     &      * xdiff / distance
      ay2 = ay2
     &      - (coulomb_force / (PROTONMASS * m2))
     &      * ydiff / distance
      az2 = az2
     &      - (coulomb_force / (PROTONMASS * m2))
     &      * zdiff / distance
      end subroutine coulomb
      
      subroutine fly_cloud(   record_interval, interps, particles, 
     &                        xxs, yys, zzs, vxxs, vyys, vzzs,
     &                        potential_maps, voltages, step_times_in,
     &                        time_steps, dimensions, is_electrode,
     &                        n_electrodes, ms, qs, din, maxdist, maxt,
     &                        x_trajs, y_trajs, z_trajs, ts, its)
     & bind(c, name = "fly_cloud")
      integer(c_int), intent(in)
     &      :: particles, interps, time_steps, n_electrodes
      real(c_double), intent(in)
     &      :: din, maxdist, maxt, record_interval
      real(c_double), dimension(particles), intent(in) 
     &      :: xxs, yys, zzs, vxxs, vyys, vzzs, ms, qs
      real(c_double), dimension(time_steps, n_electrodes), intent(in) 
     &      :: voltages
      integer(c_int), dimension(3), intent(in)
     &      :: dimensions
      integer(c_int), intent(in),
     &      dimension(dimensions(1), dimensions(2), dimensions(3))
     &      :: is_electrode
      real(c_double), 
     &      dimension(n_electrodes, dimensions(1), 
     &      dimensions(2), dimensions(3)), intent(in) 
     &      :: potential_maps
      real(c_double), dimension(time_steps), intent(in) 
     &      :: step_times_in
      real(c_double), dimension(particles, interps), intent(out) 
     &      :: x_trajs, y_trajs, z_trajs
      real(c_double), dimension(interps), intent(out) 
     &      :: ts
      real(c_double), intent(out)
     &      :: its

C     Local variables
      real(c_double), dimension(time_steps) 
     &      :: step_times
      real(c_double), dimension(particles)
     &      :: xs, ys, zs, vxs, vys, vzs, cmrs,
     &         exs, eys, ezs, exs_old, eys_old, ezs_old,
     &         axs, ays, azs, axs_old, ays_old, azs_old,
     &         tvs, tas, as, vs
      real(c_double), dimension(interps)
     &      :: temp_traj
      real(c_double) 
     &      :: t, mdist, mt, c, tv, ta, tstep, distance_sqrd, distance
     &         xdiff, ydiff, zdiff, coulomb_force, d, ex, ey, ez
      integer 
     &      :: idx, iter, alive
      logical, dimension(particles)
     &      :: dead
      real(c_double), allocatable
     &      :: vectorPath(:, :, :)
      real(c_double), allocatable
     &      :: tPath(:)

      allocate(vectorPath(3, particles, MAX_TRAJECTORY_POINTS))
      vectorPath = 0.0;
      allocate(tPath(MAX_TRAJECTORY_POINTS))
      tPath = 0.0;
      
C     Unit conversions.  For simplicity, we work with SI units within 
C           this function.
      ! positions, velocities, cmrs
      xs   = xxs  * 1.0e-3         ! mm -> m
      ys   = yys  * 1.0e-3         ! mm -> m
      zs   = zzs  * 1.0e-3         ! mm -> m
      vxs  = vxxs * 1.0e+3         ! mm/us -> m/s
      vys  = vyys * 1.0e+3         ! mm/us -> m/s
      vzs  = vzzs * 1.0e+3         ! mm/us -> m/s
      cmrs =   (ECHARGE / PROTONMASS) * (qs / ms) ! C/kg
      
      ! times
      step_times = step_times_in * 1.0e-6 ! us -> s
      
      ! other scalars
      t     = 0.0                  ! s
      d     = din * 1e-3           ! mm -> m
      mdist = maxdist * 1.0e-3     ! mm -> m
      mt    = maxt * 1.0e-6        ! us -> s

      axs = 0.0                    ! m/s^2
      ays = 0.0                    ! m/s^2
      azs = 0.0                    ! m/s^2
      exs = 0.0                    ! N/C
      eys = 0.0                    ! N/C
      ezs = 0.0                    ! N/C

C     MAIN INTEGRATION LOOP:
C     Velocity-Verlet
C     PHASE 1: Calculate/Inherit accelerations
C     PHASE 2: Adaptive timestep
C     PHASE 3: Half-timestep propagation of position
C     PHASE 4: Calculate new acceleration
C     PHASE 5: Half-timestep propagation of velocities mean acceleration
C     Subroutine for Verlet integration of ion trajectory.

C           PHASE 1: Calculate/Inherit Accelerations

      dead = .FALSE.
      call check_dead_cloud(particles, dead, xs, ys, zs, d,
     &                      dimensions, is_electrode)
      do idx = 1,particles
            if (.not. dead(idx)) then
                  call field_at(t, voltages, step_times, 
     &                    n_electrodes, time_steps,
     &                    xs(idx), ys(idx), zs(idx), d,
     &                    potential_maps, dimensions,
     &                    ex, ey, ez)
                  axs(idx) = ex * cmrs(idx)
                  ays(idx) = ey * cmrs(idx)
                  azs(idx) = ez * cmrs(idx)
                  ! Coulomb interaction, Newton's third law
                  do jdx = 1,(idx - 1)
                        call coulomb(xs(idx) - xs(jdx),            !m
     &                               ys(idx) - ys(jdx),            !m
     &                               zs(idx) - zs(jdx),            !m
     &                               qs(idx),  qs(jdx),            !e
     &                               ms(idx),  ms(jdx),            !m_p
     &                               axs(idx), ays(idx), azs(idx), !m/s2
     &                               axs(jdx), ays(jdx), azs(jdx))
                  end do
            end if
      end do

      iter = 0
      tstep = 0
      do while (iter .lt. MAX_TRAJECTORY_POINTS)
            iter = iter + 1

C           PHASE 1: Calculate/Inherit Accelerations
C           PHASE 2: Adaptive timestep
            prev_tstep = tstep
            as = SQRT((axs*axs) + (ays*ays) + (azs*azs)) + 1.0e-15
            vs = SQRT((vxs*vxs) + (vys*vys) + (vzs*vzs)) + 1.0e-15
            tvs = mdist / vs
            tas = SQRT(2.0 * mdist / as)
            tstep = MINVAL(tvs * tas / (tvs + tas))
            tstep = MIN(tstep, mt * 1.0e-6)
C             tstep = mt * 1.0e-6
            ! I record time and position before incrementing.
            ! However, I want position to take into account the state
            ! of the ion (dead/alive) so I leave that to a bit later
            tPath(iter) = t
            t = t + tstep

            do idx = 1,particles
                  if (.not. dead(idx)) then
                      ! Record the results
                      ! I do this first so the first data point == ICs
                      vectorPath(1, idx, iter) = xs(idx)
                      vectorPath(2, idx, iter) = ys(idx)
                      vectorPath(3, idx, iter) = zs(idx)
            ! This is from the previous step, to avoid a second loop:
C           PHASE 5: Half-timestep propagation of velocities,
                      vxs(idx) =   vxs(idx) 
     &                     + prev_tstep * (axs(idx) + axs_old(idx)) / 2.
                      vys(idx) =   vys(idx) 
     &                     + prev_tstep * (ays(idx) + ays_old(idx)) / 2.
                      vzs(idx) =   vzs(idx) 
     &                     + prev_tstep * (azs(idx) + azs_old(idx)) / 2.
C           PHASE 3: Half-timestep propagation of position
                      xs(idx) = xs(idx) + (vxs(idx) * tstep)
     &                          + ((axs(idx) * (tstep * tstep) / 2.0))
                      ys(idx) = ys(idx) + (vys(idx) * tstep)
     &                          + ((ays(idx) * (tstep * tstep) / 2.0))
                      zs(idx) = zs(idx) + (vzs(idx) * tstep)
     &                          + ((azs(idx) * (tstep * tstep) / 2.0))
                  end if
            end do

            ! Check if the simulation should be finished
            if (t.ge.mt) then
                  its = dble(iter) - 0.5
                  goto 2
            end if

            ! Check for casualties
            call check_dead_cloud(particles, dead, xs, ys, zs, d,
     &                            dimensions, is_electrode)

            ! Since some are left, we need to calculate the fields for:
C           PHASE 4: Calculate new acceleration
            alive = 1!0
            do idx = 1,particles
                  if (.not. dead(idx)) then
                        alive = alive + 1
                        axs_old(idx) = axs(idx)
                        ays_old(idx) = ays(idx)
                        azs_old(idx) = azs(idx)
C                         call field_at(t, voltages, step_times, 
C      &                          n_electrodes, time_steps,
C      &                          xs(idx), ys(idx), zs(idx), d,
C      &                          potential_maps, dimensions,
C      &                          ex, ey, ez)
                        ex = 0.0;
                        ey = 0.0;
                        ez = 0.0;
                        axs(idx) = ex * cmrs(idx)
                        ays(idx) = ey * cmrs(idx)
                        azs(idx) = ez * cmrs(idx)
                        ! Coulomb interaction, Newton's third law
                        do jdx = 1,(idx - 1)
                              call coulomb(xs(idx) - xs(jdx), 
     &                                     ys(idx) - ys(jdx), 
     &                                     zs(idx) - zs(jdx), 
     &                                     qs(idx), qs(jdx), 
     &                                     ms(idx), ms(jdx), 
     &                                     axs(idx), ays(idx), azs(idx),
     &                                     axs(jdx), ays(jdx), azs(jdx))
                        end do
                  end if
            end do

            ! If nobody is left, end the simulation
            if (alive .eq. 0) then
                  its = dble(iter)
                  goto 2
            end if
      end do
      its = dble(iter)
 2    continue

      tPath = tPath * 1.0e+6
      vectorPath = vectorPath * 1.0e+3

      iter = MAX(2, iter)
      call linspace_fortran(real(0.0, c_double), 
     &      tPath(iter), interps, ts)
      do i = 1,particles
            call interpolate_monotonic_1d(iter, interps, 
     &           tPath(1:iter), vectorPath(1, i, 1:iter), ts, temp_traj)
            x_trajs(i, :) = temp_traj

            call interpolate_monotonic_1d(iter, interps, 
     &           tPath(1:iter), vectorPath(2, i, 1:iter), ts, temp_traj)
            y_trajs(i, :) = temp_traj

            call interpolate_monotonic_1d(iter, interps, 
     &           tPath(1:iter), vectorPath(3, i, 1:iter), ts, temp_traj)
            z_trajs(i, :) = temp_traj
      end do
C       ts = tPath(1:interps)
C       do i = 1,particles
C             do j = 1,interps
C                   x_trajs(i, j) = vectorPath(1, i, j)
C                   y_trajs(i, j) = vectorPath(2, i, j)
C                   z_trajs(i, j) = vectorPath(3, i, j)
C             end do
C       end do
      end subroutine fly_cloud





      subroutine fly_aqs(amp_scales, a_res, off_scales, o_res, reps, 
     &                        xx, yy, zz, vxx, vyy, vzz,
     &                        rf_potential, endcap_potential,
     &                        voltages, step_times_in,
     &                        time_steps, dimensions, is_electrode,
     &                        n_electrodes, m, q, din, maxdist, maxt,
     &                        lifetimes)
     & bind(c, name = "fly_aqs")
C     Variable declarations:
C     Dummy variables:
      integer(c_int), intent(in) :: reps, a_res, o_res
      real(c_double), dimension(a_res), intent(in) :: amp_scales
      real(c_double), dimension(o_res), intent(in) :: off_scales
      real(c_double), intent(in) 
     &      :: xx, yy, zz, vxx, vyy, vzz
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
     &      :: rf_potential, endcap_potential
      real(c_double), dimension(time_steps), intent(in) 
     &      :: step_times_in
      real(c_double), dimension(a_res, o_res), intent(out) :: lifetimes

C     Local variables:
      real(c_double), dimension(MAX_TRAJECTORY_POINTS)
     &      :: x, y, z, t, ex, ey, ez
      real(c_double) :: its, increment
      real(c_double), allocatable :: potential_maps(:,:,:,:)
      integer :: i, j, k

C     Run the trajectory simulations
      lifetimes = 0.0
      allocate(potential_maps(n_electrodes, 
     &      dimensions(1), dimensions(2), dimensions(3)))

C       !$OMP PARALLEL DO DEFAULT(SHARED) 
C      &!$ PRIVATE(i, j, k, x, y, z, t, ex, ey, ez, its, potential_maps) 
C      &!$ SCHEDULE(STATIC, 16)
      do i = 1,a_res
            do j = 1,o_res
                  do k = 1,reps
                        potential_maps = amp_scales(i) * rf_potential
     &                        + endcap_potential
                        potential_maps(2,:,:,:)=potential_maps(2,:,:,:)
     &                        + rf_potential(1,:,:,:) * off_scales(j)
                        call integrate_trajectory(xx,yy,zz,vxx,vyy,vzz,
     &                        potential_maps, 
     &                        voltages, step_times_in,
     &                        time_steps, dimensions, is_electrode,
     &                        n_electrodes, m, q, din, maxdist, maxt,
     &                        x, y, z, t, ex, ey, ez, its)
C                         !$OMP ATOMIC
                        lifetimes(i, j) = lifetimes(i, j) 
     &                     + (t(int(its)) / dble(reps))
C                         !$OMP END ATOMIC
                  end do
            end do
      end do
C       !$OMP END PARALLEL DO

      end subroutine

      end module
