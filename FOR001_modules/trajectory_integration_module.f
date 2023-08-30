C     sudo gfortran -Wall -O3 -fPIC -cpp -c ./trajectory_integration_module.f
C                                      -> option "-cpp" enables C-style preprocessing for macros.
C     sudo gfortran -cpp -fc-prototypes -fsyntax-only ./trajectory_integration_module.f > ./trajectory_integration_module.h

#define MAX_TRAJECTORY_POINTS 131072
C     ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE C FILE ALSO

      module trajectory_integration
      use iso_c_binding, only: c_int, c_double
      implicit none

      contains
C     Function for linear interpolation on a 4x4x4 lattice
      function lininterpolate3D(matrix_in, xin, yin, zin, d_grid) 
     & bind(c, name = "lininterpolate3D")
      real(c_double), intent(in) :: xin, yin, zin, d_grid
      real(c_double), dimension(64), intent(in) :: matrix_in
      real(c_double), dimension(4, 4, 4) :: matrix
      real(c_double) :: lininterpolate3D
      integer(c_int) :: x_grid, y_grid, z_grid
      real(c_double) :: x_rel, y_rel, z_rel
      real(c_double) :: x, y, z
      real(c_double) :: a, b, c, d, e, f, g, h
      real(c_double) :: cg, ad, eh, bf
      real(c_double) :: fx, xp

C     Reshape the input matrix
      matrix = RESHAPE(matrix_in, (/4, 4, 4/))

C     Rescale to match the lattice
      x = xin / d_grid;
      y = yin / d_grid;
      z = zin / d_grid;

C     Find the index which represents the coordinates
      x_grid = FLOOR(x)
      y_grid = FLOOR(y)
      z_grid = FLOOR(z)

C     Find the coordinate relative to the calculated index
      x_rel = x - x_grid
      y_rel = y - y_grid
      z_rel = z - z_grid

C       x_grid = x_grid + 1
C       y_grid = y_grid + 1
C       z_grid = z_grid + 1

C     Store 8 of the nearest gridpoints to interpolate between
      a = matrix(x_grid, y_grid, z_grid)
      b = matrix(x_grid + 1, y_grid, z_grid)
      c = matrix(x_grid, y_grid + 1, z_grid)
      d = matrix(x_grid, y_grid, z_grid + 1)
      e = matrix(x_grid + 1, y_grid + 1, z_grid)
      f = matrix(x_grid + 1, y_grid, z_grid + 1)
      g = matrix(x_grid, y_grid + 1, z_grid + 1)
      h = matrix(x_grid + 1, y_grid + 1, z_grid + 1)

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
      cz = FLOOR(z / d)     ! Center coordinates
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
      tx = (2 + (x / d) - cx) * d
      ty = (2 + (y / d) - cy) * d
      tz = (2 + (z / d) - cz) * d
      call potEfunc(potential, tx, ty, tz, d, ex, ey, ez)
      end subroutine


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
      real(c_double), dimension(time_steps, n_electrodes)
     &      , intent(in) :: voltages
      integer(c_int), intent(in) :: time_steps, n_electrodes
      integer(c_int), dimension(3), intent(in) :: dimensions
      integer(c_int), intent(in),
     &      dimension(dimensions(1), dimensions(2), dimensions(3))
     &      :: is_electrode
      real(c_double), 
     &      dimension(n_electrodes, dimensions(1), 
     &      dimensions(2), dimensions(3)), 
     &      intent(in) :: potential_maps
      real(c_double), dimension(time_steps), intent(in) :: step_times_in
      real(c_double), dimension(time_steps) :: step_times
      real(c_double) :: x, y, z, vx, vy, vz, t, mdist, mt, cmr,
     &      a, v, tv, ta, tstep, ex, ey, ez, ex_new, ey_new, ez_new, d
      real(c_double), dimension(MAX_TRAJECTORY_POINTS), intent(out)
     &      :: x_traj, y_traj, z_traj, ts, exs, eys, ezs
      real(c_double), intent(out) :: its
      integer :: idx, iter
      real(c_double), dimension(n_electrodes) :: interpolated_voltages


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
      cmr = ((1.602176634e-19) * q) / ((1.660539067e-27) * m);

      call field_at(t, voltages, step_times, 
     &       n_electrodes, time_steps,
     &       x, y, z, d, potential_maps, dimensions,
     &       ex, ey, ez)    

C     Main loop of integration
      iter = 0
      do while ((t < mt) .and. (iter < MAX_TRAJECTORY_POINTS))
            iter = iter + 1

C           Record current state before it gets modified by the 
C                 integration step
            x_traj(iter) = x * 1.0e+3
            y_traj(iter) = y * 1.0e+3
            z_traj(iter) = z* 1.0e+3
            exs(iter)    = ex * 1.0e-3
            eys(iter)    = ey * 1.0e-3
            ezs(iter)    = ez * 1.0e-3
            ts(iter)     = t  * 1.0e+6

C           Integration timestep calculation
            a = SQRT(ex*ex + ey*ey + ez*ez) * cmr
            v = SQRT(vx*vx + vy*vy + vz*vz)

            a = a + 1.0e-15
            v = v + 1.0e-15
            tv = mdist / v;
            ta = SQRT(2.0 * mdist / a)
            tstep = tv * ta / (tv + ta)
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
            if ((x < 2 * d) .or. (y < 2 * d) .or. (z < 2 * d)
     &            .or. (x > (dimensions(1) - 2) * d)
     &            .or. (y > (dimensions(2) - 2) * d)
     &            .or. (z > (dimensions(3) - 2) * d)) then
                  t = maxt
C             else if (is_electrode(NINT(x/d), NINT(y/d), NINT(z/d)) == 1)
C      &       then
C                   t = maxt
            end if

            vx = vx + tstep * (ex_new + ex) * cmr / 2
            vy = vy + tstep * (ey_new + ey) * cmr / 2
            vz = vz + tstep * (ez_new + ez) * cmr / 2
            ex = ex_new
            ey = ey_new
            ez = ez_new

      end do
      its = dble(iter)
      end subroutine
      end module