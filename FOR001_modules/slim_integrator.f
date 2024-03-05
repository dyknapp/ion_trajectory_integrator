#define DEBUG 1
#define PI 3.14159265358979323846
#define ECHARGE 1.602176634e-19
#define VACEPSILON 8.8541878128e-12
#define PROTONMASS 1.67262192369e-27
#define MAX_TRAJECTORY_POINTS 1048576
C     ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE OTHER FILES

      module slim_integrators
      use iso_c_binding, only: c_int, c_double
      contains

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

C     Function for linear interpolation on a 4x4 lattice
      function lininterpolate2D_unit_grid(matrix, x, y)
      real(c_double), intent(in) :: x, y
      real(c_double), dimension(4, 4), intent(in) :: matrix
      real(c_double) :: lininterpolate2D_unit_grid
      real(c_double) :: x_rel, y_rel
      real(c_double) :: a, b, c, d
      real(c_double) :: fx, xp
      integer :: x_grid, y_grid

C     Grid coordinates
      x_grid = FLOOR(x)
      y_grid = FLOOR(y)

C     Find the coordinate relative to the calculated index
      x_rel = x - DBLE(x_grid)
      y_rel = y - DBLE(y_grid)

C     Store 8 of the nearest gridpoints to interpolate between
      a = matrix(x_grid    , y_grid + 1)
      b = matrix(x_grid + 1, y_grid + 1)
      c = matrix(x_grid    , y_grid    )
      d = matrix(x_grid + 1, y_grid    )

C     interplate matrix along faces of the evaluated grid
      fx  = y_rel * a + (1 - y_rel) * c;
      xp  = y_rel * b + (1 - y_rel) * d;

      lininterpolate2D_unit_grid = (x_rel * xp) + ((1 - x_rel) * fx)
      end function

      subroutine potEfunc_unit_grid(potential, x, y, Ex, Ey)
      real(c_double), intent(in) :: x, y
      real(c_double), dimension(4, 4), intent(in) :: potential
      real(c_double), intent(out) :: Ex, Ey

      Ex = -(lininterpolate2D_unit_grid(potential, x + 0.5, y      ) 
     &     - linInterpolate2D_unit_grid(potential, x - 0.5, y      )
     &      )
      Ey = -(linInterpolate2D_unit_grid(potential, x      , y + 0.5) 
     &     - linInterpolate2D_unit_grid(potential, x      , y - 0.5)
     &      )
      end subroutine

C     Subroutine for computing the electric field at a point
      subroutine field_at2D_t_radial(r, z, d, potential_map, dimensions,
     &      n_electrodes, voltages, er, ez)
C     Variable declarations:
      real(c_double), intent(in) :: r, z, d
      integer(c_int), dimension(2), intent(in) :: dimensions
      integer(c_int), intent(in) :: n_electrodes
      real(c_double), 
     &      dimension(n_electrodes, dimensions(1), dimensions(2)), 
     &      intent(in) :: potential_map
      real(c_double), dimension(n_electrodes), intent(in) :: voltages
      real(c_double), intent(out) :: er, ez
      real(c_double), dimension(4, 4) :: potential
      real(c_double) :: tr, tz
      integer :: cr, cz, ir, iz, idx

C     Initialize potential array
      potential = 0.0

C     Cut out a 4x4 area of the potential maps and compute linear sum
      cr = FLOOR(ABS(r) / d)
      cz = FLOOR(    z  / d) ! Center coordinates
      ! Array slicing:
      do n = 1,n_electrodes
            do ir = 1,4
                  do iz = 1,4
C       if (DEBUG .eq. 1) then
C             open (10, file='output_file.txt', position="append",
C      &                  status='unknown')
C             write(10, *) "potential_map"
C             write(10, *) n, 1 + ABS(cr + ir - 2), cz + iz - 2
C             close(10)
C       end if
                        potential(ir, iz) = potential(ir, iz)
     &                        + (voltages(n)
     &                        * potential_map(n,
     &                                        1 + ABS(cr + ir - 2), 
     &                                        cz + iz - 2))
                        ! Adding 1 to cr here corrects the offset from
                        !     indexing starting at 1 instead of zero.
                  end do
            end do
      end do     
      !potential = potential_map((cx-1):(cx+2), (cy-1):(cy+2))

C     Correct the coordinates' offset and compute the E field
      tr = 2.0 * d + ABS(r) - DBLE(cr) * d
      tz = 2.0 * d +     z  - DBLE(cz) * d
      call potEfunc_unit_grid(potential, tr/d, tz/d, er, ez)
      er = er / d
      ez = ez / d
      end subroutine

      function how_dead(dimensions, is_electrode, r, z, d)
      integer(c_int), dimension(2), intent(in) :: dimensions
      integer(c_int), dimension(dimensions(1), dimensions(2)), 
     &      intent(in):: is_electrode
      real(c_double), intent(in) :: r, z, d
      integer(c_int) :: how_dead
      how_dead = 0
      if        (z < 2 * d                  ) then
            how_dead = 10
            return
      elseif    (r > (dimensions(1) - 2) * d) then
            how_dead = 20
            return

      elseif    (z > (dimensions(2) - 2) * d) then
            how_dead = 30
            return
      else
            if (is_electrode(ABS(NINT(r/d)) + 1, NINT(z/d)).eq.1) then
                  how_dead = 2
                  return
            end if
      end if
      end function

      subroutine cylindrical_3D_spaced(position, velocity, sample_dist,
     &            is_electrode, potential_maps, voltages, voltage_lines,
     &            dimensions, n_electrodes, m, q, din, maxdist, maxt,
     &            trajectory, death, its, datas)
     & bind(c, name='cylindrical_3D_spaced')
C     Variable declarations:
C     Dummy variables:
      real(c_double), dimension(3), intent(in) :: position, velocity
      real(c_double), intent(in) :: m, q, din, sample_dist
      real(c_double), intent(in) :: maxdist, maxt
      integer(c_int), intent(in) :: n_electrodes, voltage_lines
      real(c_double), dimension(3, voltage_lines), intent(in) 
     &      :: voltages
      integer(c_int), dimension(2), intent(in) :: dimensions
      integer(c_int), dimension(dimensions(1), dimensions(2)), 
     &      intent(in):: is_electrode
      real(c_double),
     &      dimension(n_electrodes, dimensions(1), dimensions(2)), 
     &      intent(in) :: potential_maps
      real(c_double), dimension(1024, 4), intent(out) :: trajectory
      integer(c_int), intent(out) :: its, datas, death

C     Local variables
      real(c_double) :: t, mdist, mt, cmr,
     &      a, v, tv, ta, tstep, er, ez, d, r
      real(c_double), allocatable :: 
     &      pos(:), vel(:), last_rec(:), accel(:), accel_n(:), volts(:)
      integer(c_int) :: idx, iter, sample_interval, data_points, v_idx
      logical :: dead

      sample_interval = MAX_TRAJECTORY_POINTS / 1024

      allocate(pos(3))
      allocate(vel(3))
      allocate(last_rec(3))
      allocate(accel(3))
      allocate(accel_n(2))
      allocate(volts(n_electrodes))
      
C     Unit conversions.  For simplicity, we work with SI units within 
C           this function.
      pos = position * 1.0e-3 ! mm -> m
      last_rec = pos
      vel = velocity * 1.0e+3 ! mm/us -> m/s
      
      t = 0.0;
      d = din * 1.0e-3                          ! mm -> m
      mdist = maxdist * 1.0e-3
      mt = maxt * 1.0e-6                        ! us -> s

      trajectory = PI

C     Charge-to-mass ratio
      cmr = ((1.602176634e-19) * q) / ((1.660539067e-27) * m);


C     Main loop of integration
      iter = 0

C     Initial voltages
      volts = 0.0
      v_idx = 1
      if (v_idx .le. voltage_lines) then
            do while (voltages(v_idx, 2) .eq. 0.0)
                  volts(NINT(voltages(v_idx, 1))) 
     &                  = voltages(v_idx, 3)
                  v_idx = v_idx + 1;
            end do
      end if

C     Check if particle is alive
      r = NORM2(pos(1:2))
      er = 0.0
      ez = 0.0
      death = how_dead(dimensions, is_electrode, r, pos(3), d)
      if (death.gt.0) then
            t = mt;
      else
            call field_at2D_t_radial(r, pos(3), d, 
     &            potential_maps, dimensions,
     &            n_electrodes, volts,
     &            er, ez)
            if (r .gt. 0.0) then
                  accel = (/er * pos(1), er * pos(2), ez * r/) * (cmr/r)
            else
                  accel = (/er * 0.0, er * 0.0, ez * cmr/)
            end if
      end if
      data_points = 0
      x_traj = 0.0
      y_traj = 0.0
      ts     = 0.0
      do while ((t.lt.mt) .and. (data_points.lt.1024))
            iter = iter + 1

C           Record current state before it gets modified by the 
C                 integration step
            if (NORM2(pos - last_rec).ge.(sample_dist*1.0e-3)) then
                  data_points = data_points + 1
                  last_rec = pos
                  trajectory(data_points, 1:3) = pos * 1.0e+3
                  trajectory(data_points, 4)   = t   * 1.0e+6
            end if
            
C           Integration timestep calculation
            a = NORM2(accel)
            v = NORM2(vel)

            a = a + 1.0e-15
            v = v + 1.0e-15
            tv = mdist / v
            ta = SQRT(2.0 * mdist / a)
            tstep = tv * ta / (tv + ta)
            tstep = MIN(tstep, mt * 1.0e-6)
            t = t + tstep

C           Integration step
            pos = pos + tstep*vel + tstep*tstep*accel/2.

C           Check if particle is alive
            r = NORM2(pos(1:2))
            death = how_dead(dimensions,is_electrode,r,pos(3),d)
            if (death.gt.0) then
                  t = mt;
            else
C                 Update voltages
                  if (v_idx .le. voltage_lines) then
                        do while (voltages(v_idx, 2).le.t)
                              volts(NINT(voltages(v_idx, 1)))
     &                              = voltages(v_idx, 3)
                              v_idx = v_idx + 1;
                        end do
                  end if
C                 Calculate new fields
                  call field_at2D_t_radial(r, pos(3), d, 
     &                  potential_maps, dimensions,
     &                  n_electrodes, volts,
     &                  er, ez)
                  if (r .gt. 0.0) then
                        accel_n = 
     &                    (/er * pos(1), er * pos(2), ez * r/) * cmr / r
                  else
                        accel_n = (/er * 0.0, er * 0.0, ez * cmr/)
                  end if
                  vel = vel + tstep*(accel + accel_n)/2.

                  accel = accel_n
            end if
      end do

      its = iter
      datas = data_points
      deallocate(pos)
      deallocate(vel)
      deallocate(last_rec)
      deallocate(accel)
      deallocate(accel_n)
      end subroutine

      end module
