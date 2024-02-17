#define MAX_TRAJECTORY_POINTS 1048576
C     ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE OTHER FILES

      module ion_optics
      use iso_c_binding, only: c_int, c_double
      contains

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


C     Calculate the electric field by computing the gradient of a point 
C           of a given 4x4x4 potential lattice.
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
      subroutine field_at2D(x, y, d, potential_map, dimensions,
     &      ex, ey)
C     Variable declarations:
      real(c_double), intent(in) :: x, y, d
      integer(c_int), dimension(2), intent(in) :: dimensions
      real(c_double), 
     &      dimension(dimensions(1), dimensions(2)), 
     &      intent(in) :: potential_map
      real(c_double), intent(out) :: ex, ey
      real(c_double), dimension(4, 4) :: potential
      real(c_double) :: tx, ty
      integer :: cx, cy, ix, iy, idx

C     Initialize potential array
      potential = 0.0

C     Cut out a 4x4 area of the potential maps and compute linear sum
      cx = FLOOR(x / d)
      cy = FLOOR(y / d)            ! Center coordinates
      ! Array slicing:
      do ix = 1,4
            do iy = 1,4
                  potential(ix, iy) =
     &                  potential_map(cx + ix - 2, cy + iy - 2)
            end do
      end do
      !potential = potential_map((cx-1):(cx+2), (cy-1):(cy+2))

C     Correct the coordinates' offset and compute the E field
      tx = 2.0 * d + x - DBLE(cx) * d
      ty = 2.0 * d + y - DBLE(cy) * d
      call potEfunc_unit_grid(potential, tx/d, ty/d, ex, ey)
      ex = ex / d
      ey = ey / d
      end subroutine


      function is_dead(dimensions, x, y, d)
      integer(c_int), dimension(2), intent(in) :: dimensions
      real(c_double), intent(in) :: x, y, d
      logical :: is_dead
            if        ((x < 2 * d                  ) 
     &            .or. (y < 2 * d                  )
     &            .or. (x > (dimensions(1) - 2) * d)
     &            .or. (y > (dimensions(2) - 2) * d)) then
                  is_dead = .true.
            else
                  is_dead = .false.
            end if
      end function

      function how_dead(dimensions, is_electrode, x, y, d)
      integer(c_int), dimension(2), intent(in) :: dimensions
      integer(c_int), dimension(dimensions(1), dimensions(2)), 
     &      intent(in):: is_electrode
      real(c_double), intent(in) :: x, y, d
      integer(c_int) :: how_dead
      how_dead = 0
      if        ((x < 2 * d                  ) 
     &      .or. (y < 2 * d                  )
     &      .or. (x > (dimensions(1) - 2) * d)
     &      .or. (y > (dimensions(2) - 2) * d)) then
            how_dead = 1
            return
      else
            if (is_electrode(NINT(x/d), NINT(y/d)).eq.1) then
                  how_dead = 2
                  return
            end if
      end if
      end function

C     Subroutine for Verlet integration of ion trajectory.
      subroutine ray_optics_2D(xx, yy, vxx, vyy, is_electrode,
     &                        potential_maps, voltages, dimensions,
     &                        n_electrodes, m, q, din, maxdist, maxt,
     &                        x_traj, y_traj, ts, its)
     & bind(c, name = "ray_optics_2D")
C     Variable declarations:
C     Dummy variables:
      real(c_double), intent(in) :: xx, yy, vxx, vyy, m, q, din
      real(c_double), intent(in) :: maxdist, maxt
      integer(c_int), intent(in) :: n_electrodes
      real(c_double), dimension(n_electrodes), intent(in) :: voltages
      integer(c_int), dimension(2), intent(in) :: dimensions
      real(c_double), 
     &      dimension(n_electrodes, dimensions(1), dimensions(2)), 
     &      intent(in) :: potential_maps
      real(c_double), dimension(MAX_TRAJECTORY_POINTS), intent(out)
     &      :: x_traj, y_traj, ts
      integer(c_int), dimension(dimensions(1), dimensions(2)),
     &      intent(in) :: is_electrode
      real(c_double), intent(out) :: its

C     Local variables
      real(c_double), 
     &      dimension(dimensions(1), dimensions(2)) :: potential_map
      real(c_double) :: x, y, vx, vy, t, mdist, mt, cmr,
     &      a, v, tv, ta, tstep, ex, ey, ex_new, ey_new, d
      integer :: idx, iter
      logical :: dead

C     Calculate the potential based on the given electrode voltages
      potential_map = 0.0;

      do idx = 1, n_electrodes
            do ix = 1, dimensions(1)
                  do iy = 1, dimensions(2)
                        potential_map(ix, iy) = potential_map(ix, iy)
     &                        + potential_maps(idx, ix, iy)
     &                        * voltages(idx)
                  end do
            end do
      end do
      
C     Unit conversions.  For simplicity, we work with SI units within 
C           this function.
      x = xx * 1.0e-3
      y = yy * 1.0e-3                           ! mm -> m
      vx = vxx * 1.0e+3
      vy = vyy * 1.0e+3                         ! mm/us -> m/s
      
      t = 0.0;
      d = din * 1.0e-3                          ! mm -> m
      mdist = maxdist * 1.0e-3
      mt = maxt * 1.0e-6                        ! us -> s

C     Charge-to-mass ratio
      cmr = ((1.602176634e-19) * q) / ((1.660539067e-27) * m);


C     Main loop of integration
      iter = 0
C     Check if particle is alive
      dead = is_dead(dimensions, x, y, d)
      if (dead) then
            t = mt;
      else
            call field_at2D(x, y, d, potential_map, dimensions,
     &            ex, ey)
      end if
      
      do while ((t < mt) .and. (iter < MAX_TRAJECTORY_POINTS))
            iter = iter + 1

C           Record current state before it gets modified by the 
C                 integration step
            x_traj(iter) = x  * 1.0e+3
            y_traj(iter) = y  * 1.0e+3
            ts(iter)     = t  * 1.0e+6

C           Integration timestep calculation
            a = SQRT(ex*ex + ey*ey) * cmr
            v = SQRT(vx*vx + vy*vy)

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


C           Check if particle is alive
            dead = is_dead(dimensions, x, y, d)
            if (dead) then
                  t = mt;
            end if
            call field_at2D(x, y, d, potential_map, dimensions,
     &            ex_new, ey_new)

            vx = vx + tstep * (ex_new + ex) * cmr / 2
            vy = vy + tstep * (ey_new + ey) * cmr / 2

            ex = ex_new
            ey = ey_new
      end do
      its = dble(iter)
      end subroutine


      subroutine ray_optics_2D_1k_points(xx, yy, vxx, vyy,
     &                        potential_maps, voltages, dimensions,
     &                        n_electrodes, m, q, din, maxdist, maxt,
     &                        x_traj, y_traj, ts, its, datas)
C     Variable declarations:
C     Dummy variables:
      real(c_double), intent(in) :: xx, yy, vxx, vyy, m, q, din
      real(c_double), intent(in) :: maxdist, maxt
      integer(c_int), intent(in) :: n_electrodes
      real(c_double), dimension(n_electrodes), intent(in) :: voltages
      integer(c_int), dimension(2), intent(in) :: dimensions
      real(c_double), 
     &      dimension(n_electrodes, dimensions(1), dimensions(2)), 
     &      intent(in) :: potential_maps
      real(c_double), dimension(1024), intent(out)
     &      :: x_traj, y_traj, ts
      real(c_double), intent(out) :: its, datas

C     Local variables
      real(c_double), 
     &      dimension(dimensions(1), dimensions(2)) :: potential_map
      real(c_double) :: x, y, vx, vy, t, mdist, mt, cmr,
     &      a, v, tv, ta, tstep, ex, ey, ex_new, ey_new, d
      integer :: idx, iter, sample_interval, data_points
      logical :: dead

      sample_interval = MAX_TRAJECTORY_POINTS / 1024

C     Calculate the potential based on the given electrode voltages
      potential_map = 0.0;

      do idx = 1, n_electrodes
            do ix = 1, dimensions(1)
                  do iy = 1, dimensions(2)
                        potential_map(ix, iy) = potential_map(ix, iy)
     &                        + potential_maps(idx, ix, iy)
     &                        * voltages(idx)
                  end do
            end do
      end do
      
C     Unit conversions.  For simplicity, we work with SI units within 
C           this function.
      x = xx * 1.0e-3
      y = yy * 1.0e-3                           ! mm -> m
      vx = vxx * 1.0e+3
      vy = vyy * 1.0e+3                         ! mm/us -> m/s
      
      t = 0.0;
      d = din * 1.0e-3                          ! mm -> m
      mdist = maxdist * 1.0e-3
      mt = maxt * 1.0e-6                        ! us -> s

C     Charge-to-mass ratio
      cmr = ((1.602176634e-19) * q) / ((1.660539067e-27) * m);


C     Main loop of integration
      iter = 0
C     Check if particle is alive
      dead = is_dead(dimensions, x, y, d)
      if (dead) then
            t = mt;
      end if
      call field_at2D(x, y, d, potential_map, dimensions,
     &      ex, ey)

      data_points = 1
      x_traj = 0.0
      y_traj = 0.0
      ts     = 0.0

      do while ((t < mt) .and. (iter < MAX_TRAJECTORY_POINTS))
            iter = iter + 1

C           Record current state before it gets modified by the 
C                 integration step
            if (MOD(iter - 1, sample_interval) .EQ. 0) then
                  x_traj(data_points) = x  * 1.0e+3
                  y_traj(data_points) = y  * 1.0e+3
                  ts(data_points)     = t  * 1.0e+6
                  data_points = data_points + 1
            end if
            
C           Integration timestep calculation
            a = SQRT(ex*ex + ey*ey) * cmr
            v = SQRT(vx*vx + vy*vy)

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


C           Check if particle is alive
            dead = is_dead(dimensions, x, y, d)
            if (dead) then
                  t = mt;
            end if

            call field_at2D(x, y, d, potential_map, dimensions,
     &            ex_new, ey_new)
            vx = vx + tstep * (ex_new + ex) * cmr / 2
            vy = vy + tstep * (ey_new + ey) * cmr / 2

            ex = ex_new
            ey = ey_new
      end do
      its = dble(iter)
      datas = dble(data_points - 1)
      end subroutine

      subroutine ray_optics_ensemble(particles, xxs, yys, vxxs, vyys,
     &                        potential_maps, voltagess, dimensions,
     &                        n_electrodes, ms, qs, din, maxdist, maxt,
     &                        x_trajs, y_trajs, tss, itss, datass)
     & bind(c, name="ray_optics_ensemble")
      integer(c_int), intent(in) :: particles, n_electrodes
      real(c_double), dimension(particles), intent(in)
     &      :: xxs, yys, vxxs, vyys, ms, qs
      real(c_double), intent(in) :: maxdist, maxt, din
      integer(c_int), dimension(2), intent(in) :: dimensions
      real(c_double), 
     &      dimension(n_electrodes, dimensions(1), dimensions(2)), 
     &      intent(in) :: potential_maps
      real(c_double), dimension(particles, n_electrodes), intent(in) 
     &      :: voltagess

      real(c_double), dimension(particles, 1024), intent(out)
     &      :: x_trajs, y_trajs, tss
      real(c_double), dimension(particles), intent(out)
     &      :: itss, datass

      real(c_double), dimension(1024) :: x_traj, y_traj, ts
      real(c_double) :: its, datas

      do idx = 1, particles
            call ray_optics_2D_1k_points(xxs(idx), yys(idx), 
     &                        vxxs(idx), vyys(idx),
     &                        potential_maps, voltagess(idx, :), 
     &                        dimensions, n_electrodes, 
     &                        ms(idx), qs(idx), 
     &                        din, maxdist, maxt,
     &                        x_traj, y_traj, ts, its, datas)
            do k = 1, 1024
                  x_trajs(idx, k) = x_traj(k)
                  y_trajs(idx, k) = y_traj(k)
                  tss(idx, k)     = ts(k)
            end do
            itss(idx) = its
            datass(idx) = datas
      end do

      end subroutine



      subroutine ray_optics_spaced(position, velocity, sample_dist,
     &            is_electrode, potential_maps, voltages, dimensions,
     &            n_electrodes, m, q, din, maxdist, maxt,
     &            trajectory, death, its, datas)
     & bind(c, name='ray_optics_spaced')
C     Variable declarations:
C     Dummy variables:
      real(c_double), dimension(2), intent(in) :: position, velocity
      real(c_double), intent(in) :: m, q, din, sample_dist
      real(c_double), intent(in) :: maxdist, maxt
      integer(c_int), intent(in) :: n_electrodes
      real(c_double), dimension(n_electrodes), intent(in) :: voltages
      integer(c_int), dimension(2), intent(in) :: dimensions
      integer(c_int), dimension(dimensions(1), dimensions(2)), 
     &      intent(in):: is_electrode
      real(c_double),
     &      dimension(n_electrodes, dimensions(1), dimensions(2)), 
     &      intent(in) :: potential_maps
      real(c_double), dimension(1024, 3), intent(out) :: trajectory
      integer(c_int), intent(out) :: its, datas, death

C     Local variables
      real(c_double), 
     &      dimension(dimensions(1), dimensions(2)) :: potential_map
      real(c_double) :: t, mdist, mt, cmr,
     &      a, v, tv, ta, tstep, ex, ey, ex_new, ey_new, d
      real(c_double), dimension(2) :: pos, vel, last_rec, accel, accel_n
      integer(c_int) :: idx, iter, sample_interval, data_points
      logical :: dead

      sample_interval = MAX_TRAJECTORY_POINTS / 1024

C     Calculate the potential based on the given electrode voltages
      potential_map = 0.0;

      do idx = 1, n_electrodes
            do ix = 1, dimensions(1)
                  do iy = 1, dimensions(2)
                        potential_map(ix, iy) = potential_map(ix, iy)
     &                        + (potential_maps(idx, ix, iy)
     &                        * voltages(idx))
                  end do
            end do
      end do
      
C     Unit conversions.  For simplicity, we work with SI units within 
C           this function.
      pos = position * 1.0e-3 ! mm -> m
      last_rec = pos
      vel = velocity * 1.0e+3 ! mm/us -> m/s
      
      t = 0.0;
      d = din * 1.0e-3                          ! mm -> m
      mdist = maxdist * 1.0e-3
      mt = maxt * 1.0e-6                        ! us -> s

C     Charge-to-mass ratio
      cmr = ((1.602176634e-19) * q) / ((1.660539067e-27) * m);


C     Main loop of integration
      iter = 0
C     Check if particle is alive
      death = how_dead(dimensions, is_electrode, pos(1), pos(2), d)
      if (death.gt.0) then
            t = mt;
      else
            call field_at2D(pos(1),pos(2), d, potential_map, dimensions,
     &            ex, ey)
      end if
      accel = (/ex, ey/) * cmr
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
                  trajectory(data_points, 1:2) = pos * 1.0e+3
                  trajectory(data_points, 3)   = t   * 1.0e+6
            end if
            
C           Integration timestep calculation
            a = NORM2(accel)
            v = NORM2(vel)

            a = a + 1.0e-15
            v = v + 1.0e-15
            tv = mdist / v
            ta = SQRT(2.0 * mdist / a)
            tstep = tv * ta / (tv + ta)
            tstep = MIN(tstep, mt * 1.0e-3)
            t = t + tstep

C           Integration step
            pos = pos + tstep*vel + tstep*tstep*accel/2.

C           Check if particle is alive
            death = how_dead(dimensions,is_electrode,pos(1),pos(2),d)
            if (death.gt.0) then
                  t = mt;
            end if

C           Calculate new fields
            call field_at2D(pos(1),pos(2), d, potential_map, dimensions,
     &            ex, ey)
            accel_n = (/ex, ey/) * cmr
            vel = vel + tstep*(accel + accel_n)/2.

            accel = accel_n
      end do
      its = iter
      datas = data_points
      end subroutine



      subroutine ray_optics_spaced_ensmble(particles, 
     &            positions, velocities,
     &            sample_dist, is_electrode, potential_maps, voltages, 
     &            dimensions, n_electrodes, m, q, din, maxdist, maxt,
     &            trajectories, deaths, itss, datass)
     & bind(c, name='ray_optics_spaced_ensemble')
C     Variable declarations:
C     Dummy variables:
      integer(c_int), intent(in) :: n_electrodes, particles
      real(c_double), dimension(2, particles), intent(in) 
     &      :: positions, velocities
      real(c_double), intent(in) :: m, q, din, sample_dist
      real(c_double), intent(in) :: maxdist, maxt
      real(c_double), dimension(n_electrodes), intent(in) :: voltages
      integer(c_int), dimension(2), intent(in) :: dimensions
      integer(c_int), dimension(dimensions(1), dimensions(2)), 
     &      intent(in):: is_electrode
      real(c_double),
     &      dimension(n_electrodes, dimensions(1), dimensions(2)), 
     &      intent(in) :: potential_maps
      real(c_double), dimension(particles, 1024 * 3), intent(out) 
     &      :: trajectories
      integer(c_int), dimension(particles), intent(out) 
     &      :: itss, datass, deaths

      integer(c_int) :: idx, its, datas, death
      real(c_double), dimension(1024, 3) :: trajectory

      !$OMP PARALLEL DO
      do idx = 1,particles
            trajectory = 0.
            call ray_optics_spaced(positions(:,idx), velocities(:,idx), 
     &            sample_dist, is_electrode, potential_maps,
     &            voltages, dimensions,
     &            n_electrodes, m, q, din, maxdist, maxt,
     &            trajectory, death, its, datas)
            trajectories(idx,:)
     &            = RESHAPE(trajectory, (/1024 * 3/))
            deaths(idx) = death
            itss(idx) = its
            datass(idx) = datas
      end do
      !$OMP END PARALLEL DO
      end subroutine

      end module















