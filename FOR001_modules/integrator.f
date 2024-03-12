#define DEBUG 0
#define STINTLENGTH 65536

#define PI 3.14159265358979323846
#define ECHARGE 1.602176634e-19
#define VACEPSILON 8.8541878128e-12
#define PROTONMASS 1.67262192369e-27
#define BOLTZMANN 1.380649e-23
#define ROOMTEMP 300.0


      module polyphemus_integrator
      use iso_c_binding, only: c_int, c_double, c_double_complex
      use iso_fortran_env, only: qp => real128
      implicit none
      contains


      subroutine inf_trap_field(position, t, omega, depth, R,
     &                          e_field)
      real(c_double), dimension(3), intent(in) :: position
      real(c_double), intent(in) :: t, depth, omega, R
      real(c_double), dimension(3), intent(out) :: e_field
      e_field(1) = -depth * position(1) * COS(omega * t) / (R * R)
      e_field(2) = -depth * position(2) * COS(omega * t) / (R * R)
      e_field(3) = 0.0e0
      end subroutine

      subroutine inf_trap_field_qp(position, t, omega, depth, R,
     &                          e_field)
      real(qp), dimension(3), intent(in) :: position
      real(qp), intent(in) :: t, depth, omega, R
      real(qp), dimension(3), intent(out) :: e_field
      e_field(1) = -depth * position(1) * COS(omega * t) / (R * R)
      e_field(2) = -depth * position(2) * COS(omega * t) / (R * R)
      e_field(3) = 0.
      end subroutine

      subroutine coulomb(posdiff, q1, q2, m1, m2,
     &                   a1, a2)
      real(c_double), dimension(3), intent(in)
     &      :: posdiff
      real(c_double), intent(in)
     &      :: q1, q2, m1, m2
      real(c_double), dimension(3), intent(inout)
     &      :: a1, a2
      real(c_double) :: coulomb_force, distance
      integer :: cdx
      distance = NORM2(posdiff)
      coulomb_force =
     &        (q1 * q2 
     &        / (VACEPSILON * 4.0 * PI))
     &        / (distance * distance)
      do cdx = 1,3
            a1(cdx) =  (coulomb_force / m1) 
     &                  * posdiff(cdx) / distance
            a2(cdx) = -(coulomb_force / m2) 
     &                  * posdiff(cdx) / distance
      end do
      end subroutine coulomb

      subroutine coulomb_qp(posdiff, q1, q2, m1, m2,
     &                   a1, a2)
      real(qp), dimension(3), intent(in)
     &      :: posdiff
      real(qp), intent(in)
     &      :: q1, q2, m1, m2
      real(qp), dimension(3), intent(inout)
     &      :: a1, a2
      real(qp) :: coulomb_force, distance
      integer :: cdx
      distance = NORM2(posdiff)
      coulomb_force =
     &        (q1 * q2 
     &        / (VACEPSILON * 4.0q0 * PI))
     &        / (distance * distance)
      do cdx = 1,3
            a1(cdx) =  (coulomb_force / m1) 
     &                  * posdiff(cdx) / distance
            a2(cdx) = -(coulomb_force / m2) 
     &                  * posdiff(cdx) / distance
      end do
      end subroutine coulomb_qp

      subroutine polyphemus(
     &    particles, 
     &    positions,
     &    velocities, 
     &    ms_in, 
     &    qs_in,
     &    omega, 
     &    depth, 
     &    R, 
     &    max_t,
     &    max_dist,
     &    record_step,
     &    burst_time,
     &    trajectory, 
     &    times, 
     &    its,
     &    recorded
     &    )
     & bind(c, name="nbody")
C     Argument specifications:
C           particles        : number of ions to simulate
C           positions        : 3x[particles] array: initial ( x, y, z)
C           velocities       : 3x[particles] array: initial (vx,vy,vz)
C           ms               : 1x[particles] array: masses (amu)
C           qs               : 1x[particles] array: charges (au)
C           omega            : RF drive frequency
C           depth            : Trap depth at distance R from trap center
C           max_t            : Time for cutting off the simulation
C           max_dist         : for setting adaptive timestep travel
C           record_step      : how often to record a trajectory point

C     Variable declarations
      integer(c_int), intent(in) :: particles
      real(c_double), dimension(particles, 3), intent(in) 
     &      :: positions, velocities
      real(c_double), dimension(particles), intent(in) 
     &      :: ms_in, qs_in
      real(c_double), intent(in) :: omega, depth, R
      real(c_double), intent(in) :: max_t, max_dist
      real(c_double), intent(in) :: record_step, burst_time
      real(c_double), dimension(STINTLENGTH, particles, 6), 
     &      intent(out) :: trajectory
      real(c_double), dimension(STINTLENGTH), 
     &      intent(out) :: times
      integer(c_int), intent(out) :: its, recorded

C     Local variables
      real(c_double), allocatable 
     &      :: pos(:, :),   vel(:, :), vs(:, :), 
     &         accel(:, :), accel_new(:, :), as(:, :)
      real(c_double), dimension(3) :: field
      real(c_double), dimension(particles)
     &      :: ms, qs, cmrs
      real(c_double) :: scalar_a, scalar_v
      logical, dimension(particles) :: dead
      real(c_double) :: t, tstep, last_rec_t, tstep_candidate
      real(c_double) :: md, mt, rt
      integer(c_int) :: rec_idx, alive
      integer(c_int) :: idx, jdx, cdx

      allocate(pos(particles, 3))
      allocate(vel(particles, 3))
      allocate(accel(particles, 3))
      allocate(accel_new(particles, 3))
      allocate(as(particles, 3))
      allocate(vs(particles, 3))

      if (DEBUG .eq. 1) then
            open (10, file='output_file.txt',
     &            status='unknown')
      endif

      trajectory = 0.
      times = 0.

      md  = max_dist * 1.0e-3                         ! mm -> m
      mt  = max_t * 1.0e-6                            ! us -> s
      rt  = R * 1.0e-3                                ! mm -> m
      do idx = 1,particles
            do cdx = 1,3
                  pos(idx, cdx) = positions(idx, cdx)* 1.0e-3  ! mm -> m
                  vel(idx, cdx) = velocities(idx, cdx) * 1.0e+3! mm/us -> m/s
            end do
      end do
      mt = max_t * 1.0e-6                             ! us -> s
      md = max_dist * 1.0e-3                          ! mm -> m
      ms = ms_in * PROTONMASS                         ! au -> kg
      qs = qs_in * ECHARGE                            ! au -> C
      do idx = 1,particles
            cmrs(idx) = qs(idx) / ms(idx)             ! C/kg
      end do


C     Initial accelerations (PHASE 1)
      accel     = 0.
      accel_new = 0.
      field     = 0.
      do idx = 1,particles
            call inf_trap_field(pos(idx, :), t, 
     &                          omega, depth, rt, 
     &                          field)
            do cdx = 1,3
                  accel_new(idx, cdx) = field(cdx) * cmrs(idx)
            end do
C             ! Coulomb interaction, Newton's third law
C             do jdx = 1,(idx - 1)
C                   call coulomb(pos(idx, :) - pos(jdx, :), 
C      &                         qs(idx), qs(jdx), 
C      &                         ms(idx), ms(jdx), 
C      &                         accel_new(idx, :), 
C      &                         accel_new(jdx, :))
C             end do
      end do

      if (DEBUG.eq.1) then
            write(10,*) "pos:"
            write(10,*) pos
            write(10,*) "vel:"
            write(10,*) vel
            write(10,*) "field:"
            write(10,*) field
            write(10,*) "accel_new"
            write(10,*) accel_new
      end if

      its = 0
      rec_idx = 1
      t = 0.
      last_rec_t = 0.
C     MAIN INTEGRATION LOOP:
C     Velocity-Verlet
C     PHASE 0: Record current state
C     PHASE 1: Calculate/Inherit accelerations
C     PHASE 2: Adaptive timestep
C     PHASE 3: Half-timestep propagation of position
C     PHASE 4: Calculate new acceleration
C     PHASE 5: Half-timestep propagation of velocities mean acceleration
C     Subroutine for Verlet integration of ion trajectory.
      do while ((t < mt)
     &     .and.(rec_idx.le.STINTLENGTH))
            its = its + 1

C           PHASE 0: Record current state
            if (last_rec_t.le.t) then
                  trajectory(rec_idx, :, 1:3) = pos
                  trajectory(rec_idx, :, 4:6) = vel
                  times(rec_idx) = t
                  rec_idx = rec_idx + 1

                  if ((last_rec_t + (burst_time * 1.0e-6)).le.t) then
                        last_rec_t =
     &                      MIN(last_rec_t + (record_step * 1.0e-6), mt)
                  end if
            end if

C           PHASE 1: Calculate/Inherit accelerations
            do idx = 1,particles
                  do cdx = 1,3
                        accel(idx, cdx) = accel_new(idx, cdx)
                  end do
            end do

C           PHASE 2: Adaptive timestep
            tstep = mt / DBLE(STINTLENGTH)
            do idx = 1,particles
                  scalar_a = NORM2(accel(idx, :))
                  scalar_v = NORM2(  vel(idx, :))
                  scalar_a = scalar_a + 1.0e-15
                  scalar_v = scalar_v + 1.0e-15
                  scalar_a = SQRT(2.0 * md / scalar_a) ! tas
                  scalar_v = md / scalar_v             ! tvs
                  tstep_candidate = scalar_v * scalar_a
     &                               / (scalar_v + scalar_a)
                  if (DEBUG .eq. 1) then
                        write(10,*) "scalar_a, scalar_v, tstep"
                        write(10,*) scalar_a, scalar_v, tstep_candidate
                  endif
                  tstep = MIN(tstep, tstep_candidate)
            end do
            t = t + tstep

            ! Since some are left, we need to calculate the fields for:
C           PHASE 4: Calculate new acceleration
            do idx = 1,particles
C                 PHASE 3: Half-timestep propagation of position
                  do cdx = 1,3
                        pos(idx,cdx) = pos(idx,cdx) 
     &                        + (vel(idx,cdx) * tstep) 
     &                        + tstep * tstep * accel(idx,cdx) / 2.
                  end do
                  call inf_trap_field(pos(idx, :), t, 
     &                                omega, depth, rt, 
     &                                field)
                  accel_new(idx, :) = field * cmrs(idx)
                  ! Coulomb interaction, Newton's third law
                  do jdx = 1,(idx - 1)
                        call coulomb(pos(idx, :) - pos(jdx, :), 
     &                               qs(idx), qs(jdx), 
     &                               ms(idx), ms(jdx), 
     &                               accel_new(idx, :), 
     &                               accel_new(jdx, :))
                  end do
            end do
C           PHASE 5: Half-timestep propagation of velocities
            do idx = 1,particles
                  do cdx = 1,3
                        vel(idx, cdx) = vel(idx, cdx) 
     &                                  + tstep * (accel(idx, cdx) 
     &                                  + accel_new(idx, cdx)) / 2.
                  end do
            end do
      end do
 10   continue
      times = times * 1.0e+6
      trajectory(:, :, 1:3) = trajectory(:, :, 1:3) * 1.0e+3
      trajectory(:, :, 4:6) = trajectory(:, :, 4:6) * 1.0e-3
      recorded = rec_idx - 1

      if (DEBUG .eq. 1) then
            close (10)
      endif
      end subroutine

      subroutine calculate_accelerations_qp(
     &    particles, 
     &    cmrs,
     &    ms,
     &    qs,
     &    positions,
     &    t, 
     &    omega, 
     &    depth, 
     &    R, 
     &    accelerations
     &    )
      integer(c_int), intent(in) :: particles
      real(qp), dimension(particles, 3), intent(in) 
     &      :: positions
      real(qp), dimension(particles), intent(in) 
     &      :: cmrs, ms, qs
      real(qp), intent(in) :: t, omega, depth, R
      real(qp), dimension(particles, 3), intent(out)
     &      :: accelerations
      real(qp), dimension(3) :: field
      integer :: idx, jdx, cdx
      do idx = 1,particles
            call inf_trap_field_qp(positions(idx, :), t, 
     &                          omega, depth, R, 
     &                          field)
            do cdx = 1,3
                  accelerations(idx, cdx) = field(cdx) * cmrs(idx)
            end do
            ! Coulomb interaction, Newton's third law
            do jdx = 1,(idx - 1)
                  call coulomb_qp(positions(idx, :) - positions(jdx, :), 
     &                         qs(idx), qs(jdx), 
     &                         ms(idx), ms(jdx), 
     &                         accelerations(idx, :), 
     &                         accelerations(jdx, :))
            end do
      end do
      end subroutine

      subroutine polyphemus4(
     &    particles, 
     &    positions,
     &    velocities, 
     &    ms_in, 
     &    qs_in,
     &    omega_in, 
     &    depth_in, 
     &    R, 
     &    max_t,
     &    max_dist,
     &    record_step_in,
     &    burst_time_in,
     &    trajectory, 
     &    times, 
     &    its,
     &    recorded
     &    )
     & bind(c, name="nbody4")
C     Argument specifications:
C           particles        : number of ions to simulate
C           positions        : 3x[particles] array: initial ( x, y, z)
C           velocities       : 3x[particles] array: initial (vx,vy,vz)
C           ms               : 1x[particles] array: masses (amu)
C           qs               : 1x[particles] array: charges (au)
C           omega            : RF drive frequency
C           depth            : Trap depth at distance R from trap center
C           max_t            : Time for cutting off the simulation
C           max_dist         : for setting adaptive timestep travel
C           record_step      : how often to record a trajectory point

C     Variable declarations
      integer(c_int), intent(in) :: particles
      real(c_double), dimension(particles, 3), intent(in) 
     &      :: positions, velocities
      real(c_double), dimension(particles), intent(in) 
     &      :: ms_in, qs_in
      real(c_double), intent(in) :: omega_in, depth_in, R
      real(c_double), intent(in) :: max_t, max_dist
      real(c_double), intent(in) :: record_step_in, burst_time_in
      real(c_double), dimension(STINTLENGTH, particles, 6), 
     &      intent(out) :: trajectory
      real(c_double), dimension(STINTLENGTH), 
     &      intent(out) :: times
      integer(c_int), intent(out) :: its, recorded

C     Local variables
      real(qp), allocatable 
     &      :: pos(:, :),   vel(:, :), vs(:, :), 
     &         accel(:, :), accel_new(:, :), as(:, :)
      real(qp), dimension(3) :: field
      real(qp), dimension(particles)
     &      :: ms, qs, cmrs
      real(qp) :: scalar_a, scalar_v
      logical, dimension(particles) :: dead
      real(qp) :: t, tstep, last_rec_t, tstep_candidate
      real(qp) :: md, mt, rt, burst_time, record_step
      real(qp) :: omega, depth
      integer(c_int) :: rec_idx, alive
      integer(c_int) :: idx, jdx, cdx
      real(qp) :: w0, w1, c1, c2, c3, c4, d1, d2, d3, d4

      w0 = -1.7024143839193152680953756179429217q0
      w1 =  1.3512071919596576340476878089714608q0
      c1 = w1 / 2.0q0
      c2 = (w0 + w1) / 2.0q0
      c3 = (w0 + w1) / 2.0q0
      c4 = w1 / 2.0q0
      d1 = w1
      d2 = w0
      d3 = w1
      d4 = 0.0q0
C       c1 = 0.0
C       c2 = 1.0
C       d1 = 0.5
C       d2 = 0.5


      allocate(pos(particles, 3))
      allocate(vel(particles, 3))
      allocate(accel(particles, 3))
      allocate(accel_new(particles, 3))
      allocate(as(particles, 3))
      allocate(vs(particles, 3))

      if (DEBUG .eq. 1) then
            open (10, file='output_file.txt',
     &            status='unknown')
      endif

      trajectory = 0.
      times = 0.

      md  = max_dist * 1.0q-3                         ! mm -> m
      mt  = max_t * 1.0q-6                            ! us -> s
      rt  = R * 1.0q-3                                ! mm -> m
      pos = positions * 1.0q-3                        ! mm -> m
      vel = velocities * 1.0q+3                       ! mm/us -> m/s
      mt = max_t * 1.0q-6                             ! us -> s
      md = max_dist * 1.0q-3                          ! mm -> m
      ms = ms_in * PROTONMASS * 1.0q0                 ! au -> kg
      qs = qs_in * ECHARGE    * 1.0q0                 ! au -> C
      cmrs = qs / ms                                  ! C/kg
      omega = omega_in * 1.0q0
      depth = depth_in * 1.0q0
      burst_time = burst_time_in * 1.0q-6
      record_step = record_step_in * 1.0q-6

      its = 0
      rec_idx = 1

      t = 0.
      last_rec_t = 0.
      dead = .false.
      alive = particles

      call calculate_accelerations_qp(
     &      particles,
     &      cmrs,
     &      ms,
     &      qs,
     &      pos,
     &      t, 
     &      omega, 
     &      depth, 
     &      rt, 
     &      accel
     &      )

C     MAIN INTEGRATION LOOP:
C     4th-order Yoshida
C     PHASE 0: Record current state
C     PHASE 1: Calculate/Inherit accelerations
C     PHASE 2: Adaptive timestep
C     PHASE 3: Propagate x1, v1
C     PHASE 4: Calculate a(x1)
C     PHASE 5: Propagate x2, v2
C     PHASE 6: Calculate a(x2)
C     PHASE 7: Propagate x3, v3
C     PHASE 8: Calculate a(x3)
C     PHASE 9: Calculate new x, v

      do while ((t < mt)
     &      .and.(rec_idx.le.STINTLENGTH))
            its = its + 1

C           PHASE 0: Record current state
            if (last_rec_t.le.t) then
                  trajectory(rec_idx, :, 1:3) = pos
                  trajectory(rec_idx, :, 4:6) = vel
                  times(rec_idx) = t
                  rec_idx = rec_idx + 1

                  if ((last_rec_t + burst_time).le.t) then
                        last_rec_t =
     &                      MIN(last_rec_t + record_step, mt)
                  end if
            end if

C           PHASE 2: Adaptive timestep
            tstep = mt / DBLE(STINTLENGTH)
            do idx = 1,particles
                  scalar_a = NORM2(accel(idx, :))
                  scalar_v = NORM2(  vel(idx, :))
                  scalar_a = scalar_a + 1.0q-15
                  scalar_v = scalar_v + 1.0q-15
                  scalar_v = md / scalar_v             ! tvs
                  scalar_a = SQRT(2.0 * md / scalar_a) ! tas
                  tstep_candidate = scalar_v * scalar_a
     &                               / (scalar_v + scalar_a)
                  tstep = MIN(tstep, tstep_candidate)
            end do
            t = t + tstep

C           PHASE 3: Propagate x1, v1
            ! ORDER 1
            do idx = 1,particles
                  do cdx = 1,3
                        pos(idx,cdx) = pos(idx, cdx) 
     &                        + c1 * vel(idx,cdx)    * tstep
                  end do
            end do
            call calculate_accelerations_qp(
     &            particles,
     &            cmrs,
     &            ms,
     &            qs,
     &            pos,
     &            t, 
     &            omega, 
     &            depth, 
     &            rt, 
     &            accel
     &            )
            do idx = 1,particles
                  do cdx = 1,3
                        vel(idx, cdx) = vel(idx, cdx) 
     &                        + d1 * accel(idx, cdx) * tstep
                  end do
            end do

            ! ORDER 2
            do idx = 1,particles
                  do cdx = 1,3
                        pos(idx,cdx) = pos(idx, cdx) 
     &                        + c2 * vel(idx,cdx)    * tstep
                  end do
            end do
            call calculate_accelerations_qp(
     &            particles,
     &            cmrs,
     &            ms,
     &            qs,
     &            pos,
     &            t, 
     &            omega, 
     &            depth, 
     &            rt, 
     &            accel
     &            )
            do idx = 1,particles
                  do cdx = 1,3
                        vel(idx, cdx) = vel(idx, cdx) 
     &                        + d2 * accel(idx, cdx) * tstep
                  end do
            end do

            ! ORDER 3
            do idx = 1,particles
                  do cdx = 1,3
                        pos(idx,cdx) = pos(idx, cdx) 
     &                        + c3 * vel(idx,cdx)    * tstep
                  end do
            end do
            call calculate_accelerations_qp(
     &            particles,
     &            cmrs,
     &            ms,
     &            qs,
     &            pos,
     &            t, 
     &            omega, 
     &            depth, 
     &            rt, 
     &            accel
     &            )
            do idx = 1,particles
                  do cdx = 1,3
                        vel(idx, cdx) = vel(idx, cdx) 
     &                        + d3 * accel(idx, cdx) * tstep
                  end do
            end do

            ! ORDER 4
            do idx = 1,particles
                  do cdx = 1,3
                        pos(idx,cdx) = pos(idx, cdx) 
     &                        + c4 * vel(idx,cdx)    * tstep
                  end do
            end do
      end do
 10   continue
      times = times * 1.0e+6
      trajectory(:, :, 1:3) = trajectory(:, :, 1:3) * 1.0e+3
      trajectory(:, :, 4:6) = trajectory(:, :, 4:6) * 1.0e-3
      recorded = rec_idx - 1

      if (DEBUG .eq. 1) then
            close (10)
      endif
      end subroutine

      end module

