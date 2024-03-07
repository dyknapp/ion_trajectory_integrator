#define DEBUG 1
#define STINTLENGTH 65536

#define PI 3.14159265358979323846
#define ECHARGE 1.602176634e-19
#define VACEPSILON 8.8541878128e-12
#define PROTONMASS 1.67262192369e-27
#define BOLTZMANN 1.380649e-23
#define ROOMTEMP 300.0


      module multipole_integration
      use iso_c_binding, only: c_int, c_double, c_double_complex
      implicit none
      contains


      subroutine inf_trap_field(position, t, omega, depth, R,
     &                          e_field)
      real(c_double), dimension(3), intent(in) :: position
      real(c_double), intent(in) :: t, depth, omega, R
      real(c_double), dimension(3), intent(out) :: e_field
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
      distance = NORM2(posdiff)
      coulomb_force =
     &        (q1 * q2 * (ECHARGE**2.0) 
     &        / ( VACEPSILON * 4.0 * PI ))
     &        / (distance * distance)
      a1 =  (coulomb_force / m1) * posdiff / distance
      a2 = -(coulomb_force / m2) * posdiff / distance
      end subroutine coulomb

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
     &    its)
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
      real(c_double), dimension(STINTLENGTH, particles, 3), 
     &      intent(out) :: trajectory
      real(c_double), dimension(STINTLENGTH), 
     &      intent(out) :: times
      integer(c_int), intent(out) :: its

C     Local variables
      real(c_double), allocatable 
     &      :: pos(:, :),   vel(:, :), vs(:, :), 
     &         accel(:, :), accel_new(:, :), as(:, :)
      real(c_double), allocatable 
     &      :: ms(:), qs(:), cmrs(:), field(:), 
     &         scalar_as(:), scalar_vs(:)
      logical, allocatable :: dead(:)
      real(c_double) :: t, tstep, last_rec_t
      real(c_double) :: md, mt, rt
      integer(c_int) :: rec_idx, alive
      integer(c_int) :: idx, jdx

      allocate(pos(particles, 3))
      allocate(vel(particles, 3))
      allocate(accel(particles, 3))
      allocate(accel_new(particles, 3))

      allocate(field(3))

      allocate(as(particles, 3))
      allocate(vs(particles, 3))

      allocate(ms(particles))
      allocate(qs(particles))
      allocate(cmrs(particles))
      allocate(scalar_as(particles))
      allocate(scalar_vs(particles))

      allocate(dead(particles))

      trajectory = 0.

      md  = max_dist * 1.0e-3                         ! mm -> m
      mt  = max_t * 1.0e-6                            ! us -> s
      rt  = R * 1.0e-3                                ! mm -> m
      pos = positions * 1.0e-3                        ! mm -> m
      vel = velocities * 1.0e+3                       ! mm/us -> m/s
      mt = max_t * 1.0e-6                             ! us -> s
      md = max_dist * 1.0e-3                          ! mm -> m
      ms = ms_in * PROTONMASS                         ! au -> kg
      qs = qs_in * ECHARGE                            ! au -> C
      cmrs = qs / ms                                  ! C/kg

      its = 0
      rec_idx = 1

      t = 0.
      last_rec_t = 0.
      dead = .false.
      alive = particles

      do idx = 1,particles
            if (.not. dead(idx)) then
                  call inf_trap_field(pos, t, omega, depth, rt, field)
                  accel(idx, :) = field * cmrs(idx)
            end if
      end do
      do while (((t < mt)
     &      .and.(rec_idx.le.STINTLENGTH))
     &      .and.(alive.gt.0))
            its = its + 1

            if (last_rec_t.le.t) then
                  trajectory(rec_idx, :, 1:3) = pos
                  times(rec_idx) = t
                  rec_idx = rec_idx + 1

                  if ((last_rec_t + burst_time).gt.t) then
                        last_rec_t = last_rec_t + burst_time
                  end if
            end if

            do idx = 1,particles
                  scalar_as(idx) = NORM2(accel(idx, :))
                  scalar_vs(idx) = NORM2(  vel(idx, :))
            end do
            scalar_as = scalar_as + 1.0e-15
            scalar_vs = scalar_vs + 1.0e-15
            scalar_vs = md / scalar_vs              ! tvs
            scalar_as = SQRT(2.0 * md / scalar_as)  ! tas
            tstep = MINVAL(scalar_vs * scalar_as 
     &            / (scalar_vs + scalar_as))
            tstep = MIN(tstep, mt * 1.0e-6)
            t = t + tstep

            do idx = 1,particles
                  if (.not. dead(idx)) then
C           PHASE 5: Half-timestep propagation of velocities,
            ! This is from the previous step, to avoid a second loop
                      vel = vel + tstep * (accel + accel_new) / 2.
                      accel = accel_new
                      
C           PHASE 3: Half-timestep propagation of position
                      pos = pos + (vel * tstep) 
     &                        + tstep*tstep * accel / 2.
                  end if
            end do

            ! Since some are left, we need to calculate the fields for:
C           PHASE 4: Calculate new acceleration
            alive = 0
            do idx = 1,particles
                  if (.not. dead(idx)) then
                        alive = alive + 1

                        call inf_trap_field(pos, t, omega, depth, rt, 
     &                                          field)
                        accel_new(idx, :) = field * cmrs(idx)
                        ! Coulomb interaction, Newton's third law
                        do jdx = 1,(idx - 1)
                              call coulomb(pos(idx, :) - pos(jdx, :), 
     &                                     qs(idx), qs(jdx), 
     &                                     ms(idx), ms(jdx), 
     &                                     as(idx, :), as(jdx, :))
                        end do
                  end if
            end do            
      end do
      end subroutine

      end module

