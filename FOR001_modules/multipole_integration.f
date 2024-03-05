#define DEBUG 0
#define PI 3.14159265358979323846
#define ECHARGE 1.602176634e-19
#define VACEPSILON 8.8541878128e-12
#define PROTONMASS 1.67262192369e-27
#define BOLTZMANN 1.380649e-23
#define ROOMTEMP 300.0

      module multipole_integration
      use iso_c_binding, only: c_int, c_double, c_double_complex
      contains

      function legendre(l, m, x) bind(c, name="legendre_f")
      ! Calculate Legendre polynomial
      ! Taken from F77 numerical recipes, section 6.8
C     Input variables
      integer(c_int), intent(in) :: l, m
      real(c_double), intent(in) :: x
      real(c_double) :: legendre

C     Local variables
      integer(c_int) :: i, ll
      real(c_double) :: fact, pll, pmm, pmmp1, somx2

      pmm = 1
      if (m.gt.0) then
            somx2 = SQRT((1.-x)*(1.+x))
            fact = 1.
            do i = 1,m
                  pmm = -pmm * fact * somx2
                  fact = fact + 2.
            end do
      end if
      if (l.eq.m) then
            legendre = pmm
      else
            pmmp1 = x * (2.*m + 1.) * pmm
            if (l.eq.(m + 1)) then
                  legendre = pmmp1
            else
                  do ll = (m + 2), l
                        pll = (x*(2.*ll-1.)*pmmp1-(ll+m-1.)*pmm)/(ll-m)
                        pmm = pmmp1
                        pmmp1 = pll
                  end do
                  legendre = pll
            end if
      end if
      return
      end function

      function Ylm(l, m, theta, phi)
      ! Calculate spherical harmonic
C     Input variables
      integer(c_int), intent(in) :: l, m
      real(c_double), intent(in) :: phi, theta
      complex(c_double_complex) :: Ylm, lc, mc, phic, thetac
      real(c_double) :: fact_ratio
      integer(c_int) :: sv, ev, factor

      sv = MIN(l - m, l + m) + 1
      ev = MAX(l - m, l + m)
      fact_ratio = 1.
      if (ev.gt.sv) then
            do factor = sv,ev
                  fact_ratio = fact_ratio * DBLE(factor)
            end do
      end if
      if ((l - m).lt.(l + m)) then
            fact_ratio = 1. / fact_ratio
      end if

      lc = CMPLX(DBLE(l), 0.);
      mc = CMPLX(DBLE(m), 0.);
      phic   = CMPLX(phi,   0.);
      thetac = CMPLX(theta, 0.);

      Ylm = EXP(CMPLX(0., 1.) * mc * phic)
      Ylm = Ylm
     &    * SQRT((CMPLX(2., 0.)*lc+CMPLX(1., 0.))*CMPLX(fact_ratio,0.)
     &         / CMPLX(4.*PI, 0.))
      Ylm = Ylm * CMPLX(legendre(l, m, COS(theta)), 0.)
      return
      end function

      subroutine Ylm_wrapper(l, m, phi, theta, r, c)
     & bind(c, name="Ylm_f")
      integer(c_int), intent(in) :: l, m
      real(c_double), intent(in) :: phi, theta
      real(c_double), intent(out) :: r, c
      complex(c_double_complex) :: Y

      Y = Ylm(l, m, theta, phi)
      r = REAL(Y)
      c = AIMAG(Y)
      end subroutine

C       function multipole_force()

C       end function

C       subroutine integrate_multipole(position, velocity,
C      &                               multipoles, osci_multipoles, omega,
C      &                               n_moments, max_t,
C      &                               max_dist, averaging_window,
C      &                               pts, four_trajectory, its)
C      & bind(c, name = "integrate_multipole")
C C     Variable declarations
C       integer(c_int), intent(in) :: pts, n_moments, averaging_window
C       real(c_double), dimension(3), intent(in) :: 
C      &      position_sph, velocity_sph
C       real(c_double), dimension(n_moments) :: 
C      &      multipoles, osci_multipoles
C       real(c_double), intent(in) :: omega
C       real(c_double), intent(in) :: max_t, max_dist

C C     Output variables
C       integer(c_int), intent(out) :: its
C       real(c_double), dimension(pts, 4), intent(out) :: four_trajectory

C C     Local variables
C       real(c_double) :: cmr, md, mt
C       real(c_double) :: t
C       real(c_double), dimension(3) :: pos, vel

C C     The correct units
C       cmr = ((ECHARGE) * q) / ((PROTONMASS) * m);
C       pos = position_sph
C       vel = velocity_sph
C       md  = max_dist * 1.0e-3 ! mm -> m
C       mt  = max_t * 1.0e-6    ! us -> s

C C     Main loop
C       its = 0
C       t = 0
C       do while (t < mt)

C       end do
C       end subroutine

C       subroutine inf_trap_field(position, t, omega, depth, R,
C      &                          e_field)
C       real(c_double), dimension(3), intent(in) :: position
C       real(c_double), intent(in) :: t, depth, omega, R
C       real(c_double), dimension(3), intent(out) :: e_field
C       e_field(1) = -depth * position(1) * COS(omega * t) / (R * R)
C       e_field(2) = -depth * position(2) * COS(omega * t) / (R * R)
C       e_field(3) = 0.
C       end subroutine

C       subroutine inf_trap_bg_gas(position, velocity, m, q,
C      &              omega, depth, R, max_t, max_dist, averaging_window,
C      &              pts, four_trajectory, its, rpts)
C      & bind(c, name="inf_trap_bg_gas")
C C     Variable declarations
C       integer(c_int), intent(in) :: pts, averaging_window
C       real(c_double), dimension(3), intent(in) :: position, velocity
C       real(c_double), intent(in) :: m, q
C       real(c_double), intent(in) :: omega, depth, R
C       real(c_double), intent(in) :: max_t, max_dist
C       real(c_double), dimension(pts, 4), intent(out) :: four_trajectory
C       integer(c_int), intent(out) :: its, rpts

C C     Local variables
C       real(c_double), dimension(3) :: pos, vel, accel, field, accel_new
C       real(c_double) :: t, ta, tv, a, v, sample_interval, next_sample
C       real(c_double) :: mt, md, rt
C       integer(c_int) :: window_end, sample_idx

C       four_trajectory = 0.

C       cmr = ((ECHARGE) * q) / ((PROTONMASS) * m);
C       md  = max_dist * 1.0e-3 ! mm -> m
C       mt  = max_t * 1.0e-6    ! us -> s
C       rt  = R * 1.0e-3
C       pos = position * 1.0e-3
C       vel = velocity * 1.0e+3
C       mt = max_t * 1.0e-6
C       md = max_dist * 1.0e-3
C       sample_interval = mt / DBLE(pts)

C       its = 0
C       t = 0.
C       next_sample = 0.
C       window_end = 0
C       sample_idx = 0
C       call inf_trap_field(pos, t, omega, depth, rt, field)
C       accel = field * cmr
C       do while ((t < mt).and.(its.lt.1e+10))
C             its = its + 1

C             if (its.le.window_end) then
C                   four_trajectory(sample_idx, 1:3) = 
C      &                      four_trajectory(sample_idx, 1:3) + pos
C                   four_trajectory(sample_idx, 4) =
C      &                      four_trajectory(sample_idx, 4) + t
C                   if (its.eq.window_end) then
C                         four_trajectory(sample_idx, 1:4) = 
C      &                      four_trajectory(sample_idx, 1:4)
C      &                    / DBLE(averaging_window)
C                   end if
C             else
C                   if (next_sample.le.t) then
C                         next_sample = t + sample_interval
C                         window_end = its + averaging_window
C                         sample_idx = sample_idx + 1
C                   end if
C             end if

C             a = NORM2(accel)
C             v = NORM2(vel)
C             a = a + 1.0e-15
C             v = v + 1.0e-15
C             tv = md / v
C             ta = SQRT(2.0 * md / a)
C             tstep = tv * ta / (tv + ta)
C             tstep = MIN(tstep, mt * 1.0e-3)
C             t = t + tstep

C             pos = pos + (vel * tstep) + tstep*tstep * accel / 2.

C             call inf_trap_field(pos, t, omega, depth, rt, field)
C             accel_new = field * cmr

C             vel = vel + tstep * (accel + accel_new) / 2.
C             accel = accel_new
C       end do
C       rpts = sample_idx - 1
C       end subroutine

      subroutine inf_trap_field(position, t, omega, depth, R,
     &                          e_field)
      real(c_double), dimension(3), intent(in) :: position
      real(c_double), intent(in) :: t, depth, omega, R
      real(c_double), dimension(3), intent(out) :: e_field
      e_field(1) = -depth * position(1) * COS(omega * t) / (R * R)
      e_field(2) = -depth * position(2) * COS(omega * t) / (R * R)
      e_field(3) = 0.
      end subroutine

      function maxwell_distribution(cdf, vs)
      real(c_double), dimension(1024), intent(in) :: cdf, vs
      real(c_double) :: maxwell_distribution
      real(c_double) :: random
      integer(c_int) :: idx
      call RANDOM_NUMBER(random)
      idx = 1
      do while (random.lt.cdf(idx))
            idx = idx + 1
      end do
      maxwell_distribution = vs(idx)
      end function

      subroutine maxwell_distribution_test(room_temp, samples, v)
     & bind(c, name="maxwell_distribution_test")
      integer(c_int), intent(in) :: samples
      real(c_double), intent(in) :: room_temp
      real(c_double), intent(out) :: v

C     Local variables
      real(c_double), dimension(1024) :: cdf, vs
      real(c_double) :: dv
      integer(c_int) :: idx

C     Initialize the CDF for the Maxwell-Boltzmann
      vs(1) = 0.
      cdf(1) = 0.
      dv = 100.*SQRT(2.*BOLTZMANN*room_temp/(m*PROTONMASS))/1024.0
      do idx = 2,1024
            vs(idx) = dv + vs(idx - 1)
            cdf(idx) = cdf(idx - 1) + dv * EXP(-m * PROTONMASS 
     &            * vs(idx-1)*vs(idx-1) / (2*BOLTZMANN*room_temp))
      end do
      cdf = cdf / cdf(1024)
      end subroutine

      subroutine inf_trap_bg_gas(position, velocity, m, q, coll_rate,
     &    room_temp, omega, depth, R, max_t, max_dist, averaging_window,
     &    pts, four_trajectory, its, rpts, collisions)
     & bind(c, name="inf_trap_bg_gas")
C     Variable declarations
      integer(c_int), intent(in) :: pts, averaging_window
      real(c_double), dimension(3), intent(in) :: position, velocity
      real(c_double), intent(in) :: m, q, coll_rate, room_temp
      real(c_double), intent(in) :: omega, depth, R
      real(c_double), intent(in) :: max_t, max_dist
      real(c_double), dimension(pts, 4), intent(out) :: four_trajectory
      integer(c_int), intent(out) :: its, rpts, collisions

C     Local variables
      real(c_double), dimension(3) :: pos, vel, accel, field, accel_new,
     &                                vc, normal
      real(c_double) :: t, ta, tv, a, v, sample_interval, next_sample,
     &   mt, md, rt, last_coll, expected_free, random, coll_random
      real(c_double) :: vm, theta, phi
      integer(c_int) :: window_end, sample_idx, idx
      real(c_double), dimension(1024) :: cdf, vs

      four_trajectory = 0.

      cmr = ((ECHARGE) * q) / ((PROTONMASS) * m);
      md  = max_dist * 1.0e-3 ! mm -> m
      mt  = max_t * 1.0e-6    ! us -> s
      rt  = R * 1.0e-3
      pos = position * 1.0e-3
      vel = velocity * 1.0e+3
      mt = max_t * 1.0e-6
      md = max_dist * 1.0e-3
      expected_free = 1. / (coll_rate * 1.0e+6);
      sample_interval = mt / DBLE(pts)

C     Initialize the CDF for the Maxwell-Boltzmann
      vs(1) = 0.
      cdf(1) = 0.
      dv = 100.*SQRT(2.*BOLTZMANN*room_temp/(m*PROTONMASS))/1024.0
      do idx = 2,1024
            vs(idx) = dv + vs(idx - 1)
            cdf(idx) = cdf(idx - 1) + dv * EXP(-m * PROTONMASS 
     &            * vs(idx-1)*vs(idx-1) / (2*BOLTZMANN*room_temp))
      end do
      cdf = cdf / cdf(1024)

      its = 0
      t = 0.
      last_coll = 0.
      collisions = 0
      call RANDOM_NUMBER(coll_random)
      next_sample = 0.
      window_end = 0
      sample_idx = 0
      call inf_trap_field(pos, t, omega, depth, rt, field)
      accel = field * cmr
      do while ((t < mt).and.(its.lt.1e+10))
            its = its + 1

            if (its.le.window_end) then
                  four_trajectory(sample_idx, 1:3) = 
     &                      four_trajectory(sample_idx, 1:3) + pos
                  four_trajectory(sample_idx, 4) =
     &                      four_trajectory(sample_idx, 4) + t
                  if (its.eq.window_end) then
                        four_trajectory(sample_idx, 1:4) = 
     &                      four_trajectory(sample_idx, 1:4)
     &                    / DBLE(averaging_window)
                  end if
            else
                  if (next_sample.le.t) then
                        next_sample = t + sample_interval
                        window_end = its + averaging_window
                        sample_idx = sample_idx + 1
                  end if
            end if

            a = NORM2(accel)
            v = NORM2(vel)
            a = a + 1.0e-15
            v = v + 1.0e-15
            tv = md / v
            ta = SQRT(2.0 * md / a)
            tstep = tv * ta / (tv + ta)
            tstep = MIN(tstep, mt * 1.0e-6)
            t = t + tstep

            pos = pos + (vel * tstep) + tstep*tstep * accel / 2.

            call inf_trap_field(pos, t, omega, depth, rt, field)
            accel_new = field * cmr

            ! Should a collision occur?
            if (coll_random.ge.EXP(-(t-last_coll)/expected_free)) then
                  call RANDOM_NUMBER(coll_random)
                  last_coll = t
                  collisions = collisions + 1
                  
                  vm = 1.0e+3 !* maxwell_distribution(cdf, vs)
                  call RANDOM_NUMBER(random)
                  theta = 2 * PI * random
                  call RANDOM_NUMBER(random)
                  phi = ACOS(2. * random - 1.)
                  vc(1) = vm * COS(theta) * COS(PHI)
                  vc(2) = vm * COS(theta) * SIN(PHI)
                  vc(3) = vm * SIN(theta)

                  call RANDOM_NUMBER(random)
                  theta = 2 * PI * random
                  call RANDOM_NUMBER(random)
                  phi = ACOS(2. * random - 1.)
                  normal(1) = 1. * COS(theta) * COS(PHI)
                  normal(2) = 1. * COS(theta) * SIN(PHI)
                  normal(3) = 1. * SIN(theta)

                  vel = vel
     &                  - DOT_PRODUCT(vel - vc, normal) * normal
            end if

            vel = vel + tstep * (accel + accel_new) / 2.
            accel = accel_new
      end do
      rpts = sample_idx - 1
      end subroutine

      subroutine inf_trap_cloud(
     &    particles, 
     &    positions,
     &    velocities, 
     &    ms, 
     &    qs,
     &    omega, 
     &    depth, 
     &    R, 
     &    max_t,
     &    max_dist,
     &    record_step,
     &    four_trajectory, 
     &    its, 
     &    rpts, 
     &    collisions)
     & bind(c, name="inf_trap_cloud")

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
      integer(c_int), intent(in) :: pts, averaging_window
      real(c_double), dimension(3), intent(in) :: position, velocity
      real(c_double), intent(in) :: m, q, coll_rate, room_temp
      real(c_double), intent(in) :: omega, depth, R
      real(c_double), intent(in) :: max_t, max_dist
      real(c_double), dimension(pts, 4), intent(out) :: four_trajectory
      integer(c_int), intent(out) :: its, rpts, collisions

C     Local variables
      real(c_double), dimension(3) :: pos, vel, accel, field, accel_new,
     &                                vc, normal
      real(c_double) :: t, ta, tv, a, v, sample_interval, next_sample,
     &   mt, md, rt, last_coll, expected_free, random, coll_random
      real(c_double) :: vm, theta, phi
      integer(c_int) :: window_end, sample_idx, idx
      real(c_double), dimension(1024) :: cdf, vs

      four_trajectory = 0.

      cmr = ((ECHARGE) * q) / ((PROTONMASS) * m);
      md  = max_dist * 1.0e-3 ! mm -> m
      mt  = max_t * 1.0e-6    ! us -> s
      rt  = R * 1.0e-3
      pos = position * 1.0e-3
      vel = velocity * 1.0e+3
      mt = max_t * 1.0e-6
      md = max_dist * 1.0e-3
      expected_free = 1. / (coll_rate * 1.0e+6);
      sample_interval = mt / DBLE(pts)

C     Initialize the CDF for the Maxwell-Boltzmann
      vs(1) = 0.
      cdf(1) = 0.
      dv = 100.*SQRT(2.*BOLTZMANN*room_temp/(m*PROTONMASS))/1024.0
      do idx = 2,1024
            vs(idx) = dv + vs(idx - 1)
            cdf(idx) = cdf(idx - 1) + dv * EXP(-m * PROTONMASS 
     &            * vs(idx-1)*vs(idx-1) / (2*BOLTZMANN*room_temp))
      end do
      cdf = cdf / cdf(1024)

      its = 0
      t = 0.
      last_coll = 0.
      collisions = 0
      next_sample = 0.
      window_end = 0
      sample_idx = 0
      call inf_trap_field(pos, t, omega, depth, rt, field)
      accel = field * cmr
      do while ((t < mt).and.(its.lt.1e+10))
            its = its + 1

            if (its.le.window_end) then
                  four_trajectory(sample_idx, 1:3) = 
     &                      four_trajectory(sample_idx, 1:3) + pos
                  four_trajectory(sample_idx, 4) =
     &                      four_trajectory(sample_idx, 4) + t
                  if (its.eq.window_end) then
                        four_trajectory(sample_idx, 1:4) = 
     &                      four_trajectory(sample_idx, 1:4)
     &                    / DBLE(averaging_window)
                  end if
            else
                  if (next_sample.le.t) then
                        next_sample = t + sample_interval
                        window_end = its + averaging_window
                        sample_idx = sample_idx + 1
                  end if
            end if

            a = NORM2(accel)
            v = NORM2(vel)
            a = a + 1.0e-15
            v = v + 1.0e-15
            tv = md / v
            ta = SQRT(2.0 * md / a)
            tstep = tv * ta / (tv + ta)
            tstep = MIN(tstep, mt * 1.0e-6)
            t = t + tstep

            pos = pos + (vel * tstep) + tstep*tstep * accel / 2.

            call inf_trap_field(pos, t, omega, depth, rt, field)
            accel_new = field * cmr

            vel = vel + tstep * (accel + accel_new) / 2.
            accel = accel_new
      end do
      rpts = sample_idx - 1
      end subroutine

      end module

