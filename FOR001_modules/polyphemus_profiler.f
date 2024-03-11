#define STINTLENGTH 65536

      program polyphemus_profiler
      use iso_c_binding, only: c_int, c_double, c_double_complex
      use polyphemus_integrator

      !INPUTS
      integer(c_int) :: particles
      real(c_double), dimension(2, 3) 
     &      :: positions, velocities
      real(c_double), dimension(2) 
     &      :: ms, qs
      real(c_double) :: omega, depth, R
      real(c_double) :: end_time, max_dist
      real(c_double) :: record_step, burst_time

      ! OUTPUTS
      integer, parameter :: dp = kind(0.0d0)
      real(c_double), dimension(STINTLENGTH, 2, 3) :: trajectory
      real(c_double), dimension(STINTLENGTH) :: times
      integer(c_int) :: its, recorded

      ! LOCALS
      integer :: c1, c2, cr, cm
      real :: rate


      particles = 2
      positions =   RESHAPE((/ 0.0,  0.0096, 0.0,
     &                         0.0, -0.0096, 0.0/),
     &                      (/2, 3/), ORDER=(/2, 1/))
      velocities =  RESHAPE((/  0.0601, 0.0, 0.0,
     &                         -0.0601, 0.0, 0.0/),
     &                      (/2, 3/), ORDER=(/2, 1/))
      ms = (/1.0, -1.0/)
      qs = (/1.0, -1.0/)
      omega = 1.0
      depth = 0.0
      R     = 1.0
      end_time = 1.0
      max_dist = 1.0e-8
      record_step = 1.0e-3
      burst_time  = 0.0
      
      ! First initialize the system_clock
      CALL system_clock(count_rate=cr)
      CALL system_clock(count_max=cm)
      rate = REAL(cr)
      WRITE(*,*) "System_clock rate : ", rate
      CALL SYSTEM_CLOCK(c1)
      call polyphemus(
     &                particles, 
     &                positions,
     &                velocities, 
     &                ms,
     &                qs,
     &                omega, 
     &                depth, 
     &                R, 
     &                end_time,
     &                max_dist,
     &                record_step,
     &                burst_time,
     &                trajectory, 
     &                times, 
     &                its,
     &                recorded
     &               )
      CALL SYSTEM_CLOCK(c2)
      WRITE(*,*) "System_clock time : ", (c2 - c1)/rate
      WRITE(*,*) "its               : ", its
      WRITE(*,*) "recorded          : ", recorded
      WRITE(*,*) "its/s             : ", DBLE(its) / DBLE((c2-c1)/rate)

      end program