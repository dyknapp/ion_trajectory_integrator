#define MAX_TRAJECTORY_POINTS 1048576
C     ^^ IF YOU CHANGE IT HERE, YOU NEED TO CHANGE IT IN THE OTHER FILES

      module trajectory_integration
      use iso_c_binding, only: c_int, c_double
      use trajectory_integration

      contains

      subroutine fly_ensemble(particles, xs,ys,zs,vxs,vys,vzs,
     &                        potential_maps, voltages, step_times_in,
     &                        time_steps, dimensions, is_electrode,
     &                        n_electrodes, m, q, din, maxdist, maxt,
     &                        x_traj,y_traj,z_traj,ts,exs,eys,ezs,its)
     & bind(c, name = "fly_ensemble")
C     Local variables:
      real(c_double), dimension(MAX_TRAJECTORY_POINTS)
     &      :: x_traj, y_traj, z_traj, ts, exs, eys, ezs

C     Variable declarations:
C     Dummy variables:
      integer(c_int), intent(in) :: particles!, interps

      real(c_double), dimension(particles), intent(in) 
     &      :: xx, yy, zz, vxx, vyy, vzz, m, q, din

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

      real(c_double), dimension(particles), intent(out) :: its

      real(c_double), dimension(:, :), allocatable, intent(out)
     &      :: x_trajs, y_trajs, z_trajs, tss, exss, eyss, ezss, itss

C     Allocate the output variables
      !interps = MAX_TRAJECTORY_POINTS ! For now, no interpolation
      allocate(x_trajs(MAX_TRAJECTORY_POINTS, particles))
      allocate(y_trajs(MAX_TRAJECTORY_POINTS, particles))
      allocate(z_trajs(MAX_TRAJECTORY_POINTS, particles))
      allocate(tss(MAX_TRAJECTORY_POINTS, particles))
      allocate(exss(MAX_TRAJECTORY_POINTS, particles))
      allocate(eyss(MAX_TRAJECTORY_POINTS, particles))
      allocate(ezss(MAX_TRAJECTORY_POINTS, particles))

C     Run the trajectory simulations
      do i = 1,particles
            call integrate_trajectory(xs(i), ys(i), zs(i), 
     &                          vxs(i), vys(i), vzs(i),
     &                          potential_maps, voltages, step_times_in,
     &                          time_steps, dimensions, is_electrode,
     &                          n_electrodes, m, q, din, maxdist, maxt,
     &                          x_traj,y_traj,z_traj,ts,exs,eys,ezs,its)
            do j = 1,MAX_TRAJECTORY_POINTS
                  x_trajs(i, j) = x_traj(j)
                  y_trajs(i, j) = y_traj(j)
                  z_trajs(i, j) = z_traj(j)
                  tss(i, j)     = ts(j)
                  exss(i, j)    = exs(j)
                  eyss(i, j)    = eys(j)
                  ezss(i, j)    = ezs(j)
                  itss(j)       = its
            end
      end

      end subroutine

      end module