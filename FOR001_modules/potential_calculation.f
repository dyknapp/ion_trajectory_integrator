      module potential_calculation
      use iso_c_binding, only: c_int, c_double

      contains

C     Assumes that Dirichlet boundary conditions are applied in guess
C           -> For the cells corresponding to the boundary condition,
C               the bc_mask variable prevents them from being updated.
C               in theory, this could also be used for some kind of 
C               weighted evaluation... not sure if that is useful.

C     Neumann boundary conditions are applied along the border.
      subroutine iterate_laplace(guess, bc_mask, threshold, maxits, 
     &                           d1,d2,
     &                           its, refined)
     & bind(c, name="iterate_laplace")
C     Dummies     
      integer(c_int), intent(in) :: d1, d2
      real(c_double), dimension(d1, d2), intent(in) :: guess, bc_mask
      real(c_double), intent(in) :: threshold
      real(c_double), intent(in) :: maxits
      real(c_double), intent(out) :: its
      real(c_double), dimension(d1, d2), intent(out) :: refined
C     Locals
      integer, dimension(2) :: input_shape

      input_shape = SHAPE(guess)
      its = 0.0
      change = dble(input_shape(1) * input_shape(2))
      refined = guess
      do while ((its < maxits) .and. (change > threshold))
            its = its + 1.0
            change = 0.0
C           Neumann boundary conditions
C           Assume that edge of domain is composed of "ghost cells"
C           If you think about it for 0.5 seconds, you see that corners
C                 are irrelevant.
            do i = 2,(input_shape(1) - 1)
                  refined(i, 1) = refined(i, 2)
                  refined(i, input_shape(2)) 
     &                  = refined(i, input_shape(2)-1)
            end do
            do j = 2,(input_shape(2) - 1)
                  refined(1, j) = refined(2, j)
                  refined(input_shape(2), j) 
     &                  = refined(input_shape(2)-1, j)
            end do
C           Go one way:
            do i = 2,(input_shape(1) - 1)
                  do j = 2,(input_shape(2) - 1)
                        previous = refined(i, j)
                        refined(i, j) = (0.25 * (
     &                        + (((2*i)-1.5)/(2*i-2))*refined(i + 1, j)
     &                        + (((2*i)-2.5)/(2*i-2))*refined(i - 1, j)
     &                        + refined(i, j + 1)
     &                        + refined(i, j - 1)) * bc_mask(i, j))
     &                        + refined(i, j) * (1.0 - bc_mask(i, j))
                        change = change + ABS(refined(i, j) - previous)
                  end do
            end do
C           Go the other way:
            do ir = 2,(input_shape(1) - 1)
                  do jr = 2,(input_shape(2) - 1)
                        i = 1 + input_shape(1) - ir
                        j = 1 + input_shape(2) - jr
                        previous = refined(i, j)
                        refined(i, j) = (0.25 * (
     &                        + (((2*i)-1.5)/(2*i-2))*refined(i + 1, j)
     &                        + (((2*i)-2.5)/(2*i-2))*refined(i - 1, j)
     &                        + refined(i, j + 1)
     &                        + refined(i, j - 1)) * bc_mask(i, j))
     &                        + refined(i, j) * (1.0 - bc_mask(i, j))
                        change = change + ABS(refined(i, j) - previous)
                  end do
            end do

C           Normalize the change
            change = 0.5*change / dble(input_shape(1) * input_shape(2))
      end do
      end subroutine

C     A quick and dirty interpolation function for one equally-spaced
C           array to the next
      subroutine nn_interpolate(array,old_samples,new_samples,new_array)
     &      bind(c, name="nn_interpolate")
      integer(c_int), intent(in) :: old_samples, new_samples
      real(c_double), dimension(old_samples), intent(in) :: array
      real(c_double), dimension(new_samples), intent(out)  :: new_array
      real(c_double) :: increment
      increment = dble(old_samples - 1) / dble(new_samples - 1)
      do i = 0,new_samples - 1
            new_array(i+1) = array(1 + NINT(dble(i) * increment))
      end do
      end subroutine

C     nn_interpolate, but in 2 dimensions
      subroutine nn_interpolate2D(array, odim1, odim2, ndim1, ndim2,
     &                            new_array)
     &                            bind(c, name="nn_interpolate2D")
      integer(c_int), intent(in) :: odim1, odim2, ndim1, ndim2
      real(c_double), dimension(odim1, odim2), intent(in) :: array
      real(c_double), dimension(ndim1, ndim2), intent(out) :: new_array
      real(c_double) :: increment1, increment2
      increment1 = dble(odim1 - 1) / dble(ndim1 - 1)
      increment2 = dble(odim2 - 1) / dble(ndim2 - 1)
      do i = 0,ndim1 - 1
            do j = 0,ndim2 - 1
                  new_array(i+1, j+1) =array(1+NINT(dble(i)*increment1),
     &                                       1+NINT(dble(j)*increment2))
            end do
      end do
      end subroutine

      subroutine refined_laplace(guess, bc_mask, threshold, maxits, 
     &                           d1, d2, refinements,
     &                           refined)
     & bind(c, name="refined_laplace")
C     Dummies     
      integer(c_int), intent(in) :: d1, d2, refinements
      real(c_double), dimension(d1, d2), intent(in) :: guess, bc_mask
      real(c_double), intent(in) :: threshold
      real(c_double), intent(in) :: maxits
      real(c_double), dimension(d1, d2), intent(out) :: refined
C     Locals
      real(c_double), dimension(d1, d2) :: bcs
      real(c_double), allocatable :: temp_bcs(:, :), temp_guess(:, :),
     &                               temp2(:, :)
      integer(c_int), dimension(refinements, 2) :: step_dims
      real(c_double) :: its

C     Pull out the boundary conditions
      do i = 1,d1
            do j = 1,d2
                  bcs(i, j) = (1.0 - bc_mask(i, j)) * guess(i, j);
            end do
      end do

C     Refined Laplace solver
C     Pre-calculate the refined mesh sizes
      step_dims(refinements, 1) = d1
      step_dims(refinements, 2) = d2
      do k = 1,refinements - 1
            step_dims(refinements - k, 1) =
     &            step_dims(refinements - k + 1, 1) / 2
            step_dims(refinements - k, 2) =
     &            step_dims(refinements - k + 1, 2) / 2
      end do

C     Solve the Laplace equation
      do k = 1,refinements
C           Handling the interpolation / allocation of the guess
C           For first iteration, downsample the input guess
            if (k == 1) then
                  allocate(temp_guess(step_dims(k, 1), step_dims(k, 2)))
                  allocate(temp2(step_dims(k, 1), step_dims(k, 2)))
                  call nn_interpolate2D(guess, d1, d2, 
     &                            step_dims(k, 1), step_dims(k, 2),
     &                            temp_guess)
C           For 2nd iter. onward, upsample previous result
C           Requires extra temporary variable while reallocating
            else
                  ! Read further in code... this is already true:
                  !temp2 = temp_guess
                  deallocate(temp_guess)
                  allocate(temp_guess(step_dims(k, 1), step_dims(k, 2)))
                  call nn_interpolate2D(temp2,
     &                            step_dims(k-1, 1), step_dims(k-1, 2),
     &                            step_dims(k, 1), step_dims(k, 2),
     &                            temp_guess)
                  deallocate(temp2)
                  allocate(temp2(step_dims(k, 1), step_dims(k, 2)))
            end if
C           BC masks are easier:Downsample from input mask for each step
C           We just need to remember to deallocate after the first time.
            if (k > 1) then
                  deallocate(temp_bcs)
            end if
            allocate(temp_bcs(step_dims(k, 1), step_dims(k, 2)))
            call nn_interpolate2D(bc_mask, d1, d2, 
     &                            step_dims(k, 1), step_dims(k, 2),
     &                            temp_bcs)
C           Reapply boundary conditions to rough mesh.
            call nn_interpolate2D(bcs, d1, d2,
     &                            step_dims(k, 1), step_dims(k, 2),
     &                            temp2)
            do m = 1,step_dims(k, 1)
                  do n = 1,step_dims(k, 2)
                        temp_guess(m, n) 
     &                              = (1.0 - temp_bcs(m, n))
     &                              * temp2(m, n)
                  end do
            end do
C           Now, run the Laplace solver on the rough mesh.
            call iterate_laplace(temp_guess, temp_bcs, threshold,maxits,
     &                           step_dims(k, 1), step_dims(k, 2),
     &                           its, temp2)
            temp_guess = temp2
      end do
      refined = temp_guess
      deallocate(temp_guess)
      deallocate(temp_bcs)
      deallocate(temp2)
      end subroutine
      end module