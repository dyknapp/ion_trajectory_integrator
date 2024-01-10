      module cloud_integrator
      use iso_c_binding, only: c_int, c_double

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

      end module