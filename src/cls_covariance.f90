module cls_covariance
  implicit none 
  
  type covariance
     private
     integer :: n
     double precision, allocatable :: r_inv(:, :)
   contains
     procedure :: get_inv => covariance_get_inv
  end type covariance
 
  interface covariance
     module procedure :: init_covariance
  end interface covariance

contains
  
  !---------------------------------------------------------------------

  type(covariance) function init_covariance(n, a_gauss, delta) &
       & result(self)
    integer, intent(in) :: n
    double precision, intent(in) :: a_gauss
    double precision, intent(in) :: delta
    double precision :: r_mat(n, n), s(n), u(n, n), vt(n, n)
    double precision, allocatable :: work(:)
    double precision :: diag(n, n), tmp(n, n)
    double precision :: dummy(1, 1)
    double precision :: r, r_tmp
    integer :: i, j, lwork, ierr
    
    self%n = n
    allocate(self%r_inv(n, n))

    
    r = exp(-a_gauss**2 * delta**2)
    r_mat(:,:) = 0.d0
    do i = 1, n
       do j = 1, n
          r_tmp = r ** ((i - j) ** 2)
          !if (r_tmp < 1.d-8) then
          !   r_tmp = 0.d0
          !end if
          
          r_mat(j, i) = r_tmp
       end do
    end do
    
    ! check 1
    !do i = 1, n
    !   do j = 1, n
    !      write(999,*)i, j, r_mat(j, i)
    !  end do
    !end do
    
    ! Use Lapack
    tmp(:,:) = r_mat(:,:)
    call dgesvd('A', 'A', n, n, tmp, n, s, u, n, &
         & vt, n, dummy, -1, ierr)
    lwork = nint(dummy(1,1))
    allocate(work(lwork))
    call dgesvd('A', 'A', n, n, tmp, n, s, u, n, &
         & vt, n, work, lwork, ierr)
    if (ierr /= 0) then
       write(0, *)"ERROR; in inverting covariance matrix"
       stop
    end if
    
    diag(:,:) = 0.d0
    do i = 1, n
       if (s(i) > 1.0d-3) then
          diag(i, i) = 1.d0 / s(i)
       else
          diag(i,i) = 0.d0
       end if
    end do
    self%r_inv(:, :) = &
         & matmul(matmul(transpose(vt),diag),transpose(u))
    
    ! Check2
!    tmp = matmul(self%r_inv, r_mat)
!    do i = 1, n
!       do j = 1, n
!          if (i /= j .and. abs(tmp(j,i)) >1.0e-1 ) then
!             write(*,*)i, j, tmp(j, i)
!          else if (i ==  j .and. abs(tmp(j,i) - 1.d0) > 1.0e-1) then
!             write(*,*)i, j, tmp(j, i)
!          end if
!          write(998,*)i, j, tmp(j,i), self%r_inv(j, i)
!      end do
!    end do


    return 
  end function init_covariance

  !---------------------------------------------------------------------

  function covariance_get_inv(self) result(inv)
    class(covariance), intent(in) :: self
    double precision :: inv(self%n, self%n)
    
    inv = self%r_inv(:,:)
    
    return 
  end function covariance_get_inv
  
  !---------------------------------------------------------------------
  
end module cls_covariance
