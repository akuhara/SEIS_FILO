module mod_forward
  use cls_trans_d_model
  use cls_interpreter
  use cls_observation_recv_func
  use cls_recv_func
  use cls_observation_disper
  use cls_disper
  use cls_vmodel
  use mod_const, only: minus_infty
  use cls_covariance
  
  implicit none 
  
  
contains
  
  !---------------------------------------------------------------------
  
  subroutine forward_recv_func(tm, intpr, obs, rf, cov, log_likelihood)
    type(trans_d_model), intent(in) :: tm
    type(interpreter), intent(inout) :: intpr
    type(observation_recv_func), intent(in) :: obs
    type(recv_func), intent(inout) :: rf(:)
    type(covariance), intent(in) :: cov(:)
    double precision, intent(out) :: log_likelihood
    double precision, allocatable :: misfit(:), phi1(:)
    double precision :: phi, s
    type(vmodel) :: vm
    integer :: i, j, n

    ! Forward computation
    vm = intpr%get_vmodel(tm)
    do i = 1, obs%get_n_rf()
       call rf(i)%set_vmodel(vm)
       call rf(i)%compute()
    end do
    
  ! calc misfit
    log_likelihood = 0.d0
    do i = 1, obs%get_n_rf()
       n = obs%get_n_smp(i)
       s = obs%get_sigma(i)
       allocate(misfit(n), phi1(n))
       do  j = 1, n
          misfit(j) = obs%get_rf_data(j, i) - rf(i)%get_rf_data(j)
       end do
       phi1 = matmul(misfit, cov(i)%get_inv())
       phi = dot_product(phi1, misfit)
       log_likelihood = log_likelihood - 0.5d0 * phi / (s * s) &
            & - dble(n) * log(s)
       deallocate(misfit, phi1)
    end do
    
    
    return 
  end subroutine forward_recv_func

  !-----------------------------------------------------------------------
  
  subroutine forward_disper(tm, intpr, obs, disp, log_likelihood)
    implicit none 
    type(trans_d_model), intent(in) :: tm
    type(interpreter), intent(inout) :: intpr
    type(observation_disper), intent(in) :: obs
    type(disper), intent(inout) :: disp(:)
    double precision, intent(out) :: log_likelihood
    type(vmodel) :: vm
    integer :: i, j
    
    ! calculate synthetic dispersion curves
    vm = intpr%get_vmodel(tm)
    do j = 1, obs%get_n_disp()
       call disp(j)%set_vmodel(vm)
       call disp(j)%dispersion()
    end do
    
    ! calc misfit
    log_likelihood = 0.d0
    all_disp: do j = 1, obs%get_n_disp()
       do i = 1, obs%get_nf(j)
          if (disp(j)%get_c(i) == 0.d0 .or. &
               & disp(j)%get_u(i) == 0.d0) then
             log_likelihood = minus_infty
             exit all_disp
          end if
          if (obs%get_sig_c(i,j) > 0.d0) then
             log_likelihood = log_likelihood &
                  & - (disp(j)%get_c(i) - obs%get_c(i,j)) ** 2 &
                  & / (obs%get_sig_c(i,j) ** 2)
          end if
          if (obs%get_sig_u(i,j) > 0.d0) then
             log_likelihood = log_likelihood &
                  & - (disp(j)%get_u(i) - obs%get_u(i,j)) ** 2 &
                  & / (obs%get_sig_u(i,j) ** 2)
          end if
       end do
    end do all_disp
    log_likelihood = 0.5d0 * log_likelihood
    
    return 
  end subroutine forward_disper
  
  !---------------------------------------------------------------------
  

end module mod_forward
