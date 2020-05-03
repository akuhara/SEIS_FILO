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
  use cls_hyper_model
  
  implicit none 
  
  
  double precision, parameter, private :: log_2pi_half &
       & = 0.5d0 * log(2.d0 * acos(-1.d0)) 
  
contains
  
  !---------------------------------------------------------------------
  
  subroutine forward_recv_func(tm, hyp, intpr, obs, rf, cov, &
       & log_likelihood, is_ok)
    type(trans_d_model), intent(in) :: tm
    type(hyper_model), intent(in) :: hyp
    type(interpreter), intent(inout) :: intpr
    type(observation_recv_func), intent(in) :: obs
    type(recv_func), intent(inout) :: rf(:)
    type(covariance), intent(in) :: cov(:)
    double precision, intent(out) :: log_likelihood
    logical, intent(out) :: is_ok
    double precision, allocatable :: misfit(:), phi1(:)
    double precision :: phi, s, rms
    type(vmodel) :: vm
    integer :: i, j, n


    ! Forward computation
    call intpr%construct_vmodel(tm, vm, is_ok)
    if (.not. is_ok) then
       log_likelihood = minus_infty
       return
    end if
    !log_likelihood = 0.d0  !!!!!!!!!!!!!!!!!!!!!!
    !return !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    do i = 1, obs%get_n_rf()
       call rf(i)%set_vmodel(vm)
       call rf(i)%compute()
    end do
    !log_likelihood = 0.d0 !!!!!!!!!!!!!!!!
    !return !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! calc misfit
    log_likelihood = 0.d0
    do i = 1, obs%get_n_rf()
       n = obs%get_n_smp(i)
       s = hyp%get_x(i)
       allocate(misfit(n), phi1(n))
       rms = 0.d0
       do  j = 1, n
          misfit(j) = obs%get_rf_data(j, i) - rf(i)%get_rf_data(j)
          rms = rms + misfit(j) ** 2
       end do
       phi1 = matmul(misfit, cov(i)%get_inv())
       phi = dot_product(phi1, misfit)
       log_likelihood = log_likelihood - 0.5d0 * phi / (s * s) &
            & - n * log_2pi_half - n * log(s) &
            & - 0.5d0 * cov(i)%get_log_det_r()
       deallocate(misfit, phi1)
    end do
    
    
    return 
  end subroutine forward_recv_func

  !-----------------------------------------------------------------------
  
  subroutine forward_disper(tm, hyp, intpr, obs, disp, log_likelihood, &
       & is_ok)
    implicit none 
    type(trans_d_model), intent(in) :: tm
    type(hyper_model), intent(in) :: hyp
    type(interpreter), intent(inout) :: intpr
    type(observation_disper), intent(in) :: obs
    type(disper), intent(inout) :: disp(:)
    double precision, intent(out) :: log_likelihood
    logical, intent(out) :: is_ok
    double precision :: misfit_c, misfit_u
    type(vmodel) :: vm
    double precision :: sc, su
    integer :: i, j, nc, nu
    
    ! calculate synthetic dispersion curves
    call intpr%construct_vmodel(tm, vm, is_ok)
    if (.not. is_ok) then
       log_likelihood = minus_infty
       return
    end if
    
    do j = 1, obs%get_n_disp()
       call disp(j)%set_vmodel(vm)
       call disp(j)%dispersion()
    end do
    
    ! calc misfit
    log_likelihood = 0.d0
    all_disp: do j = 1, obs%get_n_disp()
       misfit_c = 0.d0
       misfit_u = 0.d0
       nc = 0
       nu = 0
       sc = hyp%get_x(2*j-1)
       su = hyp%get_x(2*j)
       do i = 1, obs%get_nf(j)
          if (disp(j)%get_c(i) == 0.d0 .or. &
               & disp(j)%get_u(i) == 0.d0) then
             ! Failer in root search
             log_likelihood = minus_infty
             is_ok = .false.
             return 
          end if
          if (obs%get_c_use(i,j)) then
             misfit_c = misfit_c &
                  & + (disp(j)%get_c(i) - obs%get_c(i,j)) ** 2
             nc = nc + 1
          end if
          if (obs%get_u_use(i,j)) then
             misfit_u = misfit_u &
                  & + (disp(j)%get_u(i) - obs%get_u(i,j)) ** 2 
             nu = nu + 1
          end if
       end do
       log_likelihood = log_likelihood  &
            & - 0.5d0 * misfit_c / (sc * sc) &
            & - nc * log_2pi_half - nc * log(sc)
       log_likelihood = log_likelihood  &
            & - 0.5d0 * misfit_u / (su * su) &
            & - nu * log_2pi_half - nu * log(su)
    end do all_disp

    
    return 
  end subroutine forward_disper
  
  !---------------------------------------------------------------------
  

end module mod_forward
