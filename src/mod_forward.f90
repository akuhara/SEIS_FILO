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
       & is_sphere, r_earth, simulate_prior, log_likelihood, is_ok)
    type(trans_d_model), intent(in) :: tm
    type(hyper_model), intent(in) :: hyp
    type(interpreter), intent(inout) :: intpr
    type(observation_recv_func), intent(in) :: obs
    type(recv_func), intent(inout) :: rf(:)
    type(covariance), intent(in) :: cov(:)
    logical, intent(in) :: is_sphere
    double precision, intent(in) :: r_earth
    logical, intent(in) :: simulate_prior
    double precision, intent(out) :: log_likelihood
    logical, intent(out) :: is_ok
    double precision, allocatable :: misfit(:), phi1(:)
    double precision :: phi, s
    type(vmodel) :: vm, vm_flattened
    integer :: i, j, n


    call intpr%construct_vmodel(tm, vm, is_ok)
    !call vm%display()
    
    ! Earth flattening
    if (is_sphere) then
       call vm%sphere2flat(r_earth, vm_flattened)
    end if

    
    if (.not. is_ok) then
       log_likelihood = minus_infty
       return
    end if
    
    if (simulate_prior) then
       log_likelihood = 0.d0
       return
    end if
    
    do i = 1, obs%get_n_rf()
       if (is_sphere) then
          call rf(i)%set_vmodel(vm_flattened)
       else
          call rf(i)%set_vmodel(vm)
       end if
       call rf(i)%compute()
    end do
    
    ! calc misfit
    log_likelihood = 0.d0
    do i = 1, obs%get_n_rf()
       n = obs%get_n_smp(i)
       s = hyp%get_x(i)
       allocate(misfit(n), phi1(n))
       do  j = 1, n
          misfit(j) = obs%get_rf_data(j, i) - rf(i)%get_rf_data(j)
       end do
       phi1 = matmul(misfit, cov(i)%get_inv())
       phi = dot_product(phi1, misfit)
       log_likelihood = log_likelihood - 0.5d0 * phi / (s * s) &
            & - n * log_2pi_half - n * log(s) &
            & - 0.5d0 * cov(i)%get_log_det_r()
       !log_likelihood = log_likelihood - sqrt(2.d0 * phi) &
       !     & - 0.25d0 * (n - 1) * log(phi) - n * log(s)

       deallocate(misfit, phi1)
    end do
    
    
    return 
  end subroutine forward_recv_func

  !-----------------------------------------------------------------------
  
  subroutine forward_disper(tm, hyp, intpr, obs, disp, &
       & is_sphere, r_earth, simulate_prior, log_likelihood, is_ok)
    implicit none 
    type(trans_d_model), intent(in) :: tm
    type(hyper_model), intent(in) :: hyp
    type(interpreter), intent(inout) :: intpr
    type(observation_disper), intent(in) :: obs
    type(disper), intent(inout) :: disp(:)
    logical, intent(in) :: is_sphere
    double precision, intent(in) ::r_earth
    logical, intent(in) :: simulate_prior
    double precision, intent(out) :: log_likelihood
    logical, intent(out) :: is_ok
    double precision :: misfit_c, misfit_u, misfit_hv, misfit_ra
    type(vmodel) :: vm, vm_flattened
    double precision :: sc, su, shv, sra
    integer :: i, j, nc, nu, nhv, nra
    
    call intpr%construct_vmodel(tm, vm, is_ok)
    
    ! Earth flattening
    if (is_sphere) then
       call vm%sphere2flat(r_earth, vm_flattened)
    end if
    
    if (.not. is_ok) then
       log_likelihood = minus_infty
       return
    end if
    
    if (simulate_prior) then
       log_likelihood = 0.d0
       return
    end if

    do j = 1, obs%get_n_disp()
       if (is_sphere) then
          call disp(j)%set_vmodel(vm_flattened)
       else
          call disp(j)%set_vmodel(vm)
       end if
       call disp(j)%dispersion()
    end do
    
    ! calc misfit
    log_likelihood = 0.d0
    all_disp: do j = 1, obs%get_n_disp()
       misfit_c = 0.d0
       misfit_u = 0.d0
       misfit_hv = 0.d0
       nc = 0
       nu = 0
       nhv = 0
       sc = hyp%get_x(4*j-3)
       su = hyp%get_x(4*j-2)
       shv = hyp%get_x(4*j-1)
       sra = hyp%get_x(4*j)
       do i = 1, obs%get_nx(j)
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
          if (obs%get_hv_use(i,j)) then
             misfit_hv = misfit_hv &
                  & + (disp(j)%get_hv(i) - obs%get_hv(i,j)) ** 2
             nhv = nhv + 1
          end if
          if (obs%get_ra_use(i,j)) then
             misfit_ra = misfit_ra &
                  & + (disp(j)%get_ra(i) - obs%get_ra(i,j)) ** 2
             nra = nra + 1
          end if

       end do
       log_likelihood = log_likelihood  &
            & - 0.5d0 * misfit_c / (sc * sc) &
            & - nc * log_2pi_half - nc * log(sc)
       log_likelihood = log_likelihood  &
            & - 0.5d0 * misfit_u / (su * su) &
            & - nu * log_2pi_half - nu * log(su)
       log_likelihood = log_likelihood  &
            & - 0.5d0 * misfit_hv / (shv * shv) &
            & - nhv * log_2pi_half - nhv * log(shv)
       log_likelihood = log_likelihood  &
            & - 0.5d0 * misfit_ra / (sra * sra) &
            & - nra * log_2pi_half - nra * log(sra)
       
    end do all_disp

    
    return 
  end subroutine forward_disper
  
  !---------------------------------------------------------------------
  

end module mod_forward
