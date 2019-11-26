module mod_MCMC_step
  use mod_random
  use mod_trans_d_model
  implicit none 
  
contains

  !---------------------------------------------------------------------
  
  subroutine propose_model(tm_current, tm_proposed)
    type(trans_d_model), intent(inout) :: tm_current
    type(trans_d_model), intent(out) :: tm_proposed
    integer :: nparam, itype, iparam
    logical :: is_ok
    
    nparam = tm_current%get_n_x()
    itype = int(rand_u() * (nparam + 2))
    tm_proposed = tm_current
    if (itype == 0) then
       ! Birth
       call tm_proposed%birth(is_ok)
    else if (itype == 1) then
       ! Death
       call tm_proposed%death(is_ok)
    else
       iparam = itype - 1
       call tm_proposed%perturb(iparam, is_ok)
    end if
    return
    
  end subroutine propose_model

  !---------------------------------------------------------------------
  
end module mod_MCMC_step
