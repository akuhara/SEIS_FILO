module mod_mcmc
  use mod_random
  use mod_trans_d_model
  implicit none 
  
  type mcmc 
     private
     
     type(trans_d_model) :: tm
     
   contains
     procedure :: propose_model => mcmc_propose_model
     procedure :: accept_model  => mcmc_accept_model
  end type mcmc
  
  interface mcmc
     module procedure init_mcmc
  end interface mcmc

contains

  !---------------------------------------------------------------------
  
  type(mcmc) function init_mcmc(tm)
    type(trans_d_model), intent(in) :: tm
    
    init_mcmc%tm = tm
    
    return 
  end function init_mcmc

  !---------------------------------------------------------------------
  
  subroutine mcmc_propose_model(self, tm_proposed, is_ok)
    class(mcmc), intent(inout) :: self
    type(trans_d_model), intent(out) :: tm_proposed
    logical, intent(out) :: is_ok
    integer :: nparam, itype, iparam

    
    nparam = self%tm%get_n_x()
    itype = int(rand_u() * (nparam + 2))
    tm_proposed = self%tm
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
    
  end subroutine mcmc_propose_model

  !---------------------------------------------------------------------

  subroutine mcmc_accept_model(self, tm)
    class(mcmc), intent(inout) :: self
    type(trans_d_model), intent(in) :: tm

    self%tm = tm

    return 
  end subroutine mcmc_accept_model
  
end module mod_mcmc
