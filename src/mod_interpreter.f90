module mod_interpreter
  use mod_trans_d_model
  use mod_vmodel
  implicit none 
  
  type interpreter
     private
     logical :: ocean_flag = .false.
     double precision :: ocean_thick =0.d0
   contains
     procedure :: get_vmodel => interpreter_get_vmodel
  end type interpreter
  
  interface interpreter
     module procedure init_interpreter
  end interface interpreter
  
contains

  !---------------------------------------------------------------------
  type(interpreter) function init_interpreter(ocean_flag, ocean_thick)
    logical, intent(in), optional :: ocean_flag
    double precision, intent(in), optional :: ocean_thick

    if (present(ocean_flag)) then
       init_interpreter%ocean_flag = ocean_flag
    end if

    if (present(ocean_thick)) then
       init_interpreter%ocean_thick = ocean_thick
    end if

    return 
  end function init_interpreter

  !---------------------------------------------------------------------
  
  type(vmodel) function interpreter_get_vmodel(self, tm) result(vm)
    class(interpreter), intent(inout) :: self
    type(trans_d_model), intent(inout) :: tm
    integer :: nlay

    
    vm = init_vmodel()
    if (self%ocean_flag) then
       nlay = tm%get_k() + 1
    else
       nlay = tm%get_k()
    end if
    call vm%set_nlay(nlay)
    
    return 
  end function interpreter_get_vmodel
  
  
end module mod_interpreter

