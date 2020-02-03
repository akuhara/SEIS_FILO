module mod_observation_recv_func
  implicit none
  
  type observation_recv_func
     private
     
     ! receiver function
     integer :: n_rf
     double precision, allocatable :: a_gauss(:)
     double precision, allocatable :: rayp(:)
     double precision, allocatable :: delta(:)
     double precision, allocatable :: sigma_min(:)
     double precision, allocatable :: sigma_max(:)
     double precision, allocatable :: rf_data(:,:)
     double precision, allocatable :: t_start(:)
     double precision, allocatable :: t_end(:)
     character(len=1), allocatable :: phase(:)
     
   contains
  end type observation_recv_func
  
  interface observation_recv_func
     module procedure :: init_observation_recv_func
  end interface observation_recv_func

contains

  !--------------------------------------------------------------------- 
  
  type(observation_recv_func) function &
     & init_observation_recv_func(recv_func_in) &
     & result(self)
    character(len=*), intent(in) :: recv_func_in
    integer :: io, ierr

    write(*,*)"Reading observed receiver functions from ", &
         & trim(recv_func_in)

    open(newunit = io, file = recv_func_in, status = "old", &
         & iostat = ierr)
    if (ierr /= 0) then
       write(0,*)"ERROR: cannot open ", trim(recv_func_in)
       write(0,*)"     : (init_observation_recv_func)"
       stop
    end if
    
    ! Read file
    
    
    
    return 
  end function init_observation_recv_func
  !---------------------------------------------------------------------
end module mod_observation_recv_func
