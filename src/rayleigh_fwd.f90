program main
  use mod_rayleigh
  use mod_vmodel
  use mod_param
  implicit none 
  integer :: n_arg
  character(len=200) :: param_file
  type(param) :: para
  
  ! Get parameter file name from command line argument
  n_arg = command_argument_count()
  if (n_arg /= 1) then
     write(0, *)"USAGE: rayleigh_fwd [parameter file]"
     stop
  end if
  call get_command_argument(1, param_file)

  ! Read parameter file
  para = init_param(param_file)
  
  

  
  
  

  stop
end program main
