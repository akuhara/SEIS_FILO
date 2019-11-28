program main
  use mod_param
  use mod_vmodel
  use mod_rayleigh
  implicit none 
  integer :: n_arg
  character(len=200) :: param_file
  type(param) :: para
  type(vmodel) :: vm
  
  ! Get parameter file name from command line argument
  n_arg = command_argument_count()
  if (n_arg /= 1) then
     write(0, *)"USAGE: rayleigh_fwd [parameter file]"
     stop
  end if
  call get_command_argument(1, param_file)

  ! Read parameter file
  para = init_param(param_file)
  
  ! Set velocity model
  vm = init_vmodel()
  call vm%read_file(para%get_vmodel_file())
  

  
  
  

  stop
end program main
