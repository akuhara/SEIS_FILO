program main
  use mod_param
  use mod_vmodel
  use mod_rayleigh
  implicit none 
  integer :: n_arg
  character(len=200) :: param_file
  type(param) :: para
  type(vmodel) :: vm
  type(rayleigh) :: ray
  
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
  call vm%read_file(para%get_vmod_in())
  
  ! Calculate dispersion curve
  ray = init_rayleigh(&
       & vm   = vm, &
       & fmin = para%get_fmin(), &
       & fmax = para%get_fmax(), &
       & df   = para%get_df(), &
       & cmin = para%get_cmin(), &
       & cmax = para%get_cmax(), &
       & dc   = para%get_dc(), &
       & ray_out = para%get_ray_out() &
       & )
  
  call ray%dispersion()
 
  
  
  

  stop
end program main
