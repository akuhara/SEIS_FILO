program main
  use mod_rayleigh
  use mod_vmodel
  implicit none 
  
  type(vmodel) :: vm
  type(rayleigh) :: ray
  
  vm = init_vmodel()
  call vm%set_example_ocean()
  call vm%display()

  ray = init_rayleigh(vm, fmin=0.1d0, fmax=1.d0, df=0.01d0)
  call ray%dispersion()
  stop
end program main
