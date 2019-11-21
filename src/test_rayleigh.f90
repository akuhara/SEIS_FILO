program main
  use mod_rayleigh
  use mod_vmodel
  implicit none 
  
  type(vmodel) :: vm
  type(rayleigh) :: ray
  
  vm = init_vmodel()
  call vm%set_example_land()
  call vm%display()

  ray = init_rayleigh(vm, fmin=1d0, fmax=100.d0, df=0.1d0, &
       & cmin=0.1d0, cmax=0.29d0, dc=0.001d0)
  call ray%dispersion()
  stop
end program main
