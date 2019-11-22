program main
  use mod_rayleigh
  use mod_vmodel
  implicit none 
  
  type(vmodel) :: vm
  type(rayleigh) :: ray
  
  vm = init_vmodel()
  call vm%set_example_land()
  call vm%display()

  ray = init_rayleigh(vm, fmin=80.d0, fmax=500.d0, df=0.1d0, &
       & cmin=0.15d0, cmax=0.35d0, dc=0.01d0)
  call ray%dispersion()
  stop
end program main
