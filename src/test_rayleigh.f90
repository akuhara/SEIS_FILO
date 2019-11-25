program main
  use mod_rayleigh
  use mod_vmodel
  implicit none 
  
  type(vmodel) :: vm
  type(rayleigh) :: ray
  
  vm = init_vmodel()
  call vm%set_example_ocean()
  call vm%display()

  ray = init_rayleigh(vm, fmin=0.3d0, fmax=10.0d0, df=0.01d0, &
       & cmin=1.0d0, cmax=1.4d0, dc=0.002d0)
  call ray%set_full_calculation(.false.)
  call ray%dispersion()
  stop
end program main
