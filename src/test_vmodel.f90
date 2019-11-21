program main
  use mod_vmodel
  implicit none 
  type(vmodel) :: vm
  integer :: nlay
  double precision :: vs

  vm = init_vmodel()
  call vm%set_nlay(10)
  nlay = vm%get_nlay()
  write(*,*)"nlay = ", nlay
  call vm%display()

  call vm%set_example_ocean()
  call vm%display()
  
  vs = vm%get_vs(2)
  write(*,*)vs

  
  stop
end program main
