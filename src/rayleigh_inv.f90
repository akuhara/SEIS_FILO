!-----------------------------------------------------------------------
program main
  use mod_random
  use mod_trans_d_model
  use mod_MCMC_step
  use mod_rayleigh
  use mod_interpreter
  implicit none 
  include 'mpif.h'

  integer, parameter :: n_iter = 100000
  integer, parameter :: k_min = 1, k_max = 21
  integer, parameter :: n_rx = 3
  integer, parameter :: id_vs = 1, id_vp = 2, id_z = 3
  integer, parameter :: id_uni = 1, id_gauss = 2
  double precision, parameter :: vs_min = 2.5d0, vs_max = 5.0d0
  double precision, parameter :: vp_min = 4.0d0, vp_max = 5.0d0
  double precision, parameter :: z_min = 0.d0, z_max = 50.d0
  double precision, parameter :: dev_vs = 0.05d0
  double precision, parameter :: dev_vp = 0.05d0
  double precision, parameter :: dev_z  = 0.2d0
  
  integer :: i
  type(trans_d_model) :: tm, tm2

  
  call init_random(12345678, 23456789, 34567890, 45678901)
  
  ! Set model parameter
  tm = init_trans_d_model(k_min=k_min, k_max=k_max, n_rx=n_rx)
  call tm%set_prior(id_vs, id_uni, vs_min, vs_max)
  call tm%set_prior(id_vp, id_uni, vp_min, vp_max)
  call tm%set_prior(id_z,  id_uni, z_min,  z_max )
  call tm%set_birth(id_vs, id_uni, vs_min, vs_max)
  call tm%set_birth(id_vp, id_uni, vp_min, vp_max)
  call tm%set_birth(id_z,  id_uni, z_min,  z_max )
  call tm%set_perturb(id_vs, dev_vs)
  call tm%set_perturb(id_vp, dev_vp)
  call tm%set_perturb(id_z,  dev_z)
  call tm%generate_model()

  ! Main
  do i = 1, n_iter
     call propose_model(tm, tm2)
  end do
  
  

  stop
end program main
