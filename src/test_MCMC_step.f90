program main
  use mod_random
  use mod_trans_d_model
  use mod_MCMC_step
  implicit none 
  type(trans_d_model) :: tm, tm2
  integer :: i

  call init_random(1114444411, 21412414, 34343434, 5151515)

  tm = init_trans_d_model(k_min = 1, k_max = 99, n_rx = 1)

  call tm%set_birth(1, 1, -0.01d0, 0.01d0)
  call tm%set_prior(1, 1, -100.d0, 100.d0)
  call tm%set_perturb(1, 1.d0)
  call tm%generate_model()
  call tm%display()
  do i = 1, 2000000
     call propose_model(tm, tm2)
     tm = tm2
  end do
  call tm2%display()

  stop
end program main
