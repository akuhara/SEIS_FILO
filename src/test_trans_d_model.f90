program main
  use mod_trans_d_model
  use mod_random
  implicit none 
  type(trans_d_model) :: tm
  integer :: k_max, k_min, n_rx, n_ix
  
  k_min = 1
  k_max = 80
  n_rx = 1
  n_ix = 1
  
  call init_random(114444, 222, 33333, 444564)

  tm = init_trans_d_model(k_min=k_min, k_max=k_max, n_rx=n_rx, &
       & n_ix=n_ix)

  call tm%set_birth(1, 2, -100.d0, 100.d0)
  call tm%set_birth(2, 2, -100.d0, 100.d0)
  call tm%generate_model()
  call tm%display()

  stop
end program main
  
