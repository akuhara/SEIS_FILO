program main
  use mod_trans_d_model
  use mod_random
  implicit none 
  type(trans_d_model) :: tm
  integer :: k_max, k_min, n_rgx, i, n_rux, n_igx, n_iux
  double precision :: sig_x, mean_x, x_min, x_max
  integer :: i_min, i_max
  logical :: is_ok
  k_min = 1
  k_max = 20
  n_rgx = 1
  n_rux = 0
  n_igx = 3
  n_iux = 6
  
  sig_x = 20.d0
  mean_x = 0.0d0
  x_min = -100.d0
  x_max = 100.d0
  i_min = 10
  i_max = 20
  
  call init_random(114444, 222, 33333, 444564)

  tm = init_trans_d_model(k_min=k_min, k_max=k_max, n_rgx=n_rgx, &
       & n_rux=n_rux, n_igx=n_igx, n_iux=n_iux)
  do i = 1, n_rgx
     call tm%set_x_gauss_birth(i, mean_x, sig_x)
  end do
  do i = 1, n_rux
     call tm%set_x_uni_birth(i, x_min, x_max)
  end do
  do i = 1, n_igx
     call tm%set_i_gauss_birth(i, mean_x, sig_x)
  end do
  do i = 1, n_iux
     call tm%set_i_uni_birth(i, i_min, i_max)
  end do
  call tm%generate_model()
  call tm%display()
  
  do i = 1, 100
     call tm%death(is_ok)
     if (.not. is_ok) then
        write(0, *)"ERROR"
        stop
     end if
     call tm%display()
     call tm%birth(is_ok)
     call tm%display()
     if (.not. is_ok) then
        write(0, *)"ERROR"
        stop
     end if
     write(*,*)i
  end do

  stop
end program main
  
