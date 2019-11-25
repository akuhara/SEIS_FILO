program main
  use mod_trans_d_model
  implicit none 
  type(trans_d_model) :: tm
  integer :: k_max, k_min, n_ix, n_rx
  logical :: is_ok
  k_min = 1
  k_max = 10
  n_ix = 0
  n_rx = 3

  tm = init_trans_d_model(k_min=k_min, k_max=k_max, n_ix=n_ix, &
       n_rx=n_rx)
  call tm%set_k(6)
  call tm%display()
  is_ok = .true.
  do while (is_ok)
     call tm%birth(is_ok)
     if (is_ok) then
        call tm%display()
     end if
  end do

  is_ok = .true.
  do while (is_ok)
     call tm%death(is_ok)
     if (is_ok) then
        call tm%display()
     end if

  end do
  

  stop
end program main
  
