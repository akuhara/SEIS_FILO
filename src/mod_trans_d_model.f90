module mod_trans_d_model
  implicit none 
  
  type trans_d_model
     private
     integer :: k
     integer :: k_min
     integer :: k_max
     integer :: n_rx
     integer :: n_ix
     double precision, allocatable :: rx(:, :)
     integer, allocatable :: ix(:, :)
     
   contains
     procedure :: set_k => trans_d_model_set_k
     procedure :: birth => trans_d_model_birth
     procedure :: death => trans_d_model_death
     procedure :: display => trans_d_model_display

  end type trans_d_model
  
  interface trans_d_model
     module procedure init_trans_d_model
  end interface trans_d_model
  
contains

  !---------------------------------------------------------------------
  
  type(trans_d_model) function init_trans_d_model(k_min, k_max, &
       & n_rx, n_ix)
    integer, intent(in) :: k_min, k_max, n_rx, n_ix
    
    if (k_max < k_min) then
       write(0,*)"ERROR: k_max must be >= k_min"
       write(0,*)"     : k_min =", k_min, "k_max=", k_max
       stop
    end if
    init_trans_d_model%k_min = k_min
    init_trans_d_model%k_max = k_max
    
    init_trans_d_model%n_rx = n_rx
    if (n_rx >= 0) then
       allocate(init_trans_d_model%rx(n_rx ,k_max))
    end if

    init_trans_d_model%n_ix = n_ix
    if (n_ix >= 0) then
       allocate(init_trans_d_model%ix(n_ix ,k_max))
    end if
    
    return 
  end function init_trans_d_model
  
  !---------------------------------------------------------------------
  
  subroutine trans_d_model_set_k(self, k)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: k
    
    if (k < self%k_min .or. k >= self%k_max) then
       write(0,*)"ERROR: k must be k_min <= k < k_max ", &
            & "(trans_d_model_set_k)"
       write(0,*)"     : k=", k, " k_min=", self%k_min, &
            & " k_max=", self%k_max
    end if

    self%k = k
    
  end subroutine trans_d_model_set_k
  
  
  !---------------------------------------------------------------------

  subroutine trans_d_model_birth(self, is_ok)
    class(trans_d_model), intent(inout) :: self
    logical, intent(out) :: is_ok

    if (self%k < self%k_max - 1) then
       self%k = self%k + 1
       is_ok = .true.
    else
       is_ok = .false.
    end if

    return 
  end subroutine trans_d_model_birth

!---------------------------------------------------------------------

  subroutine trans_d_model_death(self, is_ok)
    class(trans_d_model), intent(inout) :: self
    logical, intent(out) :: is_ok

    if (self%k > self%k_min) then
       self%k = self%k - 1
       is_ok = .true.
    else
       is_ok = .false.
    end if

    return 
  end subroutine trans_d_model_death

  !---------------------------------------------------------------------

  subroutine trans_d_model_display(self)
    class(trans_d_model), intent(inout) :: self
    integer :: i, j
    write(*,*)
    write(*,*)"k = ", self%k
    write(*,*)
    do i = 1, self%n_ix
       write(*,*)"------------------------------"
       write(*,*)"Integer prameter", i
       do j = 1, self%k
          write(*,*)j, self%ix(i, j)
       end do
    end do
                 
    do i = 1, self%n_rx
       write(*,*)"------------------------------"
       write(*,*)"Double precision prameter", i

       do j = 1, self%k
          write(*,*)j, self%rx(i, j)
       end do
    end do
    write(*,*)
    
    return 
  end subroutine trans_d_model_display

end module mod_trans_d_model
  
