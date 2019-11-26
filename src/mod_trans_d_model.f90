module mod_trans_d_model
  use mod_random
  implicit none 
  
  type trans_d_model
     private
     integer :: k
     integer :: k_min
     integer :: k_max
     
     integer :: n_rgx ! Number of real parameter 
                      ! with the birth from Gaussian
     integer :: n_igx ! Number of integer parameter
                      ! with the birth from Gaussian
     integer :: n_rux ! Number of real parameter
                      ! with the birth from uniform
     integer :: n_iux ! Number of integer parameter
                      ! with the birth from uniform
     double precision, allocatable :: rgx(:, :) ! real Gaussian param.
     integer, allocatable          :: igx(:, :) ! int. Gaussian param.
     double precision, allocatable :: rux(:, :) ! real Uniform param.
     integer, allocatable          :: iux(:, :) ! int. Uniform param.
     double precision, allocatable :: sig_x_birth(:), mean_x_birth(:)
     double precision, allocatable :: sig_i_birth(:), mean_i_birth(:)
     double precision, allocatable :: x_min_birth(:), x_max_birth(:)
     integer         , allocatable :: i_min_birth(:), i_max_birth(:)

   contains
     procedure :: set_k => trans_d_model_set_k
     procedure :: set_x_gauss_birth => trans_d_model_set_x_gauss_birth
     procedure :: set_i_gauss_birth => trans_d_model_set_i_gauss_birth
     procedure :: set_x_uni_birth => trans_d_model_set_x_uni_birth
     procedure :: set_i_uni_birth => trans_d_model_set_i_uni_birth
     procedure :: generate_model => trans_d_model_generate_model
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
       & n_rgx, n_igx, n_rux, n_iux)
    integer, intent(in) :: k_min, k_max
    integer, intent(in), optional :: n_rgx, n_igx, n_rux, n_iux
    logical :: is_given

    ! Get dimension
    if (k_max <= k_min) then
       write(0,*)"ERROR: k_max must be > k_min"
       write(0,*)"     : k_min =", k_min, "k_max=", k_max, &
            & "(init_trans_d_model)"
       stop
    end if
    init_trans_d_model%k_min = k_min
    init_trans_d_model%k_max = k_max

    ! Get number of model parameters
    init_trans_d_model%n_rgx = 0
    init_trans_d_model%n_igx = 0
    init_trans_d_model%n_rux = 0
    init_trans_d_model%n_iux = 0
    is_given = .false.
    if (present(n_rgx)) then
       init_trans_d_model%n_rgx = n_rgx
       allocate(init_trans_d_model%rgx(n_rgx, k_max))
       allocate(init_trans_d_model%mean_x_birth(n_rgx))
       allocate(init_trans_d_model%sig_x_birth(n_rgx))
       is_given = .true.
    end if
    if (present(n_igx)) then
       init_trans_d_model%n_rgx = n_igx
       allocate(init_trans_d_model%igx(n_igx, k_max))
       allocate(init_trans_d_model%mean_i_birth(n_igx))
       allocate(init_trans_d_model%sig_i_birth(n_igx))
       is_given = .true.
    end if
    if (present(n_rux)) then
       init_trans_d_model%n_rux = n_rux
       allocate(init_trans_d_model%rux(n_rux, k_max))
       allocate(init_trans_d_model%x_min_birth(n_rux))
       allocate(init_trans_d_model%x_max_birth(n_rux))
       is_given = .true.
    end if
    if (present(n_iux)) then
       init_trans_d_model%n_iux = n_iux
       allocate(init_trans_d_model%iux(n_iux, k_max))
       allocate(init_trans_d_model%i_min_birth(n_iux))
       allocate(init_trans_d_model%i_max_birth(n_iux))
       is_given = .true.
    end if
    if (.not. is_given) then
       write(0,*)"ERROR: number of model parameters must be given"
       write(0,*)"     : (init_trans_d_model)"
       stop
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
    
    return 
  end subroutine trans_d_model_set_k

  !---------------------------------------------------------------------
  
  subroutine trans_d_model_set_x_gauss_birth(self, i, &
       & mean_x_birth, sig_x_birth)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: sig_x_birth, mean_x_birth
    
    if(i < 1 .or. i > self%k_max) then
       write(0,*) "ERROR: out of range (trans_d_model)"
       write(0,*) "     : i=", i
       stop
    end if
    self%sig_x_birth(i) = sig_x_birth
    self%mean_x_birth(i) = mean_x_birth
    
    return 
  end subroutine trans_d_model_set_x_gauss_birth

  !---------------------------------------------------------------------
  
  subroutine trans_d_model_set_i_gauss_birth(self, i, &
       & mean_i_birth, sig_i_birth)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: sig_i_birth,  mean_i_birth
    
    if(i < 1 .or. i > self%k_max) then
       write(0,*) "ERROR: out of range (trans_d_model)"
       write(0,*) "     : i=", i
       stop
    end if
    self%sig_i_birth(i) = sig_i_birth
    self%mean_i_birth(i) = mean_i_birth
    
    return 
  end subroutine trans_d_model_set_i_gauss_birth

  !---------------------------------------------------------------------
  
  subroutine trans_d_model_set_x_uni_birth(self, i, &
       & x_min_birth, x_max_birth)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: x_min_birth, x_max_birth
    
    if(i < 1 .or. i > self%k_max) then
       write(0,*) "ERROR: out of range (trans_d_model)"
       write(0,*) "     : i=", i
       stop
    end if
    if (x_min_birth >= x_max_birth) then
       write(0,*)"ERROR: x_min_birth must be < x_max_birth"
       write(0,*)"     : x_min_birth = ", x_min_birth
       write(0,*)"     : x_max_birth = ", x_max_birth
       write(0,*)"     : (set_x_uni_birth)"
       stop
    end if
    self%x_min_birth(i) = x_min_birth
    self%x_max_birth(i) = x_max_birth
    
    return 
  end subroutine trans_d_model_set_x_uni_birth
  
  !---------------------------------------------------------------------

  subroutine trans_d_model_set_i_uni_birth(self, i, &
       & i_min_birth, i_max_birth)
    class(trans_d_model), intent(inout) :: self
    integer, intent(in) :: i
    integer, intent(in) :: i_min_birth, i_max_birth
    
    if(i < 1 .or. i > self%k_max) then
       write(0,*) "ERROR: out of range (trans_d_model)"
       write(0,*) "     : i=", i
       stop
    end if
    if (i_min_birth >= i_max_birth) then
       write(0,*)"ERROR: i_min_birth must be < i_max_birth"
       write(0,*)"     : i_min_birth = ", i_min_birth
       write(0,*)"     : i_max_birth = ", i_max_birth
       write(0,*)"     : (set_i_uni_birth)"
       stop
    end if
    self%i_min_birth(i) = i_min_birth
    self%i_max_birth(i) = i_max_birth
    
    return 
  end subroutine trans_d_model_set_i_uni_birth
  
  !---------------------------------------------------------------------

  subroutine trans_d_model_generate_model(self)
    class(trans_d_model), intent(inout) :: self
    integer :: i, k
    logical :: is_ok
    
    ! select k
    k = self%k_min + int(rand_u() * (self%k_max - self%k_min))
    
    ! generate model parameters
    do i = 1, k
       call self%birth(is_ok)
       if (.not. is_ok) then
          write(0,*)"ERROR: something wrong occurs (generate_model)"
          stop
       end if
    end do

    return 
  end subroutine trans_d_model_generate_model
  
  !---------------------------------------------------------------------

  subroutine trans_d_model_birth(self, is_ok)
    class(trans_d_model), intent(inout) :: self
    logical, intent(out) :: is_ok
    integer :: i

    if (self%k < self%k_max - 1) then
       self%k = self%k + 1
       is_ok = .true.
    else
       is_ok = .false.
       return
    end if

    do i = 1, self%n_rgx
       self%rgx(i,  self%k) = &
            & rand_g() * self%sig_x_birth(i) + &
            & self%mean_x_birth(i)
    end do
    do i = 1, self%n_igx
       self%igx(i,  self%k) =                  &
            & nint(                            &
            & rand_g() * self%sig_i_birth(i) + &
            & self%mean_i_birth(i)             &
            & )
    end do
    do i = 1, self%n_rux
       self%rux(i,  self%k) = &
            & rand_u() *                                    &
            & (self%x_max_birth(i) - self%x_min_birth(i)) + &
            & self%x_min_birth(i)
    end do
    do i = 1, self%n_iux
       self%iux(i,  self%k) = &
            & int(rand_u() *                                 &
            & (self%i_max_birth(i) - self%i_min_birth(i))) + &
            & self%i_min_birth(i)
    end do
    
    return 
  end subroutine trans_d_model_birth

!---------------------------------------------------------------------

  subroutine trans_d_model_death(self, is_ok)
    class(trans_d_model), intent(inout) :: self
    logical, intent(out) :: is_ok
    integer :: i, j, k_target

    k_target = self%k_min + int(rand_u()*((self%k + 1) - self%k_min))

    
    if (self%k > self%k_min) then
       self%k = self%k - 1
       is_ok = .true.
    else
       is_ok = .false.
       return
    end if
    

    
    do i = k_target, self%k
       do j = 1, self%n_rgx
          self%rgx(j,  i) = self%rgx(j, i + 1) 
       end do
       do j = 1, self%n_igx
          self%igx(j,  i) = self%igx(j, i + 1) 
       end do
       do j = 1, self%n_rux
          self%rux(j,  i) = self%rux(j, i + 1) 
       end do
       do j = 1, self%n_iux
          self%iux(j,  i) = self%iux(j, i + 1) 
       end do
    end do

    return 
  end subroutine trans_d_model_death

  !---------------------------------------------------------------------

  subroutine trans_d_model_display(self)
    class(trans_d_model), intent(inout) :: self
    integer :: i, j
    write(*,*)
    write(*,*)"k = ", self%k
    
    do i = 1, self%n_rgx
       write(*,*)"------------------------------"
       write(*,*)"Real parameter with birth from Guassian", i
       do j = 1, self%k
          write(*,*)j, self%rgx(i, j)
       end do
    end do

    do i = 1, self%n_igx
       write(*,*)"------------------------------"
       write(*,*)"Integer prameter with birth from Gaussian", i
       do j = 1, self%k
          write(*,*)j, self%igx(i, j)
       end do
    end do

    do i = 1, self%n_rux
       write(*,*)"------------------------------"
       write(*,*)"Real parameter with birth from Uniform", i
       do j = 1, self%k
          write(*,*)j, self%rux(i, j)
       end do
    end do

    do i = 1, self%n_iux
       write(*,*)"------------------------------"
       write(*,*)"Integer prameter with birth from Uniform", i
       do j = 1, self%k
          write(*,*)j, self%iux(i, j)
       end do
    end do
                 
    
    write(*,*)
    
    return 
  end subroutine trans_d_model_display

end module mod_trans_d_model
  
