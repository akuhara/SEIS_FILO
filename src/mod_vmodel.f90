module mod_vmodel
  implicit none 
  
  type vmodel
     private
     integer :: nlay
     double precision, allocatable :: vp(:), vs(:), rho(:), h(:)
   
   contains
     ! setter
     procedure :: set_nlay => vmodel_set_nlay
     procedure :: set_vp => vmodel_set_vp
     procedure :: set_vs => vmodel_set_vs
     procedure :: set_rho => vmodel_set_rho
     procedure :: set_h => vmodel_set_h
     ! getter
     procedure :: get_nlay => vmodel_get_nlay
     procedure :: get_vp => vmodel_get_vp
     procedure :: get_vs => vmodel_get_vs
     procedure :: get_rho => vmodel_get_rho
     procedure :: get_h => vmodel_get_h
     ! utilities
     procedure :: set_example_ocean => vmodel_set_example_ocean
     procedure :: set_example_land => vmodel_set_example_land
     procedure :: display => vmodel_display
     procedure :: vp2rho_brocher => vmodel_vp2rho_brocher
     procedure :: vp2vs_brocher => vmodel_vp2vs_brocher
     procedure :: vs2vp_brocher => vmodel_vs2vp_brocher

  end type vmodel
  
  interface vmodel
     module procedure init_vmodel
  end interface vmodel
  
contains
  
  !---------------------------------------------------------------------

  type(vmodel) function init_vmodel()
    
    init_vmodel%nlay = -999
    
    return 
  end function init_vmodel
  
  !---------------------------------------------------------------------
  
  subroutine vmodel_set_nlay(self, nlay)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: nlay

    ! Reset if already set
    if (self%nlay > 0) then
       deallocate(self%vp, self%vs, self%rho, self%h)
    end if
    
    self%nlay = nlay
    allocate(self%vp(nlay), self%vs(nlay), self%rho(nlay), self%h(nlay))
    self%vp(1:nlay) = -999.d0
    self%vs(1:nlay) = -999.d0
    self%rho(1:nlay) = -999.d0
    self%h(1:nlay) = -999.d0
    
    return 
  end subroutine vmodel_set_nlay

  !---------------------------------------------------------------------

  subroutine vmodel_set_vp(self, i, vp)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: vp
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (set_vp)"
       write(0,*) "     : i=", i
       stop
    end if
    if (vp < 0) then
       write(0,*) "ERROR: invalid Vp (set_vp)"
       write(0,*) "     : Vp=", vp
       stop
    end if
    self%vp(i) = vp

    return 
  end subroutine vmodel_set_vp
  
  !---------------------------------------------------------------------

  subroutine vmodel_set_vs(self, i, vs)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: vs
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (set_vs)"
       write(0,*) "     : i=", i
       stop
    end if
    if (vs < 0 .and. i /= 1) then
       write(0,*) "ERROR: invalid Vs (set_vs)"
       write(0,*) "     : Vs=", vs
       stop
    end if
    self%vs(i) = vs

    return 
  end subroutine vmodel_set_vs
  
  !---------------------------------------------------------------------
  
  subroutine vmodel_set_rho(self, i, rho)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: rho
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (set_rho)"
       write(0,*) "     : i=", i
       stop
    end if
    if (rho < 0) then
       write(0,*) "ERROR: invalid density (set_rho)"
       write(0,*) "     : rho=", rho
       stop
    end if

    self%rho(i) = rho

    return 
  end subroutine vmodel_set_rho
  
  !---------------------------------------------------------------------

  subroutine vmodel_set_h(self, i, h)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision, intent(in) :: h
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (set_h)"
       write(0,*) "     : i=", i
       stop
    end if
    if (h < 0 .and. i /= self%nlay) then
       write(0,*) "ERROR: invalid thickness (set_h)"
       write(0,*) "     : h=", h
       stop
    end if
    self%h(i) = h

    return 
  end subroutine vmodel_set_h

  !---------------------------------------------------------------------
  
  integer function vmodel_get_nlay(self) result(nlay)
    class(vmodel), intent(inout) :: self
    
    nlay = self%nlay
    
    return 
  end function vmodel_get_nlay

  !---------------------------------------------------------------------
  
  double precision function vmodel_get_vp(self, i) result(vp)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (get_vp)"
       write(0,*) "     : i=", i
       stop
    end if
    vp = self%vp(i)
    
    return 
  end function vmodel_get_vp
  
  !---------------------------------------------------------------------
  
  double precision function vmodel_get_vs(self, i) result(vs)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (get_vs)"
       write(0,*) "     : i=", i
       stop
    end if
    vs = self%vs(i)
    
    return 
  end function vmodel_get_vs
  
  !---------------------------------------------------------------------

  double precision function vmodel_get_rho(self, i) result(rho)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (get_rho)"
       write(0,*) "     : i=", i
       stop
    end if
    rho = self%rho(i)
    
    return 
  end function vmodel_get_rho

  !---------------------------------------------------------------------
  
  double precision function vmodel_get_h(self, i) result(h)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (get_h)"
       write(0,*) "     : i=", i
       stop
    end if
    h = self%h(i)
    
    return 
  end function vmodel_get_h
  
  !---------------------------------------------------------------------

  subroutine vmodel_display(self)
    class(vmodel), intent(inout) :: self
    integer :: i

    print *
    do i = 1, self%nlay
       print *, i, self%vp(i), self%vs(i), self%rho(i), self%h(i)
    end do
    print *
    
    return 
  end subroutine vmodel_display
  
  !---------------------------------------------------------------------
  
  subroutine vmodel_set_example_ocean(self)
    class(vmodel), intent(inout) :: self
    
    call self%set_nlay(3)
    
    ! Ocean
    call self%set_vp(1, 1.5d0)
    call self%set_vs(1, -1.d0)
    call self%set_rho(1, 1.d0)
    call self%set_h(1, 3.d0)

    ! Sediment
    call self%set_vp(2, 1.6d0)
    call self%vp2vs_brocher(2)
    call self%vp2rho_brocher(2)
    call self%set_h(2, 1.d0)

    ! Basement
    call self%set_vp(3, 3.d0)
    call self%vp2vs_brocher(3)
    call self%vp2rho_brocher(3)
    call self%set_h(3, 100.d0)

    return 
  end subroutine vmodel_set_example_ocean
  
  !---------------------------------------------------------------------

  subroutine vmodel_set_example_land(self)
    class(vmodel), intent(inout) :: self
    
    call self%set_nlay(3)
    
    ! Ocean
    call self%set_vp(1, 0.255d0)
    call self%set_vs(1, 0.150d0)
    call self%set_rho(1, 1.00d0)
    call self%set_h(1, 0.002d0)

    ! Sediment
    call self%set_vp(2, 0.340d0)
    call self%set_vs(2, 0.200d0)
    call self%set_rho(2, 1.00d0)
    call self%set_h(2, 0.002d0)

    ! Basement
    call self%set_vp(3, 0.510d0)
    call self%set_vs(3, 0.300d0)
    call self%set_rho(3, 1.00d0)
    call self%set_h(3, 100.d0)

    return 
  end subroutine vmodel_set_example_land
  
  !---------------------------------------------------------------------
  
  subroutine vmodel_vp2rho_brocher(self, i)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision :: a1, a2, a3, a4, a5
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (vp2rho_brocher)"
       write(0,*) "     : i=", i
       stop
    end if
    
    a1 = self%vp(i)
    a2 = a1 * a1
    a3 = a2 * a1
    a4 = a3 * a1
    a5 = a4 * a1

    self%rho(i) = &
         & 1.6612d0 * a1 - 0.4721d0 * a2 + 0.0671d0 * a3 - &
         & 0.0043d0 * a4 + 0.000106d0 * a5
    
    return 
  end subroutine vmodel_vp2rho_brocher
    
  !---------------------------------------------------------------------

  subroutine vmodel_vp2vs_brocher(self, i)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision :: a1, a2, a3, a4
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (vs2rho_brocher)"
       write(0,*) "     : i=", i
       stop
    end if
    
    a1 = self%vp(i)
    a2 = a1 * a1
    a3 = a2 * a1
    a4 = a3 * a1

    self%vs(i) = &
         & 0.7858d0 - 1.2344d0 * a1 + 0.7949d0 * a2 - &
         & 0.1238d0 * a3 + 0.0064d0 * a4
    
    return 
  end subroutine vmodel_vp2vs_brocher
    
  !---------------------------------------------------------------------

  subroutine vmodel_vs2vp_brocher(self, i)
    class(vmodel), intent(inout) :: self
    integer, intent(in) :: i
    double precision :: a1, a2, a3, a4
    
    if(i < 1 .or. i > self%nlay) then
       write(0,*) "ERROR: out of range (vs2rho_brocher)"
       write(0,*) "     : i=", i
       stop
    end if
    
    a1 = self%vs(i)
    a2 = a1 * a1
    a3 = a2 * a1
    a4 = a3 * a1

    self%vp(i) = &
         & 0.9409d0 + 2.0947d0 * a1 - 0.8206d0 * a2 + &
         & 0.2683d0 * a3 - 0.0251d0 * a4
    
    return 
  end subroutine vmodel_vs2vp_brocher
    
  !---------------------------------------------------------------------

end module mod_vmodel
