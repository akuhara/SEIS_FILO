!=======================================================================
!   SEIS_FILO: 
!   SEISmological tools for Flat Isotropic Layered structure in the Ocean
!   Copyright (C) 2019 Takeshi Akuhara
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!
!   Contact information
!
!   Email  : akuhara @ eri. u-tokyo. ac. jp 
!   Address: Earthquake Research Institute, The Univesity of Tokyo   
!           1-1-1, Yayoi, Bunkyo-ku, Tokyo 113-0032, Japan
!
!=======================================================================
module cls_proposal
  implicit none
  
  type proposal
     private
     logical :: solve_vp
     logical :: solve_rf_sig
     logical :: solve_disper_sig
     integer :: n_rf
     integer :: n_disp
    
     integer :: n_proposal
     integer :: i_birth = -999
     integer :: i_death = -999
     integer :: i_depth = -999
     integer :: i_vs    = -999
     integer :: i_vp    = -999
     integer :: i_rf_sig = -999
     integer :: i_disper_sig = -999
     
     character(16), allocatable :: label(:)
     logical :: verb = .false.

   contains
     procedure :: get_n_proposal => proposal_get_n_proposal
     procedure :: get_i_birth => proposal_get_i_birth
     procedure :: get_i_death => proposal_get_i_death
     procedure :: get_i_depth => proposal_get_i_depth
     procedure :: get_i_vs => proposal_get_i_vs
     procedure :: get_i_vp => proposal_get_i_vp
     procedure :: get_i_rf_sig => proposal_get_i_rf_sig
     procedure :: get_i_disper_sig => proposal_get_i_disper_sig
     procedure :: get_label => proposal_get_label
     
  end type proposal

  interface proposal
     module procedure :: init_proposal
  end interface proposal

contains
  
  !---------------------------------------------------------------------
  
  type(proposal) function init_proposal(solve_vp, solve_rf_sig, &
       & solve_disper_sig, n_rf, n_disp, verb) & 
       & result(self)
    logical, intent(in) :: solve_vp, solve_rf_sig, solve_disper_sig
    integer, intent(in) :: n_rf, n_disp
    logical, intent(in), optional :: verb
    integer :: i

    if (present(verb)) then
       self%verb = verb
    end if
    if (self%verb) then
       write(*,'(a)')"<< Initialize proposals >>"
    end if

    self%n_proposal = 1 ! Birth, Death, Depth, & Vs
    self%i_depth = self%n_proposal
    self%n_proposal = self%n_proposal + 1
    self%i_vs    = self%n_proposal
    
    self%solve_vp = solve_vp
    if (self%solve_vp) then
       self%n_proposal = self%n_proposal + 1
       self%i_vp = self%n_proposal
    end if
    
    self%n_proposal = self%n_proposal + 1
    self%i_birth = self%n_proposal

    self%n_proposal = self%n_proposal + 1
    self%i_death = self%n_proposal

    self%solve_rf_sig = solve_rf_sig
    self%n_rf = n_rf
    if (self%solve_rf_sig .and. self%n_rf > 0) then
       self%n_proposal = self%n_proposal + 1
       self%i_rf_sig = self%n_proposal
    end if

    self%solve_disper_sig = solve_disper_sig
    self%n_disp = n_disp
    if (self%solve_disper_sig .and. self%n_disp > 0) then
       self%n_proposal = self%n_proposal + 1
       self%i_disper_sig = self%n_proposal
    end if

    ! Make label
    allocate(self%label(self%n_proposal))
    do i = 1, self%n_proposal
       if (i == self%i_birth) then
          self%label(i) = "Birth"
       else if (i == self%i_death) then
          self%label(i) = "Death"
       else if (i == self%i_depth) then
          self%label(i) = "Depth"
       else if (i == self%i_vs) then
          self%label(i) = "Vs"
       else if (i == self%i_vp) then
          self%label(i) = "Vp"
       else if (i == self%i_rf_sig) then
          self%label(i) = "RF sigma"
       else if (i == self%i_disper_sig) then
          self%label(i) = "Dispersion sigma"
       end if
    end do

    if (self%verb) then
       write(*,'(a)')" Proposal type"
       do i = 1, self%n_proposal
          write(*,'(I3,4x,A)')i, self%label(i)
       end do
       write(*,*)
    end if
    
    return 
  end function init_proposal
  
  !---------------------------------------------------------------------
  
  integer function proposal_get_n_proposal(self) result(n_proposal)
    class(proposal), intent(in) :: self

    n_proposal = self%n_proposal
    
    return 
  end function proposal_get_n_proposal
    
  !---------------------------------------------------------------------
  
  integer function proposal_get_i_birth(self) result(i_birth)
    class(proposal), intent(in) :: self

    i_birth = self%i_birth
    
    return 
  end function proposal_get_i_birth
    
  !---------------------------------------------------------------------

  integer function proposal_get_i_death(self) result(i_death)
    class(proposal), intent(in) :: self

    i_death = self%i_death
    
    return 
  end function proposal_get_i_death
    
  !---------------------------------------------------------------------

  integer function proposal_get_i_depth(self) result(i_depth)
    class(proposal), intent(in) :: self

    i_depth = self%i_depth
    
    return 
  end function proposal_get_i_depth
    
  !---------------------------------------------------------------------

  integer function proposal_get_i_vs(self) result(i_vs)
    class(proposal), intent(in) :: self

    i_vs = self%i_vs
    
    return 
  end function proposal_get_i_vs
    
  !---------------------------------------------------------------------
  
  integer function proposal_get_i_vp(self) result(i_vp)
    class(proposal), intent(in) :: self

    i_vp = self%i_vp
    
    return 
  end function proposal_get_i_vp
    
  !---------------------------------------------------------------------

  integer function proposal_get_i_rf_sig(self) result(i_rf_sig)
    class(proposal), intent(in) :: self

    i_rf_sig = self%i_rf_sig
    
    return 
  end function proposal_get_i_rf_sig
    
  !---------------------------------------------------------------------

  integer function proposal_get_i_disper_sig(self) result(i_disper_sig)
    class(proposal), intent(in) :: self

    i_disper_sig = self%i_disper_sig
    
    return 
  end function proposal_get_i_disper_sig
    
  !---------------------------------------------------------------------

  character(16) function proposal_get_label(self, i) &
       & result(label)
    class(proposal), intent(in) :: self
    integer, intent(in) :: i
    
    label = self%label(i)

    return 
  end function proposal_get_label

  !---------------------------------------------------------------------
  
end module cls_proposal
