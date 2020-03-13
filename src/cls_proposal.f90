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
     
   contains
  end type proposal

  interface proposal
     module procedure :: init_proposal
  end interface proposal

contains
  
  !---------------------------------------------------------------------
  
  type(proposal) function init_proposal() & 
       & result(self)

    return 
  end function init_proposal
  
  !---------------------------------------------------------------------
  
end module cls_proposal
