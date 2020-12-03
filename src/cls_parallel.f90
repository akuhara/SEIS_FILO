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

module cls_parallel
  use cls_mcmc
  implicit none 
  include 'mpif.h'

  type parallel
     private
     integer :: n_proc 
     integer :: rank   
     integer :: n_chain
     type(mcmc), allocatable :: mc(:) 
     !type(trans_d_model), allocatable :: tm(:) 
     !type(hyper_model), allocatable :: hyp_rf(:)
     !type(hyper_model), allocatable :: hyp_disp(:)
     logical :: verb = .false.
   contains
     !procedure :: set_tm => parallel_set_tm 
     !procedure :: get_tm => parallel_get_tm
     !procedure :: set_hyp_rf => parallel_set_hyp_rf
     !procedure :: set_hyp_disp => parallel_set_hyp_disp
     !procedure :: get_hyp_rf => parallel_get_hyp_rf
     !procedure :: get_hyp_disp => parallel_get_hyp_disp
     procedure :: set_mc => parallel_set_mc
     procedure :: get_mc => parallel_get_mc

     procedure :: get_rank => parallel_get_rank
     procedure :: swap_temperature => parallel_swap_temperature
     procedure :: select_pair => parallel_select_pair
     procedure :: output_history => parallel_output_history
     procedure :: output_proposal => parallel_output_proposal

  end type parallel

  interface parallel
     module procedure init_parallel
  end interface parallel

  private pack_pair_info, unpack_pair_info, judge_swap
  private pack_mc_info, unpack_mc_info
  
contains
  
  !---------------------------------------------------------------------
  
  type(parallel) function init_parallel(n_proc, rank, n_chain, verb) &
       & result(self)
    integer, intent(in) :: n_proc, rank, n_chain
    logical, intent(in), optional :: verb

    if (present(verb)) then
       self%verb = verb
    end if

    if (self%verb) then
       write(*,'(A)')"<< Initialize prallel MCMC >>"
    end if

    self%n_proc = n_proc
    self%rank = rank
    self%n_chain = n_chain
    
    !allocate(self%tm(n_chain))
    !allocate(self%hyp_rf(n_chain))
    !allocate(self%hyp_disp(n_chain))
    allocate(self%mc(n_chain))

    if (self%verb) then
       write(*,'(A,I5,A,I5,A,I5)')"# of total MCMC chains = ", &
            & n_proc, " proc. * ", n_chain, " chains/proc. = ", &
            & n_chain * n_proc
       write(*,*)
    end if
    
    return 
  end function init_parallel
  
 !!---------------------------------------------------------------------
 !
 !subroutine parallel_set_tm(self, i, tm)
 !  class(parallel), intent(inout) :: self 
 !  integer, intent(in) :: i 
 !  type(trans_d_model), intent(in) :: tm 
 !  !integer :: ierr
 !  !if (i < 0 .or. i > self%n_chain) then
 !  !   write(0, *)"ERROR: in valid i (parallel_set_tm)"
 !  !   call mpi_finalize(ierr)
 !  !   stop
 !  !end if
 !  self%tm(i) = tm
 !  
 !  return 
 !end subroutine parallel_set_tm
 !
 !!---------------------------------------------------------------------
 !
 !type(trans_d_model) function parallel_get_tm(self, i) result(tm)
 !  class(parallel), intent(inout) :: self
 !  integer, intent(in) :: i
 !  !integer :: ierr
 !  !if (i < 0 .or. i > self%n_chain) then
 !  !   write(0, *)"ERROR: in valid i (parallel_get_tm)"
 !  !   call mpi_finalize(ierr)
 !  !   stop
 !  !end if
 !
 !  tm = self%tm(i)
 !
 !  return 
 !end function parallel_get_tm
 !
 !!---------------------------------------------------------------------
  
  !subroutine parallel_set_hyp_rf(self, i, hyp_rf)
  !  class(parallel), intent(inout) :: self 
  !  integer, intent(in) :: i 
  !  type(hyper_model), intent(in) :: hyp_rf
  !  
  !  self%hyp_rf(i) = hyp_rf
  !  
  !  return 
  !end subroutine parallel_set_hyp_rf
  !
  !---------------------------------------------------------------------
  !
  !subroutine parallel_set_hyp_disp(self, i, hyp_disp)
  !  class(parallel), intent(inout) :: self 
  !  integer, intent(in) :: i 
  !  type(hyper_model), intent(in) :: hyp_disp
  !  
  !  self%hyp_disp(i) = hyp_disp
  !  
  !  return 
  !end subroutine parallel_set_hyp_disp
  !
  !
  !---------------------------------------------------------------------
  !
  !type(hyper_model) function parallel_get_hyp_rf(self, i) result(hyp_rf)
  !  class(parallel), intent(inout) :: self
  !  integer, intent(in) :: i
  !  
  !  hyp_rf = self%hyp_rf(i)
  !
  !  return 
  !end function parallel_get_hyp_rf
  !
  !---------------------------------------------------------------------
  !
  !type(hyper_model) function parallel_get_hyp_disp(self, i) &
  !     & result(hyp_disp)
  !  class(parallel), intent(inout) :: self
  !  integer, intent(in) :: i
  !  
  !  hyp_disp = self%hyp_disp(i)
  !
  !  return 
  !end function parallel_get_hyp_disp

  !---------------------------------------------------------------------

  subroutine parallel_set_mc(self, i, mc)
    class(parallel), intent(inout) :: self
    integer, intent(in) :: i
    type(mcmc), intent(in) :: mc
    !integer :: ierr
    !if (i < 0 .or. i > self%n_chain) then
    !   write(0, *)"ERROR: in valid i (parallel_set_mc)"
    !   call mpi_finalize(ierr)
    !   stop
    !end if
    self%mc(i) = mc
    
    return 
  end subroutine parallel_set_mc

  !---------------------------------------------------------------------

  type(mcmc) function parallel_get_mc(self, i) result(mc)
    class(parallel), intent(inout) :: self
    integer, intent(in) :: i
    !integer :: ierr
    !if (i < 0 .or. i > self%n_chain) then
    !   write(0, *)"ERROR: in valid i (parallel_get_mc)"
    !   call mpi_finalize(ierr) 
    !   stop
    !end if

    mc = self%mc(i)

    return 
  end function parallel_get_mc
    
  !---------------------------------------------------------------------

  integer function parallel_get_rank(self) result(rank)
    class(parallel), intent(in) :: self

    rank = self%rank
    
    return 
  end function parallel_get_rank

  !---------------------------------------------------------------------

  subroutine parallel_swap_temperature(self, verb)
    class(parallel), intent(inout) :: self
    logical, intent(in), optional :: verb
    integer :: ipack(4), ierr, status(MPI_STATUS_SIZE)
    integer :: rank1, rank2, chain1, chain2
    double precision :: temp1, temp2, l1, l2, rpack(2)
    type(mcmc) :: mc1, mc2
    logical :: is_accepted

    if (self%rank == 0) then
       ipack = self%select_pair()
    end if
    call mpi_bcast(ipack, 4, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
       write(0,*)"ERROR: while MPI_BCAST"
       call mpi_finalize(ierr)
       stop
    end if
    call unpack_pair_info(ipack=ipack, rank1=rank1, rank2=rank2, &
         & chain1=chain1, chain2=chain2)
    
    if (self%rank == rank1 .and. self%rank == rank2) then
       mc1   = self%get_mc(chain1)
       temp1 = mc1%get_temp()
       l1    = mc1%get_log_likelihood()
       mc2   = self%get_mc(chain2)
       temp2 = mc2%get_temp()
       l2    = mc2%get_log_likelihood()
       
       call judge_swap(temp1, temp2, l1, l2, is_accepted)
       
       if (is_accepted) then
          call mc2%set_temp(temp1)
          call mc1%set_temp(temp2)
          call self%set_mc(chain1, mc1)
          call self%set_mc(chain2, mc2)
       end if
       
       !if (present(verb) .and. verb) then
       !   write(*,*)"----"
       !   write(*,*)"Rank1    :", rank1
       !   write(*,*)"Chain1   :", chain1
       !   write(*,*)"Temp1    : ", temp1
       !   write(*,*)"Rank2    :", rank2 
       !   write(*,*)"Chain2   :", chain2
       !   write(*,*)"Temp2    : ", temp2
       !   write(*,*)"Accepted :", is_accepted
       !   write(*,*)"----"
       !end if
       
    else if (self%rank == rank1) then
       mc1   = self%get_mc(chain1)
       temp1 = mc1%get_temp()
       l1    = mc1%get_log_likelihood()
       call mpi_recv(rpack, 2, MPI_DOUBLE_PRECISION, rank2, 1111, &
            & MPI_COMM_WORLD, status, ierr)
       if (ierr /= MPI_SUCCESS) then
          write(0,*)"ERROR: while MPI_RECV rpack 1"
          call mpi_finalize(ierr)
          stop
       end if
       call unpack_mc_info(rpack, temp=temp2, likelihood=l2)

       call judge_swap(temp1, temp2, l1, l2, is_accepted)
       
       if (is_accepted) then
          rpack = pack_mc_info(temp=temp1, likelihood=-999.d0)
          call mc1%set_temp(temp2)
          call self%set_mc(chain1, mc1)
       else
          rpack = pack_mc_info(temp=temp2, likelihood=-999.d0)
       end if
       call mpi_send(rpack, 2, MPI_DOUBLE_PRECISION, rank2, 2222, &
            & MPI_COMM_WORLD, status, ierr)
       if (ierr /= MPI_SUCCESS) then
          write(0,*)"ERROR: while MPI_SEND rpack 1"
          call mpi_finalize(ierr)
          stop
       end if
       !if (present(verb) .and. verb) then
       !   write(*,*)"----"
       !   write(*,*)"Rank1    :", rank1
       !   write(*,*)"Chain1   :", chain1
       !   write(*,*)"Temp1    : ", temp1
       !   write(*,*)"Rank2    :", rank2 
       !   write(*,*)"Chain2   :", chain2
       !   write(*,*)"Temp2    : ", temp2
       !   write(*,*)"Accepted :", is_accepted
       !   write(*,*)"----"
       !end if
    else if (self%rank == rank2) then
       ! Sender
       mc2   = self%get_mc(chain2)
       temp2 = mc2%get_temp()
       l2    = mc2%get_log_likelihood()
       rpack = pack_mc_info(temp=temp2, likelihood=l2)
       call mpi_send(rpack, 2, MPI_DOUBLE_PRECISION, rank1, 1111, &
            & MPI_COMM_WORLD, status, ierr)
       if (ierr /= MPI_SUCCESS) then
          write(0,*)"ERROR: while MPI_SEND rpack 2"
          call mpi_finalize(ierr)
          stop
       end if
       call mpi_recv(rpack, 2, MPI_DOUBLE_PRECISION, rank1, 2222, &
            & MPI_COMM_WORLD, status, ierr)
       if (ierr /= MPI_SUCCESS) then
          write(0,*)"ERROR: while MPI_RECV rpack 2"
          call mpi_finalize(ierr)
          stop
       end if
       call unpack_mc_info(rpack, temp=temp2, likelihood=l2)
       call mc2%set_temp(temp2)
       call self%set_mc(chain2, mc2)
    end if
    
    return 
  end subroutine parallel_swap_temperature
    
  !---------------------------------------------------------------------

  function parallel_select_pair(self) result(ipack)
    class(parallel), intent(in) :: self
    integer :: ipack(4)
    integer :: rank1, rank2, chain1, chain2
    integer :: i1, i2
    
    i1 = int(rand_u() * self%n_proc * self%n_chain)
    do 
       i2 = int(rand_u() * self%n_proc * self%n_chain)
       if (i1 /= i2) exit
    end do
    rank1 = int(i1 / self%n_chain)
    rank2 = int(i2 / self%n_chain)
    chain1 = mod(i1, self%n_chain) + 1
    chain2 = mod(i2, self%n_chain) + 1
    
    ipack = pack_pair_info(rank1=rank1, rank2=rank2, &
         & chain1=chain1, chain2=chain2)
    
    return 
  end function parallel_select_pair

  !---------------------------------------------------------------------

  subroutine parallel_output_history(self, filename, mode)
    class(parallel), intent(in) :: self
    character(*), intent(in) :: filename
    character(*), intent(in) :: mode
    integer :: i, ierr, io
    type(mcmc) :: mc
    integer :: icol, n_all, n_iter
    double precision, allocatable :: hist_all(:,:), hist_all2(:,:)
    character(50) :: fmt
    write(*,*)self%n_chain,  self%n_proc
    n_all = self%n_chain * self%n_proc
    mc = self%mc(1)
    n_iter = mc%get_n_iter()
    allocate(hist_all(n_iter, n_all), hist_all2(n_iter, n_all))

    ! Gather information within the same node
    do i = 1, self%n_chain
       mc = self%mc(i)
       if (mode == "l") then
          hist_all(1:n_iter, i) = mc%get_likelihood_saved()
       else if (mode == "t") then
          hist_all(1:n_iter, i) = mc%get_temp_saved() 
       end if
    end do
    
    ! MPI gather
    do i = 1, self%n_chain
       icol = (i - 1) * self%n_proc + 1
       call mpi_gather(hist_all(1:n_iter, i), n_iter, &
            & MPI_DOUBLE_PRECISION, &
            & hist_all2(1:n_iter, icol), n_iter, &
            & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end do
    
    ! Output
    if (self%rank == 0) then
       open(newunit = io, file = filename, status = "unknown", &
            & iostat = ierr)
       if (ierr /= 0) then
          write(0,*)"ERROR: cannot create ", trim(filename)
          call mpi_finalize(ierr)
          stop
       end if
       do i = 1, n_iter
          write(fmt, '("(",I0,"E16.4)")')n_all
          write(io, trim(fmt))hist_all2(i, 1:n_all)
       end do
       close(io)
    end if

    return 
  end subroutine parallel_output_history
  
  !---------------------------------------------------------------------

  subroutine parallel_output_proposal(self, filename, label)
    class(parallel), intent(inout) :: self
    character(*), intent(in) :: filename, label(:)
    type(mcmc) :: mc
    integer, allocatable :: n_accept_all(:), n_propose_all(:)
    integer, allocatable :: n_accept_sum(:), n_propose_sum(:)
    integer :: io, i, n, ierr
   
    n = size(label)
    allocate(n_accept_all(n), n_propose_all(n))
    allocate(n_accept_sum(n), n_propose_sum(n))
    n_accept_all(:) = 0
    n_propose_all(:) = 0

    ! Gather within the same node
    do i = 1, self%n_chain
       mc = self%get_mc(i)
       n_accept_all = n_accept_all + mc%get_n_accept()
       n_propose_all = n_propose_all + mc%get_n_propose()
    end do
    
    call mpi_reduce(n_accept_all, n_accept_sum, n, MPI_INTEGER4, MPI_SUM, &
         & 0, MPI_COMM_WORLD, ierr)
    call mpi_reduce(n_propose_all, n_propose_sum, n, MPI_INTEGER4, MPI_SUM, &
         & 0, MPI_COMM_WORLD, ierr)
    ! Output
    if (self%rank == 0) then
       open(newunit = io, file = filename, status = "unknown", &
            & iostat = ierr)
       do i = 1, n
          write(io, '(A,2I10)')'"' // label(i) // '"', &
               & n_propose_sum(i), n_accept_sum(i)
       end do
       close(io)
    end if
    
    return 
  end subroutine parallel_output_proposal
  
  !---------------------------------------------------------------------

  subroutine judge_swap(temp1, temp2, l1, l2, is_accepted)
    double precision, intent(in) :: temp1, temp2, l1, l2
    logical, intent(out) :: is_accepted
    double precision :: del_s
    double precision :: r
    double precision, parameter :: eps = epsilon(1.d0)
    
    del_s = (l2 - l1) * (1.d0 / temp1 - 1.d0 / temp2)
    is_accepted = .false.
    r = rand_u()
    if (r >= eps) then
       if(log(r) <= del_s) then
          is_accepted = .true.
       end if
    end if
    
    return 
  end subroutine judge_swap
  
  !---------------------------------------------------------------------

  function pack_pair_info(rank1, rank2, chain1, chain2) result(ipack)
    integer, intent(in) :: rank1, rank2, chain1, chain2
    integer :: ipack(4)
    
    ipack(1) = rank1
    ipack(2) = rank2
    ipack(3) = chain1
    ipack(4) = chain2
    
    return 
  end function pack_pair_info

  !---------------------------------------------------------------------

  subroutine unpack_pair_info(ipack, rank1, rank2, chain1, chain2)
    integer, intent(in) ::ipack(4)
    integer, intent(out) :: rank1, rank2, chain1, chain2
    
    rank1 = ipack(1)
    rank2 = ipack(2)
    chain1 = ipack(3)
    chain2 = ipack(4)
    
    return 
  end subroutine unpack_pair_info
  
  !---------------------------------------------------------------------

  function pack_mc_info(temp, likelihood) result(rpack)
    double precision, intent(in) :: temp, likelihood
    double precision :: rpack(2)
    
    rpack(1) = temp
    rpack(2) = likelihood

    return 
  end function pack_mc_info

  !---------------------------------------------------------------------

  subroutine unpack_mc_info(rpack, temp, likelihood)
    double precision, intent(in) :: rpack(2)
    double precision, intent(out) :: temp, likelihood
    
    temp = rpack(1)
    likelihood = rpack(2)

    return 
  end subroutine unpack_mc_info

  !---------------------------------------------------------------------


end module cls_parallel
