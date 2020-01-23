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
!> Module to perform parallel tempering MCMC
!> @author 
!> Takeshi Akuhara
module mod_parallel
  use mod_trans_d_model
  use mod_mcmc
  implicit none 
  include 'mpif.h'

  !---------------------------------------------------------------------
  ! Parallel
  !---------------------------------------------------------------------
  !> @brief
  !! Class to control between-chain interaction 
  !! (i.e., parallel tempering)
  !! @param n_proc  # of MCMC chains 
  !! @param rank    Process ID
  !! @param n_cahin # of MCMC chains per process
  !! @param mc      A struct that contains MCMC chain property
  !!                (see mod_mcmc.f90 for details)
  !! @param tm      A structure that contains 
  !!                transdimensional model parameters and related 
  !!                setting (see mod_trans_d_model.f90 for details)  
  !---------------------------------------------------------------------
  type parallel
     private
     integer :: n_proc 
     integer :: rank   
     integer :: n_chain
     type(mcmc), allocatable :: mc(:) 
     type(trans_d_model), &
          & allocatable :: tm(:) 
   contains
     procedure :: set_tm => parallel_set_tm 
     procedure :: get_tm => parallel_get_tm
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
  ! init_parallel
  !---------------------------------------------------------------------
  !> @brief
  !! Constructor of class 'parallel', which allocates memory for 'tm(:)'
  !! and 'mc(:)'.
  !! @param[in] n_proc   # of processes
  !! @param[in] rank     Process ID
  !! @param[in] n_chains # of MCMC chain per process
  !! @return Object
  !---------------------------------------------------------------------
  type(parallel) function init_parallel(n_proc, rank, n_chain) &
       & result(self)
    integer, intent(in) :: n_proc, rank, n_chain
    self%n_proc = n_proc
    self%rank = rank
    self%n_chain = n_chain
    
    allocate(self%tm(n_chain))
    allocate(self%mc(n_chain))
    
    return 
  end function init_parallel
  
  !---------------------------------------------------------------------
  ! parallel_set_tm 
  !---------------------------------------------------------------------
  !> @brief
  !! Set trandimensional model parameters for 'i'th MCMC chain
  !! @param[inout] self Object
  !! @param[in]    i    MCMC chain ID
  !! @param[in]    tm   A structure that contains 
  !!                    transdimensional model parameters and related 
  !!                    settings
  !!                    (see mod_trans_d_model.f90 for details)  
  !---------------------------------------------------------------------
  subroutine parallel_set_tm(self, i, tm)
    class(parallel), intent(inout) :: self 
    integer, intent(in) :: i 
    type(trans_d_model), intent(in) :: tm 
    
    if (i < 0 .or. i > self%n_chain) then
       write(0, *)"ERROR: in valid i (parallel_set_tm)"
       stop
    end if
    self%tm(i) = tm
    
    return 
  end subroutine parallel_set_tm

  !---------------------------------------------------------------------
  ! parallel_get_tm
  !---------------------------------------------------------------------
  !> @brief
  !! Get trandimensional model parameters of 'i'th MCMC chain
  !! @param[in]    self Object
  !! @param[in]    i    MCMC chain ID
  !!                    
  !! @return       A structure that contains 
  !!               transdimensional model parameters and related 
  !!               settings (see mod_trans_d_model.f90 for details)  
  !---------------------------------------------------------------------
  type(trans_d_model) function parallel_get_tm(self, i) result(tm)
    class(parallel), intent(inout) :: self
    integer, intent(in) :: i
    if (i < 0 .or. i > self%n_chain) then
       write(0, *)"ERROR: in valid i (parallel_get_tm)"
       stop
    end if

    tm = self%tm(i)

    return 
  end function parallel_get_tm

  !---------------------------------------------------------------------
  ! parallel_set_mc
  !---------------------------------------------------------------------
  !> @brief
  !! Set MCMC properties of 'i'th chain
  !! @param[inout] self Object
  !! @param[in]    i    MCMC chain ID
  !! @param[in]    mc   A structure that contains MCMC properties
  !!                    (see mod_mcmc.f90 for details)
  !---------------------------------------------------------------------
  subroutine parallel_set_mc(self, i, mc)
    class(parallel), intent(inout) :: self
    integer, intent(in) :: i
    type(mcmc), intent(in) :: mc
    
    if (i < 0 .or. i > self%n_chain) then
       write(0, *)"ERROR: in valid i (parallel_set_mc)"
       stop
    end if
    self%mc(i) = mc
    
    return 
  end subroutine parallel_set_mc

  !---------------------------------------------------------------------
  ! parallel_get_mc
  !---------------------------------------------------------------------
  !> @brief
  !! Get MCMC properties of 'i'th chain
  !! @param[in] self Object
  !! @param[in] i    MCMC chain ID
  !! @return         A structure that contains MCMC properties
  !!                 (see mod_mcmc.f90 for details)
  !---------------------------------------------------------------------
  type(mcmc) function parallel_get_mc(self, i) result(mc)
    class(parallel), intent(inout) :: self
    integer, intent(in) :: i
    if (i < 0 .or. i > self%n_chain) then
       write(0, *)"ERROR: in valid i (parallel_get_mc)"
       stop
    end if

    mc = self%mc(i)

    return 
  end function parallel_get_mc
    
  !---------------------------------------------------------------------
  ! parallel_get_rank
  !---------------------------------------------------------------------
  !> @brief
  !! Get process ID
  !! @param[in] self Object
  !! @return    process ID     
  !---------------------------------------------------------------------
  integer function parallel_get_rank(self) result(rank)
    class(parallel), intent(in) :: self

    rank = self%rank
    
    return 
  end function parallel_get_rank

  !---------------------------------------------------------------------
  ! parallel_swap_temperature
  !---------------------------------------------------------------------
  !> @brief
  !! Judge whether a temperature swap should occucr or not. 
  !! If the swap proposal is accepted, 
  !! this routine also perform the swap
  !! @param[inout] self Object
  !! @param[in]    verb Verbose flag      
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
       
       if (present(verb) .and. verb) then
          write(*,*)"----"
          write(*,*)"Rank1    :", rank1
          write(*,*)"Chain1   :", chain1
          write(*,*)"Temp1    : ", temp1
          write(*,*)"Rank2    :", rank2 
          write(*,*)"Chain2   :", chain2
          write(*,*)"Temp2    : ", temp2
          write(*,*)"Accepted :", is_accepted
          write(*,*)"----"
       end if
       
    else if (self%rank == rank1) then
       mc1   = self%get_mc(chain1)
       temp1 = mc1%get_temp()
       l1    = mc1%get_log_likelihood()
       call mpi_recv(rpack, 2, MPI_DOUBLE_PRECISION, rank2, 1111, &
            & MPI_COMM_WORLD, status, ierr)
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
       if (present(verb) .and. verb) then
          write(*,*)"----"
          write(*,*)"Rank1    :", rank1
          write(*,*)"Chain1   :", chain1
          write(*,*)"Temp1    : ", temp1
          write(*,*)"Rank2    :", rank2 
          write(*,*)"Chain2   :", chain2
          write(*,*)"Temp2    : ", temp2
          write(*,*)"Accepted :", is_accepted
          write(*,*)"----"
       end if
    else if (self%rank == rank2) then
       ! Sender
       mc2   = self%get_mc(chain2)
       temp2 = mc2%get_temp()
       l2    = mc2%get_log_likelihood()
       rpack = pack_mc_info(temp=temp2, likelihood=l2)
       call mpi_send(rpack, 2, MPI_DOUBLE_PRECISION, rank1, 1111, &
            & MPI_COMM_WORLD, status, ierr)
       call mpi_recv(rpack, 2, MPI_DOUBLE_PRECISION, rank1, 2222, &
            & MPI_COMM_WORLD, status, ierr)
       call unpack_mc_info(rpack, temp=temp2, likelihood=l2)
       call mc2%set_temp(temp2)
       call self%set_mc(chain2, mc2)
    end if
    
    return 
  end subroutine parallel_swap_temperature
    
  !---------------------------------------------------------------------
  ! parallel_select_pair
  !---------------------------------------------------------------------
  !> @brief
  !! Select pair of MCMC chains randomly
  !! @param[in] self Object
  !! @return    Process ID 1, Process ID 2, Chain ID 1, Chain ID 2
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
  ! parallel_output_history
  !---------------------------------------------------------------------
  !> @brief
  !! Output likelihood/temperature history of MCMC chains
  !! @param[in] self     Object
  !! @param[in] filename Output file name
  !! @param[in] mode     Output mode (l: likelihood, t: temperature)
  !---------------------------------------------------------------------
  subroutine parallel_output_history(self, filename, mode)
    class(parallel), intent(in) :: self
    character(*), intent(in) :: filename
    character(*), intent(in) :: mode
    integer :: i, ierr, io
    type(mcmc) :: mc
    integer :: icol, n_all, n_iter
    double precision, allocatable :: hist_all(:,:)

    
    n_all = self%n_chain * self%n_proc
    mc = self%mc(1)
    n_iter = mc%get_n_iter()
    allocate(hist_all(n_iter, n_all))

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
            & hist_all(1:n_iter, icol), n_iter, &
            & MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    end do
    
    ! Output
    if (self%rank == 0) then
       open(newunit = io, file = filename, status = "unknown", &
            & iostat = ierr)
       do i = 1, n_iter
          write(io, *)hist_all(i, 1:n_all)
       end do
       close(io)
    end if

    return 
  end subroutine parallel_output_history
  
  !---------------------------------------------------------------------
  ! parallel_output_proposal
  !---------------------------------------------------------------------
  !> @brief
  !! Output # of proposed and accepted models for each proposal type
  !! @param[inout] self     Object
  !! @param[in]    filename Output file name
  !---------------------------------------------------------------------
  subroutine parallel_output_proposal(self, filename)
    class(parallel), intent(inout) :: self
    character(*), intent(in) :: filename
    type(mcmc) :: mc
    type(trans_d_model) :: tm
    integer, allocatable :: n_accept_all(:), n_propose_all(:)
    integer, allocatable :: n_accept_sum(:), n_propose_sum(:)
    integer :: io, i, nparam, n, ierr
   
    mc = self%get_mc(1)
    tm = mc%get_tm()
    nparam = tm%get_n_x()
    n = nparam + 2
    allocate(n_accept_all(n), n_propose_all(n))
    allocate(n_accept_sum(n), n_propose_sum(n))
    n_accept_all = 0
    n_propose_all = 0

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
          write(io, *)i, n_propose_sum(i), n_accept_sum(i)
       end do
       close(io)
    end if
    
    return 
  end subroutine parallel_output_proposal
  
  !---------------------------------------------------------------------
  ! judge_swap
  !---------------------------------------------------------------------
  !> @brief
  !! Judge wether swap proposal is accepted or not
  !! @param[in]  temp1       Temperature of 1st MCMC chain
  !! @param[in]  temp2       Temperature of 2nd MCMC chain
  !! @param[in]  l1          Log-likelihood of 1st MCMC chain
  !! @param[in]  l2          Log-likelihood of 2nd MCMC chain
  !! @param[out] is_accepted Whether accepted or not
  !---------------------------------------------------------------------
  subroutine judge_swap(temp1, temp2, l1, l2, is_accepted)
    double precision, intent(in) :: temp1, temp2, l1, l2
    logical, intent(out) :: is_accepted
    double precision :: del_s
    double precision :: r
    
    del_s = (l2 - l1) * (1.d0 / temp1 - 1.d0 / temp2)
    is_accepted = .false.
    r = log(rand_u())

    if(r <= del_s) then
       is_accepted = .true.
    end if
    
    return 
  end subroutine judge_swap
  
  !---------------------------------------------------------------------
  ! pack_pair_info
  !---------------------------------------------------------------------
  !> @brief
  !! Pack MCMC chain ID info. for sending it via MPI communication
  !! @param[in]  rank1  Process ID 1
  !! @param[in]  rank2  Process ID 2
  !! @param[in]  chain1 Chain ID 1
  !! @param[in]  chain2 Chain ID 2 
  !! @return            Integer array to be passed to mpi_send
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
  ! unpack_pair_info
  !---------------------------------------------------------------------
  !> @brief
  !! Unpack MCMC chain ID info. that is received via MPI communication
  !! @param[in]  ipack  Integer array that is received via mpi_recv
  !! @param[out] rank1  Process ID 1 
  !! @param[out] rank2  Process ID 2
  !! @param[out] chain1 Chain ID 1
  !! @param[out] chain2 Chain ID 2
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
  ! pack_mc_info
  !---------------------------------------------------------------------
  !> @brief
  !! Pack MCMC chain info. for sending it via MPI communication
  !! @param[in]  temp        Temperature of MCMC chain
  !! @param[out] likelihood  Log-likelihood of MCMC chain
  !! @return                 Double precision array to be sent 
  !!                         via mpi_send
  !---------------------------------------------------------------------
  function pack_mc_info(temp, likelihood) result(rpack)
    double precision, intent(in) :: temp, likelihood
    double precision :: rpack(2)
    
    rpack(1) = temp
    rpack(2) = likelihood

    return 
  end function pack_mc_info

  !---------------------------------------------------------------------
  ! unpack_mc_info
  !---------------------------------------------------------------------
  !> @brief
  !! Unpack MCMC chain info. that is recieved via MPI communication
  !! @param[in]  rpack       Double precision array that is received by
  !!                          mpi_recv
  !! @param[out] temp        Temperature of MCMC chain
  !! @param[out] likelihood  Log-likelihood of MCMC chain
  !---------------------------------------------------------------------
  subroutine unpack_mc_info(rpack, temp, likelihood)
    double precision, intent(in) :: rpack(2)
    double precision, intent(out) :: temp, likelihood
    
    temp = rpack(1)
    likelihood = rpack(2)

    return 
  end subroutine unpack_mc_info

  !---------------------------------------------------------------------


end module mod_parallel
