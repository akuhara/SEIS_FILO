#-----------------------------------------------------------------------
# Complier  
#-----------------------------------------------------------------------

# Example 1. GNU fortran for debug
FC     = mpif90
FFLAGS = -g -Wall -pedantic -fbounds-check -O0 \
         -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace \
	 -pg

# Example 2. GNU fortran for fast computation
#FC     = mpif90
#FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check

# Example 3. Intel compmiler for fast computation
#FC = mpiifort
#FFLAGS = -assume byterecl -lmpi -xAVX -O3 

# Example 4. EIC
#FC      = mpif90
#FFLAGS  = -assume byterecl -traceback -g -check bounds -fpe0 -lmpi 
#MKLROOT = /home/app/intel/compilers_and_libraries_2019.5.281/linux/mkl


#-----------------------------------------------------------------------
# Libararies
#-----------------------------------------------------------------------

# Example 1-2. Free library
FFTW   = -I/usr/local/include -lfftw3
LAPACK = -llapack -lblas

# Example 3. Intel library
#FFTW = -I/opt/intel/mkl/include/fftw -mkl=parallel
#LAPACK =

# Example 4. EIC
#FFTW    = -I$(MKLROOT)/include/fftw -mkl=parallel
#LAPACK  =  # leave this empty


#-----------------------------------------------------------------------
# DONOT CHANGE BELOW
#-----------------------------------------------------------------------
COV =
BINDIR = ../bin
DISP_FWD = $(BINDIR)/disper_fwd
OBJS_DISP_FWD = disper_fwd.o cls_vmodel.o cls_disper.o cls_param.o \
	       cls_line_text.o mod_random.o mod_mpi.o


RF_FWD = $(BINDIR)/recv_func_fwd
OBJS_RF_FWD = recv_func_fwd.o cls_vmodel.o cls_param.o \
	 cls_line_text.o cls_signal_process.o \
	 cls_recv_func.o mod_random.o cls_covariance.o mod_mpi.o

JOINT_INV = $(BINDIR)/joint_inv
OBJS_JOINT_INV = joint_inv.o cls_parallel.o mod_random.o \
         cls_trans_d_model.o cls_mcmc.o cls_recv_func.o \
	 mod_const.o cls_interpreter.o cls_observation_recv_func.o \
	 cls_param.o cls_vmodel.o cls_signal_process.o mod_sort.o \
         cls_line_text.o cls_covariance.o cls_observation_disper.o \
	 cls_disper.o mod_output.o mod_forward.o cls_proposal.o \
	 cls_hyper_model.o mod_mpi.o

TARGET = $(DISP_FWD) $(RF_FWD) $(JOINT_INV) 

.PHONY: all
all: $(TARGET)



disper_fwd.o: cls_disper.mod cls_vmodel.mod cls_param.mod \
	        cls_recv_func.mod mod_random.mod mod_mpi.mod
recv_func_fwd.o: cls_param.mod cls_vmodel.mod cls_signal_process.mod \
                 cls_recv_func.mod mod_random.mod cls_covariance.mod \
                 mod_mpi.mod
joint_inv.o:  cls_parallel.mod mod_random.mod cls_trans_d_model.mod \
              cls_mcmc.mod mod_const.mod cls_interpreter.mod \
              cls_param.mod cls_observation_recv_func.mod \
              cls_recv_func.mod cls_covariance.mod cls_disper.mod \
              cls_observation_disper.mod mod_output.mod mod_forward.mod \
              cls_proposal.mod cls_hyper_model.mod mod_mpi.mod
cls_hyper_model.o: mod_random.mod
cls_recv_func.o: cls_vmodel.mod cls_signal_process.mod mod_random.mod \
	         cls_covariance.mod
cls_disper.o: cls_vmodel.mod mod_random.mod
cls_trans_d_model.o: mod_random.mod
cls_mcmc.o: mod_random.mod cls_trans_d_model.mod cls_hyper_model.mod \
            cls_proposal.mod
cls_interpreter.o: cls_trans_d_model.mod cls_vmodel.mod mod_sort.mod \
                   mod_const.mod cls_hyper_model.mod cls_line_text.mod
cls_parallel.o: cls_mcmc.mod mod_mpi.mod
cls_param.o: cls_line_text.mod
cls_observation_recv_func.o: cls_line_text.mod
cls_observation_disper.o: cls_line_text.mod
cls_vmodel.o: cls_line_text.mod
mod_forward.o: cls_trans_d_model.mod cls_interpreter.mod \
               cls_observation_recv_func.mod cls_recv_func.mod \
               cls_observation_disper.mod cls_disper.mod \
               cls_vmodel.mod mod_const.mod cls_covariance.mod \
	       cls_hyper_model.mod
mod_output.o: cls_parallel.mod 


$(DISP_FWD): $(OBJS_DISP_FWD)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $(COV) $^ -o $@

$(RF_FWD): $(OBJS_RF_FWD)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $(COV) $^ -o $@ $(FFTW) $(LAPACK)

$(JOINT_INV): $(OBJS_JOINT_INV)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $(COV) $^ -o $@ $(FFTW) $(LAPACK)

clean:
	rm -f *.mod $(BINDIR)/* *.o

%.o: %.f90
	$(FC) $(FFLAGS) $(COV) -c $< $(FFTW) $(LAPACK) -o $*.o 
%.mod: %.f90 %.o
	@:
