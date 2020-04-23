#-----------------------------------------------------------------------
# Complier  
#-----------------------------------------------------------------------

# Example 1. GNU fortran for debug
#FC     = mpifort
#FFLAGS = -pg -g -Wall -pedantic -std=f2008 -fbounds-check -O0 -Wuninitialized 
         -ffpe-trap=invalid,zero,overflow -fbacktrace -fno-range-check 


# Example 2. GNU fortran for fast computation
FC     = mpifort
FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check

# Example 3. Intel compiler for debug
#FC      = mpif90
#FFLAGS  = -assume byterecl -traceback -g -check bounds -fpe0 -lmpi 
#MKLROOT = /home/app/intel/compilers_and_libraries_2019.5.281/linux/mkl


#-----------------------------------------------------------------------
# Libararies
#-----------------------------------------------------------------------

# Example A. General case
FFTW   = -I/usr/local/include -lfftw3
LAPACK = -llapack -lblas

# Example B. EIC
#FFTW    = -I$(MKLROOT)/include/fftw -mkl=parallel
#LAPACK  =  # leave this empty


#-----------------------------------------------------------------------
# DONOT CHANGE BELOW
#-----------------------------------------------------------------------

BINDIR = ../bin
DISP_FWD = $(BINDIR)/disper_fwd
OBJS_DISP_FWD = disper_fwd.o cls_vmodel.o cls_disper.o cls_param.o \
	       cls_line_text.o


RF_FWD = $(BINDIR)/recv_func_fwd
OBJS_RF_FWD = recv_func_fwd.o cls_vmodel.o cls_param.o \
	 cls_line_text.o cls_signal_process.o \
	 cls_recv_func.o

JOINT_INV = $(BINDIR)/joint_inv
OBJS_JOINT_INV = joint_inv.o cls_parallel.o mod_random.o \
         cls_trans_d_model.o cls_mcmc.o cls_recv_func.o \
	 mod_const.o cls_interpreter.o cls_observation_recv_func.o \
	 cls_param.o cls_vmodel.o cls_signal_process.o mod_sort.o \
         cls_line_text.o cls_covariance.o cls_observation_disper.o \
	 cls_disper.o mod_output.o mod_forward.o cls_proposal.o \
	 cls_hyper_model.o

TARGET = $(DISP_FWD) $(RF_FWD) $(JOINT_INV)

.PHONY: all
all: $(TARGET)


disper_fwd.o: cls_disper.mod cls_vmodel.mod cls_param.mod \
	        cls_recv_func.mod
recv_func_fwd.o: cls_param.mod cls_vmodel.mod cls_signal_process.mod \
                 cls_recv_func.mod
joint_inv.o:  cls_parallel.mod mod_random.mod cls_trans_d_model.mod \
              cls_mcmc.mod mod_const.mod cls_interpreter.mod \
              cls_param.mod cls_observation_recv_func.mod \
              cls_recv_func.mod cls_covariance.mod cls_disper.mod \
              cls_observation_disper.mod mod_output.mod mod_forward.mod \
              cls_proposal.mod cls_hyper_model.mod 
cls_hyper_model.o: mod_random.mod
cls_recv_func.o: cls_vmodel.mod cls_signal_process.mod
cls_disper.o: cls_vmodel.mod
cls_trans_d_model.o: mod_random.mod
cls_mcmc.o: mod_random.mod cls_trans_d_model.mod cls_hyper_model.mod \
            cls_proposal.mod
cls_interpreter.o: cls_trans_d_model.mod cls_vmodel.mod mod_sort.mod \
                   mod_const.mod cls_hyper_model.mod
cls_parallel.o: cls_mcmc.mod 
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
	$(FC) $(FFLAGS) $^ -o $@

$(RF_FWD): $(OBJS_RF_FWD)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $^ -o $@ $(FFTW)

$(JOINT_INV): $(OBJS_JOINT_INV)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $^ -o $@ $(FFTW) $(LAPACK)

clean:
	rm -f *.mod $(BINDIR)/* *.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $< $(FFTW) $(LAPACK) -o $*.o 
%.mod: %.f90 %.o
	@:
