
FC = mpifort
#-----------------------------------------------------------------------
# For debug
FFLAGS = -g -Wall -pedantic -std=f2008 -fbounds-check -O0 -Wuninitialized \
            -ffpe-trap=invalid,zero,overflow -fbacktrace -fno-range-check
# For fast computation
#FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check

# FFTW Library
FFTW = -I/usr/local/include -lfftw3

# LAPACK Library
LAPACK = -llapack -lblas


#-----------------------------------------------------------------------
# DONOT CHANGE BELOW
#-----------------------------------------------------------------------

RAY_FWD = ../bin/rayleigh_fwd
OBJS_RAY_FWD = rayleigh_fwd.o cls_vmodel.o cls_rayleigh.o cls_param.o \
	       cls_line_text.o

RAY_INV = ../bin/rayleigh_inv	
OBJS_RAY_INV   = rayleigh_inv.o cls_mcmc.o cls_trans_d_model.o \
	 cls_rayleigh.o cls_vmodel.o mod_random.o cls_interpreter.o \
	 mod_sort.o mod_const.o cls_observation_disper.o \
	 cls_parallel.o cls_param.o cls_line_text.o mod_output.o


RF_FWD = ../bin/recv_func_fwd
OBJS_RF_FWD = recv_func_fwd.o cls_vmodel.o cls_param.o \
	 cls_line_text.o cls_signal_process.o \
	 cls_recv_func.o

RF_INV = ../bin/recv_func_inv
OBJS_RF_INV = recv_func_inv.o cls_parallel.o mod_random.o \
         cls_trans_d_model.o cls_mcmc.o cls_recv_func.o \
	 mod_const.o cls_interpreter.o cls_observation_recv_func.o \
	 cls_param.o cls_vmodel.o cls_signal_process.o mod_sort.o \
         cls_line_text.o cls_covariance.o mod_output.o

JOINT_INV = ../bin/joint_inv
OBJS_JOINT_INV = joint_inv.o cls_parallel.o mod_random.o \
         cls_trans_d_model.o cls_mcmc.o cls_recv_func.o \
	 mod_const.o cls_interpreter.o cls_observation_recv_func.o \
	 cls_param.o cls_vmodel.o cls_signal_process.o mod_sort.o \
         cls_line_text.o cls_covariance.o cls_observation_disper.o \
	 cls_rayleigh.o mod_output.o  

TARGET = $(RAY_FWD) $(RAY_INV) $(RF_FWD) $(RF_INV) $(JOINT_INV)

.PHONY: all
all: $(TARGET)


rayleigh_inv.o: mod_random.mod cls_trans_d_model.mod cls_mcmc.mod \
	        cls_rayleigh.mod cls_interpreter.mod mod_const.mod \
		cls_observation_disper.mod cls_parallel.mod \
                cls_param.mod mod_output.mod
rayleigh_fwd.o: cls_rayleigh.mod cls_vmodel.mod cls_param.mod \
	        cls_recv_func.mod
recv_func_fwd.o: cls_param.mod cls_vmodel.mod cls_signal_process.mod \
                 cls_recv_func.mod
recv_func_inv.o: cls_parallel.mod mod_random.mod cls_trans_d_model.mod \
                 cls_mcmc.mod mod_const.mod cls_interpreter.mod \
                 cls_param.mod cls_observation_recv_func.mod \
                 cls_recv_func.mod cls_covariance.mod mod_output.mod
joint_inv.o:  cls_parallel.mod mod_random.mod cls_trans_d_model.mod \
              cls_mcmc.mod mod_const.mod cls_interpreter.mod \
              cls_param.mod cls_observation_recv_func.mod \
              cls_recv_func.mod cls_covariance.mod cls_rayleigh.mod \
              cls_observation_disper.mod mod_output.mod

cls_recv_func.o: cls_vmodel.mod cls_signal_process.mod
cls_rayleigh.o: cls_vmodel.mod
cls_trans_d_model.o: mod_random.mod
cls_mcmc.o: mod_random.mod cls_trans_d_model.mod 
cls_interpreter.o: cls_trans_d_model.mod cls_vmodel.mod mod_sort.mod \
                   mod_const.mod
cls_parallel.o: cls_mcmc.mod cls_trans_d_model.mod
cls_param.o: cls_line_text.mod
cls_observation_recv_func.o: cls_line_text.mod

$(RAY_FWD): $(OBJS_RAY_FWD)
	$(FC) $(FFLAGS) $^ -o $@

$(RAY_INV): $(OBJS_RAY_INV)
	$(FC) $(FFLAGS) $^ -o $@

$(RF_FWD): $(OBJS_RF_FWD)
	$(FC) $(FFLAGS) $^ -o $@ $(FFTW)

$(RF_INV): $(OBJS_RF_INV)
	$(FC) $(FFLAGS) $^ -o $@ $(FFTW) $(LAPACK)

$(JOINT_INV): $(OBJS_JOINT_INV)
	$(FC) $(FFLAGS) $^ -o $@ $(FFTW) $(LAPACK)

clean:
	rm -f *.mod ../bin/* *.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $< $(FFTW) $(LAPACK) -o $*.o 
%.mod: %.f90 %.o
	@:
