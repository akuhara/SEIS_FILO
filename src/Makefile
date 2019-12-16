FC = mpifort
#FFLAGS = -g -Wall -pedantic -std=f2008 -fbounds-check -O0 -Wuninitialized \
            -ffpe-trap=invalid,zero,overflow -fbacktrace -fno-range-check
FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check

RAY_FWD = ../bin/rayleigh_fwd
OBJS_RAY_FWD = rayleigh_fwd.o mod_vmodel.o mod_rayleigh.o mod_param.o \
	       mod_line_text.o

RAY_INV = ../bin/rayleigh_inv	
OBJS_RAY_INV   = rayleigh_inv.o mod_mcmc.o mod_trans_d_model.o \
	 mod_rayleigh.o mod_vmodel.o mod_random.o mod_interpreter.o \
	 mod_sort.o mod_const.o mod_observation.o mod_parallel.o \
	 mod_param.o mod_line_text.o



TARGET = $(RAY_INV) $(RAY_FWD)

.PHONY: all
all: $(TARGET)

rayleigh_inv.o: mod_random.mod mod_trans_d_model.mod mod_mcmc.mod \
	        mod_rayleigh.mod mod_interpreter.mod mod_const.mod \
		mod_observation.mod mod_parallel.mod mod_param.mod
rayleigh_fwd.o: mod_rayleigh.mod mod_vmodel.mod mod_param.mod
mod_rayleigh.o: mod_vmodel.mod
mod_trans_d_model.o: mod_random.mod
mod_mcmc.o: mod_random.mod mod_trans_d_model.mod 
mod_interpreter.o: mod_trans_d_model.mod mod_vmodel.mod mod_sort.mod \
                   mod_const.mod
mod_parallel.o: mod_mcmc.mod mod_trans_d_model.mod
mod_param.o: mod_line_text.mod

$(RAY_FWD): $(OBJS_RAY_FWD)
	$(FC) $(FFLAGS) $^ -o $@

$(RAY_INV): $(OBJS_RAY_INV)
	$(FC) $(FFLAGS) $^ -o $@


clean:
	rm -f *.mod ../bin/* *.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $*.o 
%.mod: %.f90 %.o
	@:
