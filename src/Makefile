FC = mpifort
FFLAGS = -g -Wall -pedantic -std=f2008 -fbounds-check -O0 -Wuninitialized \
            -ffpe-trap=invalid,zero,overflow -fbacktrace -fno-range-check

OBJS   = rayleigh_inv.o mod_mcmc.o mod_trans_d_model.o \
	 mod_rayleigh.o mod_vmodel.o mod_random.o mod_interpreter.o \
	 mod_sort.o mod_const.o
TARGET = ../bin/rayleigh_inv


.PHONY: all
all: $(TARGET)

rayleigh_inv.o: mod_random.mod mod_trans_d_model.mod mod_mcmc.mod \
	        mod_rayleigh.mod mod_interpreter.mod mod_const.mod
mod_rayleigh.o: mod_vmodel.mod
mod_trans_d_model.o: mod_random.mod
mod_mcmc.o: mod_random.mod mod_trans_d_model.mod 
mod_interpreter.o: mod_trans_d_model.mod mod_vmodel.mod mod_sort.mod \
                   mod_const.mod

$(TARGET): $(OBJS)
	$(FC) $(FFLAGS) $^ -o $@

clean:
	rm -f *.mod ../bin/* *.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $*.o 
%.mod: %.f90 %.o
	@:
