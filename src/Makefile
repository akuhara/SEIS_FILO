FC = gfortran
FFLAGS = -pg -Wall -pedantic -std=f2008 -fbounds-check -O0 -Wuninitialized \
            -ffpe-trap=invalid,zero,overflow -fbacktrace -fno-range-check



TARGET_TEST_VMODEL = ../bin/test_vmodel
OBJS_TEST_VMODEL   = test_vmodel.o mod_vmodel.o

TARGET_TEST_RAYLEIGH = ../bin/test_rayleigh
OBJS_TEST_RAYLEIGH   = test_rayleigh.o mod_vmodel.o mod_rayleigh.o


#TARGET = $(TARGET_TEST_VMODEL) $(TARGET_TEST_RAYLEIGH)
TARGET = $(TARGET_TEST_RAYLEIGH) $(TARGET_TEST_VMODEL) 



.PHONY: all
all: $(TARGET)

test_vmodel.o: mod_vmodel.mod
test_rayleigh.o: mod_vmodel.mod mod_rayleigh.mod
mod_rayleigh.o: mod_vmodel.mod


$(TARGET_TEST_VMODEL): $(OBJS_TEST_VMODEL)
	$(FC) $(FFLAGS) $^ -o $@

$(TARGET_TEST_RAYLEIGH): $(OBJS_TEST_RAYLEIGH)
	$(FC) $(FFLAGS) $^ -o $@

clean:
	rm -f *.mod ../bin/* *.o

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $*.o 
%.mod: %.f90 %.o
	@:
