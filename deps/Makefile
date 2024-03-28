.SUFFIXS:
FC	    = gfortran
FFLAGS	= -Wall -std=f95 -Wextra -fPIC -Wimplicit-interface
LDFLAGS = -shared

TARGET  = libdftatom.so
SOURCES = $(wildcard *.f90)
OBJECTS = $(patsubst %.f90, %.o, $(SOURCES))

# Explicit rules for object files to ensure correct compilation order
types.o: types.f90
constants.o: types.o
utils.o: utils.f90 types.o
states.o: types.o utils.o
energies.o: constants.o types.o utils.o
mesh.o: mesh.f90 types.o utils.o
ode1d.o: ode1d.f90 types.o utils.o
rpoisson.o: constants.o ode1d.o types.o utils.o
mixings.o: constants.o dft_data.o ode1d.o types.o utils.o
rschroed.o: rschroed.f90 types.o utils.o mesh.o ode1d.o
rdirac.o: rdirac.f90 types.o utils.o mesh.o ode1d.o
reigen.o: reigen.f90 types.o utils.o ode1d.o rschroed.o rdirac.o
dft.o: constants.o dft_data.o energies.o mixings.o ode1d.o reigen.o \
	rpoisson.o types.o utils.o
dft_data.o: types.o
dftatom.o: dft.o drivers.o energies.o mesh.o ode1d.o reigen.o rpoisson.o \
	states.o types.o utils.o
drivers.o: constants.o dft.o dft_data.o energies.o mesh.o mixings.o ode1d.o \
	states.o types.o
wrapper.o: wrapper.f90 mesh.o utils.o ode1d.o types.o rschroed.o rdirac.o drivers.o

%.o: %.f90
	$(FC) $(FLAGS) -o $@ -c $<

$(TARGET): $(OBJECTS)
	$(FC) $(FLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

.PHONY: clean all

all: $(TARGET)

clean:
	-rm -rf *.o *.so *.mod
