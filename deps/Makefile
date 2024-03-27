.SUFFIXS:
FC	    = gfortran
FFLAGS	= -Wall -std=f95 -Wextra -fPIC -Wimplicit-interface
LDFLAGS = -shared

TARGET  = libdftatom.so
SOURCES = $(wildcard *.f90)
OBJECTS = $(patsubst %.f90, %.o, $(SOURCES))

# Explicit rules for object files to ensure correct compilation order
types.o: types.f90
	$(FC) $(FLAGS) -o $@ -c $<

utils.o: utils.f90 types.o
	$(FC) $(FLAGS) -o $@ -c $<

mesh.o: mesh.f90 types.o utils.o
	$(FC) $(FLAGS) -o $@ -c $<

ode1d.o: ode1d.f90 types.o utils.o
	$(FC) $(FLAGS) -o $@ -c $<

rschroed.o: rschroed.f90 types.o utils.o mesh.o ode1d.o
	$(FC) $(FLAGS) -o $@ -c $<

rdirac.o: rdirac.f90 types.o utils.o mesh.o ode1d.o
	$(FC) $(FLAGS) -o $@ -c $<

reigen.o: reigen.f90 types.o utils.o ode1d.o rschroed.o rdirac.o
	$(FC) $(FLAGS) -o $@ -c $<

wrapper.o: wrapper.f90 mesh.o utils.o ode1d.o types.o rschroed.o rdirac.o
	$(FC) $(FLAGS) -o $@ -c $<

%.o: %.f90
	$(FC) $(FLAGS) -o $@ -c $<

$(TARGET): $(OBJECTS)
	$(FC) $(FLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

.PHONY: clean all

all: $(TARGET)

clean:
	-rm -rf *.o *.so *.mod
