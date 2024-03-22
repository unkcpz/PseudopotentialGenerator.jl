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

ode1d.o: ode1d.f90 types.o utils.o
	$(FC) $(FLAGS) -o $@ -c $<

%.o: %.f90
	$(FC) $(FLAGS) -o $@ -c $<

$(TARGET): $(OBJECTS)
	$(FC) $(FLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

.PHONY: clean all

all: $(TARGET)

clean:
	-rm -rf *.o *.so *.mod
