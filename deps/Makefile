.SUFFIXS:
FC	    = gfortran
FFLAGS	= -Wall -std=f95 -Wextra -fPIC -Wimplicit-interface
LDFLAGS = -shared

TARGET  = libdftatom.so
SOURCES = $(wildcard *.f90)
OBJECTS = $(patsubst %.f90, %.o, $(SOURCES))

%.o: %.f90
	$(FC) $(FLAGS) -o $@ -c $<

$(TARGET): $(OBJECTS)
	$(FC) $(FLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

.PHONY: clean all

all: $(TARGET)

clean:
	-rm -rf *.o *.so *.mod
