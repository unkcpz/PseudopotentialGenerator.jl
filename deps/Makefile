.SUFFIXS:
FC	=	gfortran
FFLAGS	=	-Wall -Wextra -fPIC -shared
LDFLAGS      = -shared

TARGET = liboncv.so
SOURCES = aeo.f90 lschfb.f90 tfapot.f90 hartree.f90
OBJECTS=${SOURCES:.f90=.o}

%.o: %.f90
	$(FC) $(FLAGS) -o $@ -c $<

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(FC) $(FLAGS) $(OBJECTS) -o $@ $(LDFLAGS)

.PHONY: clean
clean:
	-rm -rf *.o *.so
