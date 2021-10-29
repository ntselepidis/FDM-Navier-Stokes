FC=nvfortran
DEBUG?=0

ifeq ($(DEBUG),1)
	FFLAGS=-Wall -fdefault-real-8 -O0 -g -fbacktrace -fcheck=all -finit-real=snan -ffpe-trap=invalid,zero,overflow
else
	FFLAGS=-Wall -r8 -O3 -Minfo=all -stdpar=gpu -gpu=cc70
endif

# Object files
OBJ=m_bc.o m_fdm.o m_stream.o m_mg.o m_sim.o main.o

all: main.x

main.x : $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o main.x

%.o : %.f90
	@echo "compiling $<"
	$(FC) -c $(FFLAGS) -o $@ $<

# Dependencies
# m_bc.o:
# m_fdm.o:
# m_stream.o:
m_mg.o: m_bc.o
m_sim.o: m_bc.o m_fdm.o m_mg.o m_stream.o
main.o: m_sim.o

# Rules
.PHONY: clean

clean:
	rm -f *.o *.mod *.x *.bin
