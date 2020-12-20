FC=gfortran
FFLAGS=-Wall -fdefault-real-8 -O3
#FFLAGS=-g -Wall -fdefault-real-8 -fbacktrace -fcheck=all -finit-real=snan -ffpe-trap=invalid,zero,overflow -O0
EXE=main.x

# Object files
OBJ=m_bc.o m_fdm.o m_stream.o m_mg.o m_sim.o main.o

exe : $(OBJ)
	@$(FC) $(FFLAGS) $(OBJ) -o $(EXE)

%.o : %.f90
	@echo "compiling $<"
	@$(FC) -c $(FFLAGS) -o $@ $<

# Rules
.PHONY: clean

clean:
	rm -f *.o *.mod *.x *.bin
