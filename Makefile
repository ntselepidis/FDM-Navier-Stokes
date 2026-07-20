FC    ?= nvfortran
DEBUG ?= 0

ifneq (,$(findstring nvfortran,$(FC)))
  ifeq ($(DEBUG),1)
    FFLAGS = -i4 -r8 -O0 -g -cpp -DNVHPC
  else
    FFLAGS = -i4 -r8 -O3 -acc=gpu -stdpar=gpu -gpu=ccnative,mem:unified -cuda \
             -Minline=reshape -Minfo=accel,inline -cpp -DNVHPC
  endif
  LDFLAGS = -cudaforlibs -cudalib=curand,nvtx3
else
  ifeq ($(DEBUG),1)
    FFLAGS = -O0 -g -fdefault-real-8 -fbacktrace \
             -fcheck=all -finit-real=snan -ffpe-trap=invalid,zero,overflow -Wall -cpp
  else
    FFLAGS = -O3 -fdefault-real-8 -Wall -cpp
  endif
  LDFLAGS =
endif

# Object files
OBJ = m_bc.o m_fdm.o m_stream.o m_nvtx.o m_mg.o m_sim.o main.o

all: main.x

main.x: $(OBJ)
	$(FC) $(FFLAGS) $(OBJ) -o main.x $(LDFLAGS)

%.o: %.f90
	@echo "compiling $<"
	$(FC) -c $(FFLAGS) -o $@ $<

# Dependencies
m_bc.o: m_nvtx.o
m_mg.o: m_bc.o m_nvtx.o
m_sim.o: m_bc.o m_fdm.o m_mg.o m_stream.o m_nvtx.o
main.o: m_sim.o m_nvtx.o

.PHONY: clean
clean:
	rm -f *.o *.mod *.x *.bin
