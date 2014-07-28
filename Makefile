PROG=nmaker

HEADERS=data.h domain.h driver.h fillcurve.h parameter.h partition.h partmesh.h ppkernel.h queue.h scheduler.h snapshot.h subcuboid.h subtree.h traversal.h
OBJS_CPU=domain.o driver.o fillcurve.o parameter.o partition.o partmesh.o ppkernel.o queue.o scheduler.o snapshot.o subcuboid.o subtree.o traversal.o
OBJS_MIC=domain.om driver.om fillcurve.om parameter.om partition.om partmesh.om ppkernel.om queue.om scheduler.om snapshot.om subcuboid.om subtree.om traversal.om

OPENMP=-openmp
OMP_NUM_THREADS := 4
KMP_AFFINITY := compact
NUM_PROCS := 4

CC=mpiicc
CFLAGS=-xHost -O3 $(OPENMP)
INC=-I/opt/intel/mkl/include/fftw

LD=mpiicc
LDFLAGS=-L/opt/intel/mkl/lib/intel64 -lfftw2x_cdft_SINGLE_lp64 -lmkl_cdft_core -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential

run: $(PROG)
	export OMP_NUM_THREADS=$(OMP_NUM_THREADS)
	export KMP_AFFINITY=$(KMP_AFFINITY)
	mpirun -np $(NUM_PROCS) ./$(PROG)

n: $(PROG)

$(PROG): $(OBJS_CPU)
	$(LD) -o $(PROG) $(OBJS_CPU) $(LDFLAGS) $(OPENMP)

.PHONY:clean
clean:
	rm -rf *.o $(PROG)

.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $<
