PROG=nmaker

HEADERS=data.h domain.h driver.h fillcurve.h parameter.h partition.h partmesh.h ppkernel.h queue.h scheduler.h snapshot.h subcuboid.h subtree.h traversal.h
OBJS_CPU=domain.o driver.o fillcurve.o parameter.o partition.o partmesh.o ppkernel.o queue.o scheduler.o snapshot.o subcuboid.o subtree.o traversal.o
OBJS_MIC=domain.om driver.om fillcurve.om parameter.om partition.om partmesh.om ppkernel.om queue.om scheduler.om snapshot.om subcuboid.om subtree.om traversal.om

OPENMP=-openmp
OMP_NUM_THREADS ?= 4
KMP_AFFINITY ?= compact
NUM_PROCS ?= 4

CC=mpiicc
CFLAGS=-xHost -O3 $(OPENMP) -g -traceback
INC=-I/opt/intel/mkl/include/fftw

LD=mpiicc
LDFLAGS=-L/opt/intel/mkl/lib/intel64 -lfftw2x_cdft_SINGLE_lp64 -lmkl_cdft_core -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential

run: $(PROG)
	OMP_NUM_THREADS=$(OMP_NUM_THREADS) \
	KMP_AFFINITY=$(KMP_AFFINITY) \
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

domain.o: domain.c domain.h data.h
driver.o: driver.c driver.h
dtime.o: dtime.c dtime.h
fillcurve.o: fillcurve.c fillcurve.h
parameter.o: parameter.c parameter.h
partition.o: partition.c partition.h domain.h data.h
partmesh.o: partmesh.c partmesh.h data.h
pp9.o: pp9.c ppkernel.h
ppkernel.o: ppkernel.c ppkernel.h
queue.o: queue.c queue.h data.h
scheduler.o: scheduler.c scheduler.h data.h
snapshot.o: snapshot.c snapshot.h data.h
subcuboid.o: subcuboid.c subcuboid.h data.h
subtree.o: subtree.c subtree.h domain.h data.h
traversal.o: traversal.c traversal.h data.h
