PROG=nmaker

HEADERS=data.h domain.h driver.h fillcurve.h offload.h parameter.h partition.h partmesh.h ppkernel.h queue.h scheduler.h snapshot.h subcuboid.h subtree.h traversal.h
OBJS_CPU=domain.o driver.o fillcurve.o parameter.o partition.o partmesh.o ppkernel.o queue.o scheduler.o snapshot.o subcuboid.o subtree.o traversal.o
OBJS_MIC=domain.om driver.om fillcurve.om parameter.om partition.om partmesh.om ppkernel.om queue.om scheduler.om snapshot.om subcuboid.om subtree.om traversal.om

OPENMP=-openmp
OMP_NUM_THREADS ?= 4
KMP_AFFINITY ?= compact
NUM_PROCS ?= 4

CC=mpiicc
MICCC=icc
CFLAGS=-xHost -fimf-domain-exclusion=15 -O3 $(OPENMP) -g -traceback -Bdynamic -fno-omit-frame-pointer 
CFLAGSMIC=-mmic -fimf-domain-exclusion=15 -O3 $(OPENMP) -g -traceback -Bdynamic -fno-omit-frame-pointer
INC=-I/opt/intel/mkl/include/fftw

LD=mpiicc
LDFLAGS=-L/opt/intel/mkl/lib/intel64 -lfftw2x_cdft_SINGLE_lp64 -lmkl_cdft_core -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -traceback -Bdynamic -fno-omit-frame-pointer

run: $(PROG)
	OMP_NUM_THREADS=$(OMP_NUM_THREADS) \
	KMP_AFFINITY=$(KMP_AFFINITY) \
	mpirun -np $(NUM_PROCS) ./$(PROG)

n: $(PROG)

$(PROG): $(OBJS_CPU)
	$(LD) -o $(PROG) $(OBJS_CPU) $(LDFLAGS) $(OPENMP)

pp: ppkernel.o
	$(CC) -o pp9 pp9.c ppkernel.o $(CFLAGS)

ppmic: ppkernel.om
	$(MICCC) -o pp9.mic pp9.c ppkernel.om $(CFLAGSMIC)

.PHONY:clean
clean:
	rm -rf *.o *.om $(PROG) pp9 pp9.mic

.SUFFIXES: .c .o .om

.c.o:
	$(CC) $(CFLAGS) $(INC) -c $<

.c.om:
	$(MICCC) $(CFLAGSMIC) $(INC) -o $*.om -c $<

domain.o: domain.c domain.h data.h fillcurve.h offload.h snapshot.h parameter.h
driver.o: driver.c driver.h parameter.h partition.h domain.h data.h fillcurve.h offload.h snapshot.h partmesh.h
dtime.o: dtime.c dtime.h
fillcurve.o: fillcurve.c fillcurve.h
offload.o: offload.c
parameter.o: parameter.c parameter.h
partition.o: partition.c partition.h domain.h data.h fillcurve.h offload.h snapshot.h parameter.h
partmesh.o: partmesh.c partmesh.h data.h fillcurve.h offload.h parameter.h domain.h snapshot.h
pp9.o: pp9.c ppkernel.h offload.h
ppkernel.o: ppkernel.c ppkernel.h offload.h
queue.o: queue.c queue.h data.h fillcurve.h offload.h domain.h snapshot.h subtree.h parameter.h subcuboid.h ppkernel.h
scheduler.o: scheduler.c scheduler.h data.h fillcurve.h offload.h domain.h snapshot.h parameter.h
snapshot.o: snapshot.c snapshot.h data.h fillcurve.h offload.h
subcuboid.o: subcuboid.c subcuboid.h data.h fillcurve.h offload.h parameter.h domain.h snapshot.h
subtree.o: subtree.c subtree.h domain.h data.h fillcurve.h offload.h snapshot.h parameter.h dtime.h
traversal.o: traversal.c traversal.h data.h fillcurve.h offload.h domain.h snapshot.h subtree.h parameter.h subcuboid.h ppkernel.h queue.h

pp9.om: pp9.c ppkernel.h offload.h
ppkernel.om: ppkernel.c ppkernel.h offload.h
