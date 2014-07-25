#mpicc *.c  -pthread -lrfftw_mpi -lfftw_mpi -lrfftw -lfftw -lm
mpiicc *.c -openmp -pthread -I/opt/intel/mkl/include/fftw -L/opt/intel/mkl/lib/intel64 -lfftw2x_cdft_SINGLE_lp64 -lmkl_cdft_core  -lmkl_blacs_intelmpi_lp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential
export OMP_NUM_THREADS=2
export KMP_AFFINITY=compact
mpirun -np 4 ./a.out
