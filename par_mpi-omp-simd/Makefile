CC = gcc

CFLAGS = -fopenmp -g -O3 -Wall -Iinc

BIN = list_res.o pathracer_mpi-omp-simd

EXEC= pathracer_mpi-omp-simd

CMD_EXEC=mpirun -n 4 -hostfile hostfile --map-by node $(EXEC)

LIB= -lm -mavx2

all : $(BIN)

list_res.o : list_res.c
	$(CC) -c $(CFLAGS) $<

pathracer_mpi-omp-simd : pathracer_mpi-omp-simd.c
	mpicc -o $@ $(CFLAGS) $< $(LIB)

exec:
	(export OMP_NUM_THREADS=2; $(CMD_EXEC))

clean : 
	rm -rf $(BIN) *~ *.o
