CC = gcc

CFLAGS = -fopenmp -g -O3 -Wall -Iinc

BIN = list_res.o pathracer_mpi-omp

EXEC= pathracer_mpi-omp

CMD_EXEC=mpirun -n 8 -hostfile hostfile --map-by node $(EXEC)

LIB= -lm

all : $(BIN)

list_res.o : list_res.c
	$(CC) -c $(CFLAGS) $<

pathracer_mpi-omp : pathracer_mpi-omp.c
	mpicc -o $@ $(CFLAGS) $< $(LIB)

exec:
	(export OMP_NUM_THREADS=2; $(CMD_EXEC))

clean : 
	rm -rf $(BIN) *~ *.o
