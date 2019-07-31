CC = gcc

CFLAGS = -fopenmp -g -O3 -Wall -Iinc

BIN= pathtracer_omp

CMD_EXEC=./pathtracer_omp

LIB= -lm

all : $(BIN)

pathtracer_omp : pathtracer_omp.c
	$(CC) -o $@ $(CFLAGS) $< $(LIB)

exec:
	(export OMP_NUM_THREADS=4; $(CMD_EXEC))

clean : 
	rm -rf $(BIN) *~ *.o
