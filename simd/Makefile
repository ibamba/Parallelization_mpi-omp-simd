CC = gcc

CFLAGS = -fopenmp -g -O3 -Wall -Iinc

BIN = pathracer_simd

EXEC= pathracer_simd

CMD_EXEC=./pathracer_simd

LIB= -lm -mavx2

all : $(BIN)

pathracer_simd : pathracer_simd.c
	mpicc -o $@ $(CFLAGS) $< $(LIB)

exec:
	$(CMD_EXEC)

clean : 
	rm -rf $(BIN) *~ *.o



