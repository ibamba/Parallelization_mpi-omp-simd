CC = gcc

CFLAGS = -g -O3 -Wall

FILES=list_res.o pathracer_mpi pathracer_mpi_statique

BIN= pathracer_mpi pathracer_mpi_statique

EXEC= pathracer_mpi

CMD_EXEC=mpirun -n 32 -hostfile hostfile --map-by node $(EXEC)

LIB= -lm

all : $(FILES)

list_res.o : list_res.c
	$(CC) -c $(CFLAGS) $<

pathracer_mpi : pathracer_mpi.c
	mpicc -o $@ $(CFLAGS) $< $(LIB)

pathracer_mpi_statique : pathracer_mpi_statique.c
	mpicc -o $@ $(CFLAGS) $< $(LIB)

exec:
	$(CMD_EXEC)

clean : 
	rm -rf $(BIN) *~ *.o
