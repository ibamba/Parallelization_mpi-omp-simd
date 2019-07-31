DIR= par_mpi/ par_mpi-omp/ par_mpi-omp-simd/ par_omp/ sequential/ simd/

all : make

make : 
	for i in $(DIR); do \
	cd $$i; make; cd ..; \
	done

clean : 
	for i in $(DIR); do \
	cd $$i; make clean; cd ..; \
	done

