CC=gcc
MCC=mpicc
CFLAGS=-ansi -Wall -lm -g -o
OMPFLAGS=-fopenmp

all: serial omp mpi

serial:
	$(CC) $(CFLAGS) lcs-serial lcs-serial.c

omp:
	$(CC) $(OMPFLAGS) $(CFLAGS) lcs-omp lcs-omp.c

mpi:
	$(MCC) $(OMPFLAGS) $(CFLAGS) lcs-mpi lcs-mpi.c
	cp lcs-mpi /tmp/lcs-mpi-13

clean:
	rm -rf *.o lcs-serial lcs-omp lcs-mpi
