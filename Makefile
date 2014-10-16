CC=g++
CFLAGS=-ansi -Wall -o

all: serial omp mpi

serial:
	$(CC) $(CFLAGS) lcs-serial lcs-serial.cpp

omp:
	$(CC) $(CFLAGS) lcs-omp lcs-omp.cpp

mpi:
	$(CC) $(CFLAGS) lcs-mpi lcs-mpi.cpp

clean:
	rm -rf *.o lcs-serial lcs-omp lcs-mpi
