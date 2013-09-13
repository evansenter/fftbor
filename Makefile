# Makefile for FFTbor2D

CFLAGS  = -c -O3 -fopenmp
LDFLAGS = -lfftw3 -L. -lgomp 
BINDIR  = /usr/local/bin # Change this to the BINDIR
CC      = g++

FFTbor2D : partition.o misc.o main.o
	$(CC) partition.o misc.o main.o $(LDFLAGS) -lRNA_1.8.5 -o FFTbor2D
	
main.o : main.cpp partition.h
	$(CC) $(CFLAGS) main.cpp

misc.o : misc.cpp misc.h
	$(CC) $(CFLAGS) misc.cpp

partition.o: partition.cpp partition.h params.h energy_par.h energy_const.h
	$(CC) $(CFLAGS) partition.cpp

clean:
	rm -f *.o FFTbor2D

install:
	cp FFTbor2D $(BINDIR)
	cp rna_turner_1.8.5.par $(BINDIR)
