# Makefile for FFTbor2D

CFLAGS  = -c -O3 -fopenmp
LDFLAGS = -lfftw3 -L. -lgomp -lRNA_2.1.2
BINDIR  = /usr/local/bin # Change this to the BINDIR
CC      = g++

FFTbor2D_212 : partition.o misc.o main.o
	$(CC) partition.o misc.o main.o $(LDFLAGS) -o FFTbor2D_212
	
main.o : main.cpp partition.h
	$(CC) $(CFLAGS) main.cpp

misc.o : misc.cpp misc.h
	$(CC) $(CFLAGS) misc.cpp

partition.o: partition.cpp partition.h params.h energy_par.h energy_const.h
	$(CC) $(CFLAGS) partition.cpp

clean:
	rm -f *.o FFTbor2D_212

install:
	cp FFTbor2D_212 $(BINDIR)
	cp rna_turner_2.1.2.par $(BINDIR)
