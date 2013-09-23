# Makefile for FFTbor2D

CCFLAGS = -c -O3 -Wall -W
LDFLAGS = -lfftw3 -L. -lRNA_2.1.2 -L/usr/local/Cellar/lapack/3.4.2/lib -llapack
BINDIR  = /usr/local/bin # Change this to the BINDIR
CC      = g++

FFTbor2D : partition.o misc.o main.o
	$(CC) partition.o misc.o main.o $(LDFLAGS) -o FFTbor2D
	
main.o : main.cpp partition.h
	$(CC) $(CCFLAGS) main.cpp

misc.o : misc.cpp misc.h
	$(CC) $(CCFLAGS) misc.cpp

partition.o: partition.cpp partition.h params.h energy_par.h energy_const.h
	$(CC) $(CCFLAGS) partition.cpp

clean:
	rm -f *.o FFTbor2D

install:
	cp FFTbor2D $(BINDIR)
	cp rna_turner_2.1.2.par $(BINDIR)
	
