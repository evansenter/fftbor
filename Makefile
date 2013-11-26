# Makefile for FFTbor

CFLAGS  = -c -O3
LDFLAGS = -L. -lm -lfftw3 -lgomp -lRNA -o
BINDIR  = /usr/local/bin # Change this to the BINDIR
CC      = g++

FFTbor : delta.o misc.o main.o
	 $(CC) -g delta.o misc.o main.o $(LDFLAGS) FFTbor
	
main.o : main.cpp delta.h
	   $(CC) -Wall -W $(CFLAGS) main.cpp

misc.o : misc.cpp misc.h
	 $(CC) -Wall -W -g $(CFLAGS) misc.cpp

delta.o: delta.cpp delta.h params.h energy_par.h energy_const.h
	$(CC) -Wall -W -g $(CFLAGS) delta.cpp

clean:
	rm -f *.o FFTbor

install:
	cp FFTbor $(BINDIR)
	cp rna_turner2004.par $(BINDIR)
