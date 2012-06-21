# Makefile for FFTbor

CFLAGS  = -c -O3
LDFLAGS = -lfftw3 -lm -L.
BINDIR  = /usr/local/bin # Change this to the BINDIR
CC      = g++

RNAbor : delta.o misc.o main.o
	 $(CC) -g delta.o misc.o main.o $(LDFLAGS) -lRNA -o RNAbor
	
main.o : main.c delta.h
	   $(CC) -Wall -W $(CFLAGS) main.c

misc.o : misc.c misc.h
	 $(CC) -Wall -W -g $(CFLAGS) misc.c

delta.o: delta.c delta.h
	$(CC) -Wall -W -g $(CFLAGS) delta.c

clean:
	rm -f *.o RNAbor

install:
	cp RNAbor $(BINDIR)
	cp energy.par $(BINDIR)
