#Makefile for computing number of delta neighbours


CFLAGS = -c -O3 # -O0
#LDFLAGS = -lm -L/usr/src/ViennaRNA-1.6.1/lib #Change this to where libRNA.a is
#LDFLAGS = -lm -L/usr/local/lib #Change this to where libRNA.a is
LDFLAGS = -lm -L.
PG = #-pg
BINDIR = /usr/local/bin                      #Change this to BINDIR
CC = gcc
#VAR = -DCOMPUTEMFE #Comment out this line if you don't want to compute MFE^delta

RNAbor : delta.o misc.o main.o
	 $(CC) $(VAR) $(PG) $(LDFLAGS) -g delta.o misc.o main.o -lRNA -o RNAbor
main.o : main.c delta.h
	   $(CC) $(VAR) -Wall -W $(PG) $(CFLAGS) $(H) main.c

misc.o : misc.c misc.h
	 $(CC) $(VAR) -Wall -W $(PG) -g $(CFLAGS) $(H) misc.c

delta.o: delta.c delta.h
	$(CC) $(VAR) -Wall -W $(PG) -g $(CFLAGS) $(H) delta.c

clean:
	rm -f *.o RNAbor

install:
	cp RNAbor  $(BINDIR)
	cp energy.par  $(BINDIR)
