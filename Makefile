# Makefile for FFTbor2D

CCFLAGS           = -c -ansi -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q -I ../../h
LDFLAGS           = -L . -L ../../lib -lfftw3 -lgomp -llapack -lgslcblas -lgsl -lRNA -lmfpt -lspectral -o
BINDIR            = ~/bin # Change this to the BINDIR
CC                = g++
GCC_VERSION      := $(shell expr `$(CC) -dumpversion`)
CC_MAJ_VERSION   := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 1` \* 10000)
CC_MIN_VERSION   := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 2` \* 100)
CC_PATCH_VERSION := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 3`)
GCC_NUM_VERSION  := $(shell expr $(CC_MAJ_VERSION) \+ $(CC_MIN_VERSION) \+ $(CC_PATCH_VERSION))
GCC_GTEQ_4.6.0   := $(shell expr $(GCC_NUM_VERSION) \>= 40600)

ifeq "$(GCC_GTEQ_4.6.0)" "1"
	CCFLAGS += -Ofast -march=native
else
	CCFLAGS += -O3
endif

FFTbor2D.out: partition.o misc.o main.o
	$(CC) partition.o misc.o main.o $(LDFLAGS) FFTbor2D.out
	
main.o: main.cpp ../../h/partition.h
	$(CC) $(CCFLAGS) main.cpp

misc.o: misc.cpp ../../h/misc.h
	$(CC) $(CCFLAGS) misc.cpp

partition.o: partition.cpp ../../h/partition.h ../../h/libmfpt_header.h ../../lib/libmfpt.a ../../h/libspectral_header.h ../../lib/libspectral.a ../../h/misc.h ../../h/energy_par.h ../../h/params.h
	$(CC) $(CCFLAGS) partition.cpp
  
../../lib/libmfpt.a:
	cd ../mfpt; make
	
../../lib/libspectral.a:
	cd ../spectral; make

clean:
	rm -f *.o FFTbor2D.out

install: FFTbor2D.out
	cp FFTbor2D.out $(BINDIR)/FFTbor2D
	cp ../../misc/rna_turner1999.par ../../misc/rna_turner2004.par $(BINDIR)
	
