# Makefile for FFTbor2D

CCFLAGS           = -c -ansi -pedantic -fopenmp -funroll-loops -Wall -Wextra
LDFLAGS           = -L. -lfftw3 -lgomp -llapack -lRNA -o
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

FFTbor2D : partition.o energy_grid_mfpt.o misc.o main.o
	$(CC) partition.o energy_grid_mfpt.o misc.o main.o $(LDFLAGS) FFTbor2D
	
main.o : main.cpp partition.h
	$(CC) $(CCFLAGS) main.cpp

misc.o : misc.cpp misc.h
	$(CC) $(CCFLAGS) misc.cpp
	
energy_grid_mfpt.o : energy_grid_mfpt.c energy_grid_mfpt.h
	$(CC) $(CCFLAGS) energy_grid_mfpt.c

partition.o: partition.cpp partition.h energy_grid_mfpt.h misc.h energy_par.h params.h
	$(CC) $(CCFLAGS) partition.cpp

clean:
	rm -f *.o FFTbor2D

install:
	cp FFTbor2D $(BINDIR)
	cp rna_turner1999.par $(BINDIR)
	
