# Makefile for FFTbor2D

CCFLAGS           = -c -ansi -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q -I $(H)/
LDFLAGS           = -L . -L $(LIB)/ -lfftw3 -lgomp -llapack -lgslcblas -lgsl -lRNA -lmfpt -lspectral -o
BINDIR            = ~/bin # Change this to the BINDIR
CC                = g++
GCC_VERSION      := $(shell expr `$(CC) -dumpversion`)
CC_MAJ_VERSION   := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 1` \* 10000)
CC_MIN_VERSION   := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 2` \* 100)
CC_PATCH_VERSION := $(shell expr `echo $(GCC_VERSION) | cut -d . -f 3`)
GCC_NUM_VERSION  := $(shell expr $(CC_MAJ_VERSION) \+ $(CC_MIN_VERSION) \+ $(CC_PATCH_VERSION))
GCC_GTEQ_4.6.0   := $(shell expr $(GCC_NUM_VERSION) \>= 40600)
GCC_GTEQ_4.9.0   := $(shell expr $(GCC_NUM_VERSION) \>= 40900)
LIB              := ../../lib
H                := ../../h

ifeq "$(GCC_GTEQ_4.6.0)" "1"
	CCFLAGS += -Ofast -march=native
else
	CCFLAGS += -O3
endif

ifeq "$(GCC_GTEQ_4.9.0)" "1"
	CCFLAGS += -fdiagnostics-color=always
endif

FFTbor2D.out: fftbor2d_functions.o rna_misc_functions.o fftbor2d_params.o fftbor2d.o
	$(CC) fftbor2d_functions.o rna_misc_functions.o fftbor2d_params.o fftbor2d.o $(LDFLAGS) FFTbor2D.out
	
fftbor2d.o: fftbor2d.cpp $(H)/fftbor2d_functions.h $(H)/fftbor2d_params.h
	$(CC) $(CCFLAGS) fftbor2d.cpp

fftbor2d_params.o: fftbor2d_params.cpp $(H)/fftbor2d_params.h $(H)/rna_misc_functions.h $(H)/vienna_data_structures.h
	$(CC) $(CCFLAGS) fftbor2d_params.cpp
  
rna_misc_functions.o: rna_misc_functions.cpp $(H)/rna_misc_functions.h
	$(CC) $(CCFLAGS) rna_misc_functions.cpp

fftbor2d_functions.o: fftbor2d_functions.cpp $(H)/fftbor2d_functions.h $(H)/libmfpt_header.h $(LIB)/libmfpt.a $(H)/libspectral_header.h $(LIB)/libspectral.a $(H)/rna_misc_functions.h $(H)/energy_par.h $(H)/vienna_functions.h
	$(CC) $(CCFLAGS) fftbor2d_functions.cpp
  
$(LIB)/libmfpt.a:
	cd ../mfpt; make
	
$(LIB)/libspectral.a:
	cd ../spectral; make

clean:
	rm -f *.o FFTbor2D.out

install: FFTbor2D.out
	cp FFTbor2D.out $(BINDIR)/FFTbor2D
	cp ../../misc/rna_turner1999.par ../../misc/rna_turner2004.par $(BINDIR)
	
