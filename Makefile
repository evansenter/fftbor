# Makefile for FFTbor2D

CCFLAGS           = -c -std=c++11 -pedantic -fopenmp -funroll-loops -Wall -Wextra -Wa,-q -I $(HEADER) -I $(SHARED_HEADER)
LDFLAGS           = -L . -L $(LIB)/ -L /usr/local/include -lfftw3 -lgomp -llapack -lgslcblas -lgsl -lRNA -lmfpt -lspectral -o
BINDIR           = ~/bin
LIBDIR           = ~/lib
CC               = g++
LIB              = ../../lib
SHARED_HEADER    = ../../h
HEADER           = h
CODE             = cpp
GCC_VERSION      = $(shell expr `$(CC) -dumpversion`)
CC_MAJ_VERSION   = $(shell expr `echo $(GCC_VERSION) | cut -d . -f 1` \* 10000)
CC_MIN_VERSION   = $(shell expr `echo $(GCC_VERSION) | cut -d . -f 2` \* 100)
CC_PATCH_VERSION = $(shell expr `echo $(GCC_VERSION) | cut -d . -f 3`)
GCC_NUM_VERSION  = $(shell expr $(CC_MAJ_VERSION) \+ $(CC_MIN_VERSION) \+ $(CC_PATCH_VERSION))
GCC_GTEQ_4.6.0   = $(shell expr $(GCC_NUM_VERSION) \>= 40600)
GCC_GTEQ_4.9.0   = $(shell expr $(GCC_NUM_VERSION) \>= 40900)

ifeq "$(GCC_GTEQ_4.6.0)" "1"
	CCFLAGS += -Ofast -march=native
else
	CCFLAGS += -O3
endif

ifeq "$(GCC_GTEQ_4.9.0)" "1"
	CCFLAGS += -fdiagnostics-color=always
endif

FFTbor2D.out: $(CODE)/fftbor2d_functions.o $(CODE)/fftbor2d_initializers.o $(CODE)/fftbor2d_params.o $(CODE)/fftbor2d.o $(LIB)/libmfpt.a $(LIB)/libspectral.a
	$(CC) $(CODE)/fftbor2d_functions.o $(CODE)/fftbor2d_initializers.o $(CODE)/fftbor2d_params.o $(CODE)/fftbor2d.o $(LDFLAGS) FFTbor2D.out
	
$(CODE)/fftbor2d.o: $(CODE)/fftbor2d.cpp $(HEADER)/functions.h $(HEADER)/params.h
	$(CC) $(CCFLAGS) $(CODE)/fftbor2d.cpp -o $(CODE)/fftbor2d.o

$(CODE)/fftbor2d_params.o: $(CODE)/fftbor2d_params.cpp $(HEADER)/params.h
	$(CC) $(CCFLAGS) $(CODE)/fftbor2d_params.cpp -o $(CODE)/fftbor2d_params.o

$(CODE)/fftbor2d_initializers.o: $(CODE)/fftbor2d_initializers.cpp $(HEADER)/initializers.h $(HEADER)/functions.h
	$(CC) $(CCFLAGS) $(CODE)/fftbor2d_initializers.cpp -o $(CODE)/fftbor2d_initializers.o

$(CODE)/fftbor2d_functions.o: $(CODE)/fftbor2d_functions.cpp $(HEADER)/functions.h $(HEADER)/initializers.h $(HEADER)/params.h $(SHARED_HEADER)/shared/libmfpt_header.h $(SHARED_HEADER)/shared/libspectral_header.h
	$(CC) $(CCFLAGS) $(CODE)/fftbor2d_functions.cpp -o $(CODE)/fftbor2d_functions.o
  
$(LIB)/libmfpt.a:
	cd ../mfpt; make
	
$(LIB)/libspectral.a:
	cd ../spectral; make

clean:
	rm -f $(CODE)/*.o FFTbor2D.out

install: FFTbor2D.out
	cp FFTbor2D.out $(BINDIR)/FFTbor2D
	cp ../../misc/rna_turner1999.par ../../misc/rna_turner2004.par $(BINDIR)
	
