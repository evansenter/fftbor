// Global variables that were previously provided by ViennaRNA library
// and global variables used throughout the FFTbor algorithm

#include "energy_const.h"
#include "params.h"
#include "memory_types.h"

// FFTbor algorithm global variables
int PF = 0;
int N = 0;
int PRECISION = 4;
int WINDOW_SIZE = 0;
int MIN_WINDOW_SIZE = 0;
char* ENERGY = nullptr;
fftbor::ParamPtr P = nullptr;

// Temperature in Celsius (default 37)
double temperature = 37.0;

// Maximum Ninio correction for asymmetric interior loops
int MAX_NINIO = 300;

// Temperature at which parameters were measured
double Tmeasure = 37.0;

// Parameter for logarithmic loop energy extrapolation
double lxc37 = 107.856;

// These arrays are declared in energy_par.h but we don't actually use them
// since our parameter_parser.cpp loads everything from .par files
// They are only here to satisfy the linker for legacy compatibility

int stack37[NBPAIRS+1][NBPAIRS+1] = {0};
int enthalpies[NBPAIRS+1][NBPAIRS+1] = {0};
int entropies[NBPAIRS+1][NBPAIRS+1] = {0};

int hairpin37[31] = {0};
int bulge37[31] = {0};
int internal_loop37[31] = {0};
int internal2_energy = 0;

int old_mismatch_37[NBPAIRS+1][5][5] = {{{0}}};
int mismatchI37[NBPAIRS+1][5][5] = {{{0}}};
int mismatchH37[NBPAIRS+1][5][5] = {{{0}}};
int mismatchM37[NBPAIRS+1][5][5] = {{{0}}};
int mism_H[NBPAIRS+1][5][5] = {{{0}}};

int dangle5_37[NBPAIRS+1][5] = {{0}};
int dangle3_37[NBPAIRS+1][5] = {{0}};
int dangle3_H[NBPAIRS+1][5] = {{0}};
int dangle5_H[NBPAIRS+1][5] = {{0}};

int int11_37[NBPAIRS+1][NBPAIRS+1][5][5] = {{{{0}}}};
int int11_H[NBPAIRS+1][NBPAIRS+1][5][5] = {{{{0}}}};

int int21_37[NBPAIRS+1][NBPAIRS+1][5][5][5] = {{{{{0}}}}};
int int21_H[NBPAIRS+1][NBPAIRS+1][5][5][5] = {{{{{0}}}}};

int int22_37[NBPAIRS+1][NBPAIRS+1][5][5][5][5] = {{{{{{0}}}}}};
int int22_H[NBPAIRS+1][NBPAIRS+1][5][5][5][5] = {{{{{{0}}}}}};

int ML_BASE37 = 0;
int ML_closing37 = 340;
int ML_intern37 = 40;

int F_ninio37[5] = {0, 0, 0, 60, 60};

int TerminalAU = 50;
int DuplexInit = 410;

char Tetraloops[1401] = "";
int TETRA_ENERGY37[200] = {0};
int TETRA_ENTH37 = 0;

char Triloops[241] = "";
int Triloop_E37[40] = {0};
