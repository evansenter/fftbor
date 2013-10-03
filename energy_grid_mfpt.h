#ifndef ENERGY_GRID_MFPT_H
#define ENERGY_GRID_MFPT_H

double** convertEnergyGridToTransitionMatrix(double*, int);
double computeMFPT(int*, int*, double**, int, int);
void inverse(double*, int);

#endif
