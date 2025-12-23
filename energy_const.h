/*

   energy constants, formerly defined in
             energy_par.h

   customized for use with RNAedit by
   S.Kopp, IMB-Jena, Germany, Mar 1996
   
*/

#ifndef _ENERGY_CONST_H
#define _ENERGY_CONST_H 1

// Energy calculation constants
constexpr double GASCONST = 1.98717;  // in [cal/K]
constexpr double K0 = 273.15;
constexpr int INF = 1000000;
constexpr int FORBIDDEN = 9999;
constexpr int BONUS = 10000;
constexpr int NBPAIRS = 7;
constexpr int TURN = 3;
constexpr int MAXLOOP = 30;

#endif
