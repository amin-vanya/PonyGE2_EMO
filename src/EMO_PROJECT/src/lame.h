#ifndef _LS_H
#define _LS_H

#include "common.h"

void EMO_LS_range(EMO_MOP *mop);
double EMO_LAME_natmin(double r);
double EMO_LS_ed2(double r);
void EMO_LS_difficulty(EMO_MOP *mop, int v);
void EMO_LS_lame(EMO_MOP *mop, double *f, double *x);
void EMO_LS_mirror(EMO_MOP *mop, double *f, double *x);

#endif
