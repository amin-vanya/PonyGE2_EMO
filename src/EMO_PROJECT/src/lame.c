
#include <stdlib.h>
#include <math.h>

#include "lame.h"

#define LS_A 0.051373
#define LS_B 0.0253235

void EMO_LS_range(EMO_MOP *mop) {
  int i;

  for(i = 0; i < mop->nobj - 1; i++) {
    mop->xmin[i] = 0.0;
    mop->xmax[i] = M_PI / 2.0;
  }

  for(i = mop->nobj - 1; i < mop->nvar; i++) {
    mop->xmin[i] = 0.0;
    mop->xmax[i] = 1.0;
  }
}

// multimodal test problem with equidistant local Pareto fronts
double EMO_LS_natmin(double r) {
  return LS_B + r - LS_A + 0.5 + 0.5 * cos(2.0 * M_PI * (r - LS_A) + M_PI);  
}

// many-to-one mappings
double EMO_LS_ed2(double r) {
  return sin(r / M_PI);
}

void EMO_LS_difficulty(EMO_MOP *mop, int v) {

  switch(v) {
    case 0:  mop->g = NULL;
             printf("Lame Superspheres: unimodal\n");
             break;
    case 1:  mop->g = &EMO_LS_natmin;
             printf("Lame Superspheres: multimodal\n");
             break;
    case 2:  mop->g = &EMO_LS_ed2;
             printf("Lame Superspheres: many-to-one mappings\n");
             break;
    default: printf("Error, unknown value of %d for the lame_difficulty in the configuration file\n", v);
             exit(1);
  }
}

// Lame Superspheres
void EMO_LS_lame(EMO_MOP *mop, double *f, double *x) {
  int i, j, m;
  double g, r;

  m = mop->nobj;

  f[0] = cos(x[0]);
  f[0] = pow(f[0] * f[0], 1.0 / mop->gamma);

  // fj, j \in {1,..m-2} 
  for(j = m-2; j > 0; j--) {
    f[j] = 1.0;

    for(i = j-1; i > -1; i--)
      f[j] *= sin(x[i]);

    f[j] *= cos(x[j]);
    f[j] = pow(f[j] * f[j], 1.0 / mop->gamma);
  }

  f[m-1] = 1.0;

  for(i = m-2; i > -1; i--)
    f[m-1] *= sin(x[i]);

  f[m-1] = pow(f[m-1] * f[m-1], 1.0 / mop->gamma);

  r = 0.0;

  for(i = m-1; i < mop->nvar; i++)
    r += x[i] * x[i];

  if(mop->g == NULL)  // unimodal
    g = r;
  else
    g = mop->g(r);

  for(i = 0; i < m; i++)
    f[i] *= (1.0 + g);
}

// Lame Superspheres
void EMO_LS_mirror(EMO_MOP *mop, double *f, double *x) {
  int i, j, m;
  double g, r;

  m = mop->nobj;

  f[0] = cos(x[0]);
  f[0] = pow(f[0] * f[0], 1.0 / mop->gamma);

  // fj, j \in {1,..m-2} 
  for(j = m-2; j > 0; j--) {
    f[j] = 1.0;

    for(i = j-1; i > -1; i--)
      f[j] *= sin(x[i]);

    f[j] *= cos(x[j]);
    f[j] = pow(f[j] * f[j], 1.0 / mop->gamma);
  }

  f[m-1] = 1.0;

  for(i = m-2; i > -1; i--)
    f[m-1] *= sin(x[i]);

  f[m-1] = pow(f[m-1] * f[m-1], 1.0 / mop->gamma);

  r = 0.0;

  for(i = m-1; i < mop->nvar; i++)
    r += x[i] * x[i];

  if(mop->g == NULL)  // unimodal
    g = r;
  else
    g = mop->g(r);

  for(i = 0; i < m; i++)
    f[i] = 1.0 - f[i] / (1.0 + g);
}

#undef LS_A
#undef LS_B


