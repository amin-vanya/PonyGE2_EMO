
/* Value path as density estimator in SMS-EMOA */

#ifndef _MOVAP_
#define _MOVAP_

#include "dominance.h"
#include "vpath.h"
#include "common.h"
#include "param.h"
#include "debug.h"
#include "stop.h"
#include "plot.h"

typedef struct {
  EMO_NDSort nd;
  EMO_VPath vpath;
  EMO_List lst;
  int xrs;          /* x-axis resolution of vpath */
  int *filter;      /* auxiliar for enabling/disabling individuals */
  double *tmp;      /* temporary array for storing an individual */
  double *min;      /* Reference points */
  double *max0;
  double *ideal;
  double *nadir;
  double **sort;    /* temporary array for sorting population */
  double *norm;     /* normalized objective functions */
  double *vnorm;    /* norm value for each individual */
  int ssize;        /* temporal population size (for memory free) */
} EMO_MOVAP;

void EMO_MOVAP_alloc(EMO_MOVAP *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_MOVAP_free(EMO_MOVAP *alg);
void EMO_MOVAP_run(EMO_MOVAP *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

#endif

