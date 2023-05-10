
#ifndef _SMSEMOA_
#define _SMSEMOA_

#include "dominance.h"
#include "hv_iwfg.h"
#include "hv_wfg.h"
#include "common.h"
#include "param.h"
#include "debug.h"
#include "stop.h"
#include "plot.h"

typedef struct {
  EMO_NDSort nd;
  EMO_HV hv;
  EMO_IWFG iwfg;
  double *max, *min;
  int *filter;
  double *chv;

  /* mpirun */
  double **sort;       // temporary array for sorting population
  int ssize;           // temporal population size (for memory free)
  EMO_List lst1, lst2; // temporary lists
  int iwfg_flag;       // incremental IWFG algorithm
  double sfactor;      // reference point for the hypervolume indicator
} EMO_SMSEMOA;

void EMO_SMSEMOA_alloc(EMO_SMSEMOA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_SMSEMOA_free(EMO_SMSEMOA *alg);
void EMO_SMSEMOA_run(EMO_SMSEMOA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_SMSEMOA_prun(EMO_SMSEMOA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

#endif

