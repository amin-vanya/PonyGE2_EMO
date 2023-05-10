
#ifndef _MOMBI3_
#define _MOMBI3_

#include "refpoint.h"
#include "utility.h"
#include "random.h"
#include "common.h"
#include "debug.h"
#include "param.h"
#include "stop.h"
#include "list.h"
#include "plot.h"

typedef struct {
  EMO_List lst1, lst2; // temporary lists
  EMO_Refpoint ref;    // Reference points
  double *norm;        // normalized objective functions
  double *rank;        // rank of individuals
  double *rank2;       // rank of individuals
  double *l2;          // l2 norm
  double **sort;       // temporary array for sorting population
  int ssize;           // temporal population size (for memory free)
  double *tmp;         // temporary array, borrar RHG
  int wsize;           // number of weight vectors
  double *W;           // weight vectors
  char **H;            // hyper-heuristic
  char **T;            // target directions
  double *min;         // ideal point
  double *max;         // nadir point
  double *ideal;       // ideal point
  double *new_min;
  double *new_max;
  double *hist;
  int *update;
  EMO_Utility utl;
  int *filter0;        // selected individuals
  int *filter;         // selected individuals
  int *filter2;        // selected individuals
  int max_hist;        /* parameters of the algorithm */
  double epsilon;
  double alpha;
  char *wfile;
  int dm;
  int (*fcomp)(const void *, const void *);
  void (*fnorm)(double *, double *, int *, int, double *, double *, int);
} EMO_MOMBI3;

void EMO_MOMBI3_alloc(EMO_MOMBI3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_MOMBI3_free(EMO_MOMBI3 *alg);
void EMO_MOMBI3_run(EMO_MOMBI3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_MOMBI3_prun(EMO_MOMBI3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

#endif

