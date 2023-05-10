
#ifndef _HYPE_H
#define _HYPE_H

#include "dominance.h"
#include "common.h"
#include "random.h"
#include "plot.h"
#include "param.h"

typedef struct {
  EMO_NDSort nd;
  int *filter;        /* selected individuals */
  double *min;            /* ideal point */
  double *max;            /* nadir point */
  int *parent1, *parent2; /* candidates for reproduction */
  int *seedp;             /* enumeration of parent population */
  EMO_List lst1, lst2;    /* temporary lists */
 
  int samples;
  double bound;

  double *tmp;            //arreglo temporal para copiar obj del ultimo frente
  double *rho;            //double rho[ pop->mu+1 ];   >>hypeIndicator
  int *hitstat;           //int hitstat[ popsize ];  >>hypeSampling
  double *sample;         //double sample[ dim ];    >>hypeSampling
  double *val;            // double val[ pop->size];  >>Parameter of hypeIndicator()
  double *boundsVec;       //double boundsVec[ nobj ];  >>hypeExact
  int *indices;           //int indices[ popsize ];   >>hypeExact
  int *pvec;              //int pvec[pnts];           >>hypeExactRecursive
  double *p;              //double p[pnts*dim];       >>hypeExactRecursive
  int *beg; 	          //int beg[MAX_LEVELS];      >>rearrangeIndicesByColumn
  int *end;               //int end[MAX_LEVELS];      >>rearrangeIndicesByColumn
  double *ref;             //double ref[rows];         >>rearrangeIndicesByColumn
} EMO_HYPE;

void EMO_HYPE_alloc(EMO_HYPE *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_HYPE_free(EMO_HYPE *alg);
void EMO_HYPE_run(EMO_HYPE *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

#endif

