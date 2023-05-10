
#ifndef _SPEA2_H
#define _SPEA2_H

#include "common.h"
#include "param.h"
#include "list.h"

typedef struct {
  int k;                  /**/
  int niche;
  int wsize;              /* number of weight vectors = population size */
  double *W;              /* weight vectors */
  int *parent1, *parent2; /* candidates for reproduction */
  double *fit;            /* fitness */ 
  double *S;
  double *F;
  EMO_List *lnn;          /* array of lists for storing neighbors */
  EMO_List copy;
  double **sort;          /* temporary array for sorting */
  double **dist;
  int *filter;            /* selected individuals */
  int *seedp;             /* enumeration of parent population */
  EMO_List lst1, lst2;    /* temporary lists */ 
  int ssize;
} EMO_SPEA2;

void EMO_SPEA2_alloc(EMO_SPEA2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_SPEA2_free(EMO_SPEA2 *alg);
void EMO_SPEA2_run(EMO_SPEA2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

#endif

