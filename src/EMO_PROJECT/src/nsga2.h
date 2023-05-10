/* Author: Elmer Jesús Morales Orozco
 *         Universidad de la Sierra Sur,
 *         Miahuatlán, Oaxaca.
 *   
 * August 2015
 */


#ifndef _NSGA2_H
#define _NSGA2_H

#include "dominance.h"
#include "common.h"
#include "random.h"
#include "plot.h"
#include "param.h"

typedef struct {
  EMO_NDSort nd;
  EMO_List lst1, lst2;    /* temporary lists */
  int *filter;            /* selected individuals */
  double *cd;             /* crowding distance */
  double *min;            /* ideal point */
  double *max;            /* nadir point */
  int *parent1, *parent2; /* candidates for reproduction */
  int *seedp;             /* enumeration of parent population */ 
  double **sort;          /* temporary array for sorting */
  int ssize;
} EMO_NSGA2;

void EMO_NSGA2_alloc(EMO_NSGA2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_NSGA2_free(EMO_NSGA2 *alg);
void EMO_NSGA2_run(EMO_NSGA2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

#endif

