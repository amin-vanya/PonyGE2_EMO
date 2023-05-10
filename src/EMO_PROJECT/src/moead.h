/* Authors: Elmer Jesús Morales Orozco
 *          Universidad de la Sierra Sur,
 *          Miahuatlán, Oaxaca.
 * 
 *          Mariano Orozco Garcia  (constraint handling)
 *          UPIITA-IPN
 *          Ciudad de Mexico, Mexico
 *   
 * August 2016
 */

#ifndef _MOEAD_H
#define _MOEAD_H

#include "utility.h"
#include "common.h"
#include "param.h"

typedef struct {
  int niche;         /* neighborhood size */
  int wsize;         /* number of weight vectors = population size */
  double *W;         /* weight vectors */
  double *min;       /* ideal point */
  double *max;       /* nadir point */
  double *diff1;     /* temporary vector for subtraction */
  double *diff2;     /* temporary vector for subtraction */
  EMO_Utility utl;   /* Utility function: tchebycheff, pbi, etc. */
  EMO_List lp;       /* list for storing population members */
  EMO_List lt1;      /* temporary list */
  EMO_List lt2;      /* temporary list */
  EMO_List *lnn;     /* array of lists for storing neighbors */
  double **sort;     /* temporary array for sorting */
  int ssize;         /* temporary variable for the population size */
  char *wfile;       /* file containing the weight vectors */
  char *utl_name;    /* scalarizing or utility function */
  int norm;          /* flag that indicates if objectives should be normalized */
  int ver;           /* MOEA/D version */
  double delta;         /* probability that parent solutions are selected from the neighborhood */
  int nr;            /* maximal number of solutions replaced by each child solution */
  void (*fdiff)(double *, double *, double *, double *, int);
} EMO_MOEAD;

void EMO_MOEAD_alloc(EMO_MOEAD *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_MOEAD_free(EMO_MOEAD *alg);
void EMO_MOEAD_run(EMO_MOEAD *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

#endif

