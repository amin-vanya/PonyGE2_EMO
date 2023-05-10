/* Authors: Elmer Jesús Morales Orozco
 *          Universidad de la Sierra Sur,
 *          Miahuatlán, Oaxaca.
 *
 *          Raquel Hernandez Gomez
 *          Cinvestav-IPN,
 *          Mexico, D.F.
 *       
 * August 2015
 */


#ifndef _IBEA_
#define _IBEA_

#define EMO_HV_IBEA  0 
#define EMO_EPS_IBEA 1 
#define EMO_R2_IBEA  2

#include "common.h"
#include "debug.h"
#include "param.h"
#include "stop.h"
#include "plot.h"
#include "utility.h" 

typedef struct EMO_IBEA {
  int indicator;
  double (*indf)(struct EMO_IBEA *, double *, double *, int); /* indicator-pointer function */
  double *indv;           /* matrix of indicator values */
  double *fit;            /* fitness */
  int *filter;            /* selected individuals */
  int *seedp;             /* enumeration of parent population */
  int *parent1, *parent2; /* candidates for reproduction */
  EMO_List lst1, lst2;    /* temporary lists */

  double kappa;           /* scaling factor */
  double rho;             /* scaling factor for the reference point (hv) */

  double *min;            /* ideal point (hv) */
  double *max;            /* nadir point (hv) */
  double *diff;           /* difference between max and min (eps, hv) */
  double *zref;           /* reference point (hv, r2) */
  double max_vol;         /* maximum volume (hv) */
  double *norm;           /* normalized objective functions (r2) */
  int wsize;              /* number of weight vectors (r2) */
  double *W;              /* weight vectors (r2) */
  EMO_Utility utl;        /* utility function (r2) */

  /* mpirun */
  double **sort;          /* temporary array for sorting population */
  int ssize;
  char *wfile;
} EMO_IBEA;

void EMO_IBEA_alloc(EMO_IBEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_IBEA_free(EMO_IBEA *alg);
void EMO_IBEA_run(EMO_IBEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_IBEA_prun(EMO_IBEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

#endif

