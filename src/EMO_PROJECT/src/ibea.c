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


/* Adaptive IBEA */

#include <stdlib.h>
#include <float.h>
#include <string.h>
#include <math.h>

#include "ibea.h"
#include "evop.h"
#include "numeric.h"
#include "dominance.h"
#include "vector.h"
#include "sort.h" // mpi
#include "io.h" 
#include "hv.h"

#define _MAXCHAR 2000

/* calculates the maximum epsilon value by which individual a must be
 * decreased in all objectives such that individual b is weakly dominated 
 */
double calcAddEpsIndicator(EMO_IBEA *alg, double *ind_a, double *ind_b, int dim){
  double temp_eps, eps = 0;
  int i;

  eps = (ind_a[0] - ind_b[0]) / alg->diff[0];

  for (i = 1; i < dim; i++) {
    temp_eps = (ind_a[i] - ind_b[i]) / alg->diff[i];

    if (temp_eps > eps)
      eps = temp_eps;
  }
  return eps;
}

double calcHVIndicator(EMO_IBEA *alg, double *ind_a, double *ind_b, int dim) {
  double va, vb, vc, vr, c;
  int i;

  va = vb = vc = vr = 1.0;

  for(i = 0; i < dim; i++) {
    va *= alg->zref[i] - ind_a[i];
    vb *= alg->zref[i] - ind_b[i];

    c = max(ind_a[i], ind_b[i]);  /* Auxiliary point */
    vc *= alg->zref[i] - c;
  }

  /* Check if a dominates b */
  if(EMO_Dominance_strict(ind_a, ind_b, dim) == 1)
    return (vb - va) / alg->max_vol;
  else
    return (vb - vc) / alg->max_vol;
}

/* Calculates the hypervolume of that portion of the objective space
 * that is dominated by individual a but not by individual b
 * (correct version, pointed out by Dimo Brockhoff in EMO 2015)
 * hv-IBEA (recursive version)
 */
double calcHypervolume_rec(EMO_IBEA *alg, double *ind_a, double *ind_b, int dim) {
  double a, b, r, max;
  double volume = 0;

  /* Reference point for the hypervolume indicator at (1.1,...,1.1),
   * objectives must be normalized 
   */
  //r = alg->rho * (alg->max[dim-1] - alg->min[dim-1]);
  r = alg->rho * alg->diff[dim-1];
  max = alg->min[dim-1] + r;

  a = ind_a[dim-1];

  if(ind_b == NULL)
    b = max;
  else
    b = ind_b[dim-1];

  if(dim == 1) {
    if(a < b)
      volume = (b - a) / r;
    else
      volume = 0;
  }
  else {
    if(a < b) {
      volume = calcHypervolume_rec(alg, ind_a, NULL, dim - 1) * (b - a) / r;
      volume += calcHypervolume_rec(alg, ind_a, ind_b, dim - 1) * (max - b) / r;
    }
    else {
      volume = calcHypervolume_rec(alg, ind_a, ind_b, dim - 1) * (max - a) / r; 
    }
  }
  return volume;
}

double calcHypervolume(EMO_IBEA *alg, double *ind_a, double *ind_b, int dim) {
  double v, w;

  if(EMO_Dominance_strict(ind_a, ind_b, dim) == 1) // a dominates b 
    v = -calcHypervolume_rec(alg, ind_a, ind_b, dim);
  else 
    v = calcHypervolume_rec(alg, ind_b, ind_a, dim);

  w = calcHVIndicator(alg, ind_a, ind_b, dim);

  printf(" v vs w , %f %f\n", v, w);

  return v;
}

double calcR2Indicator(EMO_IBEA *alg, double *ind_a, double *ind_b, int dim) {
  double va, vb, sab = 0, sa = 0;
  int i;

  for (i = 0; i < alg->wsize; i++) {
    va = alg->utl.uf(&alg->utl, &alg->W[i * dim], ind_a);
    vb = alg->utl.uf(&alg->utl, &alg->W[i * dim], ind_b);

    sa += va / (double) alg->wsize;

    if(va < vb)
      sab += va / (double) alg->wsize;
    else
      sab += vb / (double) alg->wsize ;
  }

  return sa - sab;
}


void apply_ref_point(EMO_IBEA *alg, EMO_Population *pop, int dim, int size){
  double d;
  int i, j;
 
  d = EMO_dmax(NULL, alg->diff, NULL, dim);

  for(i = 0; i < dim; i++)
    alg->zref[i] = alg->min[i] - 2.0 * d + alg->diff[i];
  
  for(i=0; i < size; i++) {
    j = i * dim;
    EMO_vdiff(&alg->norm[j], &pop->obj[j], alg->zref, dim);
  }
}


void fitness_assignment(EMO_IBEA *alg, EMO_Population *pop, EMO_MOP *mop, double *factor, int size) {
  double *obj = NULL;
  int i, j, k, o1, o2;
  double d, e;

  /* upper and lower bounds */
  EMO_maxminBound(alg->max, alg->min, pop->obj, NULL, size, mop->nobj);

  /* difference between upper and lower bounds */
  EMO_vdiff(alg->diff,  alg->max,  alg->min, mop->nobj);

  switch(alg->indicator) {
    case EMO_HV_IBEA:

      /* calculate the reference point and the maximum volume to achieve */
      alg->max_vol = 1;

      for(i = 0; i < mop->nobj; i++) {
        alg->zref[i] = alg->rho * alg->diff[i] + alg->min[i];
        alg->max_vol *= (alg->zref[i] - alg->min[i]);
      }
      obj = pop->obj;

      break;

    case EMO_EPS_IBEA:
      obj = pop->obj;
      break;

    case EMO_R2_IBEA:
      apply_ref_point(alg, pop, mop->nobj, size);
      obj = alg->norm;
  }

  *factor = -DBL_MAX;

  for(i = 0; i < size; i++) {
    o1 = i * mop->nobj;

    /* Calculates indicator values and determines the maximum absolute indicator value */
    for(j = 0; j < size; j++) {

      if(i != j) {
        o2 = j * mop->nobj;
        k = i * pop->size + j;

        alg->indv[k] = alg->indf(alg, obj+o1, obj+o2, mop->nobj);
        d = fabs(alg->indv[k]);

        if(d > *factor) *factor = d;
      }
    }
  }

  if(alg->indicator == EMO_R2_IBEA)
    e = alg->kappa;
  else
    e = (*factor) * alg->kappa;

  for(i = 0; i < size; i++) {
    alg->fit[i] = 0;

    for(j = 0; j < size; j++) {
      if(i != j) {
        k = j * pop->size + i;

        alg->fit[i] -= exp(-alg->indv[k] / e);
      }
    }
  }
}

/* Load specific parameters for the algorithm */
void IBEA_load_param(EMO_IBEA *alg, EMO_Param *param, int nvar) {

  if(!EMO_Param_get_int(param, &param->mu, "psize")) {
    printf("Error, psize is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_double(param, &param->Pc, "pc")) {
    printf("Error, pc is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_double(param, &param->Pm, "pm")) {
    printf("Error, pm is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_double(param, &param->Nc, "nc")) {
    printf("Error, nc is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_double(param, &param->Nm, "nm")) {
    printf("Error, nm is not defined in the configuration file.\n");
    exit(1);
  }

  param->Pm = (param->Pm == -1)? 1.0 / (double) nvar : param->Pm ;
  EMO_Debug_printf(param->dbg, "pm updated %f", param->Pm);

  if(strstr(param->subalgorithm, "R2") != NULL) {
    if(!EMO_Param_get_double(param, &alg->kappa, "r2ibea_kappa")) {
      printf("Error, r2ibea_kappa is not defined in the configuration file.\n");
      exit(1);
    }
  }
  else {
    if(!EMO_Param_get_double(param, &alg->kappa, "ibea_kappa")) {
      printf("Error, ibea_kappa is not defined in the configuration file.\n");
      exit(1);
    }
  }
 
  if(strstr(param->subalgorithm, "HV") != NULL)
    alg->indicator = EMO_HV_IBEA;
  else if(strstr(param->subalgorithm, "EPS") != NULL)
    alg->indicator = EMO_EPS_IBEA;
  else if(strstr(param->subalgorithm, "R2") != NULL)
    alg->indicator = EMO_R2_IBEA;
  else if(!EMO_Param_get_int(param, &alg->indicator, "indicator")) {
    printf("Error, indicator is not defined in the configuration file.\n");
    exit(1);
  } 

  if(alg->indicator == EMO_R2_IBEA) {
    if(!EMO_Param_get_char(param, alg->wfile, "wfile")) {
      printf("Error, wfile is not defined in the configuration file.\n");
      exit(1);
    } 
  }
  else if(alg->indicator == EMO_HV_IBEA) {
    if(!EMO_Param_get_double(param, &alg->rho, "hvibea_rho")) {
      printf("Error, hvibea_rho is not defined in the configuration file.\n");
      exit(1);
    }
  }
}

void EMO_IBEA_alloc(EMO_IBEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i;

  printf("IBEA\n");

  if((alg->wfile = (char *) malloc(sizeof(char) * _MAXCHAR)) == NULL) {
    printf("Error, not enough memory in IBEA.\n");
    exit(1);
  }

  IBEA_load_param(alg, param, mop->nvar);

  if(param->mu % 2 != 0) {
    printf("Warning, adjusting population size to be an even number (%d vs %d).\n", param->mu, param->mu+1);
   param->mu++;
  }

  EMO_Population_alloc(pop, mop, param->mu, param->mu);

  switch(alg->indicator) {
    case EMO_HV_IBEA:
      alg->indf = calcHVIndicator;
      break;

    case EMO_EPS_IBEA:
      alg->indf = calcAddEpsIndicator;
      break;

    case EMO_R2_IBEA:
      alg->indf = calcR2Indicator;
      break;

    default:
      printf("Error, indicador no definido.\n");
      exit(1);
  }

  if((alg->indv = (double *) malloc(sizeof(double) * pop->size * pop->size)) == NULL) {
    printf("Error, not enough memory in ibea.\n");
    exit(1);
  }

  if((alg->fit = (double *) malloc(sizeof(double) * pop->size)) == NULL) {
    printf("Error, not enough memory in ibea.\n");
    exit(1);
  }

  if((alg->filter = (int *) calloc(sizeof(int), pop->size)) == NULL) {
    printf("Error, not enough memory in ibea.\n");
    exit(1);
  }

  if((alg->seedp = (int *) malloc(sizeof(int) * pop->mu)) == NULL) {
    printf("Error, not enough memory in ibea.\n");
    exit(1);
  }

  for(i = 0; i < pop->mu; i++)
    alg->seedp[i] = i;

  if((alg->parent1 = (int *) malloc(sizeof(int) * pop->mu / 2)) == NULL) {
    printf("Error, not enough memory in ibea.\n");
    exit(1);
  }

  if((alg->parent2 = (int *) malloc(sizeof(int) * pop->mu / 2)) == NULL) {
    printf("Error, not enough memory in ibea.\n");
    exit(1);
  }

  EMO_List_alloc(&alg->lst1, pop->mu);
  EMO_List_alloc(&alg->lst2, pop->mu);

  if((alg->min = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in ibea.\n");
    exit(1);
  }

  if((alg->max = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in ibea.\n");
    exit(1);
  }

  if((alg->diff = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in ibea.\n");
    exit(1);
  }

  if((alg->zref = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in ibea.\n");
    exit(1);
  }
 
  if((alg->norm = (double *) malloc(sizeof(double) * pop->size * mop->nobj)) == NULL) {
    printf("Error, not enough memory in ibea.\n");
    exit(1);
  }

  if(alg->indicator == EMO_R2_IBEA) {
    alg->wsize = 0;
    alg->W = EMO_File_read(NULL, &alg->wsize, &mop->nobj, alg->wfile, 0); 
    EMO_Utility_alloc(&alg->utl, param, mop->nobj, "chebyshev");
  }

  #ifdef EMO_MPI
    if((alg->sort = (double **) malloc(sizeof(double *) * pop->size)) == NULL) {
      printf("Error, not enough memory in ibea.\n");
      exit(1);
    }

    for(i = 0; i < pop->size; i++) {
      if((alg->sort[i] = (double *) malloc(sizeof(double) * 2)) == NULL) {
        printf("Error, not enough memory in ibea.\n");
        exit(1);
      }
    }

    alg->ssize = pop->size;
  #endif
}


void EMO_IBEA_free(EMO_IBEA *alg) {

  #ifdef EMO_MPI
    int i;
  #endif
 
  free(alg->wfile);
  free(alg->indv);
  free(alg->fit);
  free(alg->filter);
  free(alg->seedp);
  free(alg->parent1);
  free(alg->parent2);
  EMO_List_free(&alg->lst1);
  EMO_List_free(&alg->lst2);
  free(alg->min);
  free(alg->max);
  free(alg->diff);
  free(alg->zref);
  free(alg->norm);

 if(alg->indicator == EMO_R2_IBEA) {
    free(alg->W);
    EMO_Utility_free(&alg->utl);
  }
 
  #ifdef EMO_MPI
    for(i = alg->ssize - 1; i > -1; i--)
      free(alg->sort[i]);

    free(alg->sort);
  #endif
}


void EMO_IBEA_run(EMO_IBEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int p1, p2, v1, v2, i, j, idx;
  double c;

  fitness_assignment(alg, pop, mop, &c, pop->mu);
    
  EMO_Debug_printf(param->dbg, "Run IBEA");

  while(!EMO_Stop_end(param->stop)) {

    /* Mating selection */
    if(alg->indicator == EMO_R2_IBEA) {
      tournament_selection_r2(alg->parent1, param->rand, alg->seedp, alg->indv, pop->mu, pop->size);
      tournament_selection_r2(alg->parent2, param->rand, alg->seedp, alg->indv, pop->mu, pop->size);
    }
    else {
      tournament_selection_fitness(alg->parent1, param->rand, alg->seedp, alg->fit, pop->mu);
      tournament_selection_fitness(alg->parent2, param->rand, alg->seedp, alg->fit, pop->mu);
    }

   /* Variation */
    j = 0;

    for(i = 0; i < pop->lambda; i+=2) {

      p1 = alg->parent1[j];
      p2 = alg->parent2[j++];
      p1 *= mop->nvar;
      p2 *= mop->nvar;

      v1 = pop->mu + i;
      v2 = pop->mu + i + 1;
      v1 *= mop->nvar;
      v2 *= mop->nvar;

      EMO_crossSBX(pop->var+v1, pop->var+v2, pop->var+p1, pop->var+p2, param->rand, mop, param->Pc, param->Nc);
      EMO_mutatePolynom(pop->var+v1, param->rand, mop, param->Pm, param->Nm);
      EMO_mutatePolynom(pop->var+v2, param->rand, mop, param->Pm, param->Nm);
      EMO_Population_evaluate(pop, mop, pop->mu + i, 2);
    }
 
    fitness_assignment(alg, pop, mop, &c, pop->size);

    if(alg->indicator == EMO_R2_IBEA)
      c = 1.0;

    /* Environmental selection */
    for(i = 0; i < pop->size; i++)
      alg->filter[i] = 1;

    for(i = 0; i < pop->lambda; i++) {
      /* Find the individual with the smallest fitness value */
      EMO_dmin(&idx, alg->fit, alg->filter, pop->size);

      alg->filter[idx] = 0;  /* Ignore its itself value */

      /* Update fitness values */
      for(j = 0; j < pop->size; j++) {
        if(alg->filter[j]) {
          alg->fit[j] += exp(-alg->indv[idx * pop->size + j] / (c * alg->kappa));
        }
      }
    }

    EMO_Population_survive(pop, NULL, alg->fit, mop, &alg->lst1, &alg->lst2, alg->filter);
    EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 0);
  }
}

#undef _MAXCHAR

