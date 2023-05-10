#define MAX_LEVELS 300

#include <stdlib.h>
#include <stdlib.h>
#include <unistd.h>
#include <float.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "random.h"
#include "hype.h"
#include "numeric.h"
#include "vector.h"
#include "evop.h"
#include "io.h"

void EMO_HYPE_load_param(EMO_HYPE *alg, EMO_Param *param, int nvar) {

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
  printf("pm updated %f\n", param->Pm);

  if(!EMO_Param_get_int(param, &alg->samples, "hype_samples")) {
    printf("Error, hype_samples is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_double(param, &alg->bound, "hype_bound")) {
    printf("Error, hype_bound is not defined in the configuration file.\n");
    exit(1);
  }
}

void EMO_HYPE_alloc(EMO_HYPE *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i;

  EMO_HYPE_load_param(alg, param, mop->nvar);

  if(param->mu % 2 != 0) {
    printf("Warning, adjusting population size to be an even number (%d vs %d).\n", param->mu, param->mu+1);
   param->mu++;
  }

  EMO_Population_alloc(pop, mop, param->mu, param->mu);
  EMO_NDSort_alloc(&alg->nd, pop->size);
  EMO_List_alloc(&alg->lst1, pop->mu);
  EMO_List_alloc(&alg->lst2, pop->mu);

  if((alg->filter = (int *) calloc(sizeof(int), pop->size)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->min = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->max = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->parent1 = (int *) malloc(sizeof(int) * pop->mu / 2)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->parent2 = (int *) malloc(sizeof(int) * pop->mu / 2)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->seedp = (int *) malloc(sizeof(int) * pop->mu)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  for(i = 0; i < pop->mu; i++)
    alg->seedp[i] = i;


  if((alg->tmp = (double *) malloc(sizeof(double) * (pop->size * mop->nobj))) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->val = (double *) malloc(sizeof(double) * pop->size)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->rho = (double *) malloc(sizeof(double) * (pop->size + 1) )) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->hitstat = (int *) malloc(sizeof(int) * pop->size)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->sample = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->boundsVec = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->indices = (int *) malloc(sizeof(int) * pop->size)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->pvec = (int *) malloc(sizeof(int) * pop->size)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->p = (double *) malloc(sizeof(double) * (pop->size * mop->nobj) )) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->beg = (int *) malloc(sizeof(int) * MAX_LEVELS)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->end = (int *) malloc(sizeof(int) * MAX_LEVELS)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }

  if((alg->ref = (double *) malloc(sizeof(double) * pop->size)) == NULL) {
    printf("Error, not enough memory in HYPE.\n");
    exit(1);
  }
}

void EMO_HYPE_free(EMO_HYPE *alg) {

  EMO_NDSort_free(&alg->nd);
  EMO_List_free(&alg->lst1);
  EMO_List_free(&alg->lst2);
  free(alg->filter);
  free(alg->min);
  free(alg->max);
  free(alg->parent1);
  free(alg->parent2);
  free(alg->seedp);

  free(alg->tmp);
  free(alg->rho);
  free(alg->hitstat);
  free(alg->sample);
  free(alg->val);
  free(alg->boundsVec);
  free(alg->indices);
  free(alg->pvec);
  free(alg->p);
  free(alg->beg);
  free(alg->end);
  free(alg->ref);

}

/**
 * Internal function used by hypeExact
 */
void rearrangeIndicesByColumn(EMO_HYPE *alg, double *mat, int rows, int columns, int col, int *ind) {
  int i = 0, L, R, swap;
  double pref, pind;
  
  for(i = 0; i < rows; i++)
    alg->ref[i] = mat[col + ind[i]*columns];

  i = 0;

  alg->beg[0] = 0; 
  alg->end[0] = rows;

  while(i >= 0) {
    L = alg->beg[i];
    R = alg->end[i]-1;

    if( L < R ) {
      pref = alg->ref[L];
      pind = ind[L];

      while(L < R) {
        while(alg->ref[R] >= pref && L < R)
          R--;
        if(L < R) {
          alg->ref[L] = alg->ref[R];
          ind[L++] = ind[R];
        }
        while(alg->ref[L] <= pref && L < R)
          L++;
        if(L < R) {
          alg->ref[R] = alg->ref[L];
          ind[R--] = ind[L];
        }
      }

      alg->ref[L] = pref; 
      ind[L] = pind;
      alg->beg[i+1] = L+1; 
      alg->end[i+1] = alg->end[i];
      alg->end[i++] = L;

      if(alg->end[i] - alg->beg[i] > alg->end[i-1] - alg->beg[i-1]) {
        swap = alg->beg[i]; 
        alg->beg[i] = alg->beg[i-1]; 
        alg->beg[i-1] = swap;
        swap = alg->end[i]; 
        alg->end[i] = alg->end[i-1]; 
        alg->end[i-1] = swap;
      }
    }
    else {
      i--;
    }
  }
}

/**
 * Internal function used by hypeExact
 */
void hypeExactRecursive(EMO_HYPE *alg, double *input_p, int pnts, int dim, int nrOfPnts, int actDim, double *bounds, int *input_pvec, double *fitness, double* rho, int param_k ) {
  double *tmpfit, extrusion;
  int i, j;
  //int pvec[pnts];
  //double p[pnts*dim];

  for(i = 0; i < pnts; i++) {
    fitness[i] = 0;
    alg->pvec[i] = input_pvec[i];
  }

  for(i = 0; i < pnts*dim; i++)
    alg->p[i] = input_p[i];

  rearrangeIndicesByColumn(alg, alg->p, nrOfPnts, dim, actDim, alg->pvec);
  
  for(i = 0; i < nrOfPnts; i++) {
    if( i < nrOfPnts - 1 )
      extrusion = alg->p[(alg->pvec[i+1])*dim + actDim] - alg->p[alg->pvec[i]*dim + actDim];
    else
      extrusion = bounds[actDim] - alg->p[alg->pvec[i]*dim + actDim];

    if(actDim == 0) {
      if(i+1 <= param_k)
        for(j = 0; j <= i; j++)
          fitness[alg->pvec[j]] = fitness[alg->pvec[j]] + extrusion*rho[i+1];
    }
    else if(extrusion > 0) {
      //double tmpfit[ pnts ];
      tmpfit = (double*)malloc(sizeof(double ) * pnts );
			
      hypeExactRecursive(alg, alg->p, pnts, dim, i+1, actDim-1, bounds, alg->pvec, tmpfit, alg->rho, param_k);
      for(j = 0; j < pnts; j++)
        fitness[j] += extrusion*tmpfit[j];
      free(tmpfit);
    }
  }
}

/**
 * Calculating the hypeIndicator
 * \f[ \sum_{i=1}^k \left( \prod_{j=1}^{i-1} \frac{k-j}{|P|-j} \right) \frac{ Leb( H_i(a) ) }{ i } \f]
 */
void hypeExact(EMO_HYPE *alg, int popsize, double lowerbound, double upperbound, int param_k, double *points, int dim) {
  int i;

  for(i = 0; i < dim; i++)
    alg->boundsVec[i] = alg->bound;

  for(i = 0; i < popsize; i++  )
    alg->indices[i] = i;

  /** Recursively calculate the indicator values */
   hypeExactRecursive(alg, points, popsize, dim, popsize, dim-1, alg->boundsVec, alg->indices, alg->val, alg->rho, param_k);
}

int weaklyDominates( double *point1, double *point2, int nobj) {
  int better;
  int i = 0;
  better = 1;

  while( i < nobj && better ) {
    better = point1[i] <= point2[i];
    i++;
  }
  return better;
}

/**
 * Sampling the hypeIndicator
 * \f[ \sum_{i=1}^k \left( \prod_{j=1}^{i-1} \frac{k-j}{|P|-j} \right) \frac{ Leb( H_i(a) ) }{ i } \f]
 *
 * @param[out] val vector of all indicators
 * @param[in] popsize size of the population \f$ |P| \f$
 * @param[in] lowerbound scalar denoting the lower vertex of the sampling box
 * @param[in] upperbound scalar denoting the upper vertex of the sampling box
 * @param[in] samples the total number of samples
 * @param[in] param_k the variable \f$ k \f$
 * @param[in] points matrix of all objective values dim*popsize entries
 * @param[in] rho weight coefficients
 * @pre popsize >= 0 && lowerbound <= upperbound && param_k >= 1 &&
 * 		param_k <= popsize
 */
void hypeSampling(EMO_HYPE *alg, EMO_Rand *rnd, int popsize, double lowerbound, double upperbound, int param_k, double *points, int dim) {
  int i, s, k, domCount;

  assert(popsize >= 0);
  assert(lowerbound <= upperbound);
  assert(param_k >= 1);
  assert(param_k <= popsize);

  for(s = 0; s < alg->samples; s++) {
    for(k = 0; k < dim; k++)
      alg->sample[k] =  EMO_Rand_real1(rnd, lowerbound, upperbound);

    domCount = 0;
    for(i = 0; i < popsize; i++) {
      if(weaklyDominates(points + (i * dim), alg->sample, dim)) {
        domCount++;
        if(domCount > param_k)
          break;
        alg->hitstat[i] = 1;
      }
      else {
        alg->hitstat[i] = 0;
      }
    }
    if(domCount > 0 && domCount <= param_k) {
      for(i = 0; i < popsize; i++)
        if(alg->hitstat[i] == 1)
          alg->val[i] += alg->rho[domCount];
    }
  }
  for(i = 0; i < popsize; i++)
    alg->val[i] = alg->val[i] * pow( (upperbound-lowerbound), dim ) / (double)alg->samples;
}

/**
 * Determine the hypeIndicator
 * \f[ \sum_{i=1}^k \left( \prod_{j=1}^{i-1} \frac{k-j}{|P|-j} \right) \frac{ Leb( H_i(a) ) }{ i } \f]
 *
 * if samples < 0, then do exact calculation, else sample the indicator
 *
 * @param[out] val vector of all indicator values
 * @param[in] popsize size of the population \f$ |P| \f$
 * @param[in] lowerbound scalar denoting the lower vertex of the sampling box
 * @param[in] upperbound scalar denoting the upper vertex of the sampling box
 * @param[in] samples the total number of samples or, if negative, flag
 * 		that exact calculation should be used.
 * @param[in] param_k the variable \f$ k \f$
 * @param[in] points matrix of all objective values dim*popsize entries
 * @param[in] rho weight coefficients
 */

void hypeIndicator(EMO_HYPE *alg, EMO_Rand *rnd, int popsize, double lowerbound, double upperbound, int param_k, double *points, int dim) {
  int i,j;

  /** Set alpha */
  alg->rho[0] = 0;

  for(i = 1; i <= param_k; i++) {
    alg->rho[i] = 1.0 / (double)i;

    for(j = 1; j <= i-1; j++)
      alg->rho[i] *= (double)(param_k - j ) / (double)(popsize - j);
  }

  for(i = 0; i < popsize; i++)
    alg->val[i] = 0.0;

  if(alg->samples < 0)
    hypeExact(alg, popsize, lowerbound, upperbound, param_k, points, dim);
  else
    hypeSampling(alg, rnd, popsize, lowerbound, upperbound, param_k, points, dim);
}


void getObjectivesArray( double *tmp, double *obj, int nobj, EMO_List *list) {
  int i, idx;

  for(i = 0; i < list->size; i++) {
    EMO_List_get(list, &idx, i);
    memcpy(tmp + i * nobj, obj + idx * nobj, sizeof(double) * nobj);
  }
}

void tournament_selection_rank_fit(int *ind, EMO_Rand *rand, int *seedp, int *rank, double *fit, int n) {
 int p1, p2, i, j = 0;

  EMO_Rand_shuffle(rand, seedp, n);

  for(i = 0; i < n; i+=2) {
    p1 = seedp[i];
    p2 = seedp[i+1];

    if(rank[p1] < rank[p2]) {
      ind[j++] = p1;
    }
    else if(rank[p2] < rank[p1]) {
      ind[j++] = p2;
    }
    else if(fit[p1] > fit[p2]) {
        ind[j++]=p1;
    }
    else if(fit[p2] > fit[p1]) {
      ind[j++]=p2;
    }
    else if(EMO_Rand_flip(rand, 0.5))
      ind[j++] = p1;
    else
      ind[j++] = p2;
  }
}

void EMO_HYPE_run(EMO_HYPE *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int p1, p2, v1, v2, i, j, idx, cont, lastfront;  //o1, o2,

  while(!EMO_Stop_end(param->stop)) {
   /* Offspring generation */
    j = 0 ;

    /* Non-dominated sorting algorithm */
    EMO_NDSort_run(&alg->nd, pop->obj, mop->nobj, NULL, NULL, pop->mu);

    EMO_maxminBound(alg->max, alg->min, pop->obj, NULL, pop->mu, mop->nobj);

    EMO_normalize(alg->tmp, pop->obj, NULL, pop->size, alg->min, alg->max, mop->nobj);
    hypeIndicator(alg, param->rand, pop->size, 0.0, alg->bound, pop->mu, alg->tmp, mop->nobj); //pop->obj

    tournament_selection_rank_fit(alg->parent1, param->rand, alg->seedp, alg->nd.rank, alg->val, pop->mu);
    tournament_selection_rank_fit(alg->parent2, param->rand, alg->seedp, alg->nd.rank, alg->val, pop->mu);

    for(i = 0; i < pop->lambda; i+=2) {
      p1 = alg->parent1[j];
      p2 = alg->parent2[j++];

      p1 *= mop->nvar;
      p2 *= mop->nvar;

      v1 = pop->mu + i;
      v2 = pop->mu + i + 1;
      v1 *= mop->nvar;
      v2 *= mop->nvar;
     
     /* Generate an offspring by variation operators */
      EMO_crossSBX(pop->var+v1, pop->var+v2, pop->var+p1, pop->var+p2, param->rand, mop, param->Pc, param->Nc);
      EMO_mutatePolynom(pop->var+v1, param->rand, mop, param->Pm, param->Nm);
      EMO_mutatePolynom(pop->var+v2, param->rand, mop, param->Pm, param->Nm);
      EMO_Population_evaluate(pop, mop, pop->mu + i, 2);
    }

    /* Non-dominated sorting algorithm */
    EMO_NDSort_run(&alg->nd, pop->obj, mop->nobj, NULL, NULL, pop->size);  

    j = cont = 0;

    do {
      cont += alg->nd.front[j++].size;
    } while(cont < pop->mu);

    lastfront = j;

    /* Select solutions from the last front */
    memset(alg->filter, 0, sizeof(int) * pop->size);

    if(cont > pop->mu) {
      lastfront--;

      EMO_maxminBound(alg->max, alg->min, pop->obj, NULL, pop->size, mop->nobj);

      getObjectivesArray(alg->tmp, pop->obj, mop->nobj, &alg->nd.front[lastfront]); 
      EMO_normalize(alg->tmp, alg->tmp, NULL, alg->nd.front[lastfront].size, alg->min, alg->max, mop->nobj);

      for(i = 0; i < alg->nd.front[lastfront].size; i++) {
        EMO_List_get(&(alg->nd.front[lastfront]), &idx, i);
        alg->filter[idx] = 1;
      }
    }

    while(cont > pop->mu) {
      hypeIndicator(alg, param->rand, alg->nd.front[lastfront].size, 0.0, alg->bound, cont - pop->mu, alg->tmp, mop->nobj);

      EMO_dmin(&idx, alg->val, NULL, alg->nd.front[lastfront].size);

      alg->val[idx] = DBL_MAX;
      EMO_List_get(&(alg->nd.front[lastfront]), &idx, idx);

      alg->filter[idx] = 0;
      EMO_List_remove(&alg->nd.front[lastfront], idx);

      getObjectivesArray(alg->tmp, pop->obj, mop->nobj, &alg->nd.front[lastfront]); 
      EMO_normalize(alg->tmp, alg->tmp, NULL, alg->nd.front[lastfront].size, alg->min, alg->max, mop->nobj);
      cont--;
    }

    for(i = 0; i < pop->size; i++) {
      if(alg->nd.rank[i] < lastfront)
        alg->filter[i] = 1;
    }

    /* Reduce population */
    EMO_Population_survive(pop, alg->nd.rank, alg->val, mop, &alg->lst1, &alg->lst2, alg->filter);
    EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 0);
  }
}

#undef MAX_LEVELS

