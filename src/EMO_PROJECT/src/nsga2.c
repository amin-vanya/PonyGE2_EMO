/* Author: Elmer Jesús Morales Orozco
 *         Universidad de la Sierra Sur,
 *         Miahuatlán, Oaxaca.
 *
 * August 2015
 */

#include <stdlib.h>
#include <string.h>

#include "nsga2.h"
#include "random.h"
#include "evop.h"
#include "niche.h"
#include "numeric.h"

#include "vector.h"

/* Load specific parameters for the algorithm */
void NSGA2_load_param(EMO_Param *param, int nvar) {
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
}

void EMO_NSGA2_alloc(EMO_NSGA2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i;

  NSGA2_load_param(param, mop->nvar);

  EMO_Population_alloc(pop, mop, param->mu, param->mu);

  EMO_NDSort_alloc(&alg->nd, pop->size);
  EMO_List_alloc(&alg->lst1, pop->mu);
  EMO_List_alloc(&alg->lst2, pop->mu);

  if((alg->filter = (int *) calloc(sizeof(int), pop->size)) == NULL) {
    printf("Error, not enough memory in NSGA2.\n");
    exit(1);
  }

  if((alg->cd = (double *) calloc(sizeof(double), pop->size)) == NULL) {
    printf("Error, not enough memory in NSGA2.\n");
    exit(1);
  }

  if((alg->min = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in NSGA2.\n");
    exit(1);
  }

  if((alg->max = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in NSGA2.\n");
    exit(1);
  }


  if((alg->parent1 = (int *) malloc(sizeof(int) * pop->mu / 2)) == NULL) {
    printf("Error, not enough memory in NSGA2.\n");
    exit(1);
  }

  if((alg->parent2 = (int *) malloc(sizeof(int) * pop->mu / 2)) == NULL) {
    printf("Error, not enough memory in NSGA2.\n");
    exit(1);
  }

  if((alg->seedp = (int *) malloc(sizeof(int) * pop->mu)) == NULL) {
    printf("Error, not enough memory in NSGA2.\n");
    exit(1);
  }

  for(i = 0; i < pop->mu; i++)
    alg->seedp[i] = i;

  if((alg->sort = (double **) malloc(sizeof(double *) * pop->size)) == NULL) {
    printf("Error, not enough memory in NSGA2.\n");
    exit(1);
  }

  for(i = 0; i < pop->size; i++) {
    if((alg->sort[i] = (double *) malloc(sizeof(double) * 2)) == NULL) {
      printf("Error, not enough memory in NSGA2.\n");
      exit(1);
    }
  }  

  alg->ssize = pop->size;
}


void EMO_NSGA2_free(EMO_NSGA2 *alg){
  int i;

  EMO_NDSort_free(&alg->nd);
  EMO_List_free(&alg->lst1);
  EMO_List_free(&alg->lst2);
  free(alg->filter);
  free(alg->cd);
  free(alg->min);
  free(alg->max);
  free(alg->parent1);
  free(alg->parent2);
  free(alg->seedp);

  for(i = 0; i < alg->ssize; i++)
    free(alg->sort[i]);

  free(alg->sort);
}


void EMO_NSGA2_run(EMO_NSGA2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int p1, p2, v1, v2, i, j, idx, cont, lastfront;
  int a;

  EMO_maxminBound(alg->max, alg->min, pop->obj, NULL, pop->mu, mop->nobj);
    
  /* Non-dominated sorting algorithm */
  if(mop->ncon == 0)
    EMO_NDSort_run(&alg->nd, pop->obj, mop->nobj, NULL, NULL, pop->mu);
  else
    EMO_NDSort_run2(&alg->nd, pop->obj, mop->nobj, pop->con, mop->ncon, NULL, pop->mu);

  for(a = 0; a < alg->nd.nfront; a++) {
    EMO_crowdingDistance2(alg->cd, alg->sort, pop->obj, &alg->nd.front[a], alg->max, alg->min, mop->nobj);
  }

  while(!EMO_Stop_end(param->stop)) {
    /* Offspring generation */
    j = 0 ;
    tournament_selection_rank_crowding(alg->parent1, param->rand, alg->seedp, alg->nd.rank, alg->cd, pop->mu);
    tournament_selection_rank_crowding(alg->parent2, param->rand, alg->seedp, alg->nd.rank, alg->cd, pop->mu);

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

    EMO_maxminBound(alg->max, alg->min, pop->obj, NULL, pop->size, mop->nobj);

    /* Non-dominated sorting algorithm */
    if(mop->ncon == 0)
      EMO_NDSort_run(&alg->nd, pop->obj, mop->nobj, NULL, NULL, pop->size);
    else
      EMO_NDSort_run2(&alg->nd, pop->obj, mop->nobj, pop->con, mop->ncon, NULL, pop->size);

    j = cont = 0;

    do {
      EMO_crowdingDistance2(alg->cd, alg->sort, pop->obj, &alg->nd.front[j], alg->max, alg->min, mop->nobj);
      cont += alg->nd.front[j++].size;
    } while(cont < pop->mu);

    lastfront = j - 1;

    memset(alg->filter, 0, sizeof(int) * pop->size);

    /* Select solutions from the last front */
    for(i = 0; i < alg->nd.front[lastfront].size; i++) {
      EMO_List_get(&alg->nd.front[lastfront], &idx, i);
      alg->filter[idx] = 1;
    }

    while(cont > pop->mu) {
      EMO_dmin(&idx, alg->cd, alg->filter, pop->size);
      alg->filter[idx] = 0;
      cont--;
    }

    /* Select solutions from the first fronts */
    for(i = 0; i < pop->size; i++) {
      if(alg->nd.rank[i] < lastfront) {
        alg->filter[i] = 1;
      }
    }

    /* Reduce population */
    EMO_Population_survive(pop, alg->nd.rank, alg->cd, mop, &alg->lst1, &alg->lst2, alg->filter);
    EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 0);
  }
}

