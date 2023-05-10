
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "smsemoa.h"
#include "numeric.h"
#include "vector.h"
#include "evop.h"
#include "sort.h"
#include "hv.h"

#include "io.h"

/* Load specific parameters for the algorithm */
void SMSEMOA_load_param(EMO_SMSEMOA *alg, EMO_Param *param, int nvar) {
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

  if(!EMO_Param_get_double(param, &alg->sfactor, "smsemoa_sfactor")) {
    printf("Error, smsemoa_sfactor is not defined in the configuration file.\n");
    exit(1);
  }

  if(alg->sfactor == 0) {  // random number
    alg->sfactor = EMO_Rand_prob3(param->rand);
    printf("Warning, sfactor in smsemoa.c is set randomly from the interval (0,1): %f\n", alg->sfactor);
  }

  if(!EMO_Param_get_int(param, &alg->iwfg_flag, "smsemoa_iwfg")) {
    printf("Error, smsemoa_iwfg is not defined in the configuration file.\n");
    exit(1);
  }

  param->Pm = (param->Pm == -1)? 1.0 / (double) nvar : param->Pm ;
  EMO_Debug_printf(param->dbg, "pm updated %f", param->Pm);
}

// in in in out
void EMO_SMSEMOA_alloc(EMO_SMSEMOA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i, n;

  printf("smsemoa\n");

  SMSEMOA_load_param(alg, param, mop->nvar);

  #ifdef EMO_MPI
  EMO_Population_alloc(pop, mop, param->mu, param->mu);
  #else
  EMO_Population_alloc(pop, mop, param->mu, 1);
  #endif

  n = param->mu + 1;

  EMO_NDSort_alloc(&alg->nd, n);

  if(alg->iwfg_flag)
    EMO_IWFG_alloc(&alg->iwfg, n, mop->nobj);
  else
    EMO_HV_alloc(&alg->hv, n, mop->nobj);

  if((alg->max = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in smsemoa\n");
    exit(1);
  }

  if((alg->min = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in smsemoa\n");
    exit(1);
   }

  if((alg->filter = (int *) calloc(sizeof(int), n)) == NULL) {
    printf("Error, not enough memory in smsemoa\n");
    exit(1);
   }

  if((alg->chv = (double *) malloc(sizeof(double) * n)) == NULL) {
    printf("Error, not enough memory in smsemoa\n");
    exit(1);
  }

  EMO_List_alloc(&alg->lst1, pop->mu);
  EMO_List_alloc(&alg->lst2, pop->mu);

  if((alg->sort = (double **) malloc(sizeof(double *) * n)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  for(i = 0; i < n; i++) {
    if((alg->sort[i] = (double *) malloc(sizeof(double) * 2)) == NULL) {
      printf("Error, not enough memory in mombi2.\n");
      exit(1);
    }
  }  

  alg->ssize = n;
}

void EMO_SMSEMOA_free(EMO_SMSEMOA *alg) {
  int i;

  EMO_NDSort_free(&alg->nd);

  if(alg->iwfg_flag)
    EMO_IWFG_free(&alg->iwfg);
  else
    EMO_HV_free(&alg->hv);

  free(alg->max);
  free(alg->min);
  free(alg->filter);
  free(alg->chv);
  EMO_List_free(&alg->lst1);
  EMO_List_free(&alg->lst2);

  for(i = alg->ssize - 1; i > -1; i--)
    free(alg->sort[i]);

  free(alg->sort);
}

void EMO_SMSEMOA_run(EMO_SMSEMOA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i, j, x, elem, n;
  double v = 0;

  n = pop->mu + 1;
  x = pop->mu * mop->nvar;

  EMO_Debug_printf(param->dbg, "Run SMS-EMOA");

  while(!EMO_Stop_end(param->stop)) {
    /* Randomly select two parents */
    i = EMO_Rand_int1(param->rand, 0, pop->mu-1);
    while ((j = EMO_Rand_int1(param->rand, 0, pop->mu-1)) == i);

    /* Generate an offspring by variation operators */
    i *= mop->nvar;
    j *= mop->nvar;

    EMO_crossSBX(pop->var+x, pop->vdummy, pop->var+i, pop->var+j, param->rand, mop, param->Pc, param->Nc);
    EMO_mutatePolynom(pop->var+x, param->rand, mop, param->Pm, param->Nm);
    EMO_Population_evaluate(pop, mop, pop->mu, 1);

    /* Reduce population */

    EMO_NDSort_run(&alg->nd, pop->obj, mop->nobj, NULL, NULL, n);  /* Non-dominated sorting algorithm */

    i = alg->nd.nfront - 1;  /* Last front */

    if(alg->nd.front[i].size > 1) {

      /* Select solutions from the last front */
      memset(alg->filter, 0, sizeof(int) * n);

      for(j = 0; j < alg->nd.front[i].size; j++) {
        EMO_List_get(&alg->nd.front[i], &elem, j);
        alg->filter[elem] = 1;
      }

      EMO_maxminBound(alg->max, alg->min, pop->obj, NULL, n, mop->nobj);  /* Calculate reference point, last f? */
      EMO_vsum(alg->max, alg->max, alg->sfactor, mop->nobj);  /* Add sfactor to all entrance of the reference point vector */

      /* Calculate hipervolume contributions */

      if(alg->iwfg_flag) {
        j = EMO_IWFG_run(&alg->iwfg, pop->obj, alg->filter, n, alg->max, &v);
      }
      else {
        EMO_HV_contribution(&alg->hv, alg->chv, pop->obj, alg->filter, n, alg->max, mop->nobj);

        /* Find worst contribution */
        EMO_dmin(&j, alg->chv, NULL, n);
      }
    }
    else {
      EMO_List_get(&alg->nd.front[i], &j, 0);  /* One individual is in the last front */

      EMO_List_retrieve(&alg->nd.front[i], &j, 0);
      alg->nd.nfront--;  // Update fronts
    }

    if(j != n - 1)
     EMO_Population_copy(pop, NULL, alg->chv, mop, j, n - 1); /* Remove the worst individual */

    EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 0);
  }
}

