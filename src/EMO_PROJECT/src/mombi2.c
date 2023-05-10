
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "dominance.h"
#include "indicator.h"
#include "numeric.h"
#include "mombi2.h"
#include "vector.h"
#include "evop.h"
#include "stat.h"
#include "sort.h"
#include "io.h"
#include "hv.h"

#define _MAXCHAR 2000

/* Load specific parameters for the algorithm */
void EMO_MOMBI2_load_param(EMO_MOMBI2 *alg, EMO_Param *param, int nvar, int nobj) {
  int n;

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

  if(!EMO_Param_get_char(param, alg->wfile, "wfile")) {
    printf("Error, wfile is not defined in the configuration file.\n");
    exit(1);
  } 

  if(!EMO_Param_get_int(param, &alg->max_hist, "mombi2_record")) {
    printf("Error, mombi2_record is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_double(param, &alg->alpha, "mombi2_alpha")) {
    printf("Error, mombi2_alpha is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_double(param, &alg->epsilon, "mombi2_epsilon")) {
    printf("Error, mombi2_epsilon is not defined in the configuration file.\n");
    exit(1);
  }

  n = nobj;

  if(!EMO_Param_get_vector_double(param, alg->min, &n, "mombi2_refpoint")) {
    printf("Error in the definition of mombi2_refpoint, see configuration file.\n");
    exit(1);
  }

  alg->dm = (n == 0)? 0 : 1;

  if(alg->dm && nobj != n) {
    printf("Error, mismatch dimensions of refpoint and nobj (%d vs %d) in configuration file.\n", n, nobj);
    exit(1);
  }

  if(!EMO_Param_get_char(param, alg->utl_name, "mombi2_utility")) {
    strcpy(alg->utl_name, "achievement_scalarizing_function");
    printf("Warning, mombi2_utility is not defined, using %s as default.\n", alg->utl_name);
  }
}

/* Compares two individuals */
int EMO_MOMBI2_compare_min(const void **a, const void **b) {
  double *v1, *v2;
  //int i, j;

  v1 = (double *) *a;
  v2 = (double *) *b;

  if(v1[1] < v2[1]) return -1;
  else if(v1[1] > v2[1]) return 1;

  if(v1[2] < v2[2]) return -1;
  else if(v1[2] > v2[2]) return 1;

  return 0;
}

/* Compares two individuals */
int EMO_MOMBI2_compare_max(const void **a, const void **b) {
  double *v1, *v2;

  v1 = (double *) *a;
  v2 = (double *) *b;

  if(v1[1] > v2[1]) return -1;
  else if(v1[1] < v2[1]) return 1;

  if(v1[2] < v2[2]) return -1;
  else if(v1[2] > v2[2]) return 1;

  return 0;
}

/* Update reference points */
void EMO_MOMBI2_update_refpoint(EMO_MOMBI2 *alg, EMO_Param *param, double *data, int *filter, int size, int nobj) {
  static int gen = 0;
  double mu, v;
  int i, j, k;

  if(alg->dm == 0) {
    EMO_minBound(alg->new_min, data, filter, 2*size, nobj);
 
    for(i = 0; i < nobj; i++) {
      if(alg->new_min[i] < alg->min[i]) {
        alg->min[i] = alg->new_min[i];
      }
    }
  }

  EMO_maxBound(alg->new_max, data, filter, size, nobj);

  for(i = 0; i < nobj; i++) {
    k = gen % (alg->max_hist);
    alg->hist[i * alg->max_hist + k] = alg->new_max[i];
  }
 
  if(gen >= alg->max_hist - 1)  {

    for(i = 0; i < nobj; i++) { 
      mu = EMO_mean(alg->hist + i * alg->max_hist, NULL, alg->max_hist);
      v = EMO_var(alg->hist + i * alg->max_hist, NULL, mu, alg->max_hist);

      if(v > alg->alpha) {
        v = EMO_dmax(NULL, alg->new_max, NULL, nobj);
 
        for(j = 0; j < nobj; j++)
          alg->max[j] = v;
 
        break;
      }
      else {

        if(fabs(alg->max[i] - alg->min[i]) < alg->epsilon) {
          v = EMO_dmax(NULL, alg->max, NULL, nobj);
          alg->max[i] = v;
          alg->update[i] = alg->max_hist;
        }
        else if(alg->new_max[i] > alg->max[i]) {
          alg->max[i] = alg->new_max[i] + fabs(alg->new_max[i] - alg->max[i]);
          alg->update[i] = alg->max_hist;
        }
        else if(v == 0 && (alg->new_max[i] <= mu || (alg->new_max[i] - mu) > alg->epsilon) && alg->update[i] == 0) {
          v = EMO_dmax(NULL, alg->hist + i * alg->max_hist, NULL, alg->max_hist);
          alg->max[i] = (alg->max[i] + v) / 2.0;
          alg->update[i] = alg->max_hist;
        }
      }

      if(alg->update[i] > 0) alg->update[i]--;
    }
  }
  gen++;
}

void EMO_MOMBI2_alloc(EMO_MOMBI2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i;

  if((alg->utl_name = (char *) malloc(sizeof(char) * _MAXCHAR)) == NULL) {
    printf("Error, not enough memory in MOMBI2.\n");
    exit(1);
  }

  if((alg->wfile = (char *) malloc(sizeof(char) * _MAXCHAR)) == NULL) {
    printf("Error, not enough memory in MOMBI2.\n");
    exit(1);
  }

  if((alg->min = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  EMO_MOMBI2_load_param(alg, param, mop->nvar, mop->nobj);
 
  if(param->mu % 2 != 0) {
    printf("Error, population size must be even.\n");
    exit(1);
  }

  EMO_Population_alloc(pop, mop, param->mu, param->mu);
  EMO_List_alloc(&alg->lst1, pop->mu);
  EMO_List_alloc(&alg->lst2, pop->mu);

  if((alg->norm = (double *) malloc(sizeof(double) * pop->size * mop->nobj)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->rank = (double *) malloc(sizeof(double) * pop->size)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->l2 = (double *) malloc(sizeof(double) * pop->size)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->sort = (double **) malloc(sizeof(double *) * pop->size)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  for(i = 0; i < pop->size; i++) {
    if((alg->sort[i] = (double *) malloc(sizeof(double) * 3)) == NULL) { 
      printf("Error, not enough memory in mombi2.\n");
      exit(1);
    }
  }  

  alg->ssize = pop->size;

  alg->wsize = 0;
  alg->W = EMO_File_read(NULL, &alg->wsize, &mop->nobj, alg->wfile, 0);

  if((alg->max = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->new_min = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->new_max = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->hist = (double *) malloc(sizeof(double) * alg->max_hist * mop->nobj)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->update = (int *) malloc(sizeof(int) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  for(i = 0; i < mop->nobj; i++) {
    alg->update[i] = 0;
  }

  EMO_Utility_alloc(&alg->utl, param, mop->nobj, alg->utl_name); 

  if((alg->filter = (int *) calloc(sizeof(int), pop->size)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  alg->fcomp =  (int (*)(const void *, const void *))&EMO_MOMBI2_compare_min;

  if(alg->utl.inverted)  /* Inverted utility function */
    alg->fnorm = &EMO_normalize3;
  else 
    alg->fnorm = &EMO_normalize;
}

void EMO_MOMBI2_free(EMO_MOMBI2 *alg) {
  int i;

  free(alg->wfile);
  free(alg->utl_name);
  EMO_List_free(&alg->lst1);
  EMO_List_free(&alg->lst2);
  free(alg->norm);
  free(alg->rank);
  free(alg->l2);

  for(i = alg->ssize - 1; i > -1; i--)
    free(alg->sort[i]);

  free(alg->sort);

  free(alg->W);
  free(alg->min);
  free(alg->max);
  free(alg->new_min);
  free(alg->new_max);
  free(alg->hist);
  free(alg->update);

  EMO_Utility_free(&alg->utl);
  free(alg->filter);
}

/* R2 ranking algorithm of the population */
void EMO_MOMBI2_r2_ranking(EMO_MOMBI2 *alg, EMO_Utility *utl, double *data, int *filter, int size) {
  int i, j, k;

   if(utl->uf == EMO_Utility_tlpbi || utl->uf == EMO_Utility_qpbi)
    EMO_Utility_update_dstar(utl, alg->min, alg->max);

  if(filter == NULL) {
    for(j = 0; j < size; j++) {
      alg->rank[j] = DBL_MAX;
      alg->l2[j] = EMO_vnorm(data + j * utl->nobj, 2.0, utl->nobj);
    }

    for(i = 0; i < alg->wsize; i++) {
      // Calculates the individual's contribution to a weight vector
      for(j = 0; j < size; j++) {
        alg->sort[j][0] = j;
        alg->sort[j][1] = utl->uf(utl, alg->W + i * utl->nobj, data + j * utl->nobj);
        alg->sort[j][2] = alg->l2[j];
      }
      // Sorts individuals wrt. the utility value obtained in increasing order
      qsort(alg->sort, size, sizeof(alg->sort[0]), alg->fcomp);

      // Ranks individuals
      for(j = 1; j <= size; j++) {
        k = (int) alg->sort[j-1][0]; 

        if((double) j < alg->rank[k])
           alg->rank[k] = (double) j;
      }
    }
  }
  else {

    for(j = 0; j < size; j++) {
      alg->rank[j] = DBL_MAX;

      if(filter[j])
        alg->l2[j] = EMO_vnorm(data + j * utl->nobj, 2.0, utl->nobj);
    }

    for(i = 0; i < alg->wsize; i++) {
      // Calculates the individual's contribution to a weight vector
      for(j = 0; j < size; j++) {

        alg->sort[j][0] = j;

        if(filter[j]) {
          alg->sort[j][1] = utl->uf(utl, alg->W + i * utl->nobj, data + j * utl->nobj);
          alg->sort[j][2] = alg->l2[j];
        }
        else {
          alg->sort[j][1] = DBL_MAX;
          alg->sort[j][2] = DBL_MAX;
        }
      }
      // Sorts individuals wrt. the utility value obtained in increasing order
      qsort(alg->sort, size, sizeof(alg->sort[0]), alg->fcomp);

      // Ranks individuals
      for(j = 1; j <= size; j++) {
        k = (int) alg->sort[j-1][0]; 

        if((double) j < alg->rank[k])
           alg->rank[k] = (double) j;
      }
    }
  }
}


void EMO_MOMBI2_run(EMO_MOMBI2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int p1, p2, v1, v2, i, j, k, cont; //o1, o2,

  EMO_maxBound(alg->max, pop->obj, NULL, pop->mu, mop->nobj);

  for(i = 0; i < mop->nobj; i++) {
    if(alg->dm == 0)
      alg->min[i] = DBL_MAX;

    j = i * alg->max_hist;
    alg->hist[j] = alg->max[i];
  }
  
  EMO_Debug_printf(param->dbg, "run MOMBI2");

  while(!EMO_Stop_end(param->stop)) {
    /* Offspring generation */
    for(i = 0; i < pop->lambda; i+=2) {
      // Tournament selection
      /* Randomly select two parents */
      p1 = EMO_Rand_int1(param->rand, 0, pop->mu-1);
      while ((p2 = EMO_Rand_int1(param->rand, 0, pop->mu-1)) == p1);

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

    if(mop->ncon == 0) {
      EMO_MOMBI2_update_refpoint(alg, param, pop->obj, NULL, pop->mu, mop->nobj);
      alg->fnorm(alg->norm, pop->obj, NULL, pop->size, alg->min, alg->max, mop->nobj);

      EMO_MOMBI2_r2_ranking(alg, &alg->utl, alg->norm, NULL, pop->size);

      for(i = 0; i < pop->size; i++) {
        alg->sort[i][0] = (double) i; 
        alg->sort[i][1] = alg->rank[i];
      }

      qsort(alg->sort, pop->size, sizeof(alg->sort[0]), (int (*)(const void *, const void *))&EMO_compare_asc);
      memset(alg->filter, 0, sizeof(int) * pop->size);

      for(i = 0; i < pop->mu; i++) {
        j = (int) alg->sort[i][0];
        alg->filter[j] = 1; 
      }
    }
    else {  // Manejo de restricciones

      // Selecciona a todos los individuos que no violan restricciones
      memset(alg->filter, 0, sizeof(int) * pop->size);

      cont = j = 0;

      for(i = 0; i < pop->size; i++) {
        if(pop->vio[i] == 0) {
          alg->filter[i] = 1;
          cont++;
        }
        else {
          alg->sort[j][0] = (double) i;
          alg->sort[j][1] = pop->vio[i];
          alg->sort[j][2] = 0; 
          j++;
        }
      }

      // No alcanzan los individuos para la siguiente generacion, se seleccionan los que violen menos
      if(cont < pop->mu) {

        qsort(alg->sort, j, sizeof(alg->sort[0]), (int (*)(const void *, const void *))&EMO_compare_asc);
        i = 0;

        for(i = 0; i < j; i++) {
          k = (int) alg->sort[i][0];
          alg->filter[k] = 1; 
          cont++;

          if(cont == pop->mu)
            break;
        }
      }
      else if(cont > pop->mu) {  // hay muchos individuos factibles, se aplica R2

        EMO_MOMBI2_update_refpoint(alg, param, pop->obj, alg->filter, pop->mu, mop->nobj);
        alg->fnorm(alg->norm, pop->obj, alg->filter, pop->size, alg->min, alg->max, mop->nobj);
        EMO_MOMBI2_r2_ranking(alg, &alg->utl, alg->norm, alg->filter, pop->size);

        for(i = 0; i < pop->size; i++) {
          alg->sort[i][0] = (double) i; 
          alg->sort[i][1] = alg->rank[i];
        }

        qsort(alg->sort, pop->size, sizeof(alg->sort[0]), (int (*)(const void *, const void *))&EMO_compare_asc);
        memset(alg->filter, 0, sizeof(int) * pop->size);

        for(i = 0; i < pop->mu; i++) {
          j = (int) alg->sort[i][0];
          alg->filter[j] = 1; 
        }
      }
    }

    /* Reduce population */
    EMO_Population_survive(pop, NULL, alg->rank, mop, &alg->lst1, &alg->lst2, alg->filter);
    EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 0);
  } 
}

#undef _MAXCHAR

