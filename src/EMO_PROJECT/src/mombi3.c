
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "dominance.h"
#include "indicator.h"
#include "numeric.h"
#include "mombi3.h"
#include "vector.h"
#include "evop.h"
#include "stat.h"
#include "sort.h"
#include "io.h"
#include "hv.h"

#define _MAXCHAR 2000

/* Load specific parameters for the algorithm */
void EMO_MOMBI3_load_param(EMO_MOMBI3 *alg, EMO_Param *param, int nvar, int nobj) {
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
}

/* Compares two individuals */
int EMO_MOMBI3_compare_min(const void **a, const void **b) {
  double *v1, *v2;

  v1 = (double *) *a;
  v2 = (double *) *b;

  if(v1[1] < v2[1]) return -1;
  else if(v1[1] > v2[1]) return 1;

  return 0;
}

/* Update reference points */
void EMO_MOMBI3_update_refpoint(EMO_MOMBI3 *alg, EMO_Param *param, double *data, int *filter, int size, int tot, int nobj) {
  int i;

  if(alg->dm == 0) {
    EMO_minBound(alg->new_min, data, filter, tot, nobj);
 
    for(i = 0; i < nobj; i++) {
      if(alg->new_min[i] < alg->min[i]) {
        alg->min[i] = alg->new_min[i];
      }
    }
  }
}

void EMO_MOMBI3_alloc_weights(EMO_MOMBI3 *alg, EMO_MOP *mop) {
  double *tmp, *v1, *v2;
  size_t size;
  int i, j;

  alg->wsize = 0;
  alg->W = EMO_File_read(NULL, &alg->wsize, &mop->nobj, alg->wfile, 0);

  if((tmp = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in mombi3.\n");
    exit(1);
  }

  size = sizeof(double) * mop->nobj;

  for(i = 0; i < alg->wsize; i++) {
    v1 = alg->W + i * mop->nobj;

    for(j = 0; j < mop->nobj; j++) {

      if(v1[j] == 1.0) {

        v2 = alg->W + j * mop->nobj;

        if(v1 != v2) {
          memcpy(tmp, v2, size);
          memcpy(v2, v1, size);
          memcpy(v1, tmp, size);
        }
        break;
      }
      else if(v1[j] != 0.0) {
        break;
      }
    }
  }

  free(tmp);
}

void EMO_MOMBI3_alloc(EMO_MOMBI3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i;

  if((alg->wfile = (char *) malloc(sizeof(char) * _MAXCHAR)) == NULL) {
    printf("Error, not enough memory in MOMBI3.\n");
    exit(1);
  }

  if((alg->min = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  EMO_MOMBI3_load_param(alg, param, mop->nvar, mop->nobj);
 
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

  if((alg->rank2 = (double *) malloc(sizeof(double) * pop->size)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->sort = (double **) malloc(sizeof(double *) * pop->size)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  for(i = 0; i < pop->size; i++) {
    if((alg->sort[i] = (double *) malloc(sizeof(double) * 2)) == NULL) {
      printf("Error, not enough memory in mombi2.\n");
      exit(1);
    }
  }  

  alg->ssize = pop->size;


   EMO_MOMBI3_alloc_weights(alg, mop);


  if((alg->H = (char **) malloc(sizeof(char *) * pop->size)) == NULL) {
    printf("Error, not enough memory in mombi3.\n");
    exit(1);
  }

  for(i = 0; i < pop->size; i++) {
    if((alg->H[i] = (char *) malloc(sizeof(char) * _MAXCHAR)) == NULL) {
      printf("Error, not enough memory in mombi3.\n");
      exit(1);
    }
  }

  if((alg->T = (char **) malloc(sizeof(char *) * pop->size)) == NULL) {
    printf("Error, not enough memory in mombi3.\n");
    exit(1);
  }

  for(i = 0; i < pop->size; i++) {
    if((alg->T[i] = (char *) malloc(sizeof(char) * _MAXCHAR)) == NULL) {
      printf("Error, not enough memory in mombi3.\n");
      exit(1);
    }
  }

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

  EMO_Utility_alloc(&alg->utl, param, mop->nobj, "all"); 

  alg->utl.ewc_p = 100;
  alg->utl.wpo2_p = 3;
  alg->utl.wn2_p = 0.5;
  alg->utl.aasf_alpha = 1e-4;

  alg->utl.wpo_p = 3;
  alg->utl.pbi_theta = 10.0;
  alg->utl.ache_alpha = 0.01;
  alg->utl.mche_alpha = 1e-2;
  alg->utl.cs_alpha = 0.02;

  if((alg->filter0 = (int *) calloc(sizeof(int), pop->size)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->filter = (int *) calloc(sizeof(int), pop->size)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->filter2 = (int *) calloc(sizeof(int), pop->size)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  alg->fcomp =  (int (*)(const void *, const void *))&EMO_MOMBI3_compare_min;

  alg->fnorm = &EMO_normalize;
  EMO_Refpoint_alloc(&alg->ref, param, mop->nobj);
}

void EMO_MOMBI3_free(EMO_MOMBI3 *alg) {
  int i;

  free(alg->wfile);
  EMO_List_free(&alg->lst1);
  EMO_List_free(&alg->lst2);
  free(alg->norm);
  free(alg->rank);
  free(alg->rank2);

  for(i = alg->ssize - 1; i > -1; i--) {
    free(alg->sort[i]);
    free(alg->H[i]);
    free(alg->T[i]);
  }

  free(alg->sort);
  free(alg->H);
  free(alg->T);

  free(alg->W);
  free(alg->min);
  free(alg->max);
  free(alg->new_min);
  free(alg->new_max);
  free(alg->hist);
  free(alg->update);

  EMO_Utility_free(&alg->utl);
  free(alg->filter0);
  free(alg->filter);
  free(alg->filter2);

  EMO_Refpoint_free(&alg->ref);
}

// Encuentra el punto de nadir y los puntos extremos
int EMO_MOMBI3_r2_ranking2(EMO_MOMBI3 *alg, EMO_Utility *utl, double *data, int *filter, int size, int mu) {
  int i, j, k, h, w, x, y, i1, i2, r, mem[utl->nobj], warn[utl->nobj];
  double uant;
  int count, start, flag = 0, idx[utl->nobj];

  static const EMO_UtilityFunction wdic[] = { EMO_Utility_wpo2,
                                              EMO_Utility_aasf,
                                              NULL
                                            };

  EMO_Dominance_ndset(alg->filter0, data, filter, size, utl->nobj, EMO_Dominance_strict);

  memset(warn, 0, sizeof(int) * utl->nobj);

  for(w = 0; w < utl->nobj; w++) {
    mem[w] = -1;
  }

  for(j = 0; j < size; j++) {
    alg->rank[j] = DBL_MAX;
    alg->H[j][0] = '\0';
    alg->T[j][0] = '\0';
  }

  for(j = 0; j < utl->nobj; j++)
    alg->max[j] = -DBL_MAX;

  h = 0;

  while(wdic[h] != NULL) {
    for(i = 0; i < utl->nobj; i++) {

      // Calculates the individual's contribution to a weight vector
      for(j = 0; j < size; j++) {
        alg->sort[j][0] = (double) j;

        if((filter != NULL && filter[j] == 0) || alg->filter0[j] == 0) {
          alg->sort[j][1] = DBL_MAX;
          continue;
        }

        alg->sort[j][1] = wdic[h](utl, alg->W + i * utl->nobj, data + j * utl->nobj);
      }

      // Sorts individuals wrt. the utility value obtained in increasing order
      qsort(alg->sort, size, sizeof(alg->sort[0]), alg->fcomp);

      uant = alg->sort[0][1]; 
      count = start = 0;

      for(w = 1; w < size; w++) {
        k = (int) alg->sort[w][0];

        if(uant == alg->sort[w][1] && uant != DBL_MAX) {
          count++;
        }
        else {
          if(count > 0) {
            for(x = 0; x < count; x++) {
              i1 = (int) alg->sort[start + x][0];

              for(y = x + 1; y <= count; y++) {
                i2 = (int) alg->sort[start + y][0];
                            
                r = EMO_Dominance_strict(data + i1 * utl->nobj, data + i2 * utl->nobj, utl->nobj);
               
                if(r == -1) {
                  alg->sort[start + x][0] = (double) i2;
                  alg->sort[start + y][0] = (double) i1;
                }
              }
            }
          }
          uant = alg->sort[w][1]; 
          start = w;
          count = 0;
        }
      }

      k = (int) alg->sort[0][0];

      if(data[k*utl->nobj + i] > alg->max[i]) {

        // Verifica que el vector de objetivos no se haya utilizado antes como punto extremo
        for(w = 0; w < i; w++) {
          if(k == mem[w]) { // ya se habia usado
            warn[w] ++;
            warn[i] ++;
          }
        }

        alg->max[i] = data[k*utl->nobj + i];
        mem[i] = k;
      }
    }
    h++;
  }

  EMO_maxBound2(alg->new_max, idx, data, NULL, size, utl->nobj);

  flag = 0;

  for(w = 0; w < utl->nobj; w++) {   // verifica puntos extremos repetidos
    if(warn[w] > 0) {
      alg->max[w] = alg->new_max[w];
      alg->rank[idx[w]] = 0;
      flag = 1;
    }
    else {
      if(mem[w] != -1)
        alg->rank[mem[w]] = 1;
    }
  }

  if(flag == 0) {
    // Al menos el punto de nadir debe encerrar a la mitad de la poblacion
    for(i = 0; i < utl->nobj; i++) {

      for(j = 0; j < size; j++) {
        alg->sort[j][0] = (double) j;

        if(filter != NULL && filter[j] == 0) {
          alg->sort[j][1] = DBL_MAX;
          continue;
        }
     
        alg->sort[j][1] = data[j* utl->nobj + i];
      }

      // Sorts individuals wrt. the utility value obtained in increasing order
      qsort(alg->sort, size, sizeof(alg->sort[0]), alg->fcomp);

      if(alg->sort[mu-1][1] > alg->max[i]) {
        flag = 1;
        alg->max[i] = alg->new_max[i];
      }
    }
  }

  return flag;
}


/* R2 ranking algorithm of the population */
void EMO_MOMBI3_r2_ranking(EMO_MOMBI3 *alg, EMO_Utility *utl, double *data, int *filter, int size) {
  int i, j, k, h, w, x, y, i1, i2, r;
  //char str[_MAXCHAR];
  double v, uant;
  int count, start, mu;

  static const EMO_UtilityFunction fdic[] = { 
                                             EMO_Utility_ewc,
                                             EMO_Utility_aasf,
                                             EMO_Utility_wpo2,
                                             EMO_Utility_wn2,
                                             EMO_Utility_asf,
                                             EMO_Utility_che,
                                             EMO_Utility_ws,
                                             NULL
                                          };


  for(j = 0; j < size; j++) {
    alg->H[j][0] = '\0';
    alg->T[j][0] = '\0';
  }

  h = 0;
  while(fdic[h] != NULL) {

    for(i = utl->nobj; i < alg->wsize; i++) {
      mu = 0;

      for(j = 0; j < size; j++) {
        if(filter == NULL || filter[j] == 1) {
          alg->sort[mu][0] = (double) j;
          alg->sort[mu][1] = fdic[h](utl, alg->W + i * utl->nobj, data + j * utl->nobj);
          mu++;
        }
      }

      // Sorts individuals wrt. the utility value obtained in increasing order
      qsort(alg->sort, mu, sizeof(alg->sort[0]), alg->fcomp);

      // Ordena weakly Pareto y repetidos
      uant = alg->sort[0][1]; 
      count = start = 0;

      for(w = 1; w < mu; w++) {
        k = (int) alg->sort[w][0];

        if(uant == alg->sort[w][1]) {
          count++;
        }
        else {
          if(count > 0) {

            for(x = 0; x < count; x++) {
              i1 = (int) alg->sort[start + x][0];

              for(y = x + 1; y <= count; y++) {
                i2 = (int) alg->sort[start + y][0];
                            
                r = EMO_Dominance_strict(data + i1 * utl->nobj, data + i2 * utl->nobj, utl->nobj);
               
                if(r == -1) {
                  alg->sort[start + x][0] = (double) i2;
                  alg->sort[start + y][0] = (double) i1;
                }
              }
            }
          }
          uant = alg->sort[w][1]; 
          start = w;
          count = 0;
        }
      }
      
      // Ranks individuals
      for(j = 0; j < mu; j++) {
        k = (int) alg->sort[j][0];
        v = (double) (j+1);

        if(v < alg->rank[k]) {
          alg->rank[k] = v;
        }
      }
    }
    h++;
  }
}

void EMO_MOMBI3_run(EMO_MOMBI3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int p1, p2, v1, v2, i, j, k, count;
  double rcut;

  EMO_maxBound(alg->max, pop->obj, NULL, pop->mu, mop->nobj);

  for(i = 0; i < mop->nobj; i++) {
    if(alg->dm == 0)
      alg->min[i] = DBL_MAX;

    j = i * alg->max_hist;
    alg->hist[j] = alg->max[i];
  }
  
  EMO_Debug_printf(param->dbg, "run MOMBI3");

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
      // update the ideal point
      EMO_MOMBI3_update_refpoint(alg, param, pop->obj, NULL, pop->mu, pop->size, mop->nobj);
      EMO_maxBound(alg->max, pop->obj, NULL, pop->mu, mop->nobj);

      for(i = 0; i < pop->size; i++)
        alg->rank[i] = DBL_MAX;

      EMO_MOMBI3_r2_ranking2(alg, &alg->utl, pop->obj, NULL, pop->size, pop->mu);
      alg->fnorm(alg->norm, pop->obj, NULL, pop->size, alg->min, alg->max, mop->nobj);
      EMO_MOMBI3_r2_ranking(alg, &alg->utl, alg->norm, NULL, pop->size);

      for(i = 0; i < pop->size; i++) {
        alg->sort[i][0] = (double) i; 
        alg->sort[i][1] = alg->rank[i];
      }

      qsort(alg->sort, pop->size, sizeof(alg->sort[0]), alg->fcomp);

      rcut = alg->sort[pop->mu - 1][1];
      memset(alg->filter, 0, sizeof(int) * pop->size);
      count = 0;

      for(i = 0; i < pop->size; i++) {
        if(alg->rank[i] <= rcut) {
          alg->filter[i] = 1;
          count++;
        }
      }

      EMO_Indicator_senergy(alg->rank2, alg->norm, alg->filter, pop->size, mop->nobj);

      while(count > pop->mu) {
        memset(alg->filter2, 0, sizeof(int) * pop->size);
        EMO_dmax(&i, alg->rank, alg->filter, pop->size);
        rcut = alg->rank[i];

        for(i = 0; i < pop->size; i++)
          if(alg->filter[i] && alg->rank[i] == rcut)
            alg->filter2[i] = 1;

        EMO_dmax(&i, alg->rank2, alg->filter2, pop->size);
        EMO_Indicator_senergy_update(alg->rank2, alg->norm, alg->filter, pop->size, mop->nobj, i, -1);
        alg->filter[i] = 0;
        count--;
      }
    }
    else { // Manejo de restricciones

      // Selecciona a todos los individuos que no violan restricciones
      memset(alg->filter, 0, sizeof(int) * pop->size);

      count = j = 0;

      for(i = 0; i < pop->size; i++) {
        if(pop->vio[i] == 0) {
          alg->filter[i] = 1;
          count++;
        }
        else {
          alg->sort[j][0] = (double) i;
          alg->sort[j][1] = pop->vio[i];
          j++;
        }
      }

      // No alcanzan los individuos para la siguiente generacion, se seleccionan los que violen menos
      if(count < pop->mu) {

        qsort(alg->sort, j, sizeof(alg->sort[0]), (int (*)(const void *, const void *))&EMO_compare_asc);
        i = 0;

        for(i = 0; i < j; i++) {
          k = (int) alg->sort[i][0];
          alg->filter[k] = 1;
          count++;

          if(count == pop->mu) {
            break;
          }
        }
      }
      else if(count > pop->mu) {  // hay muchos individuos factibles, se aplica R2

        EMO_Refpoint_update_ideal(&alg->ref, pop->obj, alg->filter, pop->size);
        EMO_shift(alg->norm, pop->obj, alg->filter, pop->size, alg->ref.ideal, mop->nobj);
        EMO_Refpoint_update_nadir(&alg->ref, alg->norm, alg->filter, pop->size);
        EMO_normalize2(alg->norm, alg->norm, alg->filter, pop->size, alg->ref.nadir, mop->nobj); 

         for(i = 0; i < pop->size; i++)
          alg->rank[i] = DBL_MAX;

         for(i = 0; i < mop->nobj; i++)
          alg->rank[alg->ref.xtrm[i]] = 0;

        EMO_MOMBI3_r2_ranking(alg, &alg->utl, alg->norm, alg->filter, pop->size);

        for(i = 0; i < pop->size; i++) {
          alg->sort[i][0] = (double) i; 
          alg->sort[i][1] = alg->rank[i];
        }

        qsort(alg->sort, pop->size, sizeof(alg->sort[0]), alg->fcomp);
        rcut = alg->sort[pop->mu - 1][1];

        for(i = 0; i < pop->size; i++) {
          if(alg->filter[i] == 1 && alg->rank[i] > rcut) {
            alg->filter[i] = 0;
            count--;
          }
        }

        if(count > pop->mu) {
          EMO_Indicator_senergy(alg->rank2, alg->norm, alg->filter, pop->size, mop->nobj);

          while(count > pop->mu) {
            memset(alg->filter2, 0, sizeof(int) * pop->size);
            EMO_dmax(&i, alg->rank, alg->filter, pop->size);
            rcut = alg->rank[i];

            for(i = 0; i < pop->size; i++) {  // Selecciona individuos de la ultima capa
              if(alg->filter[i] == 1 && alg->rank[i] == rcut) {
                alg->filter2[i] = 1;
              }
            }

            EMO_dmax(&i, alg->rank2, alg->filter2, pop->size);
 
            EMO_Indicator_senergy_update(alg->rank2, alg->norm, alg->filter, pop->size, mop->nobj, i, -1);
            alg->filter[i] = 0;
            count--;
          } 
        }
      }
    }

    /* Reduce population */
    EMO_Population_survive(pop, NULL, alg->rank, mop, &alg->lst1, &alg->lst2, alg->filter);
    EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 0);
  }
}

#undef _MAXCHAR

