#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "movap.h"
#include "numeric.h"
#include "vector.h"
#include "list.h"
#include "evop.h"
#include "sort.h"
#include "stat.h"


/* Load specific parameters for the algorithm */
void EMO_MOVAP_load_param(EMO_MOVAP *alg, EMO_Param *param, int nvar) {
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


  if(!EMO_Param_get_int(param, &alg->xrs, "movap_xrs")) {
    printf("Error, movap_xrs is not defined in the configuration file.\n");
    exit(1);
  }

}

/* Normalize population */
int EMO_MOVAP_normalize(EMO_MOVAP *alg, double *data, int size, int nobj) {
  int i, j, k, imin, update = 0;
  static double dmin = 0;
  double *p, v, vmin;

  EMO_List_clear(&alg->lst);
  EMO_minBound(alg->min, data, NULL, size, nobj);
  vmin = EMO_dmin(NULL, alg->min, NULL, nobj);

  /* If objectives are negative, they are shifted to the origin */
  if(vmin < 0) {

    if(vmin != dmin) {
      dmin = vmin;
      update = 1;
    }

    for(j = 0; j < nobj; j++) {
      if(alg->min[j] < 0) {
        v = fabs(alg->min[j]);

        for(i = 0; i < size; i++) {
          k = i * nobj + j;
          alg->norm[k] = data[k] + v;
        }

        alg->min[j] = 0;
      }
      else {
        for(i = 0; i < size; i++) {
          k = i * nobj + j;
          alg->norm[k] = data[k];
        }
      }
    }

    p = alg->norm;
  }
  else {
    p = data;
  }

  /* Update ideal point */
  for(j = 0; j < nobj; j++) {
    if(alg->min[j] < alg->ideal[j]) {
      alg->ideal[j] = alg->min[j];
      update = 1;
    }
  }

  /* Calculate the norm of each solution */
  for(i = 0; i < size; i++) {
    alg->vnorm[i] = EMO_vnorm(p + i * nobj, 2.0, nobj);
  }
 

  /* Update nadir point looking for the individuals parallel to the axis with the lowest norm.
     Removing dominance-resistant points (outliers) */
  for(j = 0; j < nobj; j++) {

    for(i = 0; i < size; i++) {
      alg->sort[i][0] = (double) i;

     /* cos theta = dot(a,b)/(|a| |b|), a: j-axis */
      alg->sort[i][1] = p[i * nobj + j] / alg->vnorm[i]; //EMO_vnorm(p + i * nobj, 2.0, nobj);
    }

    qsort(alg->sort, size, sizeof(alg->sort[0]), (int (*)(const void *, const void *))&EMO_compare_desc);

    i = 0;
    imin = -1;
    vmin = DBL_MAX;

    do {
      k = (int) alg->sort[i][0];
      v = alg->vnorm[k];

      if(v < vmin && p[k * nobj + j] > alg->ideal[j]) {
        vmin = v;
        imin = k;
      }
      i++;

    } while(i < size && (alg->sort[0][1] - alg->sort[i][1]) <= 1e-8); //1e-12);

    if(imin != -1) {
      k = imin * nobj + j;

      EMO_List_queue(&alg->lst, imin);
      v = alg->max0[j] - p[k];

      if(alg->nadir[j] != p[k]) {
        alg->nadir[j] = p[k];
        update = 1;
      }
    }
  }

  EMO_normalize(alg->norm, p, NULL, size, alg->ideal, alg->nadir, nobj);
  return update;
}

void EMO_MOVAP_alloc(EMO_MOVAP *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i, n;

  EMO_MOVAP_load_param(alg, param, mop->nvar);

  #ifdef EMO_MPI
  EMO_Population_alloc(pop, mop, param->mu, param->mu);
  #else
  EMO_Population_alloc(pop, mop, param->mu, 1);
  #endif

  n = param->mu + 1;

  EMO_NDSort_alloc(&alg->nd, n);
  EMO_VPath_alloc(&alg->vpath, alg->xrs, param->mu, n, mop->nobj);

  if((alg->norm = (double *) malloc(sizeof(double) * n * mop->nobj)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->vnorm = (double *) malloc(sizeof(double) * n)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->tmp = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in smsemoa\n");
    exit(1);
  }

  if((alg->min = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in smsemoa\n");
    exit(1);
  }

  if((alg->max0 = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in smsemoa\n");
    exit(1);
  }

  if((alg->ideal = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->nadir = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  if((alg->filter = (int *) calloc(sizeof(int), n)) == NULL) {
    printf("Error, not enough memory in smsemoa\n");
    exit(1);
   }

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
  EMO_List_alloc(&alg->lst, mop->nobj);
}

void EMO_MOVAP_free(EMO_MOVAP *alg) {
  int i;

  EMO_NDSort_free(&alg->nd);
  EMO_VPath_free(&alg->vpath);
  free(alg->norm);
  free(alg->vnorm);
  free(alg->tmp);
  free(alg->min);
  free(alg->max0);
  free(alg->ideal);
  free(alg->nadir);
  free(alg->filter);

  for(i = alg->ssize - 1; i > -1; i--)
    free(alg->sort[i]);

  free(alg->sort);
  EMO_List_free(&alg->lst);
}

void EMO_MOVAP_run(EMO_MOVAP *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int x, i, j, k, l, elem, r, n, count, count0 = 0;
  double v, vnadir;

  n = pop->mu + 1;

  /* Initialize reference points */
  for(i = 0; i < mop->nobj; i++) {
    alg->ideal[i] = DBL_MAX;
    alg->nadir[i] = DBL_MAX;
    alg->max0[i] = DBL_MAX;
  }

  EMO_MOVAP_normalize(alg, pop->obj, pop->mu, mop->nobj);
  memcpy(alg->max0, alg->nadir, sizeof(double) * mop->nobj);

  x = pop->mu * mop->nvar;
  EMO_Debug_printf(param->dbg, "Run MOVAP");

  while(!EMO_Stop_end(param->stop)) {

    /* Randomly select two parents */
    k = EMO_Rand_int1(param->rand, 0, pop->mu-1);
    l = EMO_Rand_int1(param->rand, 0, pop->mu-1);

    if(alg->vpath.c[k] < alg->vpath.c[l])
      i = k;
    else
      i = l;

    /* Randomly select two parents */
    k = EMO_Rand_int1(param->rand, 0, pop->mu-1);
    l = EMO_Rand_int1(param->rand, 0, pop->mu-1);

    if(alg->vpath.c[k] < alg->vpath.c[l])
      j = k;
    else
      j = l;

    /* Generate an offspring by variation operators */
    i *= mop->nvar;
    j *= mop->nvar;

    EMO_crossSBX(pop->var+x, pop->vdummy, pop->var+i, pop->var+j, param->rand, mop, param->Pc, param->Nc);
    EMO_mutatePolynom(pop->var+x, param->rand, mop, param->Pm, param->Nm);
    EMO_Population_evaluate(pop, mop, pop->mu, 1);

    /* Reduce population */

    /* Normalize population and discretize value-path plot */
    r = EMO_MOVAP_normalize(alg, pop->obj, n, mop->nobj);

    /* Unselect extreme solutions*/
    for(j = 0; j < n; j++)
      alg->filter[j] = 1;

    for(j = 0; j < alg->lst.size; j++) {
      EMO_List_get(&alg->lst, &elem, j);
      alg->filter[elem] = 0;
    }

    v = EMO_dmax(&j, alg->vnorm, alg->filter, n);
    vnadir = EMO_vnorm(alg->nadir, 2.0, mop->nobj);

    if(v > vnadir) {
    }
    else {
      /* Non-dominated sorting algorithm */
      EMO_NDSort_run(&alg->nd, pop->obj, mop->nobj, NULL, NULL, n);  
 
      i = alg->nd.nfront - 1;  /* Last front */

      if(alg->nd.front[i].size == 1) {
        EMO_List_get(&alg->nd.front[i], &j, 0);  /* One individual is in the last front */
      }
      else {
        count = EMO_VPath_run(&alg->vpath, alg->norm, NULL, n, 1.0);

        if(r == 1) {
          count0 = 0;
        }

        if(r != 1 && count < count0) {
          j = pop->mu;
        }
        else {
          /* Calculate individual contributions */
          EMO_VPath_contribution(&alg->vpath, alg->norm, NULL, n, 1.0, &alg->lst);

          /* Select solutions from the last front */
          memset(alg->filter, 0, sizeof(int) * n);

          for(j = 0; j < alg->nd.front[i].size; j++) {
            EMO_List_get(&alg->nd.front[i], &elem, j);
            alg->filter[elem] = 1;
          }

          /* Find worst contribution */
          v = EMO_dmax(&j, alg->vpath.c, alg->filter, n);
        }
        count0 = count;
      }
    }

    /* Copy worst individual */
    memcpy(alg->tmp, alg->norm + j * mop->nobj, sizeof(double) * mop->nobj);

    /* Remove the worst individual */
    if(j != n - 1)
      EMO_Population_copy(pop, NULL, alg->vpath.c, mop, j, n - 1);

    EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 0);

  }
}

