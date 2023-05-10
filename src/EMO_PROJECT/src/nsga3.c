
/* Authors: Mariano Orozco Garcia
            Raquel Hernandez Gomez
 */
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include "nsga3.h"
#include "random.h"
#include "evop.h"
#include "utility.h"
#include "numeric.h"
#include "vector.h"
#include "matrix.h"
#include "io.h"

#ifdef EMO_MPI
#include "mpi.h"
#endif

#define _MAXCHAR 2000

/* Load specific parameters for the algorithm */
void NSGA3_load_param(EMO_NSGA3 *alg, EMO_Param *param, int nvar) {
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
}


void vclose(double *x, double *vset, int *filter, int size, double *y, int n, double (*dist)(double *, double *, int)) {
  int i, j, k, l, m, flag;
  double vmin, d;

  m = 0;

  for(i = 0; i < n; i++) {
    vmin = DBL_MAX;

    l = i * n;

    for(j = 0; j < size; j++) {

      if(filter == NULL || filter[j]) {
        d = dist(vset + j*n, y + l, n);

        if(d < vmin) {
          flag = 0;

          /* Check for duplicated elements */
          for(k = 0; k < i; k++) {
            if(EMO_vdist(vset + j*n, x + k*n, n) < 1e-4) {
              flag = 1;
              break;
            }
          }

          if(!flag) {
            vmin = d;
            m = j;
          }
        }
      }
    }

    if(vmin != DBL_MAX) {
      memcpy(x + l, vset + m*n, sizeof(double) * n);

      /* Non negative entries */
      for(k = 0; k < n; k++) {
        if(x[l+k] < 0)
          x[l+k] = 0;
      }
    }
  }
}

void vclose2(EMO_NSGA3 *alg, int size, int n) {
  int i, j, k, l, m = 0, flag; //o
  double vmin, d;

  for(i = 0; i < n; i++) {
    vmin = DBL_MAX;

    l = i * n;

    for(j = 0; j < size; j++) {

      if(alg->nd.rank[j] == 0 && (alg->cv == NULL || alg->cv[j] == 0.0)) {
        EMO_vorth(alg->a, alg->norm + j * n, alg->waxis + l, n);
        d = EMO_vnorm(alg->a, 2.0, n);

        if(d < vmin) {

          flag = 0;

          /* Check for duplicated or similar elements */
          for(k = 0; k < i; k++) {
            if(EMO_vdist(alg->norm + j * n, alg->xtrm + k * n, n) < 1e-4) { //== 0  //< 1e-4 
              flag = 1;
              break;
            }
          }
          
          // Near to zero vector 
          if(EMO_vnorm(alg->norm + j * n, 2.0, n) < 1e-4) {
            flag = 1;
            break;
          }

          if(!flag) {
            vmin = d;
            m = j;
          }
        }
      }
    }

    if(vmin == DBL_MAX) {
      memset(alg->xtrm + l, 0, sizeof(double) * n);
      alg->xtrm[l + i] = 1.0;
    }
    else {
      memcpy(alg->xtrm + l, alg->norm + m*n, sizeof(double) * n);
    }
  }
}

/* Normalization */
void EMO_NSGA3_normalize(EMO_NSGA3 *alg, EMO_Population *pop, int nobj) {
  int i, j, k;

  /* Calculate the ideal point */
  for(i = 0; i < pop->size; i++) {
    for(j = 0; j < nobj; j++) {

      /* Only feasible solutions are considered */
      if(alg->filter[i] && (alg->cv == NULL || alg->cv[i] == 0.0) && pop->obj[i*nobj + j] < alg->min[j])
        alg->min[j] = pop->obj[i*nobj + j];
    }
  }

  /* Translate objectives */
  for(i = 0; i < pop->size; i++) {
    if(alg->filter[i] && (alg->cv == NULL || alg->cv[i] == 0.0)) {
      for(j = 0; j < nobj; j++) {
        k = i * nobj + j;
        alg->norm[k] = pop->obj[k] - alg->min[j];
      }
    }
  }

  /* Find extreme points */
  vclose2(alg, pop->size, nobj);

  if(EMO_minverse(alg->inv, alg->xtrm, nobj, 1, 0, 0) == 0) {
    EMO_matmul(alg->a, alg->inv, alg->one, nobj, nobj, 1);

    for(i = 0; i < pop->size; i++) {
      if(alg->filter[i] && (alg->cv == NULL || alg->cv[i] == 0.0)) {

        for(j = 0; j < nobj; j++)
          alg->norm[i*nobj + j] *= fabs(alg->a[j]);
      }
    }
  }
}

// Associate each individual w.r.t. a reference point, it also stores the distance to it.
void EMO_NSGA3_associate(EMO_NSGA3 *alg, EMO_Population *pop, int nobj) {
  int i, j, idx = 0;
  double d, vmin;

  for(i = 0; i < pop->size; i++) {
    vmin = DBL_MAX;

    if(alg->filter[i] && (alg->cv == NULL || alg->cv[i] == 0.0)) {

      for(j = 0; j < alg->wsize; j++) {
        EMO_vorth(alg->a, alg->norm + i * nobj, alg->W + j * nobj, nobj);
        d = EMO_vnorm(alg->a, 2.0, nobj);

        if(d < vmin) {
          vmin = d;
          idx = j;
        }
      }

      alg->pi[i] = idx;
      alg->dist[i] = vmin;
    }
  }
}

/* Selects n individuals */
void EMO_NSGA3_niching(EMO_NSGA3 *alg, EMO_Param *param, EMO_Population *pop, int nobj, int n, int lastfront) {
  int i, j, k, w, idx = 0;
  double vmin;

  memset(alg->niche, 0, sizeof(int) * alg->wsize);

  for(i = 0; i < pop->size; i++) 
    if(alg->nd.rank[i] < lastfront)
      alg->niche[alg->pi[i]]++;

  EMO_List_clear(&alg->lst1);

  i = 0;

  while(i < n) {

    /* Look for the reference points that are isolated by the first fronts */
    if(alg->lst1.size == 0) {
      vmin = EMO_min(NULL, alg->niche, NULL, alg->wsize);

      for(j = 0; j < alg->wsize; j++)
        if(alg->niche[j] == vmin)
          EMO_List_queue(&alg->lst1, j);
    }

    /* Selects a random reference point */
    EMO_List_get(&alg->lst1, &w, EMO_Rand_int1(param->rand, 0, alg->lst1.size-1));
    EMO_List_clear(&alg->lst2);

    /* Look for the individuals in the last front that are near to the isolated reference points */
    for(j = 0; j < alg->nd.front[lastfront].size; j++) {
      EMO_List_get(&alg->nd.front[lastfront], &k, j);
      if(alg->filter[k] == 0 && alg->pi[k] == w)
        EMO_List_queue(&alg->lst2, k);
    }

    if(alg->lst2.size == 0) {  /* If there is no individual, the reference point is removed */
      alg->niche[w] = INT_MAX;
    }
    else {
      if(alg->niche[w] == 0) {   /* There is no individual in the first fronts that are associated */ 
        vmin = DBL_MAX;         /* with the reference point */

        for(j = 0; j < alg->lst2.size; j++) {
          EMO_List_get(&alg->lst2, &k, j);

          if(alg->dist[k] < vmin) {  /* Select the individual with the shortest distance */
            vmin = alg->dist[k];
            idx = k;
          }
        }
      }
      else {
        EMO_List_get(&alg->lst2, &idx, EMO_Rand_int1(param->rand, 0, alg->lst2.size-1));  /* Select a random individual */
      }

      alg->filter[idx] = 1;
      alg->niche[w]++;
      i++;
    }

    if(!EMO_List_remove(&alg->lst1, w)) {
      printf("Error al eliminar en lista.\n");
      exit(1);
    }
  }

  EMO_List_clear(&alg->lst1);
  EMO_List_clear(&alg->lst2);
}

// Constraint violation value (CV) pop[start, end)
// It returns the number of feasible solutions.
void EMO_NSGA3_calculate_cv(EMO_NSGA3 *alg, EMO_Population *pop, int start, int end, EMO_MOP *mop) {
  double *g;
  int i, j;

  if(mop->ncon > 0) {
    for(i = start; i < end; i++) {
      alg->cv[i] = 0.0;
      g = pop->con + i * mop->ncon;

      for(j = 0; j < mop->ncon; j++)
        if(g[j] < 0.0)
          alg->cv[i] += -g[j]; /* Bracket operator */
    }
  }
}

int EMO_NSGA3_count_feasible(EMO_NSGA3 *alg, int size) {
  int i, n = 0;

  if(alg->cv != NULL) {
    for(i = 0; i < size; i++)
      if(alg->cv[i] == 0.0)
        n++;
  }

  return n;
}

int EMO_tournament_selection(EMO_Rand *rnd, double *cv, int size) {
  int p1, p2;

  p1 = EMO_Rand_int2(rnd, 0, size);
  p2 = EMO_Rand_int2(rnd, 0, size);
  
  if(cv != NULL) {
    if(cv[p1] == 0.0) {
      if(cv[p2] != 0.0)
        return p1;
    }
    else {
      if(cv[p2] == 0.0)
        return p2;
      else if(cv[p1] < cv[p2])
        return p1;
      else if(cv[p1] > cv[p2])
        return p2;
    }
  }

  if(EMO_Rand_flip(rnd, 0.5))
    return p1;

  return p2;
}


void EMO_NSGA3_alloc(EMO_NSGA3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i, n;
  #ifdef EMO_MPI
  int myrank, nproc;
  #endif

  if((alg->wfile = (char *) malloc(sizeof(char) * _MAXCHAR)) == NULL) {
    printf("Error, not enough memory in NSGA3.\n");
    exit(1);
  }

  NSGA3_load_param(alg, param, mop->nvar);

  if(param->mu % 2 != 0) {
    printf("Warning, adjusting population size to be an even number (%d vs %d).\n", param->mu, param->mu+1);
    param->mu++;
  }

  EMO_Population_alloc(pop, mop, param->mu, param->mu);
  EMO_NDSort_alloc(&alg->nd, pop->size);
  EMO_List_alloc(&alg->lst1, pop->mu);
  EMO_List_alloc(&alg->lst2, pop->mu);

  if((alg->norm = (double *) calloc(sizeof(double), pop->size * mop->nobj)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->a = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->one = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  for(i = 0; i < mop->nobj; i++)
    alg->one[i] = 1.0;

  n = mop->nobj * mop->nobj;

  if((alg->axis = (double *) malloc(sizeof(double) * n)) == NULL) {  // identity
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  EMO_vaxes(alg->axis, mop->nobj);

  alg->wsize = 0;

  #ifdef EMO_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  alg->W = EMO_File_read_skip(NULL, &alg->wsize, &mop->nobj, alg->wfile, myrank, nproc-1);
  #else
  alg->W = EMO_File_read(NULL, &alg->wsize, &mop->nobj, alg->wfile, 0);
  #endif


  if((alg->waxis = (double *) malloc(sizeof(double) * n)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  vclose(alg->waxis, alg->W, NULL, alg->wsize, alg->axis, mop->nobj, EMO_vdist);

  if((alg->niche = (int *) malloc(sizeof(int) * alg->wsize)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->xtrm = (double *) malloc(sizeof(double) * n)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->inv = (double *) malloc(sizeof(double) * n)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->min = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  for(i = 0; i < mop->nobj; i++)
    alg->min[i] = DBL_MAX;

  if((alg->pi = (int *) calloc(sizeof(int), pop->size)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->dist = (double *) calloc(sizeof(double), pop->size)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  if((alg->filter = (int *) calloc(sizeof(int), pop->size)) == NULL) {
    printf("Error, not enough memory in nsga3.\n");
    exit(1);
  }

  alg->cv = NULL;

  if(mop->ncon > 0) {
    if((alg->cv = (double *) calloc(sizeof(double), pop->size)) == NULL) {
      printf("Error, not enough memory in nsga3.\n");
      exit(1);
    }
  }
}

void EMO_NSGA3_free(EMO_NSGA3 *alg) {
  free(alg->wfile);
  EMO_NDSort_free(&alg->nd);
  EMO_List_free(&alg->lst1);
  EMO_List_free(&alg->lst2);
  free(alg->norm);
  free(alg->a);
  free(alg->one);
  free(alg->axis);
  free(alg->W);
  free(alg->waxis);
  free(alg->niche);
  free(alg->xtrm);
  free(alg->inv);
  free(alg->min);
  free(alg->pi);
  free(alg->dist);
  free(alg->filter);

  if(alg->cv != NULL) free(alg->cv);
}

void EMO_NSGA3_run(EMO_NSGA3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int  p1, p2, v1, v2, i, j, cont, lastfront;

  EMO_NSGA3_calculate_cv(alg, pop, 0, pop->mu, mop);

  while(!EMO_Stop_end(param->stop)) {

    /* Offspring generation */
    for(i = 0; i < pop->lambda; i+=2) {

      if(mop->ncon > 0) {
        p1 = EMO_tournament_selection(param->rand, alg->cv, pop->mu);
        p2 = EMO_tournament_selection(param->rand, alg->cv, pop->mu);
      }
      else {
      /* Randomly select two parents */
        p1 = EMO_Rand_int1(param->rand, 0, pop->mu-1);
        while ((p2 = EMO_Rand_int1(param->rand, 0, pop->mu-1)) == p1);
      }

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

    if(mop->ncon > 0) {
      EMO_NSGA3_calculate_cv(alg, pop, pop->mu, pop->size, mop);
    }

    /* Reduce population */
    EMO_NDSort_run(&alg->nd, pop->obj, mop->nobj, alg->cv, NULL, pop->size);  /* Non-dominated sorting algorithm */

    j = cont = 0;

    do {
      cont += alg->nd.front[j++].size;
    } while(cont < pop->mu);

    lastfront = j - 1;

    /* Select all solutions to the last front */
    memset(alg->filter, 0, sizeof(int) * pop->size);

    for(i = 0; i < pop->size; i++) {
      if(alg->nd.rank[i] <= lastfront) {
        alg->filter[i] = 1;
      }
    }

    if(cont != pop->mu) {  // there are more individuals than the required
      if(mop->ncon > 0) {
        // If an individual in the last front is unfeasible, then all 
        // individuals in the last front have the same constraint violation. 
        EMO_List_get(&alg->nd.front[lastfront], &j, 0);  

        if(alg->cv[j] > 0.0) {
          i = cont - pop->mu;  // Number of individuals to be removed

          // Select random individuals
          while(i > 0) {
            j = EMO_Rand_int2(param->rand, 0, alg->nd.front[lastfront].size);

            if(EMO_List_retrieve(&alg->nd.front[lastfront], &j, j) == 0) {
              printf("Error in retrieve (size %d, index %d), lastfront %d.\n", alg->nd.front[lastfront].size, j, lastfront);
              exit(1);
            }
 
            alg->filter[j] = 0;
            i--;
          }
        }
      }

      if(alg->cv == NULL || alg->cv[j] == 0.0) {
        EMO_NSGA3_normalize(alg, pop, mop->nobj);
        EMO_NSGA3_associate(alg, pop, mop->nobj);

        // unselect the last front
        for(i = 0; i < alg->nd.front[lastfront].size; i++) {
          EMO_List_get(&alg->nd.front[lastfront], &j, i);
          alg->filter[j] = 0;
        }

        // Number of individuals to be selected: mu - (cont - size)
        EMO_NSGA3_niching(alg, param, pop, mop->nobj, pop->mu - cont + alg->nd.front[lastfront].size, lastfront);
      }
    }

    EMO_Population_survive(pop, NULL, alg->cv, mop, &alg->lst1, &alg->lst2, alg->filter); 
    EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 0);
  }
}

#undef _MAXCHAR

