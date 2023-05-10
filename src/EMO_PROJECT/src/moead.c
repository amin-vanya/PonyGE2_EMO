/* **********************************************************************************
 * moead.c   Implementation of MOEA/D, based on \cite{Zhang07b,Asif10,Jan10}.       *
 *           The latter two for handling constraints.                               *
 *                                                                                  *
 * Authors: Elmer Jesús Morales Orozco        (original version \cite{Zhang07b})    *
 *          Universidad de la Sierra Sur,                                           *
 *          Miahuatlán, Oaxaca.                                                     *
 *                                                                                  *
 *          Mariano Orozco Garcia             (constraint handling \cite{Jain14})   *
 *          UPIITA-IPN                                                              *
 *          Ciudad de Mexico, Mexico                                                *
 *                                                                                  * 
 *          Raquel Hernandez Gomez            (variant from \cite{Li09c})           *
 *          CINVESTAV-IPN                                                           *
 *                                                                                  *
 * August 2016                                                                      *
 ************************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "moead.h"
#include "dominance.h" // constraint handling
#include "numeric.h"
#include "vector.h"
#include "niche.h"
#include "evop.h"
#include "io.h"

#ifdef EMO_MPI
#include "mpi.h"
#endif
 
#define _MAXCHAR 2000

void EMO_MOEAD_diff(double *v, double *x, double *min, double *max, int nobj) {
  EMO_vdiff(v, x, min, nobj);
}

void EMO_MOEAD_diff2(double *v, double *x, double *min, double *max, int nobj) {
  EMO_vdiff(v, max, x, nobj);
}

/* Load specific parameters for the algorithm */
void EMO_MOEAD_load_param(EMO_MOEAD *alg, EMO_Param *param, EMO_MOP *mop) {

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

  param->Pm = (param->Pm == -1)? 1.0 / (double) mop->nvar : param->Pm;
  EMO_Debug_printf(param->dbg, "pm updated %f", param->Pm);

  if(!EMO_Param_get_int(param, &alg->norm, "moead_norm")) {
    printf("Error, moead_norm is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_int(param, &alg->niche, "moead_niche")) {
    printf("Error, moead_niche is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_char(param, alg->wfile, "wfile")) {
    printf("Error, wfile is not defined in the configuration file.\n");
    exit(1);
  } 

  if(!EMO_Param_get_char(param, alg->utl_name, "moead_utility")) {
    strcpy(alg->utl_name, "chebyshev");
    printf("Warning, moead_utility is not defined, using %s as default.\n", alg->utl_name);
  }

  if(!EMO_Param_get_int(param, &alg->ver, "moead_ver2007")) {
    printf("Error, moead_ver2007 is not defined in the configuration file.\n");
    exit(1);
  }

  printf("version 2007 %d\n", alg->ver);

  if(alg->ver) {  // Original version of MOEA/D \cite{Zhang07b}
    alg->delta = 1.0;
    alg->nr = alg->niche;
  }
  else {  // Variant version of MOEA/D \cite{Li09c}  
    if(!EMO_Param_get_double(param, &alg->delta, "moead_delta")) {
      printf("Error, moead_delta is not defined in the configuration file.\n");
      exit(1);
    }

    if(!EMO_Param_get_int(param, &alg->nr, "moead_nr")) {
      printf("Error, moead_nr is not defined in the configuration file.\n");
      exit(1);
    }
  }
}

void EMO_MOEAD_alloc(EMO_MOEAD *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i;
  #ifdef EMO_MPI
  int myrank, nproc;
  #endif
  
  printf("MOEA/D\n");

  if((alg->wfile = (char *) malloc(sizeof(char) * _MAXCHAR)) == NULL) {
    printf("Error, not enough memory in MOEA/D.\n");
    exit(1);
  }

  if((alg->utl_name = (char *) malloc(sizeof(char) * _MAXCHAR)) == NULL) {
    printf("Error, not enough memory in MOEA/D.\n");
    exit(1);
  }

  EMO_MOEAD_load_param(alg, param, mop);

  alg->wsize = 0;

  #ifdef EMO_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  alg->W = EMO_File_read_skip(NULL, &alg->wsize, &mop->nobj, alg->wfile, myrank, nproc-1);
  #else
  alg->W = EMO_File_read(NULL, &alg->wsize, &mop->nobj, alg->wfile, 0);
  #endif

  if(alg->niche < 0 || alg->niche > alg->wsize) {
    printf("Error, niche must be less or equal than the number of weight vectors.\n");
    exit(1);
  }

  if(alg->nr < 1 || alg->nr > alg->niche) {
    printf("Error, moead_nr should be in the range [1, moead_niche].\n");
    exit(1);
  }

  #ifdef EMO_MPI
  EMO_Population_alloc(pop, mop, alg->wsize, alg->wsize);
  #else
  EMO_Population_alloc(pop, mop, alg->wsize, 1);
  #endif
 
  if((alg->min = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in moead.\n");
    exit(1);
  }

  if((alg->max = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in moead.\n");
    exit(1);
  }

  for(i = 0; i < mop->nobj; i++) {
    alg->min[i] = DBL_MAX;
    alg->max[i] = -DBL_MAX;
  }

  if((alg->diff1 = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in moead.\n");
    exit(1);
  }

  if((alg->diff2 = (double *) malloc(sizeof(double) * mop->nobj)) == NULL) {
    printf("Error, not enough memory in moead.\n");
    exit(1);
  }

  EMO_Utility_alloc(&alg->utl, param, mop->nobj, alg->utl_name);

  EMO_List_alloc(&alg->lp, alg->wsize);
  EMO_List_alloc(&alg->lt1, alg->wsize);
  EMO_List_alloc(&alg->lt2, alg->niche);

  if((alg->lnn = (EMO_List *) malloc(sizeof(EMO_List) * alg->wsize)) == NULL) {
    printf("Error, not enogh memory in moead\n");
    exit(1);
  }

  for(i = 0; i < alg->wsize; i++)
    EMO_List_alloc(&alg->lnn[i], alg->niche);

  if((alg->sort = (double **) malloc(sizeof(double *) * alg->wsize)) == NULL) {
    printf("Error, not enough memory in moead.\n");
    exit(1);
  }

  for(i = 0; i < alg->wsize; i++) {
    if((alg->sort[i] = (double *) malloc(sizeof(double) * 2)) == NULL) {
      printf("Error, not enough memory in moead.\n");
      exit(1);
    }
  }

  alg->ssize = alg->wsize;

  /* List of all individuals considered as a parents */
  for(i = 0; i < alg->wsize; i++)
    EMO_List_queue(&alg->lp, i);

  EMO_knn(alg->lnn, alg->sort, alg->W, alg->wsize, mop->nobj, alg->niche, 1);

  if(alg->utl.inverted)  /* Inverted utility function */
    alg->fdiff = &EMO_MOEAD_diff2;
  else
    alg->fdiff = &EMO_MOEAD_diff;
}

void EMO_MOEAD_free(EMO_MOEAD *alg) {
  int i;

  free(alg->wfile);
  free(alg->utl_name);
  free(alg->W);
  free(alg->min);
  free(alg->max);
  free(alg->diff1);
  free(alg->diff2);
  EMO_Utility_free(&alg->utl);
  EMO_List_free(&alg->lp);
  EMO_List_free(&alg->lt1);
  EMO_List_free(&alg->lt2);

  for(i = alg->ssize - 1; i > -1; i--) {
    EMO_List_free(&alg->lnn[i]);
    free(alg->sort[i]);
  }

  free(alg->lnn);
  free(alg->sort);
}

double EMO_MOEAD_fitness(EMO_MOEAD *alg, double *x, int nobj, int k) {
  int i;

  alg->fdiff(alg->diff1, x, alg->min, alg->max, nobj);

  if(alg->norm) {
    for(i = 0; i < nobj; i++)
      alg->diff1[i] /= (alg->diff2[i] == 0.0)? 1.0 : alg->diff2[i];
  }

  if(*(alg->utl.uf) == EMO_Utility_tlpbi || *(alg->utl.uf) == EMO_Utility_qpbi)
    EMO_Utility_update_dstar(&alg->utl, alg->min, alg->max);

  return alg->utl.uf(&alg->utl, &alg->W[k], alg->diff1); 
}

void EMO_MOEAD_run(EMO_MOEAD *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i, j, k, l, x, y, p1, p2, c, r = 0;
  EMO_List *list, aux, *tmp;
  double v1, v2;

  /* Calculate the ideal and nadir points */
  for(i = 0; i < pop->mu; i++) {

    /* Only feasible solutions are considered (if applicable)*/
    if(mop->ncon == 0 || pop->vio[i] == 0) {
      for(j = 0; j < mop->nobj; j++) {
        v1 = pop->obj[i*mop->nobj + j];

        if(v1 < alg->min[j])
          alg->min[j] = v1;

        if(v1 > alg->max[j])
          alg->max[j] = v1;
      }
    }
  }

  x = pop->mu * mop->nvar;
  y = pop->mu * mop->nobj;

  /* Calculate the constraint violation */ 
  if(mop->ncon > 0)
    EMO_Population_constrain(pop, mop, 0, pop->mu);

  while(!EMO_Stop_end(param->stop)) {
 
    for(i = 0; i < pop->mu; i++) {
      v1 = EMO_Rand_prob3(param->rand);

      /* Probability that parents are selected from the neighbourhood */
      if(v1 < alg->delta) {
        list = alg->lnn + i;
        tmp = &alg->lt2;
      }
      else {
        list = &alg->lp;
        tmp = &alg->lt1;
      }

      p1 = EMO_Rand_int2(param->rand, 0, list->size); 
      p2 = EMO_Rand_int2(param->rand, 0, list->size);
      EMO_List_get(list, &p1, p1); 
      EMO_List_get(list, &p2, p2);
 
      p1 *= mop->nvar;
      p2 *= mop->nvar;

      EMO_crossSBX(pop->var+x, pop->vdummy, pop->var + p1, pop->var + p2, param->rand, mop, param->Pc, param->Nc);
      EMO_mutatePolynom(pop->var+x, param->rand, mop, param->Pm, param->Nm);
      EMO_Population_evaluate(pop, mop, pop->mu, 1);

      if(mop->ncon > 0) {
        EMO_Population_constrain(pop, mop, pop->mu, 1);
      }

      /* Update ideal and nadir points */
      if(mop->ncon == 0 || pop->vio[pop->mu] == 0) {
        for(j = mop->nobj - 1; j > -1; j--) {
          v1 = pop->obj[y + j];

          if(v1 < alg->min[j])
            alg->min[j] = v1;

          if(v1 > alg->max[j])
            alg->max[j] = v1;
        }
      }

      /* Calculates denominator of the normalization (nadir - ideal) */
      if(alg->norm)
        EMO_vdiff(alg->diff2, alg->max , alg->min, mop->nobj);

      c = 0;

      /* Update neighboring solutions */

      while(list->size > 0 && c < alg->nr) {  // While the list is not empty
        // Randomly pick an index j from the list 
        j = EMO_Rand_int2(param->rand, 0, list->size);   

        if(!EMO_List_retrieve(list, &l, j))
          printf("Error when retrieving element in moead.c\n");

        k = l * mop->nobj;  /* index of the weight vector */

        v1 = v2 = 0.0;

        if(mop->ncon > 0) {
          r = EMO_Dominance_constraint3(pop->obj+y, pop->obj+k, mop->nobj, pop->cv[pop->mu], pop->cv[l], NULL);
          if(r == 1)  v2 = 1;
          if(r == -1) v1 = 1;
        }

        // If both solutions are feasible, then it calculates their utility functions (if applicable)
        if(mop->ncon == 0 || r == 3) {
          v1 = EMO_MOEAD_fitness(alg, pop->obj+y, mop->nobj, k); /* child */
          v2 = EMO_MOEAD_fitness(alg, pop->obj+k, mop->nobj, k); /* neighbor */
        }

        if(v1 < v2) {
          EMO_Population_copy(pop, NULL, NULL, mop, l, pop->mu); 
          c++;
        }

        EMO_List_queue(tmp, l);
      }

      /* Returns the items to the list */

      /* Less computational effort, performs swap between lists when one has more elements */
      if(tmp->size > list->size) { 
        aux = *list;
        *list = *tmp;
        *tmp = aux;
      }

      if(tmp->size > 0) {
        if(!EMO_List_move_all(tmp, list))
          printf("Error when moving all elements from a list to a list in moead.c\n");
      }
    }

    EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 0);
  }
}

#undef _MAXCHAR

