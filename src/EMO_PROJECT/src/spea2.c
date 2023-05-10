
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "spea2.h"
#include "dominance.h"
#include "numeric.h"
#include "vector.h"
#include "niche.h"
#include "evop.h"
#include "sort.h"

/* Load specific parameters for the algorithm */
void SPEA2_load_param(EMO_SPEA2 *alg, EMO_Param *param, int nvar) {
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

void EMO_SPEA2_alloc(EMO_SPEA2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i;
  
  SPEA2_load_param(alg, param, mop->nvar);

  EMO_Population_alloc(pop, mop, param->mu , param->mu);

  EMO_List_alloc(&alg->lst1, pop->mu);
  EMO_List_alloc(&alg->lst2, pop->mu);

  if((alg->parent1 = (int *) malloc(sizeof(int) * pop->mu )) == NULL) {
    printf("Error, not enough memory in SPEA2.\n");
    exit(1);
  }

  if((alg->parent2 = (int *) malloc(sizeof(int) * pop->mu )) == NULL) {
    printf("Error, not enough memory in SPEA2.\n");
    exit(1);
  }
 
  if((alg->fit = (double *) malloc(sizeof(double) * pop->size)) == NULL) {
    printf("Error, not enough memory in SPEA2.\n");
    exit(1);
  }

  if((alg->S = (double *) malloc(sizeof(double) * pop->size)) == NULL) {
    printf("Error, not enough memory in SPEA2.\n");
    exit(1);
  }

  if((alg->F = (double *) malloc(sizeof(double) * pop->size)) == NULL) {
    printf("Error, not enough memory in SPEA2.\n");
    exit(1);
  }


  if((alg->lnn = (EMO_List *) malloc(sizeof(EMO_List) * pop->size)) == NULL) {
    printf("Error, not enogh memory in SPEA2.\n");
    exit(1);
  }

  for(i = 0; i < pop->size; i++)
    EMO_List_alloc(&alg->lnn[i], pop->size);

  EMO_List_alloc(&alg->copy, pop->size);

  if((alg->sort = (double **) malloc(sizeof(double *) * pop->size)) == NULL) {
    printf("Error, not enough memory in SPEA2.\n");
    exit(1);
  }

  for(i = 0; i < pop->size; i++) {
    if((alg->sort[i] = (double *) malloc(sizeof(double) * 2)) == NULL) {
      printf("Error, not enough memory in SPEA2.\n");
      exit(1);
    }
  }
  
  if((alg->dist = (double **) malloc(sizeof(double *) * pop->size)) == NULL) {
    printf("Error, not enough memory in SPEA2.\n");
    exit(1);
  }

  for(i = 0; i < pop->size; i++) {
    if((alg->dist[i] = (double *) malloc(sizeof(double) * pop->size)) == NULL) {
      printf("Error, not enough memory in SPEA2.\n");
      exit(1);
    }
  }

  if((alg->filter = (int *) malloc(sizeof(int) * pop->size )) == NULL) {
    printf("Error, not enough memory in SPEA2.\n");
    exit(1);
  }
 
  if((alg->seedp = (int *) malloc(sizeof(int) * pop->mu)) == NULL) {
    printf("Error, not enough memory in SPEA2.\n");
    exit(1);
  }

  for(i = 0; i < pop->mu; i++)
    alg->seedp[i] = i;

  alg->ssize = pop->size;
  alg->k = 0;
}

void EMO_SPEA2_free(EMO_SPEA2 *alg) {
  int i;

  for(i = alg->ssize - 1; i > -1; i--) {
   EMO_List_free(&alg->lnn[i]);
   free(alg->dist[i]);
   free(alg->sort[i]);
  }
  
  free(alg->lnn);
  free(alg->dist);
  free(alg->sort);
  free(alg->fit);
  free(alg->parent1);
  free(alg->parent2);  
  free(alg->S); 
  free(alg->F);
  free(alg->filter);
  free(alg->seedp);
  EMO_List_free(&alg->lst1);
  EMO_List_free(&alg->lst2); 
  EMO_List_free(&alg->copy);
}

int SPEA2_fitness_asignment(EMO_SPEA2 *alg, EMO_Population *pop, EMO_MOP *mop, int size) {
  int h, i, j, cont, elem;

  cont = 0;

  for(i=0; i < size; i++) {
    alg->S[i] = 0;
    alg->F[i] = 0;
    EMO_List_clear(&alg->lnn[i]);
  }

  for(i=0; i < size; i++){ 
    for(j=i+1; j < size; j++){
      h = EMO_Dominance_strict(pop->obj + i * mop->nobj,  pop->obj + j * mop->nobj, mop->nobj);
      if(h == 1) { 
        alg->S[i]++;
        EMO_List_queue(&alg->lnn[j] , i);
      }
      else if(h == -1) {
        alg->S[j]++;
        EMO_List_queue(&alg->lnn[i] , j);
      }
    }
  }
 
  for(i=0; i < size; i++) { 
    for(j=0; j < alg->lnn[i].size; j++) {
      EMO_List_get(&alg->lnn[i], &elem, j);
      alg->F[i] += alg->S[elem];
    }
    if(alg->F[i] == 0) cont++; 
  }

  EMO_knn2(alg->lnn, &alg->copy, alg->dist, alg->sort, pop->obj, size, mop->nobj);
  
  for(i=0; i < size; i++){
    EMO_List_get(&alg->lnn[i], &j, alg->k);
    alg->S[i] = 1.0 / (2.0 + alg->dist[i][j]);
  }

  for(i=0; i < size; i++){
    alg->F[i] = alg->F[i] + alg->S[i]; 
  }
 
  return cont;
}

/* Retrieves the k-nearest neighbor (elem) of individual p that is activated by filter */
int SPEA2_getNN(EMO_List *lnn, int *filter, int *elem, int *k, int p) {

  if(filter == NULL) { 
    while(*k < lnn[p].size) {
      EMO_List_get(&lnn[p], elem, *k); 
      (*k)++;
    }
  }
  else {
    while(*k < lnn[p].size) {
      EMO_List_get(&lnn[p], elem, *k); 
      if(filter[*elem]) return 0;
      (*k)++;
    }
  }
  *elem = -1;
  return 1;
}

int SPEA2_compare(EMO_List *lnn, double **dist, int *filter, int p1, int p2, int size) {
  int i, j, e1, e2;

  i = j = 0;

  while(i < size && j < size) {

    if(SPEA2_getNN(lnn, filter, &e1, &i, p1)) {
      return 1;
    }

    if(SPEA2_getNN(lnn, filter, &e2, &j, p2)) {
      return -1;
    }

    if(dist[p1][e1] < dist[p2][e2]) {
      return 1;
    }
    else if(dist[p1][e1] > dist[p2][e2]) {
      return -1; 
    }

    i++;
    j++;
  }

  return 0;
}

void tournament_selection_fitness3(int *ind, EMO_Rand *rand, int *seedp, double *fit, double **dist, int size, EMO_List *l,  EMO_Population *pop, EMO_MOP *mop){
 int p1, p2, r, i, j = 0;

  EMO_Rand_shuffle(rand, seedp, size);

  for(i = 0; i < size; i+=2) {
    p1 = seedp[i];
    p2 = seedp[i+1];
    
    if(fit[p1] < fit[p2]) {
      ind[j++] = p1;
    }
    else if(fit[p2] < fit[p1]) {
      ind[j++] = p2;
    }
    else {
       r = SPEA2_compare(l, dist, NULL, p1, p2, size);

       if(r == 1) {
          ind[j++] = p2;
       }
       else if(r == -1) {
          ind[j++] = p1;
       }
       else {
         if (EMO_Rand_flip(rand, 0.5))
           ind[j++] = p1;
         else
           ind[j++] = p2;
       }
    }
  }
}

int compare_vector(const void *p1, const void *p2, int start, int dim) {
  int i;

  double *x1 = (double *)p1;
  double *x2 = (double *)p2;

  for(i = start; i < dim; i++) {
    if(x1[i] < x2[i])
      return 1;
    if(x1[i] > x2[i])
      return -1;
  }

  return 0;
}


void EMO_SPEA2_run(EMO_SPEA2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int i, j, n, cont, p1, p2, v1, v2, gen=0;  //o1, o2,
  
  SPEA2_fitness_asignment( alg, pop, mop, pop->mu);
  
  while(!EMO_Stop_end(param->stop)) {
    gen++;
   
    tournament_selection_fitness3(alg->parent1, param->rand, alg->seedp, alg->F, alg->dist, pop->mu, alg->lnn, pop, mop);
    tournament_selection_fitness3(alg->parent2, param->rand, alg->seedp, alg->F, alg->dist, pop->mu, alg->lnn, pop, mop);
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

      /* Generate an offspring by variation operators */
      EMO_crossSBX(pop->var+v1, pop->var+v2, pop->var+p1, pop->var+p2, param->rand, mop, param->Pc, param->Nc);
      EMO_mutatePolynom(pop->var+v1, param->rand, mop, param->Pm, param->Nm);
      EMO_mutatePolynom(pop->var+v2, param->rand, mop, param->Pm, param->Nm);
      EMO_Population_evaluate(pop, mop, pop->mu + i, 2);
    } 
    
    cont = SPEA2_fitness_asignment(alg, pop, mop, pop->size);
   
    /* Select solutions from the first fronts */
    memset(alg->filter, 0, sizeof(int) * pop->size); 

    n = pop->size - (pop->mu - cont);

    if(cont < pop->mu) {   /* Fills the population */

      /* Selects all dominated individuals */
      for(i = 0; i < pop->size; i++)
        if(alg->F[i] >= 1.0)
          alg->filter[i] = 1;

      while(cont < n)  {
        //Find the individual with the biggest fitness (worst individual) 
        EMO_dmax(&i, alg->F, alg->filter, pop->size);
        alg->filter[i] = 0;
        cont++;
      }

      /* Selects all non-dominated individuals */
      for(i = 0; i < pop->size; i++) {
        if(alg->F[i] < 1.0)
          alg->filter[i] = 1; 
      }
    }
    else {

      /* Selects all non-dominated individuals */
      for(i = 0; i < pop->size; i++) {
        if(alg->F[i] < 1.0)
          alg->filter[i] = 1; 
      }

      if(cont > pop->mu) {
        /* Eliminates copies */ 
        while(cont > pop->mu && alg->copy.size > 0) {
          EMO_List_dequeue(&alg->copy, &j);

          if(alg->filter[j]) {
            alg->filter[j] = 0;
            cont--;
          }
        }

        /* Finds individuals with the worst distribution */
        while(cont > pop->mu) {
          i = 0;
          while(i < pop->size && alg->filter[i] == 0) i++;

          for(j=i+1; j < pop->size; j++) {
            if(alg->filter[j] == 0) continue;

            if(SPEA2_compare(alg->lnn, alg->dist, alg->filter, j, i, pop->size) == 1)
              i = j;
          }

          alg->filter[i] = 0;
          cont--;
        }
      }
    }

    EMO_Population_survive(pop, NULL, alg->F, mop, &alg->lst1, &alg->lst2, alg->filter);
    EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 0);
  }
}

