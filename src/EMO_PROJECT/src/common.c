
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "common.h"
#include "benchmark.h"
#include "io.h"

#define _MAXCHAR 2000
#define _MAXPSIZE 1000000

/* 
 * from_file: variable population
 */
void EMO_Population_alloc(EMO_Population *pop, EMO_MOP *mop, int mu, int lambda) {

  pop->size = mu + lambda;
  pop->mu = mu;
  pop->lambda = lambda;
  pop->asize = 0;

  if(pop->size < 0 || pop->size > _MAXPSIZE) {
    printf("Error, wrong value for psize %d in common.c:EMO_Population_alloc.\n", pop->size);
    exit(1);
  }

  if((pop->obj = (double *) calloc(sizeof(double), mop->nobj * pop->size)) == NULL) {
    printf("Error, not enough memory in EMO_Population_alloc 1.\n");
    exit(1);
  }

  if((pop->odummy = (double *) malloc(sizeof(double) * mop->nobj)) == NULL){
    printf("Error, not enough memory in EMO_Population_alloc 2.\n");
    exit(1);
  }

  pop->con = NULL;
  pop->vio = NULL;
  pop->cv = NULL;

  if(mop->ncon > 0) {
    if((pop->con = (double *) calloc(sizeof(double), mop->ncon * pop->size)) == NULL){
      printf("Error, not enough memory in EMO_Population_alloc 3.\n");
      exit(1);
    }

    if((pop->vio = (int *) calloc(sizeof(int), pop->size)) == NULL) {
      printf("Error, not enough memory in EMO_Population_alloc 5.\n");
      exit(1);
    }

    if((pop->cv = (double *) calloc(sizeof(double), pop->size)) == NULL) {
      printf("Error, not enough memory in EMO_Population_alloc 4.\n");
      exit(1);
    }

    if((pop->cdummy = (double *) malloc(sizeof(double) * mop->ncon)) == NULL){
      printf("Error, not enough memory in EMO_Population_alloc 6.\n");
      exit(1);
    }
  }

  if((pop->var = (double *) calloc(sizeof(double), mop->nvar * pop->size)) == NULL) {
    printf("Error, not enough memory in EMO_Population_alloc 7.\n");
    exit(1);
  }

  if((pop->vdummy = (double *) malloc(sizeof(double) * mop->nvar)) == NULL){
    printf("Error, not enough memory in EMO_Population_alloc 8.\n");
    exit(1);
  }

  if((pop->format = (char *) malloc(sizeof(char) *_MAXCHAR)) == NULL){
    printf("Error, not enough memory in EMO_Population_alloc 9.\n");
    exit(1);
  }

  strcpy(pop->format, "%e "); 
}

void EMO_Population_init(EMO_Population *pop, EMO_Rand *rand, EMO_MOP *mop) {
  int i, j, k;

  /* Random initialization */
  printf("Random inicialization of population.\n");

  printf("mu %d, var %d, obj %d\n", pop->mu, mop->nvar, mop->nobj);

  if(mop->coding == EMO_REAL) { printf("Codificacion real\n");
    for(i = 0; i < pop->mu; i++) {
      k = i * mop->nvar;

      for(j = 0; j < mop->nvar; j++)
        pop->var[k + j] = EMO_Rand_real1(rand, mop->xmin[j], mop->xmax[j]);

      EMO_Population_evaluate(pop, mop, i, 1);
    }
  }

  if(mop->coding == EMO_BINARY) { printf("Codificacion binaria\n");
    for(i = 0; i < pop->mu; i++) {
      k = i * mop->nvar;

      for(j = 0; j < mop->nvar; j++)
        pop->var[k + j] = EMO_Rand_int1(rand, (int) mop->xmin[j], (int) mop->xmax[j]);

      EMO_Population_evaluate(pop, mop, i, 1);
    }
  }
}

/* Initialization from .pos and .pof files */
int EMO_Population_init_from_file(EMO_Population *pop, EMO_MOP *mop, const char *prefix, int start) {
  int size1, size2, i, j;
  char *str;

  str = (char *) malloc(sizeof(char) * (strlen(prefix) + 5));

  sprintf(str, "%s.pos", prefix);
  printf("Reading population from file %s.\n", str);
  size1 = pop->mu;
  EMO_File_read(pop->var, &size1, &mop->nvar, str, start);
  printf("var: mu %d, size1 read %d, var %d, obj %d\n", pop->mu, size1, mop->nvar, mop->nobj);

  sprintf(str, "%s.pof", prefix);
  printf("Reading population from file %s.\n", str);
  size2 = pop->mu;
  EMO_File_read(pop->obj, &size2, &mop->nobj, str, start);
  printf("obj: mu %d, size2 read %d, var %d, obj %d\n", pop->mu, size2, mop->nvar, mop->nobj);

  if(size1 != size2) {
    printf("Error, files %s.pos and %s.pof contain different number of records (%d vs %d)\n", prefix, prefix, size1, size2);
    exit(1);
  }

  if(mop->ncon > 0) {
    sprintf(str, "%s.con", prefix);
    printf("Reading population from file %s.\n", str);
    size2 = pop->mu;
    EMO_File_read(pop->con, &size2, &mop->ncon, str, start);
    printf("obj: mu %d, size2 read %d, var %d, obj %d, ncon %d\n", pop->mu, size2, mop->nvar, mop->nobj, mop->ncon);

    if(size1 != size2) {
      printf("Error, files %s.con and %s.pof contain different number of records (%d vs %d)\n", prefix, prefix, size1, size2);
      exit(1);
    }

    /* Check violated constraints */
    for(i = 0; i < size2; i++) {
      pop->vio[start + i] = 0;

      for(j = 0; j < mop->ncon; j++)
        if(pop->con[(start + i) * mop->ncon + j] < 0.0)
          pop->vio[start + i]++;
    }
  }

  free(str);

  return size1;
}

void EMO_Population_write(EMO_Population *pop, int *filter, EMO_MOP *mop, const char *prefix, int run, unsigned long itera) {
  char file[_MAXCHAR];
  int flag;

  flag = strcmp(prefix, "stdout");

  if(flag == 0)
    strcpy(file, "stdout");
  else if(run == 0)
    sprintf(file, "%s.pos", prefix);
  else
    sprintf(file, "%s%.2d.pos", prefix, run);

  EMO_File_write(pop->var, filter, pop->mu, mop->nvar, file, pop->format, itera);

  if(flag != 0) {
    if(run == 0)
      sprintf(file, "%s.pof", prefix);
    else
      sprintf(file, "%s%.2d.pof", prefix, run);
  }

  EMO_File_write(pop->obj, filter, pop->mu, mop->nobj, file, pop->format, itera);

  if(mop->ncon > 0) {
    if(flag != 0) {
      if(run == 0)
        sprintf(file, "%s.con", prefix);
      else
        sprintf(file, "%s%.2d.con", prefix, run);
    }

    EMO_File_write(pop->con, filter, pop->mu, mop->ncon, file, pop->format, itera);
  }
}

void EMO_Population_free(EMO_Population *pop) {
  free(pop->var);
  free(pop->vdummy);
  free(pop->obj);
  free(pop->odummy);
  free(pop->format);

  if(pop->con != NULL) {
    free(pop->con);
    free(pop->vio);
    free(pop->cv);
    free(pop->cdummy);
  }
}

/* Copy individual j to the current position of individual i */
void EMO_Population_copy(EMO_Population *pop, int *fit1, double *fit2, EMO_MOP *mop, int i, int j) {  //destino i, origen j
   int x, y;

  if(i == j) return;

  // Copy variables
  x = i * mop->nvar;
  y = j * mop->nvar;
  memcpy(pop->var + x, pop->var + y, sizeof(double) * mop->nvar);

  // Copy objectives
  x = i * mop->nobj;
  y = j * mop->nobj;
  memcpy(pop->obj + x, pop->obj + y, sizeof(double) * mop->nobj);

  // Copy constraints
  if(mop->ncon > 0) {
    x = i * mop->ncon;
    y = j * mop->ncon;
    memcpy(pop->con + x, pop->con + y, sizeof(double) * mop->ncon);
    pop->vio[i] = pop->vio[j];
    pop->cv[i] = pop->cv[j];
  }
  
  if(fit1 != NULL) 
    fit1[i] = fit1[j];

  if(fit2 != NULL) 
    fit2[i] = fit2[j];
}

// copia bloques de memoria entre poblaciones
// start y size especifican a partir de donde y cuantos individuos 
// se copian de src a dest, en la posicion asize
void EMO_Population_copy2(EMO_Population *dest, EMO_Population *src, EMO_MOP *mop, int start, int size) {
  int x, y; 

  if(dest->asize + size < dest->size) {
    // Copy variables
    x = dest->asize * mop->nvar;
    y = start * mop->nvar; 
    memcpy(dest->var + x, src->var + y, sizeof(double) * mop->nvar * size);

    // Copy objectives
    x = dest->asize * mop->nobj;
    y = start * mop->nobj;
    memcpy(dest->obj + x, src->obj + y, sizeof(double) * mop->nobj * size);

    // Copy constraints
    if(mop->ncon > 0) {
      x = dest->asize * mop->ncon;
      y = start * mop->ncon;
      memcpy(dest->con + x, src->con + y, sizeof(double) * mop->ncon * size);
      memcpy(dest->vio + dest->asize, src->vio + start, sizeof(int) * size);
      memcpy(dest->cv + dest->asize, src->cv + start, sizeof(double) * size);
    }

    dest->asize += size;
  }
  else {
    printf("Error, not enough memory in EMO_Population_copy2 (%d vs %d)\n", (dest->asize + size), dest->size);
    exit(1);
  }
}

/* Interchange individual i by j in the population */
void EMO_Population_swap(EMO_Population *pop, int *fit1, double *fit2, EMO_MOP *mop, int i, int j) {
  double dummy;
  int x, y;

  if(i == j) return;

  // Copy variables
  x = i * mop->nvar;
  y = j * mop->nvar;

  memcpy(pop->vdummy,  pop->var + x, sizeof(double) * mop->nvar);
  memcpy(pop->var + x, pop->var + y, sizeof(double) * mop->nvar);
  memcpy(pop->var + y, pop->vdummy,  sizeof(double) * mop->nvar);

  // Copy objectives
  x = i * mop->nobj;
  y = j * mop->nobj;

  memcpy(pop->odummy,  pop->obj + x, sizeof(double) * mop->nobj);
  memcpy(pop->obj + x, pop->obj + y, sizeof(double) * mop->nobj);
  memcpy(pop->obj + y, pop->odummy,  sizeof(double) * mop->nobj);

  // Copy constraints
  if(mop->ncon > 0) {
    x = i * mop->ncon;
    y = j * mop->ncon;

    memcpy(pop->cdummy,  pop->con + x, sizeof(double) * mop->ncon);
    memcpy(pop->con + x, pop->con + y, sizeof(double) * mop->ncon);
    memcpy(pop->con + y, pop->cdummy,  sizeof(double) * mop->ncon);
    x = pop->vio[i];
    pop->vio[i] = pop->vio[j];
    pop->vio[j] = x;

    dummy = pop->cv[i];
    pop->cv[i] = pop->cv[j];
    pop->cv[j] = dummy;
  }

  if(fit1 != NULL) {
    x = fit1[i];
    fit1[i] = fit1[j];
    fit1[j] = x;
  }

  if(fit2 != NULL) {
    dummy = fit2[i];
    fit2[i] = fit2[j];
    fit2[j] = dummy;
  }
}

/* Keeps the selected mu individuals, indicated by filter, at the first elements of the mixed population 

   if filter == NULL:
   Replace the list of available positions by the list of missing ones for a mixed population.
   The positions of all missing elements must be in the range [mu,size) and the positions of 
   all members of available must be in [0,mu).
*/
void EMO_Population_survive(EMO_Population *pop, int *fit1, double *fit2, EMO_MOP *mop, EMO_List *missing, EMO_List *available, int *filter) {
  int i, j, n;

  if(filter != NULL) {
    EMO_List_clear(missing);
    EMO_List_clear(available);

    // Check if it is a secondary population (archive)
    n = (pop->asize == 0)? pop->size : pop->asize; 

    for(i = 0; i < n; i++) {
      if(i < pop->mu) { 
        if(filter[i] == 0)
          EMO_List_queue(available, i);
      }
      else if(filter[i] == 1) {
        EMO_List_queue(missing, i);
      }
    }
  }

  if(pop->asize == 0 && missing->size != available->size) {
    printf("Error, the cardinality of the lists available (%d) and missing (%d) must be the same in EMO_Population_survive, archive %d, filter %d\n", available->size, missing->size, (pop->asize != 0), (filter != NULL));

    EMO_List_print(available, stdout, "available");
    EMO_List_print(missing, stdout, "missing");

    exit(1);
  }

  while(missing->size != 0 && available->size != 0) {
    EMO_List_dequeue(missing, &i);
    EMO_List_dequeue(available, &j);
    EMO_Population_swap(pop, fit1, fit2, mop, i, j);
  }
}

void EMO_Population_evaluate(EMO_Population *pop, EMO_MOP *mop, int start, int size) {
  int i, j, c, o, v, end;

  end = start + size;

  if(mop->ncon == 0) {
    for(i = start; i < end; i++) {
      o = i * mop->nobj;
      v = i * mop->nvar;
      mop->f(mop, pop->obj + o, pop->var + v);
      mop->feval++;
    }
  }
  else {
    for(i = start; i < end; i++) {
      o = i * mop->nobj;
      v = i * mop->nvar;
      c = i * mop->ncon;
      mop->fc(mop, pop->obj + o, pop->con + c, pop->var + v);
      mop->feval++;

      if(pop->vio + i != NULL) {
        pop->vio[i] = 0;

        for(j = mop->ncon - 1; j > -1; j--)
          if(pop->con[c + j] < 0.0)
            pop->vio[i]++;
      }
    }
  }
}

/* Calculate the constraint violation value of NSGA3 and MOEA/D */
void EMO_Population_constrain(EMO_Population *pop, EMO_MOP *mop, int start, int size) {
  int i, j, end;
  double *g;

  end = start + size;

  if(mop->ncon > 0) {
    for(i = start; i < end; i++) {
      pop->cv[i] = 0.0;

      if(pop->vio[i] > 0) {
        g = pop->con + i * mop->ncon;

        for(j = 0; j < mop->ncon; j++)
          if(g[j] < 0.0)
            pop->cv[i] += -g[j];
      }
    }
  }
}

/* Converts the letters of str to upper case.
   Caller is responsible for deallocating src */
char *EMO_toupper(const char *src) {
  char *dest;
  int i, n;

  n = strlen(src);
  dest = (char *) malloc(sizeof(char) * (n + 1));

  if(dest == NULL) {
    printf("Error, not enough memory in EMO_toupper.\n"); 
    exit(1);
  }

  for(i = 0; i < n; i++)
    dest[i] = toupper(src[i]);

  dest[i] = '\0';
  return dest;
}

/* Converts the letters of str to upper case.
   Caller is responsible for deallocating src */
char *EMO_tolower(const char *src) {
  char *dest;
  int i, n;

  n = strlen(src);
  dest = (char *) malloc(sizeof(char) * (n + 1));

  if(dest == NULL) {
    printf("Error, not enough memory in EMO_toupper.\n"); 
    exit(1);
  }

  for(i = 0; i < n; i++)
    dest[i] = tolower(src[i]);

  dest[i] = '\0';
  return dest;
}

/* Returns the position of the first occurence of the string 
 * 'pattern' in the dictionary 'dicc'.
 * If the pattern is not found, it returns -1.
 *
 * Note: a dictionary is an array of strings, that should 
 * have NULL in the last position. */

int EMO_Dictionary_find(const char *dicc[], char *pattern) {
  int i = 0;

  while(dicc[i] != NULL) {
    if(strcmp(pattern, dicc[i]) == 0) {
      return i;
    }
    i++;
  }

  return -1;
}

void EMO_Dictionary_print(FILE *fp, const char *dicc[], const char *s) {
  int i = 0;

  if(s != NULL)
    fprintf(fp, "%s = {", s);

  while(dicc[i] != NULL) {
    if(dicc[i+1] == NULL)
      fprintf(fp, "%s", dicc[i++]);
    else
      fprintf(fp, "%s, ", dicc[i++]);
  }

  if(s != NULL)
    fprintf(fp, "}\n");
}

#undef _MAXCHAR
#undef _MAXPSIZE


