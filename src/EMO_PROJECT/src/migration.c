
#include <string.h>
#include <stdlib.h>

#include "migration.h"
#include "common.h"
#include "io.h"

void EMO_Migration_alloc(EMO_Migration *m, EMO_Rand *rnd, int max_size, int nobj) {
  int i;

  if(max_size <= 0) {
    printf("Error, invalid value of max_size %d in EMO_Migration_alloc\n", max_size);
    exit(1);
  }

  m->seq    = (int *) malloc(sizeof(int) * max_size);  // size
  m->tmp    = (int *) malloc(sizeof(int) * max_size);  // size
  m->filter = (int *) calloc(sizeof(int), max_size);   // size

  if(m->seq == NULL || m->tmp == NULL || m->filter == NULL) {
    printf("Error, not enough memory in EMO_Migration_alloc\n");
    exit(1);
  }

  for(i = 0; i < max_size; i++)  // mu
    m->seq[i] = i;

  EMO_NDSort_alloc(&m->nd, max_size);  // size
  EMO_List_alloc(&m->best, max_size);  // mu
  EMO_List_alloc(&m->worst, max_size); // mu  
  m->rnd = rnd;

  m->max_size = max_size;
  m->nobj = nobj;
}

void EMO_Migration_free(EMO_Migration *m) {
  free(m->seq);
  free(m->tmp);
  free(m->filter);
  EMO_NDSort_free(&m->nd);
  EMO_List_free(&m->best);
  EMO_List_free(&m->worst);
}

void EMO_Migration_get_type(int *type, const char *str) {
  const char *name [] = { "random",
                          "elitist_random",
                          "elitist_ranking",
                          "front",
                          "front_random",
                          "front_ranking",
                          NULL };

   const int id []  = {  EMO_MIGRATION_RANDOM,
                         EMO_MIGRATION_ELITIST_RANDOM,
                         EMO_MIGRATION_ELITIST_RANKING,
                         EMO_MIGRATION_FRONT,
                         EMO_MIGRATION_FRONT_RANDOM,
                         EMO_MIGRATION_FRONT_RANKING,
                       };

  char *migration;
  int i;

  /* Find the function's name in dictionary */
  migration = EMO_tolower(str);
  i = EMO_Dictionary_find(name, migration);

  if(i == -1) {
    printf("Error, unknown migration policy %s in EMO_Migration_get_type.\n", str);
    free(migration);
    exit(1);
  }

  *type = id[i];
  free(migration);
}

void EMO_Replacement_get_type(int *type, const char *str) {
  const char *name [] = { "random",
                          "elitist_random",
                          "elitist_ranking",
                          "elitist",
                          NULL };

   const int id []  = {  EMO_REPLACEMENT_RANDOM,
                         EMO_REPLACEMENT_ELITIST_RANDOM,
                         EMO_REPLACEMENT_ELITIST_RANKING,
                         EMO_ELITIST,
                       };

  char *replacement;
  int i;

  /* Find the function's name in dictionary */
  replacement = EMO_tolower(str);
  i = EMO_Dictionary_find(name, replacement);

  if(i == -1) {
    printf("Error, unknown replacement policy %s in EMO_Replacement_get_type.\n", str);
    free(replacement);
    exit(1);
  }

  *type = id[i];
  free(replacement);
}


/* Select nmig random individuals for migration or replacement.
 */
void EMO_Migration_random(EMO_Migration *m, EMO_List *l, int size, int nmig) {
  int i = 0;

  EMO_List_clear(l);
  EMO_Rand_shuffle(m->rnd, m->seq, m->max_size);  // m->mu

  while(i < m->max_size && l->size < nmig) { // m->mu 
    if(m->seq[i] < size)
      EMO_List_queue(l, m->seq[i]);
    i++;
  }
}

/* Migrate nmig randomly selected individuals from PFcurrent. If all nmig members 
 * cannot be found in PFcurrent, the rest is selected randomly from the remaining 
 * population.
 */
void EMO_Migration_elitist_random(EMO_Migration *m, EMO_List *l, double *data, int size, int nmig) {
  int i, j, n;

  EMO_List_clear(l);
  n = EMO_Dominance_ndset(m->filter, data, NULL, size, m->nobj, EMO_Dominance_strict);  // m->mu

  j = 0;

  if(nmig >= n) { /* Copy individuals from the non-dominated set */
    for(i = 0; i < size; i++) {
      if(m->filter[i])
        EMO_List_queue(l, i);
      else
        m->tmp[j++] = i;
    }
  }
  else { /* Mark all non-dominated set */
    for(i = 0; i < size; i++)
      if(m->filter[i])
        m->tmp[j++] = i;
  }

  /* Randomly complete the migration list */
  if(l->size != nmig) {
    EMO_Rand_shuffle(m->rnd, m->tmp, j);
    i = 0;

    while(i < j && l->size < nmig)
      EMO_List_queue(l, m->tmp[i++]);
  }
}

/* Migrate nmig randomly selected individuals from PFcurrent. If all nmig members
 * cannot be found in PFcurrent, the rest is selected from the successively 
 * ranked Pareto fronts.
 */
void EMO_Migration_elitist_ranking(EMO_Migration *m, EMO_List *l, double *data, int size, int nmig) {
  int i, j, k, n, elem;

  EMO_List_clear(l);
  EMO_NDSort_run(&m->nd, data, m->nobj, NULL, NULL, size); //mu

  j = 0;

  while(j < m->nd.nfront && l->size < nmig) {
    n = m->nd.front[j].size;

    if(nmig >= n + l->size) { /* Copy individuals from the current front */

      EMO_List_append(l, &m->nd.front[j]);

      j++;
    }
    else { /* Randomly complete the migration list from the current front */
      for(i = 0; i < n; i++) {
        EMO_List_get(&m->nd.front[j], &elem, i); 
        m->tmp[i] = elem;
      }

      EMO_Rand_shuffle(m->rnd, m->tmp, n);
      k = nmig - l->size;

      for(i = 0; i < k; i++) 
        EMO_List_queue(l, m->tmp[i]);
    }
  }
}

/* Migrate the entire locally known Pareto front (elitist non-uniform migration). */
void EMO_Migration_front(EMO_Migration *m, EMO_List *l, double *data, int size) {
  int i;

  EMO_List_clear(l);
  EMO_Dominance_ndset(m->filter, data, NULL, size, m->nobj, EMO_Dominance_strict); // m->mu

  for(i = 0; i < size; i++) {
    if(m->filter[i])
      EMO_List_queue(l, i);
  }
}

/* Migrate the entire locally known Pareto front plus nmig randomly selected
 * individuals from the remaining population. 
 */
void EMO_Migration_front_random(EMO_Migration *m, EMO_List *l, double *data, int size, int nmig) {
  int i, j, k, n;

  EMO_List_clear(l);
  n = EMO_Dominance_ndset(m->filter, data, NULL, size, m->nobj, EMO_Dominance_strict); // m->mu

  j = 0;

  for(i = 0; i < size; i++) {
    if(m->filter[i])
      EMO_List_queue(l, i);
    else
      m->tmp[j++] = i;
  }

  EMO_Rand_shuffle(m->rnd, m->tmp, j);
  k = nmig + n;
  i = 0;

  while(i < j && l->size < k)
    EMO_List_queue(l, m->tmp[i++]);
}

/* Migrate the entire locally known Pareto front plus nmig members from 
 * the successively ranked Pareto fronts.
 */
void EMO_Migration_front_ranking(EMO_Migration *m, EMO_List *l, double *data, int size, int nmig) {
  int i, j, k, h, n, elem;

  EMO_List_clear(l);
  EMO_NDSort_run(&m->nd, data, m->nobj, NULL, NULL, size); // m->mu

  j = 0;
  h = nmig + m->nd.front[j].size;

  while(j < m->nd.nfront && l->size < h) {
    n = m->nd.front[j].size;

    if(h >= n + l->size) { /* Copy individuals from the first fronts */
      EMO_List_append(l, &m->nd.front[j]);

      j++;
    }
    else { /* Randomly complete the migration list from the successive front */
      for(i = 0; i < n; i++) {
        EMO_List_get(&m->nd.front[j], &elem, i); 
        m->tmp[i] = elem;
      }

      EMO_Rand_shuffle(m->rnd, m->tmp, n);
      k = h - l->size;

      for(i = 0; i < k; i++) 
        EMO_List_queue(l, m->tmp[i]);
    }
  }
}

void EMO_Replacement_random(EMO_Migration *m, EMO_List *migrate, EMO_List *replace, double *data, int size, int nmig) {
  int i;

  EMO_List_clear(migrate);

  EMO_Migration_random(m, replace, size, nmig); //m->mu

  for(i = 0; i < nmig; i++)
    EMO_List_queue(migrate, size + i); //m->mu
}

/* Keep all members of Pcurrent determining which of the remaining population is not
 * Pareto optimal with respect to the immigrants, and then randomly replacing nmig of
 * those individuals with the immigrants. The list of individuals to be replaced may 
 * be smaller than nmig. This method does not guarantee the worst individuals are 
 * removed but does ensure PFcurrent retention. 
 * Note: migrants should be stored at data + size * nobj.
 */
void EMO_Replacement_elitist_random(EMO_Migration *m, EMO_List *migrate, EMO_List *replace, double *data, int size, int nmig) {
  int i, j, k, n, elem;

  EMO_List_clear(replace);
  EMO_List_clear(migrate);

  if(nmig > size) {
    printf("Error, nmig is bigger than the population size\n");
    exit(1);
  }

  memset(m->tmp, 0, sizeof(int) * nmig);
  memset(m->filter, 0, sizeof(int) * size);

  /* Check dominance of candidates to be removed */
  for(j = 0; j < nmig; j++) {
    for(i = 0; i < size; i++) { // m->mu
      if(EMO_Dominance_strict(data + (size + j) * m->nobj, data + i * m->nobj, m->nobj) == 1) {
        if(m->filter[i] == 0) {  /* Non-dominated solutions */
          EMO_List_queue(replace, i);
          m->filter[i] = 1;
        }
        m->tmp[j] = 1;
      }
    }
  }

  for(i = 0; i < nmig; i++)
    if(m->tmp[i] == 1)
      EMO_List_queue(migrate, i);

  n = replace->size;

  if(n > migrate->size) {  /* Candidates are randomly discarded if they are more than nmig */
    for(i = 0; i < n; i++) {
      EMO_List_get(replace, &elem, i); 
      m->tmp[i] = elem;
    }

    EMO_Rand_shuffle(m->rnd, m->tmp, n);
    k = n - nmig;

    for(i = 0; i < k; i++) 
      EMO_List_remove(replace, m->tmp[i]);
  }

  n = migrate->size;

  if(replace->size < n) {
    for(i = 0; i < n; i++) {
      EMO_List_get(migrate, &elem, i); 
      m->tmp[i] = elem;
    }

    EMO_Rand_shuffle(m->rnd, m->tmp, n);
    k = n - replace->size;

    for(i = 0; i < k; i++) 
      EMO_List_remove(migrate, m->tmp[i]);
  }
}


/* Rank all Pareto fronts and replace individuals from the 'worst' ranked front(s)
 * with the immigrants. Individuals may be randomly selected from some front(s) 
 * if the number of individuals to be replaced does not match the number of 
 * individuals represented by the front being replaced. This method provides the 
 * greatest selection pressure as the 'worst' individuals are continually removed 
 * from the population.
 * The resulting lists (replace and migrate) are of the same size and indicate the 
 * members of the population that should be interchanged. Their cardinality may be 
 * in the range [0,size).
 *
 * Note: migrants should be stored at data + size * nobj
 */
void EMO_Replacement_elitist_ranking(EMO_Migration *m, EMO_List *migrate, EMO_List *replace, double *data, int size, int nmig) {
  int i, j, k, n, elem;

  EMO_List_clear(replace);
  EMO_List_clear(migrate);
 
  EMO_NDSort_run(&m->nd, data, m->nobj, NULL, NULL, size);  // m->mu

  j = m->nd.nfront - 1;

  while(j >= 0 && replace->size < nmig) {
    n = m->nd.front[j].size;
 
    if(nmig >= n + replace->size) { /* Copy individuals from the current front */
      EMO_List_append(replace, &m->nd.front[j]);
    }
    else { /* Randomly complete the migration list from the current front */
      for(i = 0; i < n; i++) {
        EMO_List_get(&m->nd.front[j], &elem, i); 
        m->tmp[i] = elem;
      }

      EMO_Rand_shuffle(m->rnd, m->tmp, n);
      k = nmig - replace->size;

      for(i = 0; i < k; i++) 
        EMO_List_queue(replace, m->tmp[i]);
    }
    j--;
  }

  // Initialize migration list 
  if(replace->size == nmig) {
    for(i = 0; i < nmig; i++)
      EMO_List_queue(migrate, size + i);  // m->mu 
  }
  else {
    for(i = 0; i < nmig; i++)
      m->tmp[i] = i;

    EMO_Rand_shuffle(m->rnd, m->tmp, nmig);
    i = 0;

    while(migrate->size < replace->size) {
      EMO_List_queue(migrate, size + m->tmp[i]);  // m->mu
      i++;
    }
  }
}


/* Combine immigrants with the current population, rank all Pareto fronts and remove
 * individuals from the 'worst' ranked front(s).
 * The general idea is to replace the individuals of the worst front of the current 
 * population with the best migrants. When there is a tie between these two sets, 
 * individuals are randomly replaced.
 * The resulting lists (replace and migrate) are of the same size and indicate the 
 * members of the population that should be interchanged. Their cardinality may be in 
 * the range [0,size).
 *
 * Note: migrants should be stored at data + size * nobj.
 */
void EMO_Replacement_elitist(EMO_Migration *m, EMO_List *migrate, EMO_List *replace, double *data, int size, int nmig) {
  int i, j, k, sb, sw, elem;

  EMO_List_clear(replace);
  EMO_List_clear(migrate);
  EMO_List_clear(&m->best);
  EMO_List_clear(&m->worst);

  EMO_NDSort_run(&m->nd, data, m->nobj, NULL, NULL, size + nmig);  // m->mu

  i = -1;
  j = m->nd.nfront;

  do {
    while(m->best.size == 0) {
      i++;

      if(i >= m->nd.nfront) 
        break;

      for(k = m->nd.front[i].size - 1; k > -1; k--) {
        EMO_List_get(&m->nd.front[i], &elem, k);

        if(elem >= size)   // m->mu
          EMO_List_queue(&m->best, elem);
      }
    }

    while(m->worst.size == 0) {
      j--;

      if(j < 0)
        break;

      for(k = m->nd.front[j].size - 1; k > - 1; k--) {
        EMO_List_get(&m->nd.front[j], &elem, k);

        if(elem < size) // m->mu
          EMO_List_queue(&m->worst, elem);
      }
    }

    sb = m->best.size;
    sw = m->worst.size;

    if(sb != 0 && sw != 0 && i <= j) {

      if(sb > sw) {

        while(m->worst.size > 0) {
          EMO_List_dequeue(&m->worst, &elem);
          EMO_List_queue(replace, elem);

          k = EMO_Rand_int1(m->rnd, 0, m->best.size - 1);
          EMO_List_retrieve(&m->best, &elem, k);
          EMO_List_queue(migrate, elem);
        }
      }
      else if(sb < sw) {

        while(m->best.size > 0) {
          EMO_List_dequeue(&m->best, &elem);
          EMO_List_queue(migrate, elem);

          k = EMO_Rand_int1(m->rnd, 0, m->worst.size - 1);
          EMO_List_retrieve(&m->worst, &elem, k);
          EMO_List_queue(replace, elem);
        }
      }
      else {

        while(m->best.size > 0) {
          EMO_List_dequeue(&m->best, &elem);
          EMO_List_queue(migrate, elem);
          EMO_List_dequeue(&m->worst, &elem);
          EMO_List_queue(replace, elem);
        }
      }
    }
  } while(i < m->nd.nfront && j > -1 && i <= j);
}

