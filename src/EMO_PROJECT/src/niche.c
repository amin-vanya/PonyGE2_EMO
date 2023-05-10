
#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "numeric.h"
#include "vector.h"
#include "niche.h"
#include "sort.h"

#define _EPS 1e-6
#define _MAXNN 3

void EMO_Prune_alloc(EMO_Prune *p, int size, int nobj) {
  int i;

  EMO_NDSort_alloc(&p->nd, size);
  EMO_List_alloc(&p->lst, nobj);
  EMO_Refpoint_alloc(&p->ref, NULL, nobj); 

  if((p->min = (double *) malloc(sizeof(double) * nobj)) == NULL) {
    printf("Error, not enough memory in EMO_Prune_alloc\n");
    exit(1);
  }

  if((p->ideal = (double *) malloc(sizeof(double) * nobj)) == NULL) {
    printf("Error, not enough memory in EMO_Prune_alloc\n");
    exit(1);
  }

  if((p->nadir = (double *) malloc(sizeof(double) * nobj)) == NULL) {
    printf("Error, not enough memory in EMO_Prune_alloc\n");
    exit(1);
  }

  if((p->norm = (double *) malloc(sizeof(double) * size * nobj)) == NULL) {
    printf("Error, not enough memory in EMO_Prune_alloc.\n");
    exit(1);
  }

  if((p->vnorm = (double *) malloc(sizeof(double) * size)) == NULL) {
    printf("Error, not enough memory in EMO_Prune_alloc.\n");
    exit(1);
  }

  if((p->filter = (int *) calloc(sizeof(int), size)) == NULL) {
    printf("Error, not enough memory in EMO_Prune_alloc\n");
    exit(1);
   }

  if((p->tmp = (int *) calloc(sizeof(int), size)) == NULL) {
    printf("Error, not enough memory in EMO_Prune_alloc.\n");
    exit(1);
  }

  if((p->sort = (double **) malloc(sizeof(double *) * size)) == NULL) {
    printf("Error, not enough memory in mombi2.\n");
    exit(1);
  }

  for(i = 0; i < size; i++) {
    if((p->sort[i] = (double *) malloc(sizeof(double) * 2)) == NULL) {
      printf("Error, not enough memory in mombi2.\n");
      exit(1);
    }
  }  
  p->ssize = size;

  /* Initialize reference points */
  for(i = 0; i < nobj; i++) {
    p->ideal[i] = DBL_MAX;
    p->nadir[i] = DBL_MAX;
  }
}

void EMO_Prune_free(EMO_Prune *p) {
  int i;

  EMO_NDSort_free(&p->nd);
  EMO_List_free(&p->lst);
  EMO_Refpoint_free(&p->ref);
  free(p->min);
  free(p->ideal);
  free(p->nadir);
  free(p->norm);
  free(p->vnorm);
  free(p->filter);
  free(p->tmp);

  for(i = p->ssize - 1; i > -1; i--)
    free(p->sort[i]);

  free(p->sort);
}

void EMO_nicheCount(double *niche, double *data, int row, int col, double alpha, double sigma_share) {
  double d, sh;
  int i, j;
 
  for(i = 0; i < row; i++)  {
    sh = 0;

    for(j = 0; j < row; j++) {
      if (i == j) continue;

      d = EMO_vdist(data + i*col, data + j*col, col);
      //printf("niche %f\n", d);

      if(d < sigma_share)
        sh += 1.0 - pow(d / sigma_share, alpha);  /* sharing function */
    }
    niche[i] = sh;
  }
}

void EMO_crowdingDistance(double *cd, double **sort, double *data, int row, double *max, double *min, int col) {
  int i, j, k, l, m;
  double v;

  memset(cd, 0, sizeof(double) * row);

  for(i = 0; i < col; i++) {

    for(j = 0; j < row; j++) {
      sort[j][0] = (double) j;
      sort[j][1] = data[j * col + i];
    }
      
    qsort(sort, row, sizeof(sort[0]), (int (*)(const void *, const void *))&EMO_compare_asc);

    l = (int) sort[0][0];
    m = (int) sort[row-1][0];
    cd[l] = cd[m] = DBL_MAX;

    v = max[i] - min[i];

    if(v != 0.0) {
      for(j = row - 2; j > 0; j--) {
        k = (int) sort[j][0];

        if(cd[k] != DBL_MAX) {
          l = (int) sort[j+1][0];
          m = (int) sort[j-1][0];
          cd[k] += (data[l*col + i] - data[m*col+i]) / v;
        }
      }
    }
  }
}

// solo del frente que le corresponda
void EMO_crowdingDistance2(double *cd, double **sort, double *data, EMO_List *front, double *max, double *min, int col) {
  int i, j, k, l, m;
  double v;

  for(i = 0; i < col; i++) {

    for(j = 0; j < front->size; j++) {
      EMO_List_get(front, &k, j);
      if(i == 0) cd[k] = 0.0;
      sort[j][0] = (double) k;
      sort[j][1] = data[k * col + i];
    }
      
    qsort(sort, front->size, sizeof(sort[0]), (int (*)(const void *, const void *))&EMO_compare_asc);

    l = (int) sort[0][0];
    m = (int) sort[front->size-1][0];
    cd[l] = cd[m] = DBL_MAX;

    v = max[i] - min[i];

    if(v != 0.0) {
      for(j = front->size - 2; j > 0; j--) {
        k = (int) sort[j][0];

        if(cd[k] != DBL_MAX) {
          l = (int) sort[j+1][0];
          m = (int) sort[j-1][0];
          cd[k] += (data[l*col + i] - data[m*col+i]) / v;
        }
      }
    }
  }
}

// Distance from point to a weight vector (they must be normalized)
void EMO_wdist(double *wd, double *data, int size, double *W, int wsize, double *tmp, int dim) {
  double d, vmin = 0;
  int i, j;

  for(i = 0; i < size; i++) {
    vmin = DBL_MAX;

    for(j = 0; j < wsize; j++) {
      EMO_vorth(tmp, &W[j*dim], data+i*dim, dim);  // Projection
      d = EMO_vnorm(tmp, 2.0, dim);

      if(d < vmin) vmin = d; 
    }
    wd[i] = vmin;
  }
}



/* Retrieves the k-nearest neighbor (elem) of individual p that is activated by filter */
int _KNN_get(EMO_List *lnn, int *filter, int *elem, int *k, int p) {

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

/* Compares two individuals based on their distances to their neighbors */
int _KNN_compare(EMO_List *lnn, double **dist, int *filter, int p1, int p2, int size) {
  int i, j, e1, e2;

  i = j = 0;

  while(i < size && j < size) {

    if(_KNN_get(lnn, filter, &e1, &i, p1)) {
      printf("Error inesperado e1 %d i %d p1 %d p2 %d \n", e1, i, p1, p2);
      EMO_List_print(&lnn[p1], stdout, "");
      EMO_List_print(&lnn[p2], stdout, "");
      return 1;
    }

    if(_KNN_get(lnn, filter, &e2, &j, p2)) {
      printf("Error inesperado e2 %d j %d p2 %d p1 %d\n", e2, j, p2, p1);
      EMO_List_print(&lnn[p2], stdout, "");
      EMO_List_print(&lnn[p1], stdout, "");
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

void EMO_KNN_alloc(EMO_KNN *knn, int size, int dim) {
  int i;

  if((knn->lnn = (EMO_List *) malloc(sizeof(EMO_List) * size)) == NULL) {
    printf("Error, not enogh memory in KNN.\n");
    exit(1);
  }

  for(i = 0; i < size; i++)
    EMO_List_alloc(&knn->lnn[i], size);

  EMO_List_alloc(&knn->copy, size);

  if((knn->sort = (double **) malloc(sizeof(double *) * size)) == NULL) {
    printf("Error, not enough memory in KNN.\n");
    exit(1);
  }

  for(i = 0; i < size; i++) {
    if((knn->sort[i] = (double *) malloc(sizeof(double) * 2)) == NULL) { 
      printf("Error, not enough memory in KNN.\n");
      exit(1);
    }
  }

  knn->ssize = size;

  if((knn->dist = (double **) malloc(sizeof(double *) * size)) == NULL) {
    printf("Error, not enough memory in KNN.\n");
    exit(1);
  }

  for(i = 0; i < size; i++) {
    if((knn->dist[i] = (double *) malloc(sizeof(double) * size)) == NULL) {
      printf("Error, not enough memory in KNN.\n");
      exit(1);
    }
  }

  if((knn->min = (int *) malloc(sizeof(int) * dim)) == NULL) {
    printf("Error, not enough memory in KNN.\n");
    exit(1);
  }

  if((knn->max = (int *) malloc(sizeof(int) * dim)) == NULL) {
    printf("Error, not enough memory in KNN.\n");
    exit(1);
  }

  knn->dim = dim;
}

void EMO_KNN_free(EMO_KNN *knn) {
  int i;

  for(i = knn->ssize - 1; i > -1; i--) {
   EMO_List_free(&knn->lnn[i]);
   free(knn->sort[i]);
   free(knn->dist[i]);
  }

  free(knn->lnn);
  free(knn->sort);
  free(knn->dist);
  free(knn->min);
  free(knn->max);

  EMO_List_free(&knn->copy);
}

void EMO_KNN_prune(EMO_KNN *knn, int *filter, int max_size, double *data, int size) {
  int i, j, cont = 0;

  EMO_knn2(knn->lnn, &knn->copy, knn->dist, knn->sort, data, size, knn->dim);

  // Count solutions 
  for(i = 0; i < size; i++)
    if(filter[i]) cont++;

  // Eliminates copies 
  while(cont > max_size && knn->copy.size > 0) {
    EMO_List_dequeue(&knn->copy, &j);

    if(filter[j]) {
      filter[j] = 0;
      cont--;
    }
  }

  // Selects individuals having a minimum or maximum objective value (weakness of SPEA2, since it does not do that)
  EMO_findmaxminBound(knn->max, knn->min, data, filter, size, knn->dim);

  // Finds individuals with the worst distribution 
  while(cont > max_size) {
    i = 0;

    while(i < size && filter[i] == 0) i++;

    while(EMO_find(knn->max, knn->dim, i) || EMO_find(knn->min, knn->dim, i)) i++;

    for(j=i+1; j < size; j++) {
      if(filter[j] == 0) continue;
      if(EMO_find(knn->max, knn->dim, j)) continue;

      if(_KNN_compare(knn->lnn, knn->dist, filter, j, i, size) == 1)
        i = j;
    }

    filter[i] = 0;
    cont--;
  }
}

/* Computes the k-nearest neighbors of a set of vectors W */
void EMO_knn(EMO_List *l, double **sort, double *W, int size, int dim, int k, int itself) {
  int i, j, idx;
  double d;

  for(i = 0; i < size; i++)
    EMO_List_clear(&l[i]);

  for(i = 0; i < size; i++) {
    for(j = 0; j < size; j++) {
      d = EMO_vdist(W + i*dim, W + j*dim, dim);  /* euclidean distance */
      sort[j][0] = j;
      sort[j][1] = d;
    }
 
    qsort(sort, size, sizeof(sort[0]), (int (*)(const void *, const void *))&EMO_compare_asc);
    j = 0;

    if(itself) {
      while(l[i].size < k) {
        idx = (int) sort[j][0];
        EMO_List_queue(&l[i] , idx);
        j++;
      }
    }
    else {
      while(l[i].size < k) {
        if(i == j) break;
        idx = (int) sort[j][0];
        EMO_List_queue(&l[i] , idx);
        j++;
      }
    }
  }
}

/* Computes the k-nearest neighbors of a set of vectors W, returns list of neighbors for each vector, a list of clon solutions and a matrix of distances */
void EMO_knn2(EMO_List *l, EMO_List *copy, double **dist, double **sort, double *W, int size, int dim) {
  int i, j, idx;

  //O(N^2 m)
  for(i = 0; i < size; i++) {
    dist[i][i] = 0.0;

    for(j = i+1; j < size; j++) { 
      dist[i][j] = dist[j][i] = EMO_vdist(W + i*dim, W + j*dim, dim);

      if(dist[i][j] < _EPS) {
        EMO_List_queue(copy, j);
      }
    }
  }

  //O(N (N + N log N)
  for(i = 0; i < size; i++) {
    EMO_List_clear(&l[i]);

    for(j = 0; j < size; j++) {
      sort[j][0] = j;
      sort[j][1] = dist[i][j];
    }

    qsort(sort, size, sizeof(sort[0]), (int (*)(const void *, const void *))&EMO_compare_asc);

    j = 0;

    while(j < size) {
      idx = (int) sort[j][0];

      if(idx != i)
        EMO_List_queue(&l[i], idx);
      j++;
    }
  }
}

void EMO_kNN2(EMO_List *lst, double **sort, double *data, int row, int col) {
  int i, j, k, l, n;
  double v, w, y;

  n = col * 3;

  for(i = 0; i < row; i++) {
    for(j = 0; j < row; j++) {
      if(i == j) {
        sort[j][0] = 0; 
        sort[j][1] = DBL_MAX; 
      }
      else {
        sort[j][0] = (double) j;
        sort[j][1] = EMO_vdist(data+i*col, data+j*col, col);
      }
    } 

    qsort(sort, row, sizeof(sort[0]), (int (*)(const void *, const void *))&EMO_compare_asc);

    EMO_List_clear(&lst[i]);

    w = sort[0][1];
    y = 0;
    v = 0;
    j = 0;

    for(k = 1; k < n; k++) {
      w = sort[k][1] - sort[k-1][1];
      y = w - y;
      if(w > v && k != 1 && k >= col) {
        v = w;
        j = k;
      }
    }

    for(k = 0; k < j; k++) {
      l = (int) sort[k][0];
      EMO_List_queue(&lst[i], l);
    }
  }

  // Propiedad simetrica
  for(i = 0; i < row; i++) {
    for(j = 0; j < lst[i].size; j++) {
      EMO_List_get(&lst[i], &k, j);

      if(!EMO_List_seek(&lst[k], i)) 
        EMO_List_retrieve(&lst[i], NULL, j);
      
    }
  }
}


#undef _EPS
#undef _MAXNN

