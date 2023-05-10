#include <stdlib.h>
#include <math.h>
#include <stdio.h> // borrar

#include "sort.h"
#include "stat.h"

double EMO_mean(double *data, int *filter, int size) {
  double v = 0;
  int i;

  if(filter == NULL) {
    for(i = 0; i < size; i++)
      v += data[i];
  }
  else {
    for(i = 0; i < size; i++)
      if(filter[i])
        v += data[i];
  }

  return v / size;
}


/* returns the IQR
 * q: arrays of 3 elements (q1, q2, q3)
 *
 * odd: |6, 7, 15, 36, 39| 40| 41, 42, 43, 47, 49|   size = 11
 *             q1          q2          q3
 *
 * even: |7, 15, 36| 39, 40, 41|   size = 6
 *               q2  q2
 *           q1          q3
 */
double EMO_quartile(double *q, double *data, int *filter, double **sort, int size) {
  int i, n;

  if(filter == NULL) {
    n = size;

    for(i = 0; i < size; i++) {
      sort[i][0] = (double) i; 
      sort[i][1] = data[i];
    }
  }
  else {
    n = 0;

    for(i = 0; i < size; i++) {
      if(filter[i]) {
        sort[n][0] = (double) i; 
        sort[n++][1] = data[i];
      }
    }
  }

  qsort(sort, n, sizeof(sort[0]), (int (*) (const void *, const void *)) &EMO_compare_asc);

  q[0] = sort[n/4][1];
  q[2] = sort[(int) (3.0/4.0*n)][1];

  if(n % 2 == 0)  // even
    q[1] = (sort[n/2][1] + sort[n/2-1][1]) / 2;
  else  // odd
    q[1] = sort[n/2][1];

  //returns IQR
  return q[2] - q[0];
}

// representative of the sample, no average when there are odd numer of elements
double EMO_median(int *idx, double *data, int *filter, double **sort, int size) {
  int i, j, n;

  if(filter == NULL) {
    n = size;

    for(i = 0; i < n; i++) {
      sort[i][0] = (double) i; 
      sort[i][1] = data[i];
    }
  }
  else {
    n = 0;

    for(i = 0; i < n; i++) {
      sort[n][0] = (double) i; 
      sort[n++][1] = data[i];
    }
  }

  qsort(sort, n, sizeof(sort[0]), (int (*) (const void *, const void *)) &EMO_compare_asc);

  if(n % 2 == 0)
    j = sort[n/2 - 1][0];
  else
    j = sort[n/2][0];

  if(idx != NULL)
    *idx = j;

  return data[j];
}

double EMO_var(double *data, int *filter, double mean, int size) {
  double v = 0;
  int i;

  if(filter == NULL) {
    for(i = 0; i < size; i++)
      v += pow(data[i] - mean, 2.0);
  }
  else {
    for(i = 0; i < size; i++)
      if(filter[i])
        v += pow(data[i] - mean, 2.0);
  }
 
  return v / size; 
}

double EMO_std(double *data, int *filter, double mean, int size) {
  return sqrt(EMO_var(data, filter, mean, size));
}


