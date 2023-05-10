/*
 * Authors: 
 *          Elmer Jesus Morales Orozco 
 *          Adriana Menchaca Mendez
 *          Raquel Hernandez Gomez
 */

#ifndef _INDICATOR_
#define _INDICATOR_

#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "indicator.h"
#include "dominance.h"
#include "numeric.h"
#include "matrix.h"
#include "vector.h"
#include "hv.h"

double mindist(double *v, int *filter, double *f, int n, int dim) {
  double min, d;
  int i;

  min = DBL_MAX;
  d = 0;

  if(filter == NULL) {
    for(i = 0; i < n; i++) {
      d = EMO_vdist(v, f+i*dim, dim);

      if(d < min) min = d;
    }
  }
  else {
    for(i = 0; i < n; i++) {
      if(filter[i]) {
        d = EMO_vdist(v, f+i*dim, dim);

        if(d < min) min = d;
      }
    }
  }

  return min; 
}

/* v: reference vector
 * f: population
 * filter: active population
 * n: population size
 * dim: number of objectives
 */
double mindist_plus(double *v, double *f, int *filter, int n, int dim) {
  double min, d, s;
  int i, j;

  min = DBL_MAX;
  d = 0;

  if(filter == NULL) {
    for(i = 0; i < n; i++) {
      s = 0;

      for(j = 0; j < dim; j++) {
        d = max(f[i*dim + j] - v[j], 0);  /* minimization problems: a_i - zref_i */
        s += d * d;
      }

      d = sqrt(s);

      if(d < min) min = d;
    }
  }
  else {
    for(i = 0; i < n; i++) {
      if(filter[i]) {
        s = 0;

        for(j = 0; j < dim; j++) {
          d = max(f[i*dim + j] - v[j], 0);
          s += d * d;
        }

        d = sqrt(s);

        if(d < min) min = d;
      }
    }
  }

  return min; 
}

/* Average distance from the aproximation set to the
   discretization of the Pareto front */
double EMO_Indicator_gd(double *f, int *filter, int n, double *pf, int m, int dim, double p) {
  double s;
  int i;

  s = 0.0;

  if(filter == NULL) {
    for(i = 0; i < n; i++)
      s += pow(mindist(f+i*dim, NULL, pf, m, dim), p);
  }
  else {
    for(i = 0; i < n; i++)
      if(filter[i])
        s += pow(mindist(f+i*dim, NULL, pf, m, dim), p);
  }

  return pow(s / (double) n, 1.0 / p);
}

/* Average distance from the discretization of the
   Pareto front to the aproximation set */
double EMO_Indicator_igd(double *f, int *filter, int n, double *pf, int m, int dim, double p) {
  double s;
  int i;

  s = 0.0;

  for(i = 0; i < m; i++)
    s += pow(mindist(pf+i*dim, filter, f, n, dim), p);

  return pow(s / (double) m, 1.0 / p);
}

/* Average distance from the discretization of the
   Pareto front to the aproximation set */
double EMO_Indicator_igd2(double *f, int *filter, int n, double *pf, int m, int dim) {
  double s;
  int i;

  s = 0.0;

  for(i = 0; i < m; i++)
    s += mindist(pf+i*dim, filter, f, n, dim);

  return s / (double) m;
}

/* Professor Ishibuchi */
double EMO_Indicator_igd_plus(double *f, int *filter, int n, double *pf, int m, int dim) {
  double s;
  int i;

  s = 0.0;

  for(i = 0; i < m; i++) {
    s += mindist_plus(pf+i*dim, f, filter, n, dim);
  }

  return s / (double) m;
}

double EMO_Indicator_deltap(double *f, int *filter, int n, double *pf, int m, int dim, double p) {
  double igd, gd;

  igd = EMO_Indicator_igd(f, filter, n, pf, m, dim, p);
  gd = EMO_Indicator_gd(f, filter, n, pf, m, dim, p);

  if(gd > igd) return gd;

  return igd;
}

double EMO_Indicator_r2(double *data, int *filter, int size, double *W, int wsize, EMO_Utility *utl) {
  double d, vmin = 0, s = 0;
  int i, j;

  if(filter == NULL) {
    for(i = 0; i < wsize; i++) {
      vmin = DBL_MAX;

      for(j = 0; j < size; j++) {
        d = utl->uf(utl, &W[i * utl->nobj], data + j * utl->nobj);

        if(d < vmin) vmin = d; 
      }
      s += vmin;
    }
  }
  else {
    for(i = 0; i < wsize; i++) {
      vmin = DBL_MAX;

      for(j = 0; j < size; j++) {
        if(filter[j]) {
          d = utl->uf(utl, &W[i * utl->nobj], data + j * utl->nobj);

          if(d < vmin) vmin = d; 
        }
      }
      s += vmin;
    }
  }

  return s/wsize;
}

/* Compares two individuals */
int compare(const void **a, const void **b) {
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

/* R2 ranking algorithm of the population */
void EMO_Indicator_r2_ranking(double *rank, double **sort, double *norm, double *tmp, double *data, int size, double *W, int wsize, EMO_Utility *utl) {
  int i, j, k;

  for(j = 0; j < size; j++) {
    rank[j] = DBL_MAX;
    norm[j] = EMO_vnorm(data + j * utl->nobj, 2.0, utl->nobj);
  }

  for(i = 0; i < wsize; i++) {
    // Calculates the individual's contribution to a weight vector
    for(j = 0; j < size; j++) {
      sort[j][0] = j;
      sort[j][1] = utl->uf(utl, W + i * utl->nobj, data + j * utl->nobj);
      sort[j][2] = norm[j];

    }
    // Sorts individuals wrt. the utility value obtained in increasing order
    qsort(sort, size, sizeof(sort[0]), (int (*)(const void *, const void *))&compare);

    // Ranks individuals
    for(j = 1; j <= size; j++) {
      k = (int) sort[j-1][0]; 

      if((double) j < rank[k])
         rank[k] = (double) j;
    }
  }
}

int compare2(const void * a, const void * b)
{
 return ( *(int*)a - *(int*)b );
}


void valid_data(double a, double b) {
  if((a < 0 && b > 0) ||
     (a > 0 && b < 0) ||
      a == 0 || b == 0) {
      printf("Error in data a %f, b %f\n", a, b);
      exit(1);
  }
}

/* What is the smallest amount epsilon that the set A should be translated
   so that every point in B is covered?
 */
double EMO_Indicator_epsilon_multiplicative(double *a, int *filter, int na, double *b, int nb, int dim) {
  double diff, eps, eps_a, eps_obj;
  int i, j, k;

  eps = eps_a = 0;

  for(i = 0; i < nb; i++) {
    for(j = 0; j < na; j++) {
      valid_data(a[j*dim], b[i*dim]);
      eps_obj = a[j*dim] / b[i*dim];

      for(k = 1; k < dim; k++) {
        valid_data(a[j*dim + k], b[i*dim + k]);
        diff = a[j*dim + k] / b[i*dim + k];

        if(diff > eps_obj) // max k
          eps_obj = diff;
      }

      if(j == 0 || eps_obj < eps_a)  // min a
        eps_a = eps_obj;
    }

    if(i == 0 || eps_a > eps)  // max b 
      eps = eps_a;
  }
  return eps;
}

double EMO_Indicator_epsilon_additive(double *a, int *filter, int na, double *b, int nb, int dim) {
  double diff, eps, eps_a, eps_obj;
  int i, j, k;

  eps = eps_a = 0;

  for(i = 0; i < nb; i++) {
    for(j = 0; j < na; j++) {
      eps_obj = a[j*dim] - b[i*dim];

      for(k = 1; k < dim; k++) {
        diff = a[j*dim + k] - b[i*dim + k];

        if(diff > eps_obj)  // max k
          eps_obj = diff;
      }
      if(j == 0 || eps_obj < eps_a)  // min a
        eps_a = eps_obj;
    }

    if(i == 0 || eps_a > eps)  // max b
      eps = eps_a;
  }
  return eps;
}

double EMO_Indicator_maximin(double *fit, double *f, int *filter, int n, int dim) {
  double maximin, diff, min, max;
  int i, j, k;

  maximin = 0.0;

  for (i = 0; i < n; i++) {
    max = -DBL_MAX;

    for (j = 0; j < n; j++) {
      if(i != j) {
        min = f[i*dim] - f[j*dim];

        for(k = 1; k < dim; k++) {
          diff = f[i*dim + k] - f[j*dim + k];  

          if(diff < min)
            min = diff;
        }

        if(min > max)
          max = min;
      }
    }
    maximin += max;

    if(fit != NULL)
      fit[i] = max;
  }
  return maximin / (double) n; 
}

// Overall Non-dominated Vector Generation [Sch95, Vel99, Vel00]
double EMO_Indicator_onvg(double *f, int *filter, int n, int dim) {
  return (double) EMO_Dominance_ndset(NULL, f, filter, n, dim, EMO_Dominance_strict);
}

// Overall Non-dominated Vector Generation Ratio [Vel99, Vel00]
double EMO_Indicator_onvgr(double *a, int *filter, int na, double *b, int nb, int dim) {
  double v;

  v = (double) EMO_Dominance_ndset(NULL, b, filter, nb, dim, EMO_Dominance_strict);

  if(v == 0) {
    printf("Error in EMO_Indicator_onvg, division by zero. Reference set has dominated solutions.\n");
    exit(1);
  }

  return (double) EMO_Dominance_ndset(NULL, a, filter, na, dim, EMO_Dominance_strict) / v;
}

// Coverage of two sets [Zit98, Zit99, Zit99b]
double EMO_Indicator_c(double *a, int *filter, int na, double *b, int nb, int dim) {
  int i, j, cont;

  cont = 0;

  if(nb == 0) {
    printf("Error in EMO_indicator_c, set B is empty.\n");
    exit(1);
  }

  for(i = 0; i < na; i++) {
    for(j = 0; j < nb; j++) {
      if (EMO_Dominance_strict(a+i*dim,  b+j*dim, dim) == 1) {
        cont++;
        break;
      }
    }
  }

  return cont / nb; 
}

// Spacing (SP) [Schott95]
double EMO_Indicator_spacing(double *f, int *filter, int n, int dim) {
  double d, avg, vmin, v, sp;
  int i, j, k;

  avg = sp = 0;

  for(i = 0; i < n; i++) {
    vmin = DBL_MAX;

    for(j = 0; j < n; j++) {
      if(i != j) {

        for(k = 0, d = 0; k < dim; k++)
          d += fabs(f[i*dim + k] - f[j*dim + k]);

        if(d < vmin)
          vmin = d;
      }
    }
    avg += vmin;
  }

  avg /= (double) n;

  for(i = 0; i < n; i++) {
    vmin = DBL_MAX;

    for(j = 0; j < n; j++) {
      if(i != j) {

        for(k = 0, d = 0; k < dim; k++)
          d += fabs(f[i*dim + k] - f[j*dim + k]);

        if(d < vmin)
          vmin = d;
      }
    }

    v = (vmin - avg);
    v *= v;
    sp += v;
  }

  return sqrt(sp / (double) (n - 1));
}

// Discretizing Manifolds via Minimum Energy Points
// D.P. Hardin and E.B. Saff, 2004
//
// Generalized Decomposition
// Ioannis Giagkiozis, Robin C. Purshouse, and Peter J. Fleming, 2013
//
// Riesz s-energy

double EMO_Indicator_senergy(double *fit, double *f, int *filter, int n, int dim) {
  double v, s1, s2;
  int i, j;

  s2 = 0;

  if(fit != NULL) 
    memset(fit, 0, sizeof(double) * n);

  for(i = 0; i < n; i++) {
    s1 = 0;

    if(filter != NULL && filter[i] == 0) continue;

    for(j = i+1; j < n; j++) {

      if(filter != NULL && filter[j] == 0) continue;

      v = EMO_vdist(f + i*dim, f + j*dim, dim);

      if(v == 0) {
        v = 1e-4;
      }

      v = pow(v, -(dim - 1.0));

      if(fit != NULL) fit[j] += v;
      s1 += v;
      s2 += 2.0 * v;
    }

    if(fit != NULL)
      fit[i] += s1;
  }

  return s2;
}

double EMO_Indicator_senergy_update(double *fit, double *f, int *filter, int n, int dim, int idx, int action) {
  double v, s1, s2, sign;
  int j;

  if(filter != NULL && filter[idx] == 0) {
    printf("Error, filter at %d is not active (EMO_Indicator_senergy_update).\n", idx);
    exit(1);
  }

  sign = (double) action / fabs((double) action);
  s1 = s2 = 0;

  for(j = 0; j < n; j++) {

    if(idx == j) continue;

    if(filter != NULL && filter[j] == 0) continue;

    v = EMO_vdist(f + idx*dim, f + j*dim, dim);

    if(v == 0) {
      v = 1e-4;
    }

    v = sign * pow(v, -(dim - 1.0));
    if(fit != NULL) fit[j] += v;
    s1 += v;
    s2 += fit[j];
  }

  if(fit != NULL)
    fit[idx] += s1;

  return s2;
}

// Solow and Polasky Diversity
// to be maximized
// SPD returns a value between 1 and n (0 if there is an error), 
// which can be interpreted as the number of species
double EMO_Indicator_solow_polasky(double *fit, double *f, int *filter, int n, int dim) {
  double *m, *c;
  double v, theta = 10.0;
  int i, j, k, l, size;
  int *filter2;

  filter2 = (int *) malloc(sizeof(int) * n);

  for(i = 0; i < n; i++)
    filter2[i] = 1;

  size = n;

  for(i = 0; i < n; i++) {
    for(j = i+1; j < n; j++) {
      v = EMO_vdist(f + i * dim, f + j * dim, dim);

      if(v == 0 && filter2[j]) {
        filter2[j] = 0;
        size--;
      }
    } 
  }


  m = (double *) malloc(sizeof(double) * size * size);
  c = (double *) malloc(sizeof(double) * size * size);

  for(i = 0, k = 0; i < n; i++) {

    if(!filter2[i]) {
      continue;
    }

    for(j = 0, l = 0; j < n; j++) {

      if(!filter2[j]) {
        continue;
      }

      v = (i == j)? 0.0 : EMO_vdist(f + i * dim, f + j * dim, dim);
      m[k*size + l] = exp(-theta * v);
      l++;
    }
    k++;
  }

  if(EMO_minverse(c, m, size, 1, 0, 0) == 1) {
    printf("Error to calculate the inverse\n");
    free(filter2);
    free(m);
    free(c);
    return 0.0;
  }

  v = 0;

  for(i = 0; i < size; i++)
    for(j = 0; j < size; j++)
      v += c[i*size + j];

  free(filter2);
  free(m);
  free(c);

  return v;
}

#endif

