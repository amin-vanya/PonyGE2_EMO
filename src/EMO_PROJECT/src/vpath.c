
/* Value path */

#include <string.h>
#include <stdlib.h>
#include <float.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "vpath.h"
#include "niche.h"
#include "vector.h"
#include "numeric.h"
#include "dominance.h"
#include "refpoint.h"
#include "sort.h"

#define EMO_VPath_M 3
#define EMO_VPath_M2 (int) EMO_VPath_M / 2 

// mu: resolucion planeada para mu individuos, tamaño ideal
// max_mu: numero maximo de individuos (reserva memoria)
void EMO_VPath_alloc(EMO_VPath *vpath, int resolution, int mu, int max_mu, int dim) {
  int n, theta;

  if(dim < 2) {
    printf("Error, dimension must be greather or equal than two.\n");
    exit(1);
  }

  n = ((dim - 1) * dim) / 2;
  theta = resolution * mu;
  theta /= n;

  if(theta < 3) {
    printf("Error, not enough resolution.\n");
    printf("dim %d, n %d, resolution %d, mu %d, theta %d\n", dim, n, resolution, mu, theta);
    exit(1);
  }

  vpath->theta = theta;
  vpath->x = theta * n;
  vpath->y = theta;
  vpath->dim = dim;

  if((vpath->M = (int *) calloc(sizeof(int), vpath->x * vpath->y)) == NULL) {
    printf("Error, not enough memory in VPath.\n");
    exit(1);
  }

  if((vpath->c = (double *) calloc(sizeof(double), max_mu)) == NULL) {
    printf("Error, not enough memory in VPath.\n");
    exit(1);
  }
}

void EMO_VPath_free(EMO_VPath *vpath) {
  free(vpath->M);
  free(vpath->c);
}

void EMO_VPath_write(EMO_VPath *vpath, const char *str) {
  int i, j;
  FILE *fp;

  if(strcmp(str, "stdout") == 0) {
    fp = stdout;
    fflush(stdout);
  }
  else {
    if((fp = fopen(str, "w")) == NULL) {
      printf("Error to open file %s\n", str);
      exit(1);
    }
  }

  if(fprintf(fp, "# %d %d\n", vpath->y, vpath->x) == 0) {
    printf("Error to write file %s\n", str);
    exit(1);
  }
 
  for(i = vpath->y-1; i > -1; i--) {
    for(j = 0; j < vpath->x; j++) {
      if(fprintf(fp, "%.2f ", (double) vpath->M[i*vpath->x+j]) == 0) {
        printf("Error to write file %s\n", str);
        exit(1);
      }
    }
    fprintf(fp, "\n");
  }

  if(strcmp(str, "stdout") != 0) 
    fclose(fp); 
}

/* M: matrix of value path
 * data: population (objectives)
 * size: number of individuals
 * ref: reference line
 * dim: number of objectives
 * theta: [1,n] resolution (times the size of the population)
 */
int EMO_VPath_run(EMO_VPath *vpath, double *data, int *filter, int size, double z) {
  int i, j, objs, k, a, b, count, flag;
  double  v1, v2, m, t, y;

  memset(vpath->M, 0, sizeof(int) * vpath->x * vpath->y);

  for(i = 0; i < size; i++) {
    flag = 0;

    if(filter != NULL && filter[i] == 0) continue;

    for(j = 0; j < vpath->dim; j++) {
      if(data[i*vpath->dim + j] > z) {
        flag = 1;
        break;
      }
    }

    if(flag)
      continue;

    objs = vpath->dim - 1;
    count = 0;

    for(j = 0; j < objs; j++) {
      for(k = j+1; k <= objs; k++) {
        v1 = data[i * vpath->dim + j];
        v2 = data[i * vpath->dim + k];

        m = (v2 - v1) / (double) (vpath->theta - 1);

        for(b = 0; b < vpath->theta; b++) {
          y = m * ((double) b) + v1;
          t = (vpath->y - 1.0) * (1.0 - y / z);
          a = round(t);

          if(a >= 0 && a < vpath->y) { // && b >= 0 && b < vpath->x)
            vpath->M[a * vpath->x + (b + count * vpath->theta) ]++;
          }
        }
        count++;
      }
    }
  }

  count = 0;

  // arriba
  for(j = 0; j < vpath->x; j++) {
    i = 0;
 
    while(i < vpath->y && vpath->M[i*vpath->x + j] == 0) {
      vpath->M[i*vpath->x + j] = -1;
      i++;
      count++;
    }
  }
  return count;
}

// gecco16
void EMO_VPath_contribution(EMO_VPath *vpath, double *data, int *filter, int size, double z, EMO_List *lst) {
  int objs, count, cempty, cedge, cfilled, filled, flag;
  int i, j, k, a, b, bt, c, mx, my, s1, s2, e1, e2, ry;
  double  v1, v2, m, y, t, vmax, vmin;

  vmax = 0;
  vmin = DBL_MAX;

  for(i = 0; i < size; i++) {
    if(filter != NULL && filter[i] == 0) continue;
    vpath->c[i] = 0;
    flag = 0;

    for(j = 0; j < vpath->dim; j++) {
      if(data[i*vpath->dim + j] > z) {
        if(lst == NULL || (lst != NULL && EMO_List_seek(lst, i) == 0)) {
          vpath->c[i] = INT_MAX;
          flag = 1;
          break;
        }
      }

      if(data[i*vpath->dim + j] == z) {
        if(lst == NULL || (lst != NULL && EMO_List_seek(lst, i) == 1)) {
          vpath->c[i] = 0;
          flag = 1;
          break;
        }
      }
    }

    if(flag)
      continue;

    objs = vpath->dim - 1;
    count = 0;

    for(j = 0; j < objs; j++) {
      for(k = j+1; k <= objs; k++) {
        v1 = data[i * vpath->dim + j];
        v2 = data[i * vpath->dim + k];
        m = (v2 - v1) / (double) (vpath->theta - 1);

        for(bt = 0; bt < vpath->theta; bt++) {
          y = m * ((double) bt) + v1;
          t = (vpath->y - 1.0) * (1.0 - y / z);
          a = round(t);
          b = bt + count * vpath->theta;

          if(a >= 0 && a < vpath->y) {
            cedge = cempty = cfilled = filled = 0;
            ry = vpath->y - a;

            s1 = a - EMO_VPath_M2;
            e1 = a + EMO_VPath_M2;
            s2 = b - EMO_VPath_M2;
            e2 = b + EMO_VPath_M2;

            for(my = s1; my <= e1; my++) {
              for(mx = s2; mx <= e2; mx++) {

                if(mx >= 0 && mx < vpath->x && my >= 0 && my < vpath->y) {
                  c = vpath->M[my*vpath->x + mx];

                  if(c == -1)
                    cedge++;
                  else if(c == 0)
                    cempty++;
                  else
                    if(!(my == a && mx == b)) {
                      filled += c;
                      cfilled++;
                    }
                }
              }
            }

            vpath->c[i] += ry * vpath->M[a*vpath->x + b];
            if((cempty > 0 && cedge > 0) || (vpath->M[a*vpath->x + b] == 1 && cedge > 0))

              vpath->c[i] += (a + 1) * (cedge * cempty);
            else 
              vpath->c[i] -= (a + 1) * (cedge + cempty);

            if(cfilled > 0) {
              if(vpath->dim > 2 && (b % 2) == 0)
                vpath->c[i] -= ry * (cfilled);
              else
                vpath->c[i] += ry * ((double)(filled / ((double)cfilled)));
            }
          }
        }
        count++;
      }
    }

    if(vmax < vpath->c[i])
      vmax = vpath->c[i];

    if(vmin > vpath->c[i])
      vmin = vpath->c[i] - 1.0;
  }

  for(i = 0; i < size; i++) {
    if(vpath->c[i] == INT_MAX) {
      vpath->c[i] = EMO_vnorm(data + i *vpath->dim, 2.0, vpath->dim);
    }
    else if(vpath->c[i] != 0 && (vmax - vmin) != 0) {
      vpath->c[i] = (vpath->c[i] - vmin) / (vmax - vmin);
    }
  }
}

void EMO_VPath_update(EMO_VPath *vpath, double *v, double z, int val) {
  int i, j, k, a, b, bt, count, objs;
  double y, t, m;

  for(j = 0; j < vpath->dim; j++) {
    if(v[j] > z) {
      return;
    }
  }

  objs = vpath->dim - 1;
  count = 0;

  for(j = 0; j < objs; j++) {
    for(k = j+1; k <= objs; k++) {
      m = (v[k] - v[j]) / (double) (vpath->theta - 1);

      for(bt = 0; bt < vpath->theta; bt++) {
        y = m * ((double) bt) + v[j];
        t = (vpath->y - 1.0) * (1.0 - y / z);
        a = round(t);
        b = bt + count * vpath->theta;

        if(a >= 0 && a < vpath->y) {

          i = a * vpath->x + b;

          vpath->M[i] += val;

          if(val == -1 && vpath->M[i] == 0) {
            if(a == 0 || vpath->M[(a-1) * vpath->x + b] == -1) {
              vpath->M[i] = -1;

              while(++a < vpath->y && vpath->M[a * vpath->x + b] == 0)
                vpath->M[a * vpath->x + b] = -1;
            }
          }

          if(val == 1 && vpath->M[i] == 0) {
            vpath->M[i] = 1;

            while(++a < vpath->y && vpath->M[a * vpath->x + b] == -1)
              vpath->M[a * vpath->x + b] = 0;
          }
        }
      }
      count++;
    }
  }
}
/* Normalize population */
int EMO_VPath_normalize(EMO_Prune *p, double *data, int size, int nobj) {
  int i, j, k, imin, update = 0;
  static double dmin = 0;
  double *pt, v, vmin;

  EMO_List_clear(&p->lst);
  EMO_minBound(p->min, data, NULL, size, nobj);
  vmin = EMO_dmin(NULL, p->min, NULL, nobj);

  /* If objectives are negative, they are shifted to the origin */
  if(vmin < 0) {

    if(vmin != dmin) {
      dmin = vmin;
      update = 1;
    }

    for(j = 0; j < nobj; j++) {
      if(p->min[j] < 0) {
        v = fabs(p->min[j]);

        for(i = 0; i < size; i++) {
          k = i * nobj + j;
          p->norm[k] = data[k] + v;
        }

        p->min[j] = 0;
      }
      else {
        for(i = 0; i < size; i++) {
          k = i * nobj + j;
          p->norm[k] = data[k];
        }
      }
    }

    pt = p->norm;
  }
  else {
    pt = data;
  }

  /* Update ideal point */
  for(j = 0; j < nobj; j++) {
    if(p->min[j] < p->ideal[j]) {
      p->ideal[j] = p->min[j];
      update = 1;
    }
  }

  /* Calculate the norm of each solution */
  for(i = 0; i < size; i++) {
    p->vnorm[i] = EMO_vnorm(pt + i * nobj, 2.0, nobj);
  }
 

  /* Update nadir point looking for the individuals parallel to the axis with the lowest norm.
     Removing dominance-resistant points (outliers) */
  for(j = 0; j < nobj; j++) {

    for(i = 0; i < size; i++) {
      p->sort[i][0] = (double) i;
      p->sort[i][1] = pt[i * nobj + j] / p->vnorm[i];
    }

    qsort(p->sort, size, sizeof(p->sort[0]), (int (*)(const void *, const void *))&EMO_compare_desc);

    i = 0;
    imin = -1;
    vmin = DBL_MAX;

    do {
      k = (int) p->sort[i][0];

      v = p->vnorm[k];

      if(v < vmin && pt[k * nobj + j] > p->ideal[j]) {
        vmin = v;
        imin = k;
      }
      i++;

    } while(i < size && (p->sort[0][1] - p->sort[i][1]) <= 1e-8); //1e-12);

    if(imin != -1) {
      k = imin * nobj + j;

      EMO_List_queue(&p->lst, imin);

      if(p->nadir[j] != pt[k]) {
        p->nadir[j] = pt[k];
        update = 1;
      }
    }
  }

  EMO_normalize(p->norm, pt, NULL, size, p->ideal, p->nadir, nobj);
  return update;
}

void EMO_VPath_prune2(EMO_VPath *vpath, EMO_Prune *p, double *data, int size, int new_size) {
  int i, j, n, elem;

  if(size <= new_size) {
    for(i = 0; i < size; i++)
      p->filter[i] = 1; 
    return;
  }

  EMO_NDSort_run(&p->nd, data, vpath->dim, NULL, NULL, size);  

  // regresa update
  EMO_VPath_normalize(p, data, size, vpath->dim);

  n = size - new_size;    /* number of elements to be removed */
  j = p->nd.nfront - 1;   /* last front */

  /* Select all members in the population */
  for(i = 0; i < size; i++)
    p->filter[i] = 1; 
  
  while(n > 0) {
    /* Remove the last front if its size is less or equal than the number of 
     * individuals to be removed */
    if(p->nd.front[j].size <= n) {
      for(i = 0; i < p->nd.front[j].size; i++) {
        EMO_List_get(&p->nd.front[j], &elem, i); 
        p->filter[elem] = 0;
        n--;
      }
    }
    else { /* select the worst distributed individuals from the current front */

      /* Select solutions from the last front */
      memset(p->tmp, 0, sizeof(int) * size);

      for(i = 0; i < p->nd.front[j].size; i++) {
        EMO_List_get(&p->nd.front[j], &elem, i);
        p->tmp[elem] = 1;
      }

      EMO_VPath_run(vpath, p->norm, p->filter, size, 1.0);

      while(n > 0) {
        EMO_VPath_contribution(vpath, p->norm, p->filter, size, 1.0, NULL);

        /* Find worst contribution */
        EMO_dmax(&elem, vpath->c, p->tmp, size);
        p->filter[elem] = p->tmp[elem] = 0;
        n--;

        EMO_VPath_update(vpath, p->norm + elem * vpath->dim, 1.0, -1);      // elimina individuo anterior
      }
    }

    j--;
  }
}

int EMO_VPath_prune(EMO_VPath *vpath, EMO_Prune *p, double *data, int size, int new_size) {
  int i, j, t, n;
  int tot = 0;

  t = EMO_Dominance_ndset(p->filter, data, NULL, size, vpath->dim, EMO_Dominance_strict);

  EMO_Refpoint_update_ideal(&p->ref, data, NULL, size);
  EMO_Refpoint_update_nadir(&p->ref, data, p->filter, size);

  EMO_normalize(p->norm, data, NULL, size, p->ref.ideal, p->ref.nadir, vpath->dim);

  // Hay mas individuos dominados, los descarta por diversidad
  if(t <= new_size) {
    tot = t;
  }
  else {
    n = t - new_size;

    for(i = 0; i < n; i++) {
      EMO_VPath_run(vpath, p->norm, p->filter, size, 1.0);
      EMO_VPath_contribution(vpath, p->norm, p->filter, size, 1.0, NULL);
      EMO_dmax(&j, vpath->c, p->filter, size);
      
      p->filter[j] = 0;
    }

    tot = new_size;
  }
  return tot;
}

