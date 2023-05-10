#ifndef _VPATH_H
#define _VPATH_H

#include "list.h"
#include "niche.h"

typedef struct {
  int *M;  /* matrix of value path */
  double *c;
  int x, y;
  int theta;
  int dim;
} EMO_VPath;

void EMO_VPath_alloc(EMO_VPath *vpath, int resolution, int mu, int max_mu, int dim);
void EMO_VPath_free(EMO_VPath *vpath);
void EMO_VPath_write(EMO_VPath *vpath, const char *str);
int EMO_VPath_run(EMO_VPath *vpath, double *data, int *filter, int size, double z);
void EMO_VPath_update(EMO_VPath *vpath, double *v, double z, int val);
void EMO_VPath_contribution(EMO_VPath *vpath, double *data, int *filter, int size, double z, EMO_List *lst);

int EMO_VPath_prune(EMO_VPath *vpath, EMO_Prune *p, double *data, int size, int new_size);

#endif

