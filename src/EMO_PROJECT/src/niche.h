
#ifndef _NICHE_
#define _NICHE_

#include "list.h"
#include "refpoint.h"
#include "dominance.h"

typedef struct {
  EMO_List lst;
  EMO_NDSort nd;
  EMO_Refpoint ref;
  double *min;
  double *ideal;
  double *nadir;
  double **sort;    /* temporary array for sorting population */
  double *norm; 
  double *vnorm;    /* norm value for each individual */
  int *filter;
  int *tmp;
  int ssize;        /* temporal population size (for memory free) */
} EMO_Prune;

void EMO_Prune_alloc(EMO_Prune *p, int size, int nobj);
void EMO_Prune_free(EMO_Prune *p);

void EMO_nicheCount(double *niche, double *data, int row, int col, double alpha, double sigma_share);
void EMO_crowdingDistance(double *cd, double **sort, double *data, int row, double *max, double *min, int col);
void EMO_crowdingDistance2(double *cd, double **sort, double *data, EMO_List *front, double *max, double *min, int col);
void EMO_wdist(double *wd, double *data, int size, double *W, int wsize, double *tmp, int dim);


// Metodo de truncado de soluciones basado en SPEA2
typedef struct {
  EMO_List *lnn;          /* array of lists for storing neighbors */
  EMO_List copy;
  double **sort;          /* temporary array for sorting */
  double **dist;
  int *min, *max;         /* individuals that are part of the extreme points */
  int dim, ssize;
} EMO_KNN;

void EMO_KNN_alloc(EMO_KNN *knn, int size, int dim);
void EMO_KNN_free(EMO_KNN *knn);
void EMO_KNN_prune(EMO_KNN *knn, int *filter, int max_size, double *data, int size);

void EMO_knn(EMO_List *l, double **sort, double *W, int size, int dim, int k, int itself);
void EMO_knn2(EMO_List *l, EMO_List *copy, double **dist, double **sort, double *W, int size, int dim);
void EMO_kNN2(EMO_List *lst, double **sort, double *data, int row, int col);

#endif


