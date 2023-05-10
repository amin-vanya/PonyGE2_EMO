
#ifndef _STAT_H
#define _STAT_H

double EMO_mean(double *data, int *filter, int size);
double EMO_quartile(double *q, double *data, int *filter, double **sort, int size);
double EMO_median(int *idx, double *data, int *filter, double **sort, int size);
double EMO_var(double *data, int *filter, double mean, int size);
double EMO_std(double *data, int *filter, double mean, int size);

#endif

