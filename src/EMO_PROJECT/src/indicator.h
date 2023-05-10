
#include "utility.h"

double EMO_Indicator_gd(double *f, int *filter, int n, double *pf, int m, int dim, double p);
double EMO_Indicator_igd(double *f, int *filter, int n, double *pf, int m, int dim, double p);
double EMO_Indicator_igd2(double *f, int *filter, int n, double *pf, int m, int dim);
double EMO_Indicator_igd_plus(double *f, int *filter, int n, double *pf, int m, int dim);
double EMO_Indicator_deltap(double *f, int *filter, int n, double *pf, int m, int dim, double p);
double EMO_Indicator_r2(double *data, int *filter, int size, double *W, int wsize, EMO_Utility *utl);
void EMO_Indicator_r2_ranking(double *rank, double **sort, double *norm, double *tmp, double *data, int size, double *W, int wsize, EMO_Utility *utl);
double EMO_Indicator_epsilon_multiplicative(double *a, int *filter, int na, double *b, int nb, int dim);
double EMO_Indicator_epsilon_additive(double *a, int *filter, int na, double *b, int nb, int dim);
double EMO_Indicator_maximin(double *fit, double *f, int *filter, int n, int dim);
double EMO_Indicator_onvg(double *f, int *filter, int n, int dim);
double EMO_Indicator_onvgr(double *a, int *filter, int na, double *b, int nb, int dim);
double EMO_Indicator_c(double *a, int *filter, int na, double *b, int nb, int dim);
double EMO_Indicator_spacing(double *f, int *filter, int n, int dim);
double EMO_Indicator_senergy(double *fit, double *f, int *filter, int n, int dim);
double EMO_Indicator_senergy_update(double *fit, double *f, int *filter, int n, int dim, int idx, int action);
double EMO_Indicator_solow_polasky(double *fit, double *f, int *filter, int n, int dim);

