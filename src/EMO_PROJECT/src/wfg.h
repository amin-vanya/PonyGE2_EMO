#ifndef _WFG_H
#define _WFG_H

#include "common.h"

double EMO_WFG_correct_to_01(double a);
double EMO_WFG_sh_linear_1(double *x, int M);
double EMO_WFG_sh_linear_m(double *x, int M, int m);
double EMO_WFG_sh_linear_M(double *x);
void EMO_WFG_sh_linear(double *h, double *x, int M);
double EMO_WFG_sh_convex_1(double *x, int M);
double EMO_WFG_sh_convex_m(double *x, int M, int m);
double EMO_WFG_sh_convex_M(double *x);
void EMO_WFG_sh_convex(double *h, double *x, int M);
double EMO_WFG_sh_concave_1(double *x, int M);
double EMO_WFG_sh_concave_m(double *x, int M, int m);
double EMO_WFG_sh_concave_M(double *x);
void EMO_WFG_sh_concave(double *h, double *x, int M);
double EMO_WFG_sh_mixed_M(double *x, double A, double alpha);
double EMO_WFG_sh_disc_M(double *x, double A, double alpha, double beta);
double EMO_WFG_tr_b_poly(double y, double alpha);
double EMO_WFG_tr_b_flat(double y, double A, double B, double C);
double EMO_WFG_tr_b_param(double y, double u, double A, double B, double C);
double EMO_WFG_tr_s_linear(double y, double A);
double EMO_WFG_tr_s_decept(double y, double A, double B, double C);
double EMO_WFG_tr_s_multi(double y, double A, double B, double C);
double EMO_WFG_tr_r_sum(double *y, double *w, int n);
double EMO_WFG_tr_r_nonsep(double *y, double A, int n);
void EMO_WFG_wfg1_t1(EMO_MOP *mop, double *z);
void EMO_WFG_wfg1_t2(EMO_MOP *mop, double *z);
void EMO_WFG_wfg1_t3(EMO_MOP *mop, double *z);
void EMO_WFG_wfg1_t4(EMO_MOP *mop, double *z);
void EMO_WFG_wfg2_t2(EMO_MOP *mop, double *z);
void EMO_WFG_wfg2_t3(EMO_MOP *mop, double *z);
void EMO_WFG_wfg4_t1(EMO_MOP *mop, double *z);
void EMO_WFG_wfg4_t2(EMO_MOP *mop, double *z);
void EMO_WFG_wfg5_t1(EMO_MOP *mop, double *z);
void EMO_WFG_wfg6_t2(EMO_MOP *mop, double *z);
void EMO_WFG_wfg7_t1(EMO_MOP *mop, double *z);
void EMO_WFG_wfg8_t1(EMO_MOP *mop, double *z);
void EMO_WFG_wfg9_t1(EMO_MOP *mop, double *z);
void EMO_WFG_wfg9_t2(EMO_MOP *mop, double *z);
void EMO_WFG_calc_x(double *x, double *t, double *A, int M);
void EMO_WFG_calc_f(double *f, double *x, double *S, double *h, int M);
void EMO_WFG_wfg1(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg2(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg3(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg4(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg5(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg6(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg7(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg8(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg9(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg1_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg2_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg3_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg4_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg5_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg6_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg7_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg8_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg9_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_range(EMO_MOP *mop);
int EMO_WFG_alloc(EMO_MOP *mop);
void EMO_WFG_free(EMO_MOP *mop);

#endif
