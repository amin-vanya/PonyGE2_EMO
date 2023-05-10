/**************************************************************
 * utility.c    Definition of utility functions.              *
 *                                                            *
 * Computer Science Department, CINVESTAV-IPN                 *
 *                                                            *
 * Professor:   Dr. Carlos A. Coello Coello                   *
 * Author:      Raquel Hernandez Gomez                        *
 *                                                            *
 * March 2013                                                 *
 *************************************************************/
#ifndef _UTILITY_H
#define _UTILITY_H

#include "param.h"

typedef struct EMO_Utility {
  char *name;
  int nobj;
  int inverted;            /* inverted scalarizing function */
  int tlpbi_H;
  int qpbi_H;
  double wzero;            /* Small value of wi, used to avoid division by zero */
  double wcp_p;
  double wcp2_p;
  double wpo_p;
  double wpo2_p;
  double ewc_p;
  double ewc2_p;
  double wn_p;
  double wn2_p;
  double vads_p;
  double vads2_p;
  double ache_alpha;
  double ache2_alpha;
  double mche_alpha;
  double mche2_alpha;
  double aasf_alpha;
  double tlpbi_alpha;
  double qpbi_alpha;
  double gsf_alpha;
  double gsf2_alpha;
  double nsf_alpha;
  double nsf2_alpha;
  double cs_alpha;
  double cs2_alpha;
  double pbi_theta;
  double qpbi_theta;
  double tlpbi_theta1;
  double tlpbi_theta2;
  double gsf_beta;
  double gsf2_beta;
  double didass_beta;
  double didass2_beta;
  double tlpbi_dstar;
  double qpbi_dstar;
  double *pbi_v;      /* temporary vectors */
  double *tlpbi_v;
  double *qpbi_v;
  double *didass_v;
  double *didass2_v;
  double (*uf)(struct EMO_Utility *u, double *w, double *x);
} EMO_Utility;

typedef double (*EMO_UtilityFunction)(EMO_Utility *u, double *w, double *x);

void EMO_Utility_alloc(EMO_Utility *u, EMO_Param *param, int nobj, const char *str);
void EMO_Utility_free(EMO_Utility *u);

double EMO_Utility_wcp(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wcp2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wpo(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wpo2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ewc(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ewc2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wpr(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wpr2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ws(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ws2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wn(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wn2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ls(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ls2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_che(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ache(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ache2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_mche(EMO_Utility *u, double *w, double *x);
double EMO_Utility_mche2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_asf(EMO_Utility *u, double *w, double *x);
double EMO_Utility_aasf(EMO_Utility *u, double *w, double *x);
double EMO_Utility_pbi(EMO_Utility *u, double *w, double *x);
void EMO_Utility_update_dstar(EMO_Utility *u, double *zmin, double *zmax);
double EMO_Utility_tlpbi(EMO_Utility *u, double *w, double *x);
double EMO_Utility_qpbi(EMO_Utility *u, double *w, double *x);
double EMO_Utility_gsf(EMO_Utility *u, double *w, double *x);
double EMO_Utility_gsf2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_nsf(EMO_Utility *u, double *w, double *x);
double EMO_Utility_nsf2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_cs(EMO_Utility *u, double *w, double *x);
double EMO_Utility_cs2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_vads(EMO_Utility *u, double *w, double *x);
double EMO_Utility_vads2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_didass(EMO_Utility *u, double *w, double *x);
double EMO_Utility_didass2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_refsf(EMO_Utility *u, double *w, double *x);

#endif

