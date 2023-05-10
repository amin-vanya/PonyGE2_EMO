/**************************************************************
 * utility.c    Definition of utility functions.              *
 *                                                            *
 * Computer Science Department, CINVESTAV-IPN                 *
 *                                                            *
 * Professor:   Dr. Carlos A. Coello Coello                   *
 *                                                            *
 * Authors:     Miriam Pescador Rojas                         *
 *              Raquel Hernandez Gomez                        *
 *                                                            *
 * June 2016                                                  *
 *************************************************************/
// In all cases
// x[i] := f[i] - z^min[i]  or
// x[i] := (f[i] - z^min[i]) / (z^max[i] - z^min[i])  normalized

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "common.h"
#include "utility.h"
#include "vector.h"
#include "numeric.h"

const char *EMO_Utility_list[] = { "weighted_compromise_programming",
                                   "weighted_compromise_programming2",
                                   "weighted_power",
                                   "weighted_power2",
                                   "exponential_weighted_criteria",
                                   "exponential_weighted_criteria2",
                                   "weighted_product",
                                   "weighted_product2",
                                   "weighted_sum",
                                   "weighted_sum2",
                                   "weighted_norm",
                                   "weighted_norm2",
                                   "least_squares",
                                   "least_squares2",
                                   "chebyshev",
                                   "augmented_chebyshev",
                                   "augmented_chebyshev2",
                                   "modified_chebyshev",
                                   "modified_chebyshev2",
                                   "achievement_scalarizing_function",
                                   "augmented_achievement_scalarizing_function",
                                   "penalty_boundary_intersection",
                                   "two_level_penalty_boundary_intersection",
                                   "quadratic_penalty_boundary_intersection",
                                   "general_scalarizing_function",
                                   "general_scalarizing_function2",
                                   "normalized_scalarizing_function",
                                   "normalized_scalarizing_function2",
                                   "conic_scalarization",
                                   "conic_scalarization2",
                                   "vector_angle_distance_scaling",
                                   "vector_angle_distance_scaling2",
                                   "didass",
                                   "didass2",
                                   "refsf", // RHG
                                    NULL
                                 };

const EMO_UtilityFunction fdic[] = { EMO_Utility_wcp,
                                     EMO_Utility_wcp2,
                                     EMO_Utility_wpo,
                                     EMO_Utility_wpo2,
                                     EMO_Utility_ewc,
                                     EMO_Utility_ewc2,
                                     EMO_Utility_wpr,
                                     EMO_Utility_wpr2,
                                     EMO_Utility_ws,
                                     EMO_Utility_ws2,
                                     EMO_Utility_wn,
                                     EMO_Utility_wn2,
                                     EMO_Utility_ls,
                                     EMO_Utility_ls2,
                                     EMO_Utility_che,
                                     EMO_Utility_ache,
                                     EMO_Utility_ache2,
                                     EMO_Utility_mche,
                                     EMO_Utility_mche2,
                                     EMO_Utility_asf,
                                     EMO_Utility_aasf,
                                     EMO_Utility_pbi,
                                     EMO_Utility_tlpbi,
                                     EMO_Utility_qpbi,
                                     EMO_Utility_gsf,
                                     EMO_Utility_gsf2,
                                     EMO_Utility_nsf,
                                     EMO_Utility_nsf2,
                                     EMO_Utility_cs,
                                     EMO_Utility_cs2,
                                     EMO_Utility_vads,
                                     EMO_Utility_vads2,
                                     EMO_Utility_didass,
                                     EMO_Utility_didass2,
                                     EMO_Utility_refsf
                                   };

void EMO_Utility_alloc(EMO_Utility *u, EMO_Param *param, int nobj, const char *str) {
  int i, d, all = 0;
  double v;

  u->nobj = nobj;

  /* Name of function is converted to lower case */
  u->name = EMO_tolower(str);
  u->uf = NULL;
  u->pbi_v     = NULL;
  u->tlpbi_v   = NULL;
  u->qpbi_v    = NULL;
  u->didass_v  = NULL;
  u->didass2_v = NULL;

  if(strcmp(u->name, "all") == 0 || strcmp(u->name, "*") == 0) {  // Hyper-heuristic approach
    all = 1;
  }
  else {
    /* Find the function's name in dictionary */
    i = EMO_Dictionary_find(EMO_Utility_list, u->name);

    if(i == -1) {
      printf("Error, unknown %s in EMO_Utility_alloc. Available functions:\n\n", u->name);
      EMO_Dictionary_print(stdout, EMO_Utility_list, "UTILITY");
      free(u->name);
      exit(1);
    }

    u->uf = fdic[i];

    if(u->uf == NULL) {
      printf("Error, unknown utility (or scalarization) function %s.\n", u->name);
      free(u->name);
      exit(1);
    }
  }

  if(all || strcmp(u->name, "penalty_boundary_intersection") == 0) {
    if((u->pbi_v = (double *) malloc(sizeof(double) * nobj)) == NULL) {
      printf("Not enough memory in utility.c for pbi_v.\n");
      free(u->name);
      exit(1);
    }
    memset(u->pbi_v, 0, sizeof(double) * nobj);
  }

  if(all || strcmp(u->name, "two_level_penalty_boundary_intersection") == 0) {
    if((u->tlpbi_v = (double *) malloc(sizeof(double) * nobj)) == NULL) {
      printf("Not enough memory in utility.c for tlpbi_v.\n");
      free(u->name);
      exit(1);
    }
    memset(u->tlpbi_v, 0, sizeof(double) * nobj);
  }

  if(all || strcmp(u->name, "quadratic_penalty_boundary_intersection") == 0) {
    if((u->qpbi_v = (double *) malloc(sizeof(double) * nobj)) == NULL) {
      printf("Not enough memory in utility.c for qpbi_v.\n");
      free(u->name);
      exit(1);
    }
    memset(u->qpbi_v, 0, sizeof(double) * nobj);
  }

  if(all || strcmp(u->name, "didass") == 0) {
    if((u->didass_v = (double *) malloc(sizeof(double) * nobj)) == NULL) {
      printf("Not enough memory in utility.c for didass_v.\n");
      free(u->name);
      exit(1);
    }

    for(d = 0; d < nobj; d++)
      u->didass_v[d] = 0.1;
  }

  if(all || strcmp(u->name, "didass2") == 0) {
    if((u->didass2_v = (double *) malloc(sizeof(double) * nobj)) == NULL) {
      printf("Not enough memory in utility.c for didass_v.\n");
      free(u->name);
      exit(1);
    }

    for(d = 0; d < nobj; d++)
      u->didass2_v[d] = 0.1;
  }

  /* Default values */
  u->inverted = 0;
  u->wcp_p = 10;
  u->wcp2_p = 10;
  u->wpo_p = 2;
  u->wpo2_p = 2;
  u->ewc_p = 20;
  u->ewc2_p = 20;
  u->wn_p = 1.5;
  u->wn2_p = 1.5;
  u->vads_p = 200;
  u->vads2_p = 200;

  u->ache_alpha = 0.1;
  u->ache2_alpha = 0.1;
  u->mche_alpha = 1e-2;
  u->mche2_alpha = 1e-2;
  u->aasf_alpha = 1e-4;
  u->tlpbi_alpha = 0.8;
  u->qpbi_alpha = 0.2;
  u->gsf_alpha = 0.01;
  u->nsf_alpha = 0.01;
  u->nsf2_alpha = 0.01;
  u->cs_alpha = 0.02;
  u->cs2_alpha = 1e-4;

  u->pbi_theta = 5.0;
  u->tlpbi_theta1 = 0.1;
  u->tlpbi_theta2 = 10.0;
  u->qpbi_theta = 5.0;

  u->tlpbi_dstar = 1.0;
  u->qpbi_dstar = 1.0;
  u->tlpbi_H = 1.0;
  u->qpbi_H = 1.0;

  u->gsf_beta = 0.1;
  u->didass_beta = 0.1;
  u->didass2_beta = 0.1;

  u->wzero = (u->nobj == 2)? 1e-6 : 1e-2;

  if(all == 0) {
    if(param == NULL || !EMO_Param_get_int(param, &u->inverted, "utility_inverted"))
        printf("Error, param = NULL or utility_inverted is not defined in the configuration file, setting the default value %d\n", u->inverted);

    if(u->inverted)
      printf("Inverted utility: %s\n", u->name);
    else
      printf("Utility: %s\n", u->name);

    if(u->uf == EMO_Utility_wcp || u->uf == EMO_Utility_wcp2 ||
       u->uf == EMO_Utility_wpo || u->uf == EMO_Utility_wpo2 ||
       u->uf == EMO_Utility_ewc || u->uf == EMO_Utility_ewc2 ||
       u->uf == EMO_Utility_wn || u->uf == EMO_Utility_wn2 ||
       u->uf == EMO_Utility_vads || u->uf == EMO_Utility_vads2) {

      if(param == NULL || !EMO_Param_get_double(param, &v, "utility_p")) {
        printf("Warning, param = NULL or utility_p is not defined in the configuration file, setting the default value ");

        if(u->uf == EMO_Utility_wcp)
          printf("%f\n", u->wcp_p);
        else if(u->uf == EMO_Utility_wcp2)
          printf("%f\n", u->wcp2_p);
        else if(u->uf == EMO_Utility_wpo)
          printf("%f\n", u->wpo_p);
        else if(u->uf == EMO_Utility_wpo2)
          printf("%f\n", u->wpo2_p);
        else if(u->uf == EMO_Utility_ewc)
          printf("%f\n", u->ewc_p);
        else if(u->uf == EMO_Utility_ewc2)
          printf("%f\n", u->ewc2_p);
        else if(u->uf == EMO_Utility_wn)
          printf("%f\n", u->wn_p);
        else if(u->uf == EMO_Utility_wn2)
          printf("%f\n", u->wn2_p);
        else if(u->uf == EMO_Utility_vads)
          printf("%f\n", u->vads_p);
        else if(u->uf == EMO_Utility_vads2)
          printf("%f\n", u->vads2_p);
      }
      else {
        if(u->uf == EMO_Utility_wcp)
          u->wcp_p = v;
        else if(u->uf == EMO_Utility_wcp2)
          u->wcp2_p = v;
        else if(u->uf == EMO_Utility_wpo)
          u->wpo_p = v;
        else if(u->uf == EMO_Utility_wpo2)
          u->wpo2_p = v;
        else if(u->uf == EMO_Utility_ewc)
          u->ewc_p = v;
        else if(u->uf == EMO_Utility_ewc2)
          u->ewc2_p = v;
        else if(u->uf == EMO_Utility_wn)
          u->wn_p = v;
        else if(u->uf == EMO_Utility_wn2)
          u->wn2_p = v;
        else if(u->uf == EMO_Utility_vads)
          u->vads_p = v;
        else if(u->uf == EMO_Utility_vads2)
          u->vads2_p = v;
      }
    }

    if(u->uf == EMO_Utility_ache || u->uf == EMO_Utility_ache2 ||
       u->uf == EMO_Utility_mche || u->uf == EMO_Utility_mche2 ||
       u->uf == EMO_Utility_aasf ||
       u->uf == EMO_Utility_tlpbi || u->uf == EMO_Utility_qpbi ||
       u->uf == EMO_Utility_gsf || u->uf == EMO_Utility_gsf2 ||
       u->uf == EMO_Utility_nsf || u->uf == EMO_Utility_nsf2 ||
       u->uf == EMO_Utility_cs || u->uf == EMO_Utility_cs2 ) {

      if(param == NULL || !EMO_Param_get_double(param, &v, "utility_alpha")) {
        printf("Error, param = NULL or utility_alpha is not defined in the configuration file, setting the default value ");

        if(u->uf == EMO_Utility_ache)
          printf("%f\n", u->ache_alpha);
        else if(u->uf == EMO_Utility_ache2)
          printf("%f\n", u->ache2_alpha);
        else if(u->uf == EMO_Utility_mche)
          printf("%f\n", u->mche_alpha);
        else if(u->uf == EMO_Utility_mche2)
          printf("%f\n", u->mche2_alpha);
        else if(u->uf == EMO_Utility_aasf)
          printf("%f\n", u->aasf_alpha);
        else if(u->uf == EMO_Utility_tlpbi)
          printf("%f\n", u->tlpbi_alpha);
        else if(u->uf == EMO_Utility_qpbi)
          printf("%f\n", u->qpbi_alpha);
        else if(u->uf == EMO_Utility_gsf)
          printf("%f\n", u->gsf_alpha);
        else if(u->uf == EMO_Utility_gsf2)
          printf("%f\n", u->gsf2_alpha);
        else if(u->uf == EMO_Utility_nsf)
          printf("%f\n", u->nsf_alpha);
        else if(u->uf == EMO_Utility_nsf2)
          printf("%f\n", u->nsf2_alpha);
        else if(u->uf == EMO_Utility_cs)
          printf("%f\n", u->cs_alpha);
        else if(u->uf == EMO_Utility_cs2)
          printf("%f\n", u->cs2_alpha);
      }
      else {
        if(u->uf == EMO_Utility_ache)
          u->ache_alpha = v;
        else if(u->uf == EMO_Utility_ache2)
          u->ache2_alpha = v;
        else if(u->uf == EMO_Utility_mche)
          u->mche_alpha = v;
        else if(u->uf == EMO_Utility_mche2)
          u->mche2_alpha = v;
        else if(u->uf == EMO_Utility_aasf)
          u->aasf_alpha = v;
        else if(u->uf == EMO_Utility_tlpbi)
          u->tlpbi_alpha = v;
        else if(u->uf == EMO_Utility_qpbi)
          u->qpbi_alpha =v;
        else if(u->uf == EMO_Utility_gsf)
          u->gsf_alpha = v;
        else if(u->uf == EMO_Utility_gsf2)
          u->gsf2_alpha = v;
        else if(u->uf == EMO_Utility_nsf)
          u->nsf_alpha = v;
        else if(u->uf == EMO_Utility_nsf2)
          u->nsf2_alpha = v;
        else if(u->uf == EMO_Utility_cs)
          u->cs_alpha = v;
        else if(u->uf == EMO_Utility_cs2)
          u->cs2_alpha = v;
      }
    }

    if(u->uf == EMO_Utility_pbi || u->uf == EMO_Utility_qpbi) {
      if(param == NULL || !EMO_Param_get_double(param, &v, "utility_theta")) {
        printf("Error, param = NULL or utility_theta is not defined in the configuration file, setting the default value ");

        if(u->uf == EMO_Utility_pbi)
          printf("%f\n", u->pbi_theta);
        else if(u->uf == EMO_Utility_qpbi)
          printf("%f\n", u->qpbi_theta);
      }
      else {
        if(u->uf == EMO_Utility_pbi)
          u->pbi_theta = v;
        else if(u->uf == EMO_Utility_qpbi)
          u->qpbi_theta = v;
      }
    }

    if(u->uf == EMO_Utility_tlpbi) {
      if(param == NULL || !EMO_Param_get_double(param, &u->tlpbi_theta1, "utility_theta1"))
        printf("Error, param = NULL or utility_theta1 is not defined in the configuration file, setting the default value %f\n", u->tlpbi_theta1);

      if(param == NULL || !EMO_Param_get_double(param, &u->tlpbi_theta2, "utility_theta2"))
        printf("Error, param = NULL or utility_theta2 is not defined in the configuration file, setting the default value %f\n", u->tlpbi_theta2);
    }

    if(u->uf == EMO_Utility_tlpbi || u->uf == EMO_Utility_qpbi) {
      if(param == NULL || !EMO_Param_get_int(param, &d, "utility_H")) {
        printf("Error, param = NULL or utility_H is not defined in the configuration file, setting the default value ");

        if(u->uf == EMO_Utility_tlpbi)
          printf("%d\n", u->tlpbi_H);
        else if(u->uf == EMO_Utility_qpbi)
          printf("%d\n", u->qpbi_H);
      }
      else {
        if(u->uf == EMO_Utility_tlpbi)
          u->tlpbi_H = d;
        else if(u->uf == EMO_Utility_qpbi)
          u->qpbi_H = d;
      }
    }

    if(u->uf == EMO_Utility_gsf || u->uf == EMO_Utility_gsf2 ||
       u->uf == EMO_Utility_didass || u->uf == EMO_Utility_didass2) {

      if(param == NULL || !EMO_Param_get_double(param, &v, "utility_beta")) {
        printf("Error, param = NULL or utility_beta is not defined in the configuration file, setting the default value %f\n", v);

        if(u->uf == EMO_Utility_gsf)
          printf("%f\n", u->gsf_beta);
        else if(u->uf == EMO_Utility_gsf2)
          printf("%f\n", u->gsf2_beta);
        else if(u->uf == EMO_Utility_didass)
          printf("%f\n", u->didass_beta);
        else if(u->uf == EMO_Utility_didass2)
          printf("%f\n", u->didass2_beta);
      }
      else {
        if(u->uf == EMO_Utility_gsf)
          u->gsf_beta = v;
        else if(u->uf == EMO_Utility_gsf2)
          u->gsf2_beta = v;
        else if(u->uf == EMO_Utility_didass)
          u->didass_beta = v;
        else if(u->uf == EMO_Utility_didass2)
          u->didass2_beta = v;
      }

      if(u->uf == EMO_Utility_didass) {
        if(param == NULL || !EMO_Param_get_vector_double(param, u->didass_v, &nobj, "utility_gamma")) {
          printf("Error, param = NULL or utility_gamma is not defined in the configuration file, setting the default value ");
          EMO_vprint(stdout, u->didass_v, nobj, NULL);
        }
        else if(nobj != u->nobj) {
          printf("Error, mismatch dimensions of utility_gama and nobj (%d vs %d) in configuration file.\n", nobj, u->nobj);
          exit(1);
        }
      }

      if(u->uf == EMO_Utility_didass2) {
        if(param == NULL || !EMO_Param_get_vector_double(param, u->didass2_v, &nobj, "utility_gamma")) {
          printf("Error, param = NULL or utility_gamma is not defined in the configuration file, setting the default value ");
          EMO_vprint(stdout, u->didass2_v, nobj, NULL);
        }
        else if(nobj != u->nobj) {
          printf("Error, mismatch dimensions of utility_gama and nobj (%d vs %d) in configuration file.\n", nobj, u->nobj);
          exit(1);
        }
      }
    }
  }
}

void EMO_Utility_free(EMO_Utility *u) {
  free(u->name);

  if(u->pbi_v != NULL)
    free(u->pbi_v);
  if(u->tlpbi_v != NULL)
    free(u->tlpbi_v);
  if(u->qpbi_v != NULL)
    free(u->qpbi_v);
  if(u->didass_v != NULL)
    free(u->didass_v);
  if(u->didass2_v != NULL)
    free(u->didass2_v);
}

/* Weighted Compromising Programming (WCP)
 * Athan and Papalambros 1996, Lee14
 * p = {3, 5, 7, 9, 11}
 */
double EMO_Utility_wcp(EMO_Utility *u, double *w, double *x) {
  double t = 0;
  int i;

  for(i = u->nobj - 1; i > -1; i--)
    t += pow(w[i] * x[i], u->wcp_p);

  return t;
}

/* Modified utility function, where the weight vector is a divisor */
double EMO_Utility_wcp2(EMO_Utility *u, double *w, double *x) {
  double t, y;
  int i;

  t = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i]? w[i] : u->wzero;
    t += pow(x[i] / y, u->wcp2_p);
  }

  return t;
}

// Weighted t-th Power Approach
// Li[149] Figueira, Greco Ehrgott, Multiple... p. 677
double EMO_Utility_wpo(EMO_Utility *u, double *w, double *x) {
  double t = 0;
  int i;

  for(i = u->nobj - 1; i > -1; i--)
    t += w[i] * pow(x[i], u->wpo_p);

  return t;
}

/* Modified utility function, where the weight vector is a divisor */
double EMO_Utility_wpo2(EMO_Utility *u, double *w, double *x) {
  double t, y;
  int i;

  t = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i]? w[i] : u->wzero;
    t += pow(x[i], u->wpo2_p) / y;
  }

  return t;
}

/*Athan96(1996)  Marler04 */
double EMO_Utility_ewc(EMO_Utility *u, double *w, double *x) {
  double s, y;
  int i;

  s = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i] ? w[i] : u->wzero;
    s += (exp(u->ewc_p * y) - 1.0) * (exp(u->ewc_p * x[i]));
  }

  //printf ("EWC value: %lf\n", s);
  return s;
}

double EMO_Utility_ewc2(EMO_Utility *u, double *w, double *x) {
  double s, y;
  int i;

  s = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i] ? w[i] : u->wzero;
    s += (exp(u->ewc2_p / y) - 1.0) * (exp(u->ewc2_p * x[i]));
  }

  return s;
}

// Odu13
// Bridgeman, 1922
double EMO_Utility_wpr(EMO_Utility *u, double *w, double *x) {
  double t = 1.0;
  int i;

  for(i = u->nobj - 1; i > -1; i--)
    t *= pow(x[i], w[i]);

  return t;
}

/* Modified utility function, where the weight vector is a divisor */
double EMO_Utility_wpr2(EMO_Utility *u, double *w, double *x) {
  double t, y;
  int i;

  t = 1.0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i]? w[i] : u->wzero;
    t *= pow(x[i], 1.0 / y);
  }

  return t;
}

// Gass and Saaty (1995), Zadeh (1963)
// Miettinen99, p 78
// wi >= 0, sum wi = 1
// x = fi-z*
double EMO_Utility_ws(EMO_Utility *u, double *w, double *x) {
  double t = 0;
  int i;

  for(i = u->nobj - 1; i > -1; i--)
    t += w[i] * x[i];
  return t;
}

/* Modified utility function, where the weight vector is a divisor */
double EMO_Utility_ws2(EMO_Utility *u, double *w, double *x) {
  double t, y;
  int i;

  t = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i]? w[i] : u->wzero;
    t += x[i] / y;
  }
  return t;
}

/* Zeleny (1973), Miettinen99, p. 97
 * wi = 0, sum wi = 1
 * 1 <= p < infty
 * ideal point
 */
double EMO_Utility_wn(EMO_Utility *u, double *w, double *x) {
  double t = 0;
  int i;

  for(i = u->nobj - 1; i > -1; i--)
    t += w[i] * pow(fabs(x[i]), u->wn_p);

  return pow(t, 1.0 / u->wn_p);
}

/* Modified utility function, where the weight vector is a divisor */
double EMO_Utility_wn2(EMO_Utility *u, double *w, double *x) {
  double t, y;
  int i;

  t = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i]? w[i] : u->wzero;
    t += pow(fabs(x[i]), u->wn2_p) / y;
  }

  return pow(t, 1.0 / u->wn2_p);
}

// special case of weighted Lp metrics p = 2
// Miettinen99, Sundar05 (euclidean distance measure)
double EMO_Utility_ls(EMO_Utility *u, double *w, double *x) {
  double t = 0;
  int i;

  for(i = u->nobj - 1; i > -1; i--)
    t += w[i] * (x[i] * x[i]);

  return sqrt(t);
}

/* Modified utility function, where the weight vector is a divisor */
double EMO_Utility_ls2(EMO_Utility *u, double *w, double *x) {
  double t, y;
  int i;

  t = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i]? w[i] : u->wzero;
    t += (x[i] * x[i]) / y;
  }
  return sqrt(t);
}

/* Bowman (1976), Miettinen99, p. 97
 * wi = 0, sum wi = 1
 * 1 <= p < infty
 * ideal point
 */
double EMO_Utility_che(EMO_Utility *u, double *w, double *x) {
  double y, v, vmax = 0;
  int i;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i] ? w[i] : 1e-6;
    v = fabs(x[i]) * y;

    if(v > vmax)
      vmax = v;
  }
  return vmax;
}

// Augmented weighted Tchebycheff metric
// utopian reference point
// Steuer 1986, Steuer and Choo (1983)
// Miettinen99, p. 101
double EMO_Utility_ache(EMO_Utility *u, double *w, double *x) {
  double t, s, vmax;
  int i;

  s = vmax = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    t = w[i] * fabs(x[i]);
    s += fabs(x[i]);

    if(t > vmax)
      vmax = t;
  }

  return vmax + u->ache_alpha * s;
}

/* Modified utility function, where the weight vector is a divisor */
double EMO_Utility_ache2(EMO_Utility *u, double *w, double *x) {
  double t, s, y, vmax;
  int i;

  s = vmax = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i]? w[i] : u->wzero;
    t = fabs(x[i]) / y;
    s += fabs(x[i]);

    if(t > vmax)
      vmax = t;
  }

  return vmax + u->ache2_alpha * s;
}

// Modified weighted Tchebycheff metric
// utopian reference point
// Miettinen99, p. 101
double EMO_Utility_mche(EMO_Utility *u, double *w, double *x) {
  double t, s, vmax;
  int i;

  vmax = s = 0;

  for(i = u->nobj - 1; i > -1; i--)
    s += fabs(x[i]);

  for(i = u->nobj - 1; i > -1; i--) {
    t = w[i] * (fabs(x[i]) + u->mche_alpha * s);

    if(t > vmax)
      vmax = t;
  }

  return vmax;
}

/* Modified utility function, where the weight vector is a divisor */
double EMO_Utility_mche2(EMO_Utility *u, double *w, double *x) {
  double t, s, y, vmax;
  int i;

  vmax = s = 0;

  for(i = u->nobj - 1; i > -1; i--)
    s += fabs(x[i]);

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i]? w[i] : u->wzero;
    t = (fabs(x[i]) + u->mche2_alpha * s) / y;

    if(t > vmax)
      vmax = t;
  }

  return vmax;
}

double use(double y){
    return 0;
}

// Achievement scalarizing function
// Wierzbicki (1981, 1982, 1986a, b, 1977, 1980a, b)
// Jahn (1984) and Luc (1986)
// Miettinen99, p. 107
// w is some fixed positive weighting vector
double EMO_Utility_asf(EMO_Utility *u, double *w, double *x) {
  double y, v, vmax, x_sum, x_avg, x_max, x_min;
  int i;
  x_min = 1e+8;
  x_max = -1e+8;
  x_sum = x_avg = vmax = 0;
  y = 1;
  use(y); // variable 'y' may not be used in a given scalarizing function

  // Initialize additional variables
  for(i=u->nobj-1; i>-1; i--){
      x_sum += x[i];
      if(x[i] > x_max){
	         x_max = x[i];
      }
      if(x[i] < x_min){
	         x_min = x[i];
      }
  }
  x_avg = x_sum/u->nobj;

  for(i = u->nobj - 1; i > -1; i--) {

    y = w[i]? w[i] : u->wzero;
    // ASF
    //v = x[i] / y;

    // GrammaticalEvolution_ScalarizingFunction
	v = x_sum*x[i]/y-y/96.51/00.69/y;
    if(v > vmax){
        vmax = v;
    }
  }

		vmax = vmax*x_min-x_min/x_max/x_max-34.32-sqrt(x_avg);

  //EMO_vprint(stdout, w, u->nobj, "w");
  //EMO_vprint(stdout, x, u->nobj, "w");
  //printf("vmax %f, f %f\n\n", vmax, f);

  return vmax;
}

// Augmented weighted achievement function
// Wierzbicki
// Miettinen99, p. 111
// the weighting coefficients can also be dropped from the latter part s += x[i] / w[i];
// Tutum15
// Multiple ciretia decision analysis, Figueira, Greco, Ehrgott, p. 682
double EMO_Utility_aasf(EMO_Utility *u, double *w, double *x) {
  double t, s, vmax, y;
  int i;

  s = vmax = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i] ? w[i] : u->wzero;
    t = x[i] / y;
    s += x[i] / y;

    //printf(" %f\n", t);
    if(t > vmax)
      vmax = t;
  }

  //printf("alpha %f, vmax %f, sum %f, u %f \n", u->aasf_alpha, vmax, s, u->aasf_alpha * s);
  return vmax + u->aasf_alpha * s;
}

// Defined in MOEA/D's paper
double EMO_Utility_pbi(EMO_Utility *u, double *w, double *x) {
  double norm = 0, d1 = 0, d2 = 0;
  int i;

  // Normalize the weight vector (line segment)
  for(i = u->nobj - 1; i > -1; i--)
    norm  += w[i] * w[i];
  norm = sqrt(norm);

  for(i = u->nobj - 1; i > -1; i--)
    d1 += x[i] * w[i] / norm;

  d1 = fabs(d1);

  for(i = u->nobj - 1; i > -1; i--)
    u->pbi_v[i] = x[i] - d1 * w[i] / norm;

  d2 = 0;

  // Normalize b
  for(i = u->nobj - 1; i > -1; i--)
    d2  += u->pbi_v[i] * u->pbi_v[i];

  d2 = sqrt(d2);

  if(u->inverted)
    return u->pbi_theta * d2 - d1;

  return d1 + u->pbi_theta * d2;
}

void EMO_Utility_update_dstar(EMO_Utility *u, double *zmin, double *zmax) {
  double s;
  int i;

  s = 0.0;

  for(i = u->nobj - 1; i > -1; i--)
    s += zmax[i] - zmin[i];

  if(u->uf == EMO_Utility_tlpbi)
    u->tlpbi_dstar = s * u->tlpbi_alpha / ((double) (u->tlpbi_H * u->nobj));
  else //else if(u->uf == EMO_Utility_qpbi)
    u->qpbi_dstar = s * u->qpbi_alpha / ((double) (u->qpbi_H * u->nobj));
}

// Ishibuchi PPSN 2016
double EMO_Utility_tlpbi(EMO_Utility *u, double *w, double *x) {
  double norm = 0, d1 = 0, d2 = 0;
  int i;

  // Normalize the weight vector (line segment)
  for(i = u->nobj - 1; i > -1; i--)
    norm  += w[i] * w[i];
  norm = sqrt(norm);

  for(i = u->nobj - 1; i > -1; i--)
    d1 += x[i] * w[i] / norm;

  d1 = fabs(d1);

  for(i = u->nobj - 1; i > -1; i--)
    u->tlpbi_v[i] = x[i] - d1 * w[i] / norm;

  d2 = 0;

  // Normalize b
  for(i = u->nobj - 1; i > -1; i--)
    d2  += u->tlpbi_v[i] * u->tlpbi_v[i];

  d2 = sqrt(d2);


  if(d2 <= u->tlpbi_dstar)
    return d1 + u->tlpbi_theta1 * d2;
  else
    return d1 + u->tlpbi_theta1 * u->tlpbi_dstar + u->tlpbi_theta2 * (d2 - u->tlpbi_dstar);
}

double EMO_Utility_qpbi(EMO_Utility *u, double *w, double *x) {
  double norm = 0, d1 = 0, d2 = 0;
  int i;

  // Normalize the weight vector (line segment)
  for(i = u->nobj - 1; i > -1; i--)
    norm  += w[i] * w[i];
  norm = sqrt(norm);

  for(i = u->nobj - 1; i > -1; i--)
    d1 += x[i] * w[i] / norm;

  d1 = fabs(d1);

  for(i = u->nobj - 1; i > -1; i--)
    u->qpbi_v[i] = x[i] - d1 * w[i] / norm;

  d2 = 0;

  // Normalize b
  for(i = u->nobj - 1; i > -1; i--)
    d2  += u->qpbi_v[i] * u->qpbi_v[i];

  d2 = sqrt(d2);

  return d1 + u->qpbi_theta * d2 * d2 / u->qpbi_dstar;
}

// Derbel Bilel, Dimo Brockhoff 2014
double EMO_Utility_gsf(EMO_Utility *u, double *w, double *x) {
  double t, s, vmax;
  int i;

  s = vmax = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    t = w[i] * fabs(x[i]);
    s += t;

    if(t > vmax)
      vmax = t;
  }

  return u->gsf_beta * vmax + u->gsf_alpha * s;
}

double EMO_Utility_gsf2(EMO_Utility *u, double *w, double *x) {
  double t, y, s, vmax;
  int i;

  s = vmax = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i] ? w[i] : u->wzero;
    t = fabs(x[i]) / y;
    s += t;

    if(t > vmax)
      vmax = t;
  }

  return u->gsf2_beta * vmax + u->gsf2_alpha * s;
}

double EMO_Utility_nsf(EMO_Utility *u, double *w, double *x) {
  double t, s, vmax;
  int i;

  s = vmax = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    t = w[i] * fabs(x[i]);
    s += t;

    if(t > vmax)
      vmax = t;
  }

  return (1.0 - u->nsf_alpha) * vmax + u->nsf_alpha * s;
}

double EMO_Utility_nsf2(EMO_Utility *u, double *w, double *x) {
  double t, s, y, vmax;
  int i;

  s = vmax = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i] ? w[i] : u->wzero;
    t = fabs(x[i]) / y;
    s += t;

    if(t > vmax)
      vmax = t;
  }

  return (1.0 - u->nsf2_alpha) * vmax + u->nsf2_alpha * s;
}

// Conic Scalarization Method (CS)
// Kasimbeyli15
double EMO_Utility_cs(EMO_Utility *u, double *w, double *x) {
  double s1, s2;
  int i;

  s1 = s2 = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    s1 += w[i] * x[i];
    s2 += fabs(x[i]);
  }

  return s1 + u->cs_alpha * s2;
}

/* Modified utility function, where the weight vector is a divisor */
double EMO_Utility_cs2(EMO_Utility *u, double *w, double *x) {
  double s1, s2, y;
  int i;

  s1 = s2 = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i]? w[i] : u->wzero;
    s1 += x[i] / y;
    s2 += fabs(x[i]);
  }

  return s1 + u->cs2_alpha * s2;
}

// Hughes03b, u->p = 100
double EMO_Utility_vads(EMO_Utility *u, double *w, double *x) {
  double nx, nw, v;
  int i;

  nx = nw = v = 0;

  for(i = 0; i < u->nobj; i++) {
    nx += x[i] * x[i];
    nw += w[i] * w[i];
    v += x[i] * w[i];
  }

  nx = sqrt(nx);
  nw = sqrt(nw);
  v = pow(v / (nx * nw) , u->vads_p);
  return nx / (1e-3 + v);
}

// Hughes03b
double EMO_Utility_vads2(EMO_Utility *u, double *w, double *x) {
  double nx, nw, v, y;
  int i;

  nx = nw = v = 0;

  for(i = 0; i < u->nobj; i++) {
    y = w[i]? w[i] : u->wzero;
    nx += x[i] * x[i];
    nw += 1.0 / (y * y);
    v += x[i] / y;
  }

  nx = sqrt(nx);
  nw = sqrt(nw);
  v = pow(v / (nx * nw) , u->vads2_p);
  return nx / (1e-3 + v);

  //printf("v %f nx %f \n", v, nx);
  //return v;
}


// [Jean-Charles_Pomerol,__Sergio_Barba-Romero__(auth(BookZZ.org), p. 258
double EMO_Utility_didass(EMO_Utility *u, double *w, double *x) {
  double t, s1, s2, vmax;
  int i;

  vmax = s1 = s2 = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    t = w[i] * fabs(x[i]);

    if(t > vmax)
      vmax = t;

    s1 += t;
    s2 += u->didass_v[i] * fabs(x[i]);
  }

  vmax *= u->didass_beta;

  if(vmax > s1)
    return vmax + s2;

  return s1 + s2;
}

/* Modified utility function, where the weight vector is a divisor */
double EMO_Utility_didass2(EMO_Utility *u, double *w, double *x) {
  double t, y, s1, s2, vmax;
  int i;

  vmax = s1 = s2 = 0;

  for(i = u->nobj - 1; i > -1; i--) {
    y = w[i]? w[i] : u->wzero;
    t = fabs(x[i]) / y;

    if(t > vmax)
      vmax = t;

    s1 += t;
    s2 += u->didass2_v[i] * fabs(x[i]);
  }

  vmax *= u->didass2_beta;

  if(vmax > s1)
    return vmax + s2;

  return s1 + s2;
}

// for finding nadir point, points close to axis (1,0)
double EMO_Utility_refsf(EMO_Utility *u, double *w, double *x) {
  return max(EMO_Utility_aasf(u,w,x), EMO_Utility_ws(u,w,x));
}
