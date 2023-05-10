/**************************************************************
 * testf.c   Definition of test functions.                    *
 *                                                            *
 * Computer Science Department, CINVESTAV-IPN                 *
 *                                                            *
 * Professor:   Dr. Carlos A. Coello Coello                   *
 *                                                            *
 * Authors:     Mariano Orozco Garcia                         *
 *              Raquel Hernandez Gomez                        *
 *                                                            *
 * March 2013                                                 *
 *************************************************************/
/*  Reference (chapter 4):

    Carlos A. Coello Coello, Gary B. Lamont, and David A. Van Veldhuizen
    Evolutionary Algorithms for Solving Multi-Objective Problems,
    Second edition, Springer, New York
    September, 2007,
    ISBN 978-0-387-33254-3

    A Multi-Objective Optimization Problem (MOP) should be defined as:

    minimize   {f_0(x), f_1(x), ..., f_{mop->nobj-1}(x)}

    subject to g_0(x) >= 0
               g_1(x) >= 0
               ...
               g_{mop->ncon - 1}(x) >= 0

               x \in [x_min, x_max]
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "benchmark.h"
#include "numeric.h"
#include "vector.h"
#include "wfg.h"
#include "cec09.h"
#include "lame.h"
//#include "knapsack.h"

/* Dictionary of test functions (global variables) */
const char *EMO_Benchmark_list[] = {  "FON1",  "FON2",  "KUR",
                                      "LAU",   "LIS",   "MUR",
                                      "POL",   "QUA",
                                      "REN1",  "REN2",
                                      "SCH1",  "SCH2",
                                      "VIE1",  "VIE2",  "VIE3",
                                      "DEB1",  "DEB2",  "DEB3",
                                      "OKA1",  "OKA2",
                                      "BNH1",  "BNH3",  "SDD",
                                      "STZ1",  "STZ2",  "STZ3",
                                      "ZDT1",  "ZDT2",  "ZDT3",
                                      "ZDT4",  "ZDT5",  "ZDT6",
                                      "UF1",   "UF2",   "UF3",
                                      "UF4",   "UF5",   "UF6",
                                      "UF7",   "UF8",   "UF9",
                                      "UF10",
                                      "DTLZ1", "DTLZ2", "DTLZ3",
                                      "DTLZ4", "DTLZ5", "DTLZ6",
                                      "DTLZ7",
                                      "DTLZ2_CONVEX", "DTLZ3_CONVEX", "DTLZ4_CONVEX",
                                      "DTLZ1_MINUS", "DTLZ2_MINUS", "DTLZ3_MINUS",
                                      "DTLZ4_MINUS", "DTLZ5_MINUS", "DTLZ6_MINUS",
                                      "DTLZ7_MINUS",
                                      "EBN",   "LAME",  "MIRROR",
                                      "WFG1",  "WFG2",  "WFG3",
                                      "WFG4",  "WFG5",  "WFG6",
                                      "WFG7",  "WFG8",  "WFG9",
                                      "WFG1_MINUS", "WFG2_MINUS", "WFG3_MINUS",
                                      "WFG4_MINUS", "WFG5_MINUS", "WFG6_MINUS",
                                      "WFG7_MINUS", "WFG8_MINUS", "WFG9_MINUS",
                                     //"KNAPSACK",
                                      NULL };

const char *EMO_Benchmark_listc[] = { "BEL",   "BNH2",  "BNH4",
                                      "JIM",   "KITA",  "OBA",
                                      "OSY1",  "OSY2",  "SRN",
                                      "TMK",   "TNK",   "VIE4",
                                      "DTLZ8", "DTLZ9",
                                      "CF1",   "CF2",   "CF3",
                                      "CF4",   "CF5",   "CF6",
                                      "CF7",   "CF8",   "CF9",
                                      "CF10",
                                      "WATER", "CAR-SIDE_IMPACT",
                                      "2BAR_TRUSS", NULL };

//EMO_Knapsack knap;

/* Test problem FON1 (Fonseca)

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are open
   PFTrue is formed with x* = [-1, 1]
   POS connected, symmetric; PFTrue connected, concave
*/
void EMO_Benchmark_fon1(EMO_MOP *mop, double *f, double *x) {
  double s1, s2;

  s1 = - pow(x[0] - 1.0, 2.0) - pow(x[1] + 1.0, 2.0);
  s2 = - pow(x[0] + 1.0, 2.0) - pow(x[1] - 1.0, 2.0);

  f[0] = 1.0 - exp(s1);
  f[1] = 1.0 - exp(s2);
}

void EMO_Benchmark_fon1_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -1.0;
    mop->xmax[i] = 1.0;
  }
}

/* Test problem FON2, (Fonseca and Fleming, 1995)

   Real variables = n (scalable)
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [-4,4]
   PFTrue is formed with x* = [-1/sqrt(n), 1/sqrt(n)]
                         f2* = 1-exp(-(2-sqrt(-ln(1-f1*)))^2)
   POS connected, symmetric; PFTrue connected, concave
*/
void EMO_Benchmark_fon2(EMO_MOP *mop, double *f, double *x) {
  double s, s1, s2;
  int i;

  s1 = s2 = 0.0;

  s = 1.0 / sqrt((double) mop->nvar);

  for (i = 0; i < mop->nvar; i++) {
    s1 += pow(x[i] - s, 2.0);
    s2 += pow(x[i] + s, 2.0);
  }

  f[0] = 1.0 - exp(-s1);
  f[1] = 1.0 - exp(-s2);
}

void EMO_Benchmark_fon2_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -4.0;
    mop->xmax[i] = 4.0;
  }
}

/* Test problem KUR (Kursawe, 1990)

   Real variables = 3
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [-5,5]
   POS disconnected, symmetric; PFTrue disconnected
*/
void EMO_Benchmark_kur(EMO_MOP *mop, double *f, double *x) {
  double a, e1, e2;
  int i;

  a = x[1] * x[1];
  e1 = - 0.2 * sqrt(x[0] * x[0] + a);
  e2 = - 0.2 * sqrt(a + x[2] * x[2]);

  f[0] = - 10.0 * (exp(e1) + exp(e2));
  f[1] = 0.0;

  for (i = 0; i < 3; i++)
    f[1] += pow(fabs(x[i]), 0.8) + 5.0 * sin(pow(x[i], 3.0));
}

void EMO_Benchmark_kur_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -5.0;
    mop->xmax[i] = 5.0;
  }
}

/* Test problem LAU (Laumanns)

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [-50,50]
   POS disconnected; PFTrue convex
*/
void EMO_Benchmark_lau(EMO_MOP *mop, double *f, double *x) {
  f[0] = x[0] * x[0] + x[1] * x[1];
  f[1] = pow(x[0] + 2.0, 2.0) + x[1] * x[1];
}

void EMO_Benchmark_lau_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -50.0;
    mop->xmax[i] = 50.0;
  }
}

/* Test problem LIS

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [-5,10]
   POS disconnected; PFTrue disconnected and concave
*/
void EMO_Benchmark_lis(EMO_MOP *mop, double *f, double *x) {
  f[0] = pow(x[0] * x[0] + x[1] * x[1], 1.0 / 8.0);
  f[1] = pow(pow(x[0] - 0.5, 2.0) + pow(x[1] - 0.5, 2.0), 1.0 / 4.0);
}

void EMO_Benchmark_lis_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -5.0;
    mop->xmax[i] = 10.0;
  }
}

/* Test problem MUR (Murata)

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x1 is in [1,4], x2 is in [1,2]
   POS connected, PFTrue concave
*/
void EMO_Benchmark_mur(EMO_MOP *mop, double *f, double *x) {
  f[0] = 2.0 * sqrt(x[0]);
  f[1] = x[0] * (1.0 - x[1]) + 5.0;
}

void EMO_Benchmark_mur_range(EMO_MOP *mop) {
  mop->xmin[0] = 1.0;
  mop->xmax[0] = 4.0;
  mop->xmin[1] = 1.0;
  mop->xmax[1] = 2.0;
}

/* Test problem POL (Poloni, 2000)

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [-PI,PI]
   PFTrue is difficult to know
   POS disconnected; PFTrue disconnected and convex
*/
void EMO_Benchmark_pol(EMO_MOP *mop, double *f, double *x) {
  double a1, a2, b1, b2;

  a1 = 0.5 * sin(1.0) - 2.0 * cos(1.0) + sin(2.0) - 1.5 * cos(2.0);
  a2 = 1.5 * sin(1.0) - cos(1.0) + 2.0 * sin(2.0) - 0.5 * cos(2.0);
  b1 = 0.5 * sin(x[0]) - 2.0 * cos(x[0]) + sin(x[1]) - 1.5 * cos(x[1]);
  b2 = 1.5 * sin(x[0]) - cos(x[0]) + 2.0 * sin(x[1]) - 0.5 * cos(x[1]);
  f[0] = 1.0 + pow(a1 - b1, 2.0) + pow(a2 - b2, 2.0);
  f[1] = pow(x[0] + 3.0, 2.0) + pow(x[1] + 1.0, 2.0);
}

void EMO_Benchmark_pol_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -PI;
    mop->xmax[i] = PI;
  }
}

/* Test problem QUA (Quagliarella)

   Real variables >= 16
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [-5.12, 5.12]
   POS disconnected, PFTrue concave
*/
void EMO_Benchmark_qua(EMO_MOP *mop, double *f, double *x) {
  double a1, a2;
  int i;

  a1 = a2 = 0;

  for(i = 0; i < mop->nvar; i++) {
    a1 += x[i] * x[i] - 10.0 * cos(2.0 * PI * x[i]) + 10.0;
    a2 += pow(x[i] - 1.5, 2.0) - 10.0 * cos(2.0 * PI * (x[i] - 1.5)) + 10.0;
  }

  f[0] = pow(a1 / (double) mop->nvar, 1.0 / 4.0);
  f[1] = pow(a2 / (double) mop->nvar, 1.0 / 4.0);
}

void EMO_Benchmark_qua_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -5.12;
    mop->xmax[i] = 5.12;
  }
}

/* Test problem REN1 (Rendon)

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [-3, 3]
   POS connected, PFTrue convex
*/
void EMO_Benchmark_ren1(EMO_MOP *mop, double *f, double *x) {
  double a, b;

  a = x[0] * x[0];
  b = x[1] * x[1];

  f[0] = 1.0 / (a + b + 1.0);
  f[1] = a + 3.0 * b + 1.0;
}

/* Test problem REN2 (Rendon)

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [-3, 3]
   POS connected, PFTrue convex
*/
void EMO_Benchmark_ren2(EMO_MOP *mop, double *f, double *x) {
  f[0] = x[0] + x[1] + 1.0;
  f[1] = x[0] * x[0] + 2.0 * x[1] - 1.0;
}

void EMO_Benchmark_ren_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -3.0;
    mop->xmax[i] = 3.0;
  }
}

/* Test problem SCH1 (Schaffer, 1984)

   Real variables = 1
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [-A,A], A=10..10^5 (more difficulty)
   PFTrue is formed with f1* = [0,4], f2* = (sqrt(f1*) - 2)^2
   POS connected, PFTrue convex
*/
void EMO_Benchmark_sch1(EMO_MOP *mop, double *f, double *x) {
  f[0] = x[0] * x[0];
  f[1] = pow(x[0] - 2.0, 2.0);
}

void EMO_Benchmark_sch1_range(EMO_MOP *mop) {
  mop->xmin[0] = -10.0;
  mop->xmax[0] = 10.0;
}

/* Test problem SCH2 (Schaffer, 1984)

   Real variables = 1
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [-5,10]
   PFTrue is formed with x* in {[1,2] union [4,5]}
   POS disconnected; PFTrue convex are disconnected
*/
void EMO_Benchmark_sch2(EMO_MOP *mop, double *f, double *x) {

  if (x[0] <= 1.0)
    f[0] = - x[0];
  else if (x[0] <= 3.0)
    f[0] = x[0] - 2.0;
  else if (x[0] <= 4.0)
    f[0] = 4.0 - x[0];
  else f[0] = x[0] - 4.0;

  f[1] = pow(x[0] - 5.0, 2.0);
}

void EMO_Benchmark_sch2_range(EMO_MOP *mop) {
  mop->xmin[0] = -5.0;
  mop->xmax[0] = 10.0;
}

/* Test problem VIE1 (Viennet, 1996)

   Real variables = 2
   Bin variables = 0
   Objectives = 3
   Constraints = 0
   Values of x are in [-2,2]
   POS connected and symmetric, PFTrue curved surface
*/
void EMO_Benchmark_vie1(EMO_MOP *mop, double *f, double *x) {
  double a;

  a = x[0] * x[0];

  f[0] = a + pow(x[1] - 1.0, 2.0);
  f[1] = a + pow(x[1] + 1.0, 2.0) + 1.0;
  f[2] = pow(x[0] - 1.0, 2.0) + x[1] * x[1] + 2.0;
}

void EMO_Benchmark_vie1_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -2.0;
    mop->xmax[i] = 2.0;
  }
}

/* Test problem VIE2 (Viennet, 1996)

   Real variables = 2
   Bin variables = 0
   Objectives = 3
   Constraints = 0
   Values of x are in [-4,4]
   POS connected, PFTrue disconnected
*/
void EMO_Benchmark_vie2(EMO_MOP *mop, double *f, double *x) {
  f[0] = pow(x[0] - 2.0, 2.0) * 0.5 + pow(x[1] + 1.0, 2.0) / 13.0 + 3.0;
  f[1] = pow(x[0] + x[1] - 3.0, 2.0) / 36.0 + pow(-x[0] + x[1] + 2.0, 2.0) * 0.125 - 17.0;
  f[2] = pow(x[0] + 2.0 * x[1] - 1.0, 2.0) / 175.0 + pow(2.0 * x[1] - x[0], 2.0) / 17.0 - 13.0;
}

void EMO_Benchmark_vie2_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -4.0;
    mop->xmax[i] = 4.0;
  }
}

/* Test problem VIE3 (Viennet, 1996)

   Real variables = 2
   Bin variables = 0
   Objectives = 3
   Constraints = 0
   Values of x are in [-3,3]
   POS disconnected and unsymmetric, PFTrue connected
*/
void EMO_Benchmark_vie3(EMO_MOP *mop, double *f, double *x) {
  double a;

  a = x[0] * x[0] + x[1] * x[1];

  f[0] = 0.5 * a + sin(a);
  f[1] = pow(3.0 * x[0] - 2.0 * x[1] + 4.0, 2.0) * 0.125 + pow(x[0] - x[1] + 1.0, 2.0) / 27.0 + 15.0;
  f[2] = 1.0 / (a + 1.0) - 1.1 * exp(-a);
}

void EMO_Benchmark_vie3_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -3.0;
    mop->xmax[i] = 3.0;
  }
}

/* Test problem DEB1

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1]
   POS connected; PFTrue connected, concave
*/
void EMO_Benchmark_deb1(EMO_MOP *mop, double *f, double *x) {
  double g, h;

  g = 1.0 + x[1] * x[1];
  h = (x[0] <= g) ? 1.0 - pow(x[0] / g, 2.0) : 0.0;

  f[0] = x[0];
  f[1] = g * h;
}

/* Test problem DEB2

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1]
   POS disconnected, PFTrue disconnected
*/
void EMO_Benchmark_deb2(EMO_MOP *mop, double *f, double *x) {
  double g, h;

  g = 1.0 + 10.0 * x[1];
  h = 1.0 - pow(x[0] / g, 2.0) - x[0] / g * sin(12.0 * PI * x[0]);

  f[0] = x[0];
  f[1] = g * h;
}

/* Test problem DEB3

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1]
   POS connected; PFTrue connected, concave
*/
void EMO_Benchmark_deb3(EMO_MOP *mop, double *f, double *x) {
  double g, h;

  f[0] = 1.0 - exp(-4.0 * x[0]) * pow(sin(10.0 * PI * x[0]), 4.0);

  g = 1.0 + x[1] * x[1];
  h = (f[0] <= g) ? 1.0 - pow(f[0] / g, 10.0) : 0.0;

  f[1] = g * h;
}

void EMO_Benchmark_deb_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = 0.0;
    mop->xmax[i] = 1.0;
  }
}

/* Test problem OKA1

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x1 are in [6*sin(PI/12), 6*sin(PI/12)+2*PI*cos(PI/12)]
                       [1.5529142706, 7.6220052302]
             x2 are in [-2*PI*sin(PI/12), 6*cos(PI/12)]
                       [-1.6262080214 5.7955549577]
   PFTrue is at f2 = sqrt(2*PI) - sqrt(f1), f1 in [0, 2 * PI]
                x2' = 3 * cos(x1') + 3, x1' in [0, 2*PI]
   Convex Pareto-optimal front
*/

void EMO_Benchmark_oka1(EMO_MOP *mop, double *f, double *x) {
  double x0, x1;

  x0 = cos(PI/12.0) * x[0] - sin(PI/12.0) * x[1];
  x1 = sin(PI/12.0) * x[0] + cos(PI/12.0) * x[1];

  f[0] = x0;
  f[1] = sqrt(2.0 * PI) - sqrt(fabs(x0)) + 2.0 * pow(fabs(x1 - 3.0 * cos(x0) - 3.0), 1.0/3.0);
}

void EMO_Benchmark_oka1_range(EMO_MOP *mop) {
  mop->xmin[0] = 6.0 * sin(PI/12.0);
  mop->xmax[0] = mop->xmin[0] + 2 * PI * cos(PI/12.0);
  mop->xmin[1] = -2.0 * PI * sin(PI/12.0);
  mop->xmax[1] = 6.0 * cos(PI/12.0);
}

/* Test problem OKA2

   Real variables = 3
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x1 are in [-PI,PI], x2,x3 in [-5,5]
   PFTrue is at f2 = 1 - 1/(4*PI^2) * (f1 + PI)^2, f1 in [-PI, PI]
        (x1,x2,x3) = (x1, 5*cos(x1), 5*sin(x1))
   Non-convex Pareto-optimal front
*/
void EMO_Benchmark_oka2(EMO_MOP *mop, double *f, double *x) {
  f[0] = x[0];
  f[1] = 1.0 - pow(x[0] + PI, 2.0) / (4.0 * PI * PI);
  f[1] += pow(fabs(x[1] - 5.0 * cos(x[0])), 1.0/3.0);
  f[1] += pow(fabs(x[2] - 5.0 * sin(x[0])), 1.0/3.0);  // RHG, Revisar, x[2], rangos no estan definidos
}

void EMO_Benchmark_oka2_range(EMO_MOP *mop) {
  mop->xmin[0] = -PI;
  mop->xmax[0] = PI;
  mop->xmin[1] = -5.0;
  mop->xmax[1] = 5.0;
}

/* Test problem BNH1 (Binh)

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [-5,10]
   POS connected, symmetric, PFTrue connected, convex
*/
void EMO_Benchmark_bnh1(EMO_MOP *mop, double *f, double *x) {
  f[0] = x[0] * x[0] + x[1] * x[1];
  f[1] = pow(x[0] - 5.0, 2.0) + pow(x[1] - 5.0, 2.0);
}

void EMO_Benchmark_bnh1_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -5.0;
    mop->xmax[i] = 10.0;
  }
}

/* Test problem BNH3 (Binh)

   Real variables = 2
   Bin variables = 0
   Objectives = 3
   Constraints = 0
   Values of x are in [10^-6,10^6]
   PFTrue = [-1e+06, -1e-06, -2e+00]
   POS connected, PFTrue is a point
*/
void EMO_Benchmark_bnh3(EMO_MOP *mop, double *f, double *x) {
  f[0] = x[0] - pow(10.0, 6.0);
  f[1] = x[1] - 2.0 * pow(10.0, -6.0);
  f[2] = x[0] * x[1] - 2.0;
}

void EMO_Benchmark_bnh3_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = 1e-6;
    mop->xmax[i] = 1e6;
  }
}

/* Test problem SDD (Schutze, Dell'Aere, Dellnitz)

   Real variables = (scalable)
   Bin variables = 0
   Objectives = (scalable)
   Constraints = 0
   POS connected, PFTrue convex
*/
void EMO_Benchmark_sdd(EMO_MOP *mop, double *f, double *x) {
  double a, b = 0;
  int i, j;

  for(j = 0; j < mop->nobj; j++) {
    f[j] = 0;
    a = -1.0;

    for(i = 0; i < mop->nvar; i++) {
      switch(j) {
        case 0: a = 1.0;
                break;
        case 1: a = -1.0;
                break;
        default: if(i % (j-1) == 0)
                   a *= (-1.0);
                 break;
      }

      if(j == i)
        b = a;
      else
        f[j] += pow(x[i] - a, 2.0);
    }
    f[j] += pow(x[j] - b, 4.0);
  }
}

/* Test problem STZ1 (Schutze)

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   POS connected, PFTrue convex
*/
void EMO_Benchmark_stz1(EMO_MOP *mop, double *f, double *x) {
  f[0] = pow(x[0] - 1.0, 2.0) + pow(x[1], 2.0);
  f[1] = pow(x[0], 2.0) + pow(x[1] - 1.0, 2.0);
}

/* Test problem STZ2 (Schutze)

   Real variables = (scalable)
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   POS connected, PFTrue convex
*/
void EMO_Benchmark_stz2(EMO_MOP *mop, double *f, double *x) {
  int i;

  f[0] = pow(x[0] - 1.0, 4.0);

  for(i = 1; i < mop->nvar; i++)
    f[0] += pow(x[i] - 1.0, 2.0);

  f[1] = 0;

  for(i = 0; i < mop->nvar; i++)
    f[1] += pow(x[i] + 1.0, 2.0);
}

/* Test problem STZ3 (Schutze)

   Real variables = 2
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   POS connected, PFTrue mixed (convex, concave)
*/
void EMO_Benchmark_stz3(EMO_MOP *mop, double *f, double *x) {
  double a, b, c, lambda= 0.85;

  a = sqrt(1.0 + pow(x[0] + x[1], 2.0));
  b = sqrt(1.0 + pow(x[0] - x[1], 2.0));
  c = lambda * exp(-pow(x[0] - x[1], 2.0));

  f[0] = 0.5 * (a + b + x[0] - x[1]) + c;
  f[1] = 0.5 * (a + b - x[0] + x[1]) + c;
}

// verificar
void EMO_Benchmark_stz_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = -10.0;
    mop->xmax[i] = 10.0;
  }
}

/* Test problem ZDT1 (Zitzler-Deb-Thiele, 2000)

   Real variables >= 30
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1]
   PFTrue is formed with g(x) = 1
   Convex Pareto-optimal front
*/
void EMO_Benchmark_zdt1(EMO_MOP *mop, double *f, double *x) {
  double g, h;
  int i;

  f[0] = x[0];
  g = 0;

  for(i = 1; i < mop->nvar; i++)
    g += x[i];

  g = 1.0 + 9.0 * g / (mop->nvar - 1.0);
  h = 1.0 - sqrt(x[0] / g);
  f[1] = g * h;
}

/* Test problem ZDT2 (Zitzler-Deb-Thiele, 2000)

   Real variables >= 30
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1]
   PFTrue is formed with g(x) = 1
   Non-convex Pareto-optimal front
*/
void EMO_Benchmark_zdt2(EMO_MOP *mop, double *f, double *x) {
  double g, h;
  int i;

  f[0] = x[0];
  g = 0;

  for(i = 1; i < mop->nvar; i++)
    g += x[i];

  g = 1.0 + 9.0 * g / (mop->nvar - 1.0);
  h = 1.0 - pow(x[0] / g, 2.0);
  f[1] = g * h;
}

/* Test problem ZDT3 (Zitzler-Deb-Thiele, 2000)

   Real variables >= 30
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1]
   PFTrue is formed with g(x) = 1
   Pareto-optimal front disconnected, consisting of
   several non-contiguous convex parts.
*/
void EMO_Benchmark_zdt3(EMO_MOP *mop, double *f, double *x) {
  double g, h;
  int i;

  f[0] = x[0];
  g = 0;

  for(i = 1; i < mop->nvar; i++)
    g += x[i];

  g = 1.0 + 9.0 * g / (mop->nvar - 1.0);
  h = 1 - sqrt(x[0] / g) - (x[0] / g) * sin(10.0 * PI * x[0]);
  f[1] = g * h;
}

void EMO_Benchmark_zdt_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = 0.0;
    mop->xmax[i] = 1.0;
  }
}

/* Test problem ZDT4 (Zitzler-Deb-Thiele, 2000)

   Real variables = 10
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x1 in [0,1]; x2,...,xn in [-5,5]
   PFTrue is formed with g(x) = 1, PFlocal g(x) = 1.25
   21^9 local Pareto-optimal solutions (multi-frontality)
*/
void EMO_Benchmark_zdt4(EMO_MOP *mop, double *f, double *x) {
  double g, h;
  int i;

  f[0] = x[0];
  g = 0;

  for(i = 1; i < mop->nvar; i++)
    g += x[i] * x[i] - 10.0 * cos(4.0 * PI * x[i]);

  g += 1.0 + 10.0 * ((double) mop->nvar - 1.0);
  h = 1.0 - sqrt(x[0] / g);
  f[1] = g * h;
}

void EMO_Benchmark_zdt4_range(EMO_MOP *mop) {
  int i; //, n;
printf("zdt4Range\n");
  mop->xmin[0] = 0.0;
  mop->xmax[0] = 1.0;

  for(i = 1; i < mop->nvar; i++) {
    mop->xmin[i] = -5.0;
    mop->xmax[i] = 5.0;
  }
}

/* Test problem ZDT5 (Zitzler-Deb-Thiele, 2000)

   Real variables = 0
   Bin variables = 11
   Objectives = 2
   Constraints = 0
   Values of x1 is represented by a 30-bit substring,
             x2-11 are represented by 5 bits each
   PFTrue is formed with g(x) = 10, x2-x10 are all 1s.
   Convex, deceptive and multi-frontal problem.
*/
void EMO_Benchmark_zdt5(EMO_MOP *mop, double *f, double *x) {
/*void zdt5 (double *xreal, double *xbin, int **gene, double *obj, double *constr, int n, int k)

  int i, j;
  int u[11];
  int v[11];
  double f1, f2, g, h;

  for(i = 0; i < 11; i++)
    u[i] = 0.0;

  for(j = 0; j < 30; j++)
    if(gene[0][j] == 1)
      u[0]++;

  for(i = 1; i < 11; i++)
    for(j = 0; j < 4; j++)
      if(gene[i][j] == 1)
        u[i]++;

  f1 = 1.0 + u[0];

  for(i = 1; i < 11; i++)
    if(u[i] < 5)
      v[i] = 2.0 + u[i];
    else
      v[i] = 1.0;

  g = 0.0;

  for (i = 1; i < 11; i++)
    g += v[i];

  h = 1.0 / f1;
  f2 = g * h;
  f[0] = f1;
  f[1] = f2;
*/
}

/* Test problem ZDT6 (Zitzler-Deb-Thiele, 2000)

   Real variables >= 10
   Bin variables = 0
   Objectives = 2
   constraints = 0
   Values of xi in [0,1]
   PFTrue is formed with g(x) = 1 and is non-convex.
   PTrue are non-uniformly distributed along the PFTrue,
   the density of the solutions is lowest near the PFTrue
   and highest away from the front.
*/
void EMO_Benchmark_zdt6(EMO_MOP *mop, double *f, double *x) {
  double g, h;
  int i;

  f[0] = 1.0 - exp(-4.0 * x[0]) * pow(sin(6.0 * PI * x[0]), 6.0);
  g = 0;

  for(i = 1; i < mop->nvar; i++)
    g += x[i];

  //g = 1.0 + 9.0 * pow(g / 9.0, 0.25);
  g = 1.0 + 9.0 * pow(g / ((double) mop->nvar - 1.0), 0.25);
  h = 1.0 - pow(f[0] / g, 2.0);
  f[1] = g * h;
}

/* Test problem DTLZ1 (Deb-Thiele-Laumanns-Zitzler)

   Real variables = mop->nobj + k - 1, k = 5
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is linear, separable, multimodal
   POS is at x_{M}^{*} = 0, sum_{m=1}^{M} fm = 0.5
*/
void EMO_Benchmark_dtlz1(EMO_MOP *mop, double *f, double *x) {
  double g, v;
  int i, j, t, c;

  g = 0.0;
  c = mop->nvar - mop->nobj + 1;

  for(i = mop->nvar - c; i < mop->nvar; i++)
    g += pow(x[i] - 0.5, 2.0) - cos(20.0 * PI * (x[i] - 0.5));

  g = 100.0 * (c + g);
  v = 0.5 * (1.0 + g);

  for(i = 0; i < mop->nobj; i++) {
    f[i] = v;
    t = mop->nobj - i - 1;

    for(j = 0; j < t; j++)
      f[i] *= x[j];

    if(t < mop->nobj - 1)
      f[i] *= (1.0 - x[t]);
  }
}

/* Test problem DTLZ2 (Deb-Thiele-Laumanns-Zitzler)

   Real variables = NOBJ + k - 1, k = 10
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is concave
   POS is at xi = 0.5 for all xi \in x_{M}, sum_{m=1}^{M} fm^2 = 1
*/
void EMO_Benchmark_dtlz2(EMO_MOP *mop, double *f, double *x) {
  double g, v;
  int i, j, t, c;

  g = 0.0;
  c = mop->nvar - mop->nobj + 1;

  for(i = mop->nvar - c; i < mop->nvar; i++)
    g += pow(x[i] - 0.5, 2.0);

  v = 1.0 + g;

  for(i = 0; i < mop->nobj; i++) {
    f[i] = v;
    t = mop->nobj - i - 1;

    for(j = 0; j < t; j++)
      f[i] *= cos(x[j] * PI * 0.5);

    if(t < mop->nobj - 1)
      f[i] *= sin(x[t] * PI * 0.5);
  }
}

/* Test problem DTLZ3 (Deb-Thiele-Laumanns-Zitzler)

   Real variables = NOBJ + k - 1, k = 10
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is concave and multimodal
   POS is at xi = 0.5 for all xi \in x_{M}, g* = 0, sum_{m=1}^{M} fm^2 = 1
*/
void EMO_Benchmark_dtlz3(EMO_MOP *mop, double *f, double *x) {
  double g, v;
  int i, j, t, c;

  g = 0.0;
  c = mop->nvar - mop->nobj + 1;

  for(i = mop->nvar - c; i < mop->nvar; i++)
    g += pow(x[i] - 0.5, 2.0) - cos(20.0 * PI * (x[i] - 0.5));

  g = 100.0 * (c + g);
  v = 1.0 + g;

  for(i = 0; i < mop->nobj; i++) {
    f[i] = v;
    t = mop->nobj - i - 1;

    for(j = 0; j < t; j++)
      f[i] *= cos(x[j] * PI * 0.5);

    if(t < mop->nobj - 1)
      f[i] *= sin(x[t] * PI * 0.5);
  }
}

/* Test problem DTLZ4 (Deb-Thiele-Laumanns-Zitzler)

   Real variables = NOBJ + k - 1, k = 10
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is concave, separable, unimodal
   POS is at xi = 0.5 for all xi \in x_{M}, g* = 0, sum_{m=1}^{M} fm^2 = 1
*/
void EMO_Benchmark_dtlz4(EMO_MOP *mop, double *f, double *x) {
  static double alpha = 100.0;
  double g, v;
  int i, j, t, c;

  g = 0.0;
  c = mop->nvar - mop->nobj + 1;

  for(i = mop->nvar - c; i < mop->nvar; i++)
    g += pow(x[i] - 0.5, 2.0);

  v = 1.0 + g;

  for(i = 0; i < mop->nobj; i++) {
    f[i] = v;
    t = mop->nobj - i - 1;

    for(j = 0; j < t; j++)
      f[i] *= cos(pow(x[j], alpha) * PI * 0.5);

    if(t < mop->nobj - 1)
      f[i] *= sin(pow(x[t], alpha) * PI * 0.5);
  }
}

/* Test problem DTLZ5 (Deb-Thiele-Laumanns-Zitzler)

   Real variables = NOBJ + k - 1, k = 10
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is unimodal, m < 4 degenerate
   POS is at xi = 0.5 for all xi \in x_{M}, g* = 0, sum_{m=1}^{M} fm^2 = 1
*/
void EMO_Benchmark_dtlz5(EMO_MOP *mop, double *f, double *x) {
  double g, p, v;
  int i, j, t, c;

  g = 0.0;
  c = mop->nvar - mop->nobj + 1;

  for(i = mop->nvar - c; i < mop->nvar; i++)
    g += pow(x[i] - 0.5, 2.0);

  mop->t[0] = x[0] * PI * 0.5;
  p = PI / (4.0 * (1.0 + g));

  for(i = 1; i < mop->nobj-1; i++)
    mop->t[i] = p * (1.0 + 2.0 * g * x[i]);

  v = 1.0 + g;

  for(i = 0; i < mop->nobj; i++) {
    f[i] = v;
    t = mop->nobj - i - 1;

    for(j = 0; j < t; j++)
      f[i] *= cos(mop->t[j]);

    if(t < mop->nobj - 1)
      f[i] *= sin(mop->t[t]);
  }
}

/* Test problem DTLZ6 (Deb-Thiele-Laumanns-Zitzler)

   Real variables = NOBJ + k - 1, k = 10
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is unimodal, bias, many-to-one-mapping, m < 4 degenerate
   POS is at xi = 0 for all xi \in x_{M}, sum_{m=1}^{M} fm^2 = 1
*/
void EMO_Benchmark_dtlz6(EMO_MOP *mop, double *f, double *x) {
  double g, p, v;
  int i, j, t, c;

  g = 0.0;
  c = mop->nvar - mop->nobj + 1;

  for(i = mop->nvar - c; i < mop->nvar; i++)
    g += pow(x[i], 0.1);

  mop->t[0] = x[0] * PI * 0.5;
  p = PI / (4.0 * (1.0 + g));

  for(i = 1; i < mop->nobj-1; i++)
    mop->t[i] = p * (1.0 + 2.0 * g * x[i]);

  v = 1.0 + g;

  for(i = 0; i < mop->nobj; i++) {
    f[i] = v;
    t = mop->nobj - i - 1;

    for(j = 0; j < t; j++)
      f[i] *= cos(mop->t[j]);

    if(t < mop->nobj - 1)
      f[i] *= sin(mop->t[t]);
  }
}

/* Test problem DTLZ7 (Deb-Thiele-Laumanns-Zitzler)

   Real variables = NOBJ + k - 1, k = 20
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is disconnected
   POS is at xM = 0, sum_{m=1}^{M} fm^2 = 1
*/
void EMO_Benchmark_dtlz7(EMO_MOP *mop, double *f, double *x) {
  double g, h;
  int i, c;

  g = h = 0.0;
  c = mop->nvar - mop->nobj + 1;

  for(i = 0; i < mop->nobj - 1; i++)
    f[i] = x[i];

  for(i = mop->nvar - c; i < mop->nvar; i++)
    g += x[i];

  g = 1.0 + 9.0 * g / (double) c;

  for(i = 0; i < mop->nobj - 1; i++)
    h += x[i] / (1.0 + g) * (1.0 + sin(3.0 * PI * x[i]));

  h = mop->nobj - h;

  f[mop->nobj-1] = (1.0 + g) * h;
}

void EMO_Benchmark_dtlz_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = 0.0;
    mop->xmax[i] = 1.0;
  }
}

/* Convex DTLZ2 Problem
 * Defined in NSGA-III part I, \cite[p. 592]{Deb14a}
 * The Pareto-optimal surface is given as:
 * sqrt(f1) + sqrt(f2) + ... + sqrt(fm-1) + fm = 1
 */
void EMO_Benchmark_dtlz2_convex(EMO_MOP *mop, double *f, double *x) {
  int i;

  EMO_Benchmark_dtlz2(mop, f, x);

  for(i = 0; i < mop->nobj - 1; i++)
    f[i] = pow(f[i], 4.0);

  i = mop->nobj - 1;
  f[i] *= f[i];
}

/* Convex DTLZ3 Problem
 * The Pareto-optimal surface is given as:
 * sqrt(f1) + sqrt(f2) + ... + sqrt(fm-1) + fm = 1
 */
void EMO_Benchmark_dtlz3_convex(EMO_MOP *mop, double *f, double *x) {
  int i;

  EMO_Benchmark_dtlz3(mop, f, x);

  for(i = 0; i < mop->nobj - 1; i++)
    f[i] = pow(f[i], 4.0);

  i = mop->nobj - 1;
  f[i] *= f[i];
}

/* Convex DTLZ4 Problem
 * The Pareto-optimal surface is given as:
 * sqrt(f1) + sqrt(f2) + ... + sqrt(fm-1) + fm = 1
 */
void EMO_Benchmark_dtlz4_convex(EMO_MOP *mop, double *f, double *x) {
  int i;

  EMO_Benchmark_dtlz4(mop, f, x);

  for(i = 0; i < mop->nobj - 1; i++)
    f[i] = pow(f[i], 4.0);

  i = mop->nobj - 1;
  f[i] *= f[i];
}

/* Minus DTLZ1 (Ishibuchi17)

   Real variables = mop->nobj + k - 1, k = 5
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is linear, separable, multimodal
   POS is at x_{M}^{*} = 0, sum_{m=1}^{M} fm = 0.5
*/
void EMO_Benchmark_dtlz1_minus(EMO_MOP *mop, double *f, double *x) {
  int i;

  EMO_Benchmark_dtlz1(mop, f, x);

  for(i = mop->nobj - 1; i > -1; i--)
    f[i] *= -1.0;
}

/* Minus DTLZ2 (Ishibuchi17)

   Real variables = NOBJ + k - 1, k = 10
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is concave
   POS is at xi = 0.5 for all xi \in x_{M}, sum_{m=1}^{M} fm^2 = 1
*/
void EMO_Benchmark_dtlz2_minus(EMO_MOP *mop, double *f, double *x) {
  int i;

  EMO_Benchmark_dtlz2(mop, f, x);

  for(i = mop->nobj - 1; i > -1; i--)
    f[i] *= -1.0;
}

/* Minus DTLZ3 (Ishibuchi17)

   Real variables = NOBJ + k - 1, k = 10
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is concave and multimodal
   POS is at xi = 0.5 for all xi \in x_{M}, g* = 0, sum_{m=1}^{M} fm^2 = 1
*/
void EMO_Benchmark_dtlz3_minus(EMO_MOP *mop, double *f, double *x) {
  int i;

  EMO_Benchmark_dtlz3(mop, f, x);

  for(i = mop->nobj - 1; i > -1; i--)
    f[i] *= -1.0;
}

/* Minus DTLZ4 (Ishibuchi17)

   Real variables = NOBJ + k - 1, k = 10
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is concave, separable, unimodal
   POS is at xi = 0.5 for all xi \in x_{M}, g* = 0, sum_{m=1}^{M} fm^2 = 1
*/
void EMO_Benchmark_dtlz4_minus(EMO_MOP *mop, double *f, double *x) {
  int i;

  EMO_Benchmark_dtlz4(mop, f, x);

  for(i = mop->nobj - 1; i > -1; i--)
    f[i] *= -1.0;
}

/* Minus DTLZ5 (Ishibuchi17)

   Real variables = NOBJ + k - 1, k = 10
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is unimodal, m < 4 degenerate
   POS is at xi = 0.5 for all xi \in x_{M}, g* = 0, sum_{m=1}^{M} fm^2 = 1
*/
void EMO_Benchmark_dtlz5_minus(EMO_MOP *mop, double *f, double *x) {
  int i;

  EMO_Benchmark_dtlz5(mop, f, x);

  for(i = mop->nobj - 1; i > -1; i--)
    f[i] *= -1.0;
}

/* Minus DTLZ6 (Ishibuchi17)

   Real variables = NOBJ + k - 1, k = 10
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is unimodal, bias, many-to-one-mapping, m < 4 degenerate
   POS is at xi = 0 for all xi \in x_{M}, sum_{m=1}^{M} fm^2 = 1
*/
void EMO_Benchmark_dtlz6_minus(EMO_MOP *mop, double *f, double *x) {
  int i;

  EMO_Benchmark_dtlz6(mop, f, x);

  for(i = mop->nobj - 1; i > -1; i--)
    f[i] *= -1.0;
}

/* Minus DTLZ7 (Ishibuchi17)

   Real variables = NOBJ + k - 1, k = 20
   Bin variables = 0
   Objectives = (scalable)
   constraints = 0
   Values of xi in [0,1]
   PFTrue is disconnected
   POS is at xM = 0, sum_{m=1}^{M} fm^2 = 1
*/
void EMO_Benchmark_dtlz7_minus(EMO_MOP *mop, double *f, double *x) {
  int i;

  EMO_Benchmark_dtlz7(mop, f, x);

  for(i = mop->nobj - 1; i > -1; i--)
    f[i] *= -1.0;
}

/* Emmerich05, Lame super-spheres for 2 objectives
 * gamma = 1 yields linear Pareto front
 * gamma > 1 leads to convex fronts
 * gamma < 1 concave ones
 */
void EMO_Benchmark_ebn(EMO_MOP *mop, double *f, double *x) {
  double v;
  int i;

  f[0] = f[1] = 0.0;

  for(i = mop->nvar - 1; i > -1; i--) {
    f[0] += fabs(x[i]);
    f[1] += fabs(x[i] - 1.0);
  }

  v = pow((double) mop->nvar, -mop->gamma);
  f[0] = pow(f[0], mop->gamma) * v;
  f[1] = pow(f[1], mop->gamma) * v;
}

void EMO_Benchmark_ebn_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = 0.0;
    mop->xmax[i] = 1.0;
  }
}

/* Belegundu
 *
 * f: 2  objective functions
 * g: 2  constraints
 * x: 2  variables
 */
void EMO_Benchmark_bel(EMO_MOP *mop, double *f, double *g, double *x) {
  f[0] = - 2.0 * x[0] + x[1];
  f[1] = 2.0 * x[0] + x[1];

  g[0] = x[0] - x[1] + 1.0;
  g[1] = -x[0] - x[1] + 7.0;
}

void EMO_Benchmark_bel_range(EMO_MOP *mop) {
  mop->xmin[0] = 0.0;
  mop->xmax[0] = 5.0;

  mop->xmin[1] = 0.0;
  mop->xmax[1] = 3.0;
}

/* Constraint Binh 2
 *
 * f: 2  objective functions
 * g: 2  constraints
 * x: 2  variables
 */
void EMO_Benchmark_bnh2(EMO_MOP *mop, double *f, double *g, double *x) {
  f[0] = 4.0 * x[0] * x[0] + 4.0 * x[1] * x[1];
  f[1] = (x[0] - 5.0) * (x[0] - 5.0) + (x[1] - 5.0) * (x[1] - 5.0);

  g[0] = -(x[0] - 5.0) * (x[0] - 5.0) - x[1] * x[1] + 25.0;
  g[1] = (x[0] - 8.0) * (x[0] - 8.0) + (x[1] + 3.0) * (x[1] + 3.0) - 7.7;
}

void EMO_Benchmark_bnh2_range(EMO_MOP *mop) {
  mop->xmin[0] = 0.0;
  mop->xmax[0] = 5.0;

  mop->xmin[1] = 0.0;
  mop->xmax[1] = 5.0;
}

/* Constraint Binh 4
 *
 * f: 3  objective functions
 * g: 2  constraints
 * x: 2  variables
 */
void EMO_Benchmark_bnh4(EMO_MOP *mop, double *f, double *g, double *x) {
  f[0] = 1.5 - x[0] * (1.0 - x[1]);
  f[1] = 2.25 - x[0] * (1.0 - x[1] * x[1]);
  f[2] = 2.625 - x[0] * (1.0 - x[1] * x[1] * x[1]);

  g[0] = x[0] * x[0] + (x[1] - 0.5) * (x[1] - 0.5) - 9.0;
  g[1] = -(x[0] - 1.0) * (x[0] - 1.0) - (x[1] - 0.5) * (x[1] - 0.5) + 6.25;
}

void EMO_Benchmark_bnh4_range(EMO_MOP *mop) {
  mop->xmin[0] = -10.0;
  mop->xmax[0] = 10.0;

  mop->xmin[1] = -10.0;
  mop->xmax[1] = 10.0;
}

/* Jimenez
 *
 * f: 2  objective functions
 * g: 4  constraints
 * x: 2  variables
 */
void EMO_Benchmark_jim(EMO_MOP *mop, double *f, double *g, double *x) {
  f[0] = - (5.0 * x[0] + 3.0 * x[1]);
  f[1] = - (2.0 * x[0] + 8.0 * x[1]);

  g[0] = -x[0] - 4.0 * x[1] + 100.0;
  g[1] = -3.0 * x[0] - 2.0 * x[1] + 150.0;
  g[2] = -200.0 + 5.0 * x[0] + 3.0 * x[1];
  g[3] = -75.0 + 2.0 * x[0] + 8.0 * x[1];
}

void EMO_Benchmark_jim_range(EMO_MOP *mop) {
  mop->xmin[0] = 0.0;
  mop->xmin[1] = 0.0;

  mop->xmax[0] = 50;
  mop->xmax[1] = 16;
}

/* Kita
 *
 * f: 2  objective functions
 * g: 3  constraints
 * x: 2  variables
 */
void EMO_Benchmark_kita(EMO_MOP *mop, double *f, double *g, double *x) {
  f[0] = x[0] * x[0] - x[1];
  f[1] = - x[0] / 2.0 - x[1] - 1.0;

  g[0] = -x[0] / 6.0 - x[1] + 13.0 / 2.0;
  g[1] = -x[0] / 2.0 - x[1] + 15.0 / 2.0;
  g[2] = -5.0 * x[0] - x[1] + 30.0;
}

void EMO_Benchmark_kita_range(EMO_MOP *mop) {
  mop->xmin[0] = 0.0;
  mop->xmin[1] = 0.0;

  mop->xmax[0] = 3.0;
  mop->xmax[1] = 6.5;
}

/* Obayashi
 *
 * f: 2  objective functions
 * g: 1  constraints
 * x: 2  variables
 */
void EMO_Benchmark_oba(EMO_MOP *mop, double *f, double *g, double *x) {
  f[0] = x[0];
  f[1] = x[1];

  g[0] = x[0]*x[0] + x[1]*x[1] - 1.0;
}

void EMO_Benchmark_oba_range(EMO_MOP *mop) {
  mop->xmin[0] = 0.0;
  mop->xmax[0] = 1.0;

  mop->xmin[1] = 0.0;
  mop->xmax[1] = 1.0;
}

/* Osyczka
 *
 * f: 2  objective functions
 * g: 2  constraints
 * x: 2  variables
 */
void EMO_Benchmark_osy1(EMO_MOP *mop, double *f, double *g, double *x) {
  f[0] = x[0] + x[1] * x[1];
  f[1] = x[0] * x[0] + x[1];

  g[0] = 12.0 - x[0] - x[1];
  g[1] = x[0] * x[0] + 10.0 * x[0] - x[1] * x[1] + 16.0 * x[1] - 80.0;
}

void EMO_Benchmark_osy1_range(EMO_MOP *mop) {
  mop->xmin[0] = 2.0;
  mop->xmax[0] = 7.0;

  mop->xmin[1] = 5.0;
  mop->xmax[1] = 10.0;
}

/* Osyczka and Kundu (1995) function
 *
 * f: 2 objective functions
 * g: 6 constraints
 * x: 6 variables
 */
void EMO_Benchmark_osy2(EMO_MOP *mop, double *f, double *g, double *x) {

  f[0] = -25.0*(x[0]-2.0)*(x[0]-2.0) - (x[1]-2.0)*(x[1]-2.0) - (x[2]-1.0)*(x[2]-1.0) - (x[3]-4.0)*(x[3]-4.0) - (x[4]-1.0)*(x[4]-1.0);
  f[1] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] + x[3]*x[3] + x[4]*x[4] + x[5]*x[5];

  g[0] = x[0] + x[1] - 2.0;
  g[1] = 6.0 - x[0] - x[1];
  g[2] = 2.0 + x[0] - x[1];
  g[3] = 2.0 - x[0] + 3.0 * x[1];
  g[4] = 4.0 - (x[2]-3.0)*(x[2]-3.0) - x[3];
  g[5] = (x[4]-3.0)*(x[4]-3.0) + x[5] - 4.0;
}

void EMO_Benchmark_osy2_range(EMO_MOP *mop) {
  mop->xmin[0] =  0.0;
  mop->xmax[0] = 10.0;

  mop->xmin[1] =  0.0;
  mop->xmax[1] = 10.0;

  mop->xmin[2] = 1.0;
  mop->xmax[2] = 5.0;

  mop->xmin[3] = 0.0;
  mop->xmax[3] = 6.0;

  mop->xmin[4] = 1.0;
  mop->xmax[4] = 5.0;

  mop->xmin[5] =  0.0;
  mop->xmax[5] = 10.0;
}

/* Srinivas and Deb (1995) function
 *
 * f: 2 objective functions
 * g: 2 constraints
 * x: 2 variables
 *
 * Optimal solutions correspond to
 * x_1 = -2.5, x_2 \in [-14.79, 2.50]
 */
void EMO_Benchmark_srn(EMO_MOP *mop, double *f, double *g, double *x) {
  f[0] = 2.0 + (x[0]-2.0)*(x[0]-2.0) + (x[1]-1.0)*(x[1]-1.0);
  f[1] = 9.0*x[0] - (x[1]-1.0)*(x[1]-1.0);

  g[0] = 225.0 - x[0]*x[0] - x[1]*x[1];
  g[1] = 3.0*x[1] - x[0] - 10.0;
}

void EMO_Benchmark_srn_range(EMO_MOP *mop) {
  mop->xmin[0] = -20.0;
  mop->xmax[0] = 20.0;

  mop->xmin[1] = -20.0;
  mop->xmax[1] = 20.0;
}

/* Tamaki
 *
 * f: 3  objective functions
 * g: 1  constraints
 * x: 3  variables
 */
void EMO_Benchmark_tmk(EMO_MOP *mop, double *f, double *g, double *x) {
  f[0] = x[0];
  f[1] = x[1];
  f[2] = x[2];

  g[0] = x[0]*x[0] + x[1]*x[1] + x[2]*x[2] - 1.0;
}

void EMO_Benchmark_tmk_range(EMO_MOP *mop) {
  mop->xmin[0] = 0.0;
  mop->xmin[1] = 0.0;
  mop->xmin[2] = 0.0;

  mop->xmax[0] = 2.0;  // RHG
  mop->xmax[1] = 2.0;
  mop->xmax[2] = 2.0;
}

/* Tanaka  (1995)
 *
 * f: 2 objective functions
 * g: 2 constraints
 * x: 2 variables
 */
void EMO_Benchmark_tnk(EMO_MOP *mop, double *f, double *g, double *x) {
  double a, b;

  f[0] = x[0];
  f[1] = x[1];

  a = 0.1;
  b = 16.0;

  g[0] = x[0]*x[0] + x[1]*x[1] - 1.0 - a*cos(b * atan(x[1]/x[0]));
  g[1] = 0.5 - (x[0]-0.5)*(x[0]-0.5) - (x[1]-0.5)*(x[1]-0.5);
}

void EMO_Benchmark_tnk_range(EMO_MOP *mop) {
  mop->xmin[0] = 1e-16;
  mop->xmax[0] = M_PI;

  mop->xmin[1] = 1e-16;
  mop->xmax[1] = M_PI;
}

/* Constraint Viennet 4
 *
 * f: 3  objective functions
 * g: 3  constraints
 * x: 2  variables
 */
void EMO_Benchmark_vie4(EMO_MOP *mop, double *f, double *g, double *x) {
  f[0] = (x[0] - 2.0) * (x[0] - 2.0) / 2.0 + (x[1] + 1.0) * (x[1] + 1.0) / 13.0 + 3.0;
  f[1] = (x[0] + x[1] - 3.0) * (x[0] + x[1] - 3.0) / 175.0 + (2.0 * x[1] - x[0]) * (2.0 * x[1] - x[0]) / 17.0 - 13.0;
  f[2] = (3.0 * x[0] - 2.0 * x[1] + 4.0) * (3.0 * x[0] - 2.0 * x[1] + 4.0) / 8.0 + (x[0] - x[1] + 1.0) * (x[0] - x[1] + 1.0) / 27.0 + 15.0;

  g[0] = 4.0 - x[1] - 4.0 * x[0];
  g[1] = x[0] + 1.0;
  g[2] = x[1] - x[0] + 2.0;
}

void EMO_Benchmark_vie4_range(EMO_MOP *mop) {
  mop->xmin[0] = -4.0;
  mop->xmax[0] = 4.0;

  mop->xmin[1] = -4.0;
  mop->xmax[1] = 4.0;
}

/* Constraint test problem DTLZ8 (Deb-Thiele-Laumanns-Zitzler)

   Real variables = 10 * mop->nobj
   Bin variables = 0
   Objectives = (scalable)
   constraints = mop->nobj
   Values of xi in [0,1]
   PFTrue is a combination of a straight line
   and a hyper-plane.
*/
void EMO_Benchmark_dtlz8(EMO_MOP *mop, double *f, double *g, double *x) {
  double v = DBL_MAX;
  int i, j, c;

//  j | (j-1) n/m | j n/m
//------------------------
//  1     0         n/m
//  2     n/m       2 n/m
//  3     2 n/m     3 n/m
//  4     3 n/m     4 n/m

  c = (int) (mop->nvar / mop->nobj);

  for(j = 0; j < mop->nobj; j++) {
    f[j] = 0;

    for(i = j*c; i < (j+1)*c; i++)
      f[j] += x[i];

    f[j] /= (double) c;
  }

  for(i = 0; i < mop->nobj; i++)
    for(j = 0; j < mop->nobj; j++)
      if((f[i] + f[j]) < v)
        v = f[i] + f[j];

  for(j = 0; j < mop->nobj - 1; j++)
    g[j] = f[mop->nobj - 1] + 4.0 * f[j] - 1;

  g[mop->nobj - 1] = 2.0 * f[mop->nobj - 1] + v - 1.0;
}

/* Constraint test problem DTLZ9 (Deb-Thiele-Laumanns-Zitzler)

   Real variables = 10 * mop->nobj
   Bin variables = 0
   Objectives = (scalable)
   constraints = mop->nobj - 1
   Values of xi in [0,1]
   PFTrue is similar to that in DTLZ5.
   However, the density of solutions gets
   thinner towards the Pareto-optimal region.
*/
void EMO_Benchmark_dtlz9(EMO_MOP *mop, double *f, double *g, double *x) {
  int i, j, c;

  c = (int) (mop->nvar / mop->nobj);

  for(j = 0; j < mop->nobj; j++) {
    f[j] = 0;

    for(i = j*c; i < (j+1)*c; i++)
      f[j] += pow(x[i], 0.1);

    f[j] /= (double) c;
  }

  for(j = 0; j < mop->nobj - 1; j++)
    g[j] = f[mop->nobj - 1] * f[mop->nobj - 1] + f[j] * f[j] - 1.0;
}

void EMO_Benchmark_dtlz89_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = 0.0;
    mop->xmax[i] = 1.0;
  }
}

/* NSGA-III, part II; Jain & Deb 2014,
 * original described in Musselman and Talavage 1980
 * and then in Ray 2001.
 *
 * f: 5  objective functions
 * g: 7 constraints
 * x: 3  variables
 */
void EMO_Benchmark_water(EMO_MOP *mop, double *f, double *g, double *x) {
  f[0] = 106780.37 * (x[1] + x[2]) + 61704.67;
  f[1] = 3000.0 * x[0];
  f[2] = 305700.0 * 2289.0 * x[1] / pow(0.06 * 2289.0, 0.65);
  f[3] = 250.0 * 2289.0 * exp(-39.75 * x[1] + 9.9 * x[2] + 2.74);
  f[4] = 25.0 * (1.39 / (x[0]*x[1]) + 4940.0 * x[2] - 80.0);

  g[0] = 1.0 - 0.00139 / (x[0] * x[1]) - 4.94 * x[2] + 0.08;
  g[1] = 1.0 - 0.000306 / (x[0] * x[1]) - 1.082 * x[2] + 0.0986;
  g[2] = 50000.0 - 12.307 / (x[0] * x[1]) - 49408.24 * x[2] - 4051.02;
  g[3] = 16000.0 - 2.098 / (x[0] * x[1]) - 8046.33 * x[2] + 696.71;
  g[4] = 10000.0 - 2.138 / (x[0] * x[1]) - 7883.39 * x[2] + 705.04;
  g[5] = 2000.0 - 0.417 / (x[0] * x[1]) - 1721.26 * x[2] + 136.54;
  g[6] = 550.0 - 0.164 / (x[0] * x[1]) - 631.13 * x[2] + 54.48;
}

void EMO_Benchmark_water_range(EMO_MOP *mop) {
  mop->xmin[0] = 0.01;
  mop->xmax[0] = 0.45;

  mop->xmin[1] = 0.01;
  mop->xmax[1] = 0.10;

  mop->xmin[2] = 0.01;
  mop->xmax[2] = 0.10;
}

/* NSGA-III, part II; Jain & Deb 2014
 *
 * f: 3  objective functions
 * g: 10 constraints
 * x: 7  variables
 */
void EMO_Benchmark_carside_impact(EMO_MOP *mop, double *f, double *g, double *x) {
  double F, Vmbp, Vfd;

  F = 4.72 - 0.5 * x[3] - 0.19 * x[1] * x[2];
  Vmbp = 10.58 - 0.674 * x[0] * x[1] - 0.67275 * x[1];
  Vfd = 16.45 - 0.489 * x[2] * x[6] - 0.843 * x[4] * x[5];

  f[0] = 1.98 + 4.9 * x[0] + 6.67 * x[1] + 6.98 * x[2] +
         4.01 * x[3] + 1.78 * x[4] + 0.00001 * x[5] + 2.73 * x[6];
  f[1] = F;
  f[2] = 0.5 * (Vmbp + Vfd);

  g[0] = 1.0 - 1.16 + 0.3717 * x[1] * x[3] + 0.0092928 * x[2];
  g[1] = 0.32 - 0.261 + 0.0159 * x[0] * x[1] + 0.06486 * x[0] +
         0.019 * x[1] * x[6] - 0.0144 * x[2] * x[4] - 0.0154464 * x[5];
  g[2] = 0.32 - 0.214 - 0.00817 * x[4] + 0.045195 * x[0] + 0.0135168 * x[0] -
         0.03099 * x[1] * x[5] + 0.018 * x[1] * x[6] - 0.007176 * x[2] -
         0.023232 * x[2] + 0.00364 * x[4] * x[5] + 0.018 * x[1] * x[1];
  g[3] = 0.32 - 0.74 + 0.61 * x[1] + 0.031296 * x[2] + 0.031872 * x[6] -
         0.227 * x[1] * x[1];
  g[4] = 32.0 - 28.98 - 3.818 * x[2] + 4.2 * x[0] * x[1] - 1.27296 * x[5] +
         2.68065 * x[6];
  g[5] = 32.0 - 33.86 - 2.95 * x[2] + 5.057 * x[0] * x[1] + 3.795 * x[1] +
         3.4431 * x[6] - 1.45728;
  g[6] = 32.0 - 46.36 + 9.9 * x[1] + 4.4505 * x[0];
  g[7] = 4.0 - F;
  g[8] = 9.9 - Vmbp;
  g[9] = 15.7 - Vfd;
}

void EMO_Benchmark_carside_impact_range(EMO_MOP *mop) {
  mop->xmin[0] = 0.5;
  mop->xmax[0] = 1.5;

  mop->xmin[1] = 0.45;
  mop->xmax[1] = 1.35;

  mop->xmin[2] = 0.5;
  mop->xmax[2] = 1.5;

  mop->xmin[3] = 0.5;
  mop->xmax[3] = 1.5;

  mop->xmin[4] = 0.875;
  mop->xmax[4] = 2.625;

  mop->xmin[5] = 0.4;
  mop->xmax[5] = 1.2;

  mop->xmin[6] = 0.4;
  mop->xmax[6] = 1.2;
}

/* Ray 2001, original proposed in Rao
 *
 * f: 2  objective functions
 * g: 2 constraints
 * x: 2  variables
 */
void EMO_Benchmark_2bar_truss(EMO_MOP *mop, double *f, double *g, double *x) {
  const double rho = 0.283;
  const double h = 100;
  const double P = 10000;
  const double E = 3.0E07;
  const double sigma = 2.0E04;
  double u, v;

//  Amin = 1.0;
//  x0 = x[0] / h;
//  x1 = A / Amin;

  f[0] = 2.0 * rho * h * x[1] * sqrt(1.0 + x[0] * x[0]);
  f[1] = (P * h * pow(1.0 + x[0] * x[0], 1.5)  * pow(1.0 + pow(x[0], 4.0), 0.5)) /
         (2.0 * sqrt(2.0) * E * x[0] * x[0] * x[1]);

  u = pow(1.0 + x[0] * x[0], 0.5);
  v = 2.0 * sqrt(2.0) * x[0] * x[1];

  g[0] = sigma - P * (1.0 + x[0]) * u / v;
  g[1] = sigma - P * (1.0 - x[0]) * u / v;
}

void EMO_Benchmark_2bar_truss_range(EMO_MOP *mop) {
  mop->xmin[0] = 0.10;
  mop->xmax[0] = 2.25;

  mop->xmin[1] = 0.50;
  mop->xmin[1] = 2.50;
}

void EMO_Benchmark_alloc(EMO_MOP *mop, EMO_Param *param, const char *problem) {

  typedef void (*Evalc)(EMO_MOP *mop, double *, double *, double *);
  typedef void (*Eval) (EMO_MOP *mop, double *, double *);
  typedef void (*Range)(EMO_MOP *);
  int i, d, nvar, nobj, ncon, npos;
  char *str;

  /* Vector of objective function evaluation */
  const Eval fdic[] = {
                   EMO_Benchmark_fon1,  EMO_Benchmark_fon2,  EMO_Benchmark_kur,
                   EMO_Benchmark_lau,   EMO_Benchmark_lis,   EMO_Benchmark_mur,
                   EMO_Benchmark_pol,   EMO_Benchmark_qua,
                   EMO_Benchmark_ren1,  EMO_Benchmark_ren2,
                   EMO_Benchmark_sch1,  EMO_Benchmark_sch2,
                   EMO_Benchmark_vie1,  EMO_Benchmark_vie2,  EMO_Benchmark_vie3,
                   EMO_Benchmark_deb1,  EMO_Benchmark_deb2,  EMO_Benchmark_deb3,
                   EMO_Benchmark_oka1,  EMO_Benchmark_oka2,
                   EMO_Benchmark_bnh1,  EMO_Benchmark_bnh3,  EMO_Benchmark_sdd,
                   EMO_Benchmark_stz1,  EMO_Benchmark_stz2,  EMO_Benchmark_stz3,
                   EMO_Benchmark_zdt1,  EMO_Benchmark_zdt2,  EMO_Benchmark_zdt3,
                   EMO_Benchmark_zdt4,  EMO_Benchmark_zdt5,  EMO_Benchmark_zdt6,
                   EMO_Benchmark_uf1,   EMO_Benchmark_uf2,   EMO_Benchmark_uf3,
                   EMO_Benchmark_uf4,   EMO_Benchmark_uf5,   EMO_Benchmark_uf6,
                   EMO_Benchmark_uf7,   EMO_Benchmark_uf8,   EMO_Benchmark_uf9,
                   EMO_Benchmark_uf10,
                   EMO_Benchmark_dtlz1, EMO_Benchmark_dtlz2, EMO_Benchmark_dtlz3,
                   EMO_Benchmark_dtlz4, EMO_Benchmark_dtlz5, EMO_Benchmark_dtlz6,
                   EMO_Benchmark_dtlz7,
                   EMO_Benchmark_dtlz2_convex, EMO_Benchmark_dtlz3_convex,
                   EMO_Benchmark_dtlz4_convex,
                   EMO_Benchmark_dtlz1_minus, EMO_Benchmark_dtlz2_minus, EMO_Benchmark_dtlz3_minus,
                   EMO_Benchmark_dtlz4_minus, EMO_Benchmark_dtlz5_minus, EMO_Benchmark_dtlz6_minus,
                   EMO_Benchmark_dtlz7_minus,
                   EMO_Benchmark_ebn,   EMO_LS_lame,         EMO_LS_mirror,
                   EMO_WFG_wfg1,        EMO_WFG_wfg2,        EMO_WFG_wfg3,
                   EMO_WFG_wfg4,        EMO_WFG_wfg5,        EMO_WFG_wfg6,
                   EMO_WFG_wfg7,        EMO_WFG_wfg8,        EMO_WFG_wfg9,
                   EMO_WFG_wfg1_minus,  EMO_WFG_wfg2_minus,  EMO_WFG_wfg3_minus,
                   EMO_WFG_wfg4_minus,  EMO_WFG_wfg5_minus,  EMO_WFG_wfg6_minus,
                   EMO_WFG_wfg7_minus,  EMO_WFG_wfg8_minus,  EMO_WFG_wfg9_minus
                  }; //,
                   //knapsack };

  /* Functions that specify the range of variables */
  const Range range_dic[] = {
                   EMO_Benchmark_fon1_range,  EMO_Benchmark_fon2_range, EMO_Benchmark_kur_range,
                   EMO_Benchmark_lau_range,   EMO_Benchmark_lis_range,  EMO_Benchmark_mur_range,
                   EMO_Benchmark_pol_range,   EMO_Benchmark_qua_range,
                   EMO_Benchmark_ren_range,   EMO_Benchmark_ren_range,
                   EMO_Benchmark_sch1_range,  EMO_Benchmark_sch2_range,
                   EMO_Benchmark_vie1_range,  EMO_Benchmark_vie2_range, EMO_Benchmark_vie3_range,
                   EMO_Benchmark_deb_range,   EMO_Benchmark_deb_range,  EMO_Benchmark_deb_range,
                   EMO_Benchmark_oka1_range,  EMO_Benchmark_oka2_range,
                   EMO_Benchmark_bnh1_range,  EMO_Benchmark_bnh3_range, EMO_Benchmark_stz_range,
                   EMO_Benchmark_stz_range,   EMO_Benchmark_stz_range,  EMO_Benchmark_stz_range,
                   EMO_Benchmark_zdt_range,   EMO_Benchmark_zdt_range,  EMO_Benchmark_zdt_range,
                   EMO_Benchmark_zdt4_range,  EMO_Benchmark_zdt_range,  EMO_Benchmark_zdt_range,
                   EMO_Benchmark_ufa_range,   EMO_Benchmark_ufa_range,  EMO_Benchmark_ufc_range,
                   EMO_Benchmark_ufd_range,   EMO_Benchmark_ufa_range,  EMO_Benchmark_ufa_range,
                   EMO_Benchmark_ufa_range,   EMO_Benchmark_ufb_range,  EMO_Benchmark_ufb_range,
                   EMO_Benchmark_ufb_range,
                   EMO_Benchmark_dtlz_range,  EMO_Benchmark_dtlz_range, EMO_Benchmark_dtlz_range,
                   EMO_Benchmark_dtlz_range,  EMO_Benchmark_dtlz_range, EMO_Benchmark_dtlz_range,
                   EMO_Benchmark_dtlz_range,
                   EMO_Benchmark_dtlz_range,  EMO_Benchmark_dtlz_range, EMO_Benchmark_dtlz_range,
                   EMO_Benchmark_dtlz_range,  EMO_Benchmark_dtlz_range, EMO_Benchmark_dtlz_range,
                   EMO_Benchmark_dtlz_range,  EMO_Benchmark_dtlz_range, EMO_Benchmark_dtlz_range,
                   EMO_Benchmark_dtlz_range,
                   EMO_Benchmark_ebn_range,   EMO_LS_range,             EMO_LS_range,
                   EMO_WFG_range,             EMO_WFG_range,            EMO_WFG_range,
                   EMO_WFG_range,             EMO_WFG_range,            EMO_WFG_range,
                   EMO_WFG_range,             EMO_WFG_range,            EMO_WFG_range,
                   EMO_WFG_range,             EMO_WFG_range,            EMO_WFG_range,
                   EMO_WFG_range,             EMO_WFG_range,            EMO_WFG_range,
                   EMO_WFG_range,             EMO_WFG_range,            EMO_WFG_range
                  };//,

  /* Number of objectives for each test function, 0 means scalable */
  const int obj[] = {
                   2, 2, 2,    // fon1, fon2, kur
                   2, 2, 2,    // lau, lis, mur
                   2, 2,       // pol, qua
                   2, 2,       // ren1, ren2
                   2, 2,       // sch1, sch2
                   3, 3, 3,    // vie1, vie2, vie3
                   2, 2, 2,    // deb1, deb2, deb3
                   2, 2,       // oka1, oka2
                   2, 3, 0,    // bnh1, bnh3, sdd
                   2, 2, 2,    // stz1, stz2, stz3
                   2, 2, 2,    // zdt1, zdt2, zdt3
                   2, 2, 2,    // zdt4, zdt5, zdt6
                   2, 2, 2,    // uf1, uf2, uf3
                   2, 2, 2,    // uf4, uf5, uf6
                   2, 3, 3,    // uf7, uf8, uf9
                   3,          // uf10
                   0, 0, 0,    // dtlz1, dtlz2, dtlz3
                   0, 0, 0,    // dtlz4, dtlz5, dtlz6
                   0,          // dtlz7
                   0, 0, 0,    // dtlz{2,3,4}-convex
                   0, 0, 0,    // dtlz1_minus, dtlz2_minus, dtlz3_minus
                   0, 0, 0,    // dtlz4_minus, dtlz5_minus, dtlz6_minus
                   0,          // dtlz7_minus
                   2, 0, 0,    // ebn, lame, mirror
                   0, 0, 0,    // wfg1, wfg2, wfg3
                   0, 0, 0,    // wfg4, wfg5, wfg6
                   0, 0, 0,    // wfg7, wfg8, wfg9
                   0, 0, 0,    // wfg1_minus, wfg2_minus, wfg3_minus
                   0, 0, 0,    // wfg4_minus, wfg5_minus, wfg6_minus
                   0, 0, 0,    // wfg7_minus, wfg8_minus, wfg9_minus
                   0 };        // knapsack

  /* Number of variables: 0 means scalable, -1 means fixed */
  const int var[] = {
                   2, 0, 3,    // fon1, fon2, kur
                   2, 2, 2,    // lau, lis, mur
                   2, 0,       // pol, qua
                   2, 2,       // ren1, ren2
                   1, 1,       // sch1, sch2
                   2, 2, 2,    // vie1, vie2, vie3
                   2, 2, 2,    // deb1, deb2, deb3
                   2, 3,       // oka1, oka2
                   2, 2, 0,    // bnh1, bnh3, sdd
                   2, 0, 2,    // stz1, stz2, stz3
                   0, 0, 0,    // zdt1, zdt2, zdt3
                   0, 0, 0,    // zdt4, zdt5, zdt6
                   0, 0, 0,    // uf1, uf2, uf3
                   0, 0, 0,    // uf4, uf5, uf6
                   0, 0, 0,    // uf7, uf8, uf9
                   0,          // uf10
                   -1, -1, -1, // dtlz1, dtlz2, dtlz3
                   -1, -1, -1, // dtlz4, dtlz5, dtlz6
                   -1,         // dtlz7
                   -1, -1, -1, // dtlz{2,3,4}-convex
                   -1, -1, -1, // dtlz1_minus, dtlz2_minus, dtlz3_minus
                   -1, -1, -1, // dtlz4_minus, dtlz5_minus, dtlz6_minus
                   -1,         // dtlz7_minus
                   0, 0, 0,    // ebn, lame, mirror
                   0, 0, 0,    // wfg1, wfg2, wfg3
                   0, 0, 0,    // wfg4, wfg5, wfg6
                   0, 0, 0,    // wfg7, wfg8, wfg9
                   0, 0, 0,    // wfg1_minus, wfg2_minus, wfg3_minus
                   0, 0, 0,    // wfg4_minus, wfg5_minus, wfg6_minus
                   0, 0, 0,    // wfg7_minus, wfg8_minus, wfg9_minus
                   0 };        // knapsack

  /* Default number of variables (only for scalable problems): 0 means not applicable */
  /* For DTLZ and LAME test problems, entries represent the k values for NVAR = NOBJ + k - 1 */
  const int default_var[] = {
                   0, 10, 0,   // fon1, fon2, kur
                   0, 0, 0,    // lau, lis, mur
                   0, 16,      // pol, qua
                   0, 0,       // ren1, ren2
                   0, 0,       // sch1, sch2
                   0, 0, 0,    // vie1, vie2, vie3
                   0, 0, 0,    // deb1, deb2, deb3
                   0, 0,       // oka1, oka2
                   0, 0, 10,   // bnh1, bnh3, sdd
                   0, 10, 0,   // stz1, stz2, stz3
                   30, 30, 30, // zdt1, zdt2, zdt3
                   10, 11, 10, // zdt4, zdt5, zdt6
                   30, 30, 30, // uf1, uf2, uf3
                   30, 30, 30, // uf4, uf5, uf6
                   30, 30, 30, // uf7, uf8, uf9
                   30,         // uf10
                   5, 10, 10,  // dtlz1, dtlz2, dtlz3
                   10, 10, 10, // dtlz4, dtlz5, dtlz6
                   20,         // dtlz7
                   10, 10, 10, // dtlz{2,3,4}-convex
                   5, 10, 10,  // dtlz1_minus, dtlz2_minus, dtlz3_minus
                   10, 10, 10, // dtlz4_minus, dtlz5_minus, dtlz6_minus
                   20,         // dtlz7_minus
                   6, 5, 5,    // ebn, lame, mirror
                   24, 24, 24, // wfg1, wfg2, wfg3
                   24, 24, 24, // wfg4, wfg5, wfg6
                   24, 24, 24, // wfg7, wfg8, wfg9
                   24, 24, 24, // wfg1_minus, wfg2_minus, wfg3_minus
                   24, 24, 24, // wfg4_minus, wfg5_minus, wfg6_minus
                   24, 24, 24, // wfg7_minus, wfg8_minus, wfg9_minus
                   500 };      // knapsack

  /* Test functions  with constraints */

  /* Vector of objective function evaluation */
  const Evalc fdicc[] = { EMO_Benchmark_bel, EMO_Benchmark_bnh2, EMO_Benchmark_bnh4,
                          EMO_Benchmark_jim, EMO_Benchmark_kita, EMO_Benchmark_oba,
                          EMO_Benchmark_osy1, EMO_Benchmark_osy2, EMO_Benchmark_srn,
                          EMO_Benchmark_tmk, EMO_Benchmark_tnk, EMO_Benchmark_vie4,
                          EMO_Benchmark_dtlz8, EMO_Benchmark_dtlz9,
                          EMO_Benchmark_cf1, EMO_Benchmark_cf2, EMO_Benchmark_cf3,
                          EMO_Benchmark_cf4, EMO_Benchmark_cf5, EMO_Benchmark_cf6,
                          EMO_Benchmark_cf7, EMO_Benchmark_cf8, EMO_Benchmark_cf9,
                          EMO_Benchmark_cf10,
                          EMO_Benchmark_water, EMO_Benchmark_carside_impact,
                          EMO_Benchmark_2bar_truss};

  /* Functions that specify the range of variables */
  const Range range_dicc[] = { EMO_Benchmark_bel_range, EMO_Benchmark_bnh2_range, EMO_Benchmark_bnh4_range,
                               EMO_Benchmark_jim_range, EMO_Benchmark_kita_range, EMO_Benchmark_oba_range,
                               EMO_Benchmark_osy1_range, EMO_Benchmark_osy2_range, EMO_Benchmark_srn_range,
                               EMO_Benchmark_tmk_range, EMO_Benchmark_tnk_range, EMO_Benchmark_vie4_range,
                               EMO_Benchmark_dtlz89_range, EMO_Benchmark_dtlz89_range,
                               EMO_Benchmark_ufc_range, EMO_Benchmark_ufa_range, EMO_Benchmark_ufd_range,
                               EMO_Benchmark_ufd_range, EMO_Benchmark_ufd_range, EMO_Benchmark_ufd_range,
                               EMO_Benchmark_ufd_range, EMO_Benchmark_ufe_range, EMO_Benchmark_ufb_range,
                               EMO_Benchmark_ufb_range,
                               EMO_Benchmark_water_range, EMO_Benchmark_carside_impact_range,
                               EMO_Benchmark_2bar_truss_range};

  /* Number of objectives for each test function, 0 means scalable */
  const int objc[] = { 2, 2, 3, // bel,  bnh2, bnh4
                       2, 2, 2, // jim,  kita, oba
                       2, 2, 2, // osy1, osy2, srn
                       3, 2, 3, // tmk,  tnk,  vie4
                       0, 0,    // dtlz8, dtlz9
                       2, 2, 2, // cf1, cf2, cf3,
                       2, 2, 2, // cf4, cf5, cf6,
                       2, 3, 3, // cf7, cf8, cf9,
                       3,       // cf10
                       5, 3,    // water, car-side impact
                       2        // 2 bar truss
                     };

  /* Number of variables: 0 means scalable, -1 means fixed */
  const int varc[] = { 2,   2,   2, // bel,  bnh2, bnh4
                       2,   2,   2, // jim,  kita, oba
                       2,   6,   2, // osy1, osy2, srn
                       3,   2,   2, // tmk,  tnk,  vie4
                       -1, -1,      // dtlz8, dtlz9
                       0,   0,   0, // cf1, cf2, cf3,
                       0,   0,   0, // cf4, cf5, cf6,
                       0,   0,   0, // cf7, cf8, cf9,
                       0,           // cf10
                       3,   7,      // water, car-side impact
                       2            // 2 bar truss
                     };

  /* Default number of variables (only for scalable problems): 0 means not applicable */
  const int default_varc[] = { 0,   0,  0, // bel,  bnh2, bnh4
                               0,   0,  0, // jim,  kita, oba
                               0,   0,  0, // osy1, osy2, srn
                               0,   0,  0, // tmk,  tnk,  vie4
                               10,  10,    // dtlz8, dtlz9
                               10, 10, 10, // cf1, cf2, cf3,
                               10, 10, 10, // cf4, cf5, cf6,
                               10, 10, 10, // cf7, cf8, cf9,
                               10,         // cf10
                               0,   0,     // water, car-side impact
                               0           // 2 bar truss
                             };

  /* Number of constraints. (: 0 means scalable, -1 means fixed) */
  const int con[] = { 2, 2, 2, // bel,  bnh2, bnh4
                      4, 3, 1, // jim,  kita, oba
                      2, 6, 2, // osy1, osy2, srn
                      1, 2, 3, // tmk,  tnk,  vie4
                     -1, -1,   // dtlz8, dtlz9
                      1, 1, 1, // cf1, cf2, cf3,
                      1, 1, 2, // cf4, cf5, cf6,
                      2, 1, 1, // cf7, cf8, cf9,
                      1,       // cf10
                      7, 10,   // water, car-side impact
                      2        // 2 bar truss
                    };

  /* Default number of constraints (only for scalable problems): 0 means not applicable */
  const int default_con[] = { 0, 0, 0, // bel,  bnh2, bnh4
                              0, 0, 0, // jim,  kita, oba
                              0, 0, 0, // osy1, osy2, srn
                              0, 0, 0, // tmk,  tnk,  vie4
                              0, 0,    // dtlz8, dtlz9
                              0, 0, 0, // cf1, cf2, cf3,
                              0, 0, 0, // cf4, cf5, cf6,
                              0, 0, 0, // cf7, cf8, cf9,
                              0,       // cf10
                              0, 0,    // water, car-side impact
                              0        // 2 bar truss
                             };

  /* Function's name is converted to uppercase */
  mop->name = EMO_toupper(problem);

  nvar = nobj = ncon = npos = mop->gamma = 0;

  if(strncmp(mop->name, "WFG", 3) == 0) {
    if(!EMO_Param_get_int(param, &nvar, "wfg_nvar")) {
      printf("Error, wfg_nvar is not defined in the configuration file.\n");
      exit(1);
    }

    if(!EMO_Param_get_int(param, &npos, "wfg_npos")) {
      printf("Error, wfg_npos is not defined in the configuration file.\n");
      exit(1);
    }
  }
  else if(strncmp(mop->name, "EBN", 3) == 0 || strncmp(mop->name, "LAME", 4) == 0 || strncmp(mop->name, "MIRROR", 6) == 0) {
    if(!EMO_Param_get_double(param, &mop->gamma, "lame_gamma")) {
      printf("Error, lame_gamma is not defined in the configuration file.\n");
      exit(1);
    }

    if(strncmp(mop->name, "EBN", 3) != 0) {
      if(!EMO_Param_get_int(param, &d, "lame_difficulty")) {
        printf("Error, lame_difficulty is not defined in the configuration file.\n");
        exit(1);
      }

      EMO_LS_difficulty(mop, d);
    }
 }

  if(strncmp(mop->name, "WFG", 3) != 0) {
    if(!EMO_Param_get_int(param, &nvar, "nvar")) {
      printf("Error, nvar is not defined in the configuration file.\n");
      exit(1);
    }
  }

  if(!EMO_Param_get_int(param, &nobj, "nobj")) {
    printf("Error, nobj is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_int(param, &ncon, "ncon")) {
    printf("Error, ncon is not defined in the configuration file.\n");
    exit(1);
  }

  if(nvar < 0 || nobj < 0) {
    printf("Error, invalid reference to the number of variables or objectives.\n");
    free(mop->name);
    exit(1);
  }

  mop->f = NULL;
  mop->fc = NULL;

  /* Find the function's name in dictionary */
  i = EMO_Dictionary_find(EMO_Benchmark_list, mop->name);

  if(i != -1) {
    mop->f = fdic[i];

    /* Default objective values */
    if(nobj == 0 || obj[i] > 0) {
      nobj = (obj[i] <= 0)? 3 : obj[i]; /* Scalable number of objectives */

      if(nobj != obj[i] && nobj != 0 && obj[i] > 0) /* Verify that is scalable w.r.t. the number of objectives */
        EMO_Debug_printf(param->dbg, "Benchmark: warning, function %s is not scalable with respect to the number of objectives, using the default value: %d.", mop->name, nobj);
    }

    /* Default variable values */
    if(nvar == 0 || var[i] > 0) {
      nvar = (var[i] <= 0)? default_var[i] : var[i]; /* Scalable number of variables */

      if(strncmp(mop->name, "DTLZ", 4) == 0 || strncmp(mop->name, "LAME", 4) == 0) /* Fixed number of variables */
        nvar += nobj - 1;

      if(strncmp(mop->name, "WFG", 3) == 0) /* Fixed number of variables in WFG test suite */
        nvar = 24;

      if(nvar != 0 && var[i] > 0)  /* Verify that is scalable w.r.t. the number of variables */
        EMO_Debug_printf(param->dbg, "Benchmark: warning, function %s is not scalable with respect to the number of variables, using the default value: %d.", mop->name, nvar);
    }

    if(nvar != 0 && var[i] <= 0) {
      if(strncmp(mop->name, "DTLZ", 4) == 0) { /* Fixed number of variables */
        nvar = default_var[i] + nobj - 1;
        EMO_Debug_printf(param->dbg, "Benchmark: warning, fixed number of variables (%d) in %s test problem.", nvar, mop->name);
      }

      /* Number of position-related parameters in WFG test suite */
      if(strncmp(mop->name, "WFG", 3) == 0 && npos == 0) {
        npos = (nobj == 2)? 4 : 2 * nobj - 2;
        EMO_Debug_printf(param->dbg, "Benchmark: warning, number of position parameters in %s is not specified, using the default value %d.", mop->name, npos);
      }
    }
  }

  // constraints
  if (mop->f == NULL) {
    i = EMO_Dictionary_find(EMO_Benchmark_listc, mop->name);

    if(i != -1) {
      mop->fc = fdicc[i];

      /* Default objective values */
      if(nobj == 0 || objc[i] > 0) {
        nobj = (objc[i] <= 0)? 3 : objc[i]; /* Scalable number of objectives */

        if(nobj != objc[i] && nobj != 0 && objc[i] > 0) /* Verify that is scalable w.r.t. the number of objectives */
          EMO_Debug_printf(param->dbg, "Benchmark: warning, function %s is not scalable with respect to the number of objectives, using the default value: %d.", mop->name, nobj);
      }

      /* Default variable values */
      if(nvar == 0 || varc[i] > 0) {
        nvar = (varc[i] <= 0)? default_varc[i] : varc[i]; /* Scalable number of variables */

        if(strncmp(mop->name, "DTLZ8", 5) == 0 || strncmp(mop->name, "DTLZ9", 5) == 0) /* Fixed number of variables in DTLZ test suite */
          nvar *= nobj;

        if(nvar != 0 && varc[i] > 0)  /* Verify that is scalable w.r.t. the number of variables */
          EMO_Debug_printf(param->dbg, "Benchmark: warning, function %s is not scalable with respect to the number of variables, using the default value: %d.", mop->name, nvar);
      }

      /* Default constraints */
      if(ncon == 0 || con[i] > 0) {
        ncon = (con[i] <= 0)? default_con[i] : con[i]; /* Scalable number of constraints */

        if(strncmp(mop->name, "DTLZ8", 5) == 0) /* Fixed number of constraints in DTLZ test suite */
          ncon = nobj;

        if(strncmp(mop->name, "DTLZ9", 5) == 0) /* Fixed number of constraints in DTLZ test suite */
          ncon = nobj - 1;

        if(ncon != 0 && con[i] > 0)
          EMO_Debug_printf(param->dbg, "Benchmark: warning, function %s is not scalable with respect to the number of constraints, using the default value: %d.", mop->name, ncon);
      }
    }
    else {
      printf("Error, unknown function %s in EMO_Benchmark_alloc.\n", mop->name);
      free(mop->name);
      exit(1);
    }
  }

  mop->feval = 0;
  mop->nobj = nobj;
  mop->nvar = nvar;
  mop->ncon = ncon;
  mop->npos = npos;  /* DTLZ k */

  mop->xmin = (double *) calloc(sizeof(double), nvar);
  mop->xmax = (double *) calloc(sizeof(double), nvar);

  if(ncon == 0)
    range_dic[i](mop);
  else
    range_dicc[i](mop);

  if(strncmp(mop->name, "KNAPSACK", 8) == 0 || strncmp(mop->name, "ZDT5", 4) == 0)
    mop->coding = EMO_BINARY;
  else
    mop->coding = EMO_REAL;

  if(strncmp(mop->name, "WFG", 3) == 0) {
    if(EMO_WFG_alloc(mop)) {
      free(mop->name);
      exit(1);
    }
  }

  /* Memory allocation for specific functions */
  if(strncmp(mop->name, "DTLZ5", 5) == 0 || strncmp(mop->name, "DTLZ6", 5) == 0) {
    if((mop->t = (double *) malloc(sizeof(double) * (nobj-1))) == NULL) {
      printf("Not enough memory in benchmark.c\n");
      free(mop->name);
      exit(1);
      /*return NULL;*/
    }
  }

  /* Include parameters to the MOP's name */
  if(strncmp(mop->name, "EBN", 3) == 0 || strncmp(mop->name, "LAME", 4) == 0 || strncmp(mop->name, "MIRROR", 6) == 0) {

    if((str = (char *) malloc(sizeof(char) * 1000)) == NULL) {
      printf("Error, not enough memory in EMO_Benchmark_alloc.\n");
      exit(1);
    }

    if(strncmp(mop->name, "EBN", 3) == 0)
      sprintf(str, "%s_GAMMA_%.1f", mop->name, mop->gamma);
    else
      sprintf(str, "%s_GAMMA_%.1f_DEGREE_%d", mop->name, mop->gamma, d);

    free(mop->name);
    mop->name = str;
  }

  EMO_Debug_printf(param->dbg, "Benchmark:test function %s", mop->name);
  EMO_Debug_printf(param->dbg, "Benchmark:variables: %d", nvar);
  EMO_Debug_printf(param->dbg, "Benchmark:objective functions: %d", nobj);
  EMO_Debug_printf(param->dbg, "Benchmark:constraints: %d", ncon);
  EMO_Debug_printf(param->dbg, "Benchmark:position parameters: %d", npos);
  EMO_Debug_printf(param->dbg, "Benchmark:encoding: %d", mop->coding);
  EMO_Debug_printv(param->dbg, mop->xmin, mop->nvar, "Benchmark:lower bound");
  EMO_Debug_printv(param->dbg, mop->xmax, mop->nvar, "Benchmark:upper bound");
}

/* Free memory */
void EMO_Benchmark_free(EMO_MOP *mop) {

  if(strncmp(mop->name, "WFG", 3) == 0)
    EMO_WFG_free(mop);

  if(strncmp(mop->name, "DTLZ5", 5) == 0 || strncmp(mop->name, "DTLZ6", 5) == 0)
    free(mop->t);

  //if(strncmp(mop->name, "KNAPSACK", 8) == 0)
    //freeKnapsack(&knap);

  free(mop->xmin);
  free(mop->xmax);
  free(mop->name);
}
