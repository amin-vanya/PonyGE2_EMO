
#include <math.h>
#include "cec09.h"
#include "numeric.h"

#define MYSIGN(x) ((x)>0?1.0:-1.0)

/* Test problem UF1

   Real variables = 30 (scalable)
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1] x [-1,1]^(n-1)
   PFTrue is at f2 = 1 - sqrt(f1), f1 in [0,1]
                xj = sin(6*PI*x1 + j*PI/n), j = 2,...,n 0 <= x1 <= 1
   Convex Pareto-optimal front
*/
void EMO_Benchmark_uf1(EMO_MOP *mop, double *f, double *x) {
  int j, count1, count2;
  double sum1, sum2, yj;

  sum1   = sum2   = 0.0;
  count1 = count2 = 0;

  for(j = 2; j <= mop->nvar; j++) {
    yj = x[j-1] - sin(6.0 * PI * x[0] + j * PI / (double) mop->nvar);
    yj = yj * yj;

    if(j % 2 == 0) {
      sum2 += yj;
      count2++;
    } 
    else {
      sum1 += yj;
      count1++;
    }
  }

  f[0] = x[0] + 2.0 * sum1 / (double) count1;
  f[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double) count2;
}

/* Test problem UF2

   Real variables = 30 (scalable)
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1]x[-1,1]^(n-1)
   PFTrue is at f2 = 1 - sqrt(f1), f1 in [0, 1]
   xj = {(0.3*x1^2*cos(24*PI*x[0]+4*j*PI/n)+0.6*x[0])*cos(6*PI*x[0]+j*PI/n) j in J1
         (0.3*x1^2*cos(24*PI*x[0]+4*j*PI/n)+0.6*x[0])*sin(6*PI*x[0]+j*PI/n) j in J2
        0 <= x <= 1
   Convex Pareto-optimal front
*/
void EMO_Benchmark_uf2(EMO_MOP *mop, double *f, double *x) {
  int j, count1, count2;
  double sum1, sum2, yj;

  sum1   = sum2   = 0.0;
  count1 = count2 = 0;

  for(j = 2; j <= mop->nvar; j++) {
    if(j % 2 == 0) {
      yj = x[j-1]-0.3*x[0]*(x[0]*cos(24.0*PI*x[0]+4.0*j*PI/(double) mop->nvar)+2.0)*sin(6.0*PI*x[0]+j*PI/(double) mop->nvar);
      sum2 += yj*yj;
      count2++;
    }
    else {
      yj = x[j-1]-0.3*x[0]*(x[0]*cos(24.0*PI*x[0]+4.0*j*PI/(double) mop->nvar)+2.0)*cos(6.0*PI*x[0]+j*PI/(double) mop->nvar);
      sum1 += yj*yj;
      count1++;
    }
  }

  f[0] = x[0] + 2.0 * sum1 / (double)count1;
  f[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double)count2;
}

/* Test problem UF3

   Real variables = 30 (scalable)
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1]^n
   PFTrue is at f2 = 1 - sqrt(f1), f1 in [0, 1]
   xj = x1^(0.5*(1+3*(j-2) / (n-2))) j = 2,...,n 0 <= x1 <= 1
   Convex Pareto-optimal front
*/
void EMO_Benchmark_uf3(EMO_MOP *mop, double *f, double *x) {
  int j, count1, count2;
  double sum1, sum2, prod1, prod2, yj, pj;

  sum1   = sum2   = 0.0;
  count1 = count2 = 0;
  prod1  = prod2  = 1.0;

  for(j = 2; j <= mop->nvar; j++) {
    yj = x[j-1] - pow(x[0], 0.5 * (1.0 + 3.0 * (j - 2.0) / ((double) mop->nvar - 2.0)));
    pj = cos(20.0 * yj * PI / sqrt((double) j));

    if (j % 2 == 0) {
      sum2  += yj*yj;
      prod2 *= pj;
      count2++;
    } 
    else {
      sum1  += yj*yj;
      prod1 *= pj;
      count1++;
    }
  }

  f[0] = x[0] + 2.0 * (4.0 * sum1 - 2.0 * prod1 + 2.0) / ((double) count1);
  f[1] = 1.0 - sqrt(x[0]) + 2.0 * (4.0 * sum2 - 2.0 * prod2 + 2.0) / ((double) count2);
}

/* Test problem UF4

   Real variables = 30 (scalable)
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1] x [-2,2]^(n-1)
   PFTrue is at f2 = 1 - f1^2, f1 in [0, 1]
                xj = sin(6*PI*x[0] + j * PI/n), j = 2,...,n 0 <= x1 <= 1
   Concave Pareto-optimal front
*/
void EMO_Benchmark_uf4(EMO_MOP *mop, double *f, double *x) {
  int j, count1, count2;
  double sum1, sum2, yj, hj;

  sum1   = sum2   = 0.0;
  count1 = count2 = 0;

  for(j = 2; j <= mop->nvar; j++) {
    yj = x[j-1] - sin(6.0 * PI * x[0] + j * PI / (double) mop->nvar);
    hj = fabs(yj) / (1.0 + exp(2.0 * fabs(yj)));

    if (j % 2 == 0) {
      sum2 += hj;
      count2++;
    } 
    else {
      sum1 += hj;
      count1++;
    }
  }

  f[0] = x[0] + 2.0 * sum1 / (double)count1;
  f[1] = 1.0 - x[0] * x[0] + 2.0 * sum2 / (double)count2;
}
 
/* Test problem UF5

   Real variables = 30 (scalable)
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1] x [-1,1]^(n-1)
   PF has 2N+1 Pareto Optimal solutions: (i/(2N),1-i/(2N))
   Disconnected and linear Pareto-optimal front
*/
void EMO_Benchmark_uf5(EMO_MOP *mop, double *f, double *x) {
  int j, count1, count2;
  double sum1, sum2, yj, hj, N, E;

  sum1   = sum2   = 0.0;
  count1 = count2 = 0;
  N = 10.0; E = 0.1;

  for(j = 2; j <= mop->nvar; j++) {
    yj = x[j-1] - sin(6.0 * PI * x[0] + j * PI / (double) mop->nvar);
    hj = 2.0 * yj * yj - cos(4.0 * PI * yj) + 1.0;

    if (j % 2 == 0) {
      sum2  += hj;
      count2++;
    } 
    else {
      sum1  += hj;
      count1++;
    }
  }

  hj = (0.5 / N + E) * fabs(sin(2.0 * N * PI * x[0]));
  f[0] = x[0] + hj + 2.0 * sum1 / (double)count1;
  f[1] = 1.0 - x[0] + hj + 2.0 * sum2 / (double)count2;
}

/* Test problem UF6

   Real variables = 30 (scalable)
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1] x [-1,1]^(n-1)
   PFTrue is at f2 = 1 - f1, f1 in Union_{i=1}^{N} [(2i-1)/2N, 2i/2N]
   Disconnected Pareto-optimal front
*/
void EMO_Benchmark_uf6(EMO_MOP *mop, double *f, double *x) {
  int j, count1, count2;
  double sum1, sum2, prod1, prod2, yj, hj, pj, N, E;

  N = 2.0; E = 0.1;
  sum1   = sum2   = 0.0;
  count1 = count2 = 0;
  prod1  = prod2  = 1.0;

  for(j = 2; j <= mop->nvar; j++) {

    yj = x[j-1] - sin(6.0 * PI * x[0] + j * PI / (double) mop->nvar);
    pj = cos(20.0 * yj * PI / sqrt((double) j));

    if (j % 2 == 0) {
      sum2  += yj * yj;
      prod2 *= pj;
      count2++;
    }
    else {
      sum1  += yj * yj;
      prod1 *= pj;
      count1++;
    }
  }

  hj = 2.0 * (0.5 / N + E) * sin(2.0 * N * PI * x[0]);

  if(hj < 0.0) hj = 0.0;

  f[0] = x[0] + hj + 2.0 * (4.0 * sum1 - 2.0 * prod1 + 2.0) / (double)count1;
  f[1] = 1.0 - x[0] + hj + 2.0 * (4.0 * sum2 - 2.0 * prod2 + 2.0) / (double)count2;
}

/* Test problem UF7

   Real variables = 30 (scalable)
   Bin variables = 0
   Objectives = 2
   Constraints = 0
   Values of x are in [0,1] x [-1,1]^(n-1)
   PFTrue is at f2 = 1 - f1, f1 in [0,1]
                xj = sin(6*PI*x[0] + j*PI/n), j = 2,...,n, x1 in [0,1]
   Linear Pareto-optimal front
*/
void EMO_Benchmark_uf7(EMO_MOP *mop, double *f, double *x) {
  int j, count1, count2;
  double sum1, sum2, yj;

  sum1   = sum2   = 0.0;
  count1 = count2 = 0;

  for(j = 2; j <= mop->nvar; j++) {
    yj = x[j-1] - sin(6.0 * PI * x[0] + j * PI / (double) mop->nvar);

    if (j % 2 == 0) {
      sum2  += yj * yj;
      count2++;
    } 
    else {
      sum1  += yj * yj;
      count1++;
    }
  }

  yj = pow(x[0], 0.2);

  f[0] = yj + 2.0 * sum1 / (double) count1;
  f[1] = 1.0 - yj + 2.0 * sum2 / (double) count2;
}

/* Test problem UF8

   Real variables = 30 (scalable)
   Bin variables = 0
   Objectives = 3
   Constraints = 0
   Values of x are in [0,1]^2 x [-2,2]^(n-2)
   PFTrue is at f1^2 + f2^2 + f3^2 = 1, 0 <= f1,f2,f3 <= 1
   xj = 2*x2*sin(2*PI*x1 + j*PI/n), j = 3,...,n
   Concave Pareto-optimal front
*/
void EMO_Benchmark_uf8(EMO_MOP *mop, double *f, double *x) {
  int j, count1, count2, count3;
  double sum1, sum2, sum3, yj;

  sum1   = sum2   = sum3   = 0.0;
  count1 = count2 = count3 = 0;

  for(j = 3; j <= mop->nvar; j++) {
    yj = x[j-1] - 2.0 * x[1] * sin(2.0 * PI * x[0] + j * PI / (double) mop->nvar);

    if(j % 3 == 1) {
      sum1  += yj * yj;
      count1++;
    }
    else if(j % 3 == 2) {
      sum2  += yj * yj;
      count2++;
    }
    else {
      sum3  += yj * yj;
      count3++;
    }
  }

  f[0] = cos(0.5*PI*x[0]) * cos(0.5*PI*x[1]) + 2.0 * sum1 / (double)count1;
  f[1] = cos(0.5*PI*x[0]) * sin(0.5*PI*x[1]) + 2.0 * sum2 / (double)count2;
  f[2] = sin(0.5*PI*x[0]) + 2.0 * sum3 / (double)count3;
}

/* Test problem UF9

   Real variables = 30 (scalable)
   Bin variables = 0
   Objectives = 3
   Constraints = 0
   Values of x are in [0,1]^2 x [-2,2]^(n-2)
   PFTrue has two parts. The first part is
      0 <= f3 <= 1, 0 <= f1 <= 1/4 (1 - f3), f2 = 1 - f1 - f3
   and the second one is
      0 <= f3 <= 1, 3/4 (1 - f3) <= f1 <= 1, f2 = 1 - f1 - f3
   The PS also has two disconnected parts:
      x1 in [0,0.25] Union [0.75,1], 0 <= x2 <= 1,
      xj = 2*x2*sin(2*PI*x1 + j*PI/n), j = 3,...,n
   Disconnected Pareto-optimal front
*/
void EMO_Benchmark_uf9(EMO_MOP *mop, double *f, double *x) {
  int j, count1, count2, count3;
  double sum1, sum2, sum3, yj, E;

  E = 0.1;
  sum1   = sum2   = sum3   = 0.0;
  count1 = count2 = count3 = 0;

  for(j = 3; j <= mop->nvar; j++) {
    yj = x[j-1] - 2.0 * x[1] * sin(2.0 * PI * x[0] + j * PI / (double) mop->nvar);

    if(j % 3 == 1) {
      sum1  += yj * yj;
      count1++;
    } 
    else if(j % 3 == 2) {
      sum2  += yj * yj;
      count2++;
    } 
    else {
      sum3  += yj * yj;
      count3++;
    }
  }

  yj = (1.0 + E) * (1.0 - 4.0 * (2.0 * x[0] - 1.0) * (2.0 * x[0] - 1.0));

  if(yj < 0.0) yj = 0.0;

  f[0] = 0.5 * (yj + 2.0 * x[0]) * x[1] + 2.0 * sum1 / (double)count1;
  f[1] = 0.5 * (yj - 2.0 * x[0] + 2.0) * x[1] + 2.0 * sum2 / (double)count2;
  f[2] = 1.0 - x[1] + 2.0 * sum3 / (double)count3;
}

/* Test problem UF10

   Real variables = 30 (scalable)
   Bin variables = 0
   Objectives = 3
   Constraints = 0
   Values of x are in [0,1]^2 x [-2,2]^(n-2)
   PFTrue is f1^2 + f2^2 + f3^2 = 1, 0 <= f1,f2,f3 <= 1
             xj = 2*x2*sin(2*PI*x1 + j*PI/n), j = 3,...,n
   Concave Pareto-optimal front
*/
void EMO_Benchmark_uf10(EMO_MOP *mop, double *f, double *x) {
  int j, count1, count2, count3;
  double sum1, sum2, sum3, yj, hj;

  sum1   = sum2   = sum3   = 0.0;
  count1 = count2 = count3 = 0;

  for(j = 3; j <= mop->nvar; j++) {
    yj = x[j-1] - 2.0 * x[1] * sin(2.0 * PI * x[0] + j * PI / (double) mop->nvar);
    hj = 4.0 * yj * yj - cos(8.0 * PI * yj) + 1.0;

    if(j % 3 == 1) {
      sum1  += hj;
      count1++;
    }
    else if(j % 3 == 2) {
      sum2  += hj;
      count2++;
    }
    else {
      sum3  += hj;
      count3++;
    }
  }

  f[0] = cos(0.5*PI*x[0]) * cos(0.5*PI*x[1]) + 2.0 * sum1 / (double)count1;
  f[1] = cos(0.5*PI*x[0]) * sin(0.5*PI*x[1]) + 2.0 * sum2 / (double)count2;
  f[2] = sin(0.5*PI*x[0]) + 2.0 * sum3 / (double) count3;
}

/* Test problem CF1

   Real variables = 10 (scalable)
   Bin variables = 0
   Objectives = 2 
   Constraints = 1 
   Values of x are in [0,1]^n
   PFTrue consists of 2N + 1 points: 
     (i/(2N), 1-i/(2N), i = 0, 1, ..., 2N
   Linear, disconnected Pareto-optimal front
*/
void EMO_Benchmark_cf1(EMO_MOP *mop, double *f, double *g, double *x) {
  double sum1, sum2, yj, N, a;
  int j, count1, count2;

  N = 10.0; a = 1.0;
  sum1 = sum2 = 0.0;
  count1 = count2 = 0;

  for(j = 2; j <= mop->nvar; j++) {
    yj = x[j-1]-pow(x[0], 0.5 * (1.0 + 3.0 * (j - 2.0)/((double) mop->nvar - 2.0)));

    if(j % 2 == 1) {
      sum1  += yj * yj;
      count1++;
    } 
    else {
      sum2  += yj * yj;
      count2++;
    }
  }

  f[0] = x[0] + 2.0 * sum1 / (double)count1;
  f[1] = 1.0 - x[0] + 2.0 * sum2 / (double)count2;
  g[0] = f[0] + f[1] - a * fabs(sin(N * PI * (f[0] - f[1] + 1.0))) - 1.0; 
}

/* Test problem CF2

   Real variables = 10 (scalable)
   Bin variables = 0
   Objectives = 2 
   Constraints = 1 
   Values of x are in [0,1] x [-1,1]^(n-1)
   Its PF in the objective space consist of
    + an isolated Pareto optimal solution (0,1) 
      in the objective space, and
    + N disconnected parts, the i-th part is
      f2 = 1 - sqrt(f1), ((2i-1)/(2N))^2 <= f1 <= (2i/(2N))^2, i =1,...,N.
*/
void EMO_Benchmark_cf2(EMO_MOP *mop, double *f, double *g, double *x) {
  double sum1, sum2, yj, N, a, t;
  int j, count1, count2;

  N = 2.0; a = 1.0;
  sum1 = sum2 = 0.0;
  count1 = count2 = 0;

  for(j = 2; j <= mop->nvar; j++) {
    if(j % 2 == 1) {
      yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/(double) mop->nvar);
      sum1  += yj*yj;
      count1++;
    } 
    else {
      yj = x[j-1] - cos(6.0*PI*x[0] + j*PI/(double) mop->nvar);
      sum2  += yj*yj;
      count2++;
    }
  }

  f[0] = x[0] + 2.0 * sum1 / (double) count1;
  f[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double) count2;
  t = f[1] + sqrt(f[0]) - a * sin(N * PI * (sqrt(f[0]) - f[1] + 1.0)) - 1.0;
  g[0] = MYSIGN(t) * fabs(t)/(1.0 + exp(4.0 * fabs(t)));
}

/* Test problem CF3

   Real variables = 10 (scalable)
   Bin variables = 0
   Objectives = 2 
   Constraints = 1 
   Values of x are in [0,1] x [-2,2]^(n-1)
   Its PF in the objective space consist of
    + an isolated Pareto optimal solution (0,1) 
      in the objective space, and
    + N disconnected parts, the i-th part is
      f2 = 1 - f1^2, sqrt((2i-1)/(2N)) <= f1 <= sqrt(2i/(2N)), i =1,...,N.
*/
void EMO_Benchmark_cf3(EMO_MOP *mop, double *f, double *g, double *x) {
  double sum1, sum2, prod1, prod2, yj, pj, N, a;
  int j, count1, count2;

  N = 2.0; a = 1.0;
  sum1 = sum2 = 0.0;
  count1 = count2 = 0;
  prod1 = prod2 = 1.0;

  for(j = 2; j <= mop->nvar; j++) {
    yj = x[j-1] - sin(6.0 * PI * x[0] + j * PI / (double) mop->nvar);
    pj = cos(20.0 * yj * PI /sqrt((double) j));

    if(j % 2 == 0) {
      sum2 += yj * yj;
      prod2 *= pj;
      count2++;
    } 
    else {
      sum1 += yj * yj;
      prod1 *= pj;
      count1++;
    }
  }

  f[0] = x[0] + 2.0 * (4.0 * sum1 - 2.0 * prod1 + 2.0) / (double)count1;
  f[1] = 1.0 - x[0] * x[0] + 2.0 * (4.0 * sum2 - 2.0 * prod2 + 2.0) / (double)count2;
  g[0] = f[1] + f[0] * f[0] - a * sin(N * PI * (f[0] * f[0] - f[1] + 1.0)) - 1.0;
}

/* Test problem CF4

   Real variables = 10 (scalable)
   Bin variables = 0
   Objectives = 2 
   Constraints = 1 
   Values of x are in [0,1] x [-2,2]^(n-1)
   The PF in the objective space is:
    f2 := { 1 -f1          if 0   <= f1 <= 0.5
            -0.5 f1 + 3/4  if 0.5  < f1 <= 0.75
            1 - f1 + 0.125 if 0.75 < f1 <= 1 }.
*/
void EMO_Benchmark_cf4(EMO_MOP *mop, double *f, double *g, double *x) {
  double sum1, sum2, yj, t;
  int j;

  sum1 = sum2 = 0.0;

  for(j = 2; j <= mop->nvar; j++) {
    yj = x[j-1] - sin(6.0 * PI * x[0] + j * PI/(double) mop->nvar);

    if(j % 2 == 1) {
      sum1  += yj * yj;
    } 
    else {
      if(j == 2)
        sum2 += yj < 1.5 - 0.75 * sqrt(2.0) ? fabs(yj) : (0.125 + (yj - 1.0) * (yj - 1.0));
      else
        sum2 += yj * yj;
    }
  }

  f[0] = x[0] + sum1;
  f[1] = 1.0 - x[0] + sum2;
  t = x[1] - sin(6.0 * x[0] * PI + 2.0 * PI/(double) mop->nvar) - 0.5 * x[0] + 0.25;
  g[0] = MYSIGN(t) * fabs(t) / (1.0 + exp(4.0 * fabs(t)));
}

/* Test problem CF5

   Real variables = 10 (scalable)
   Bin variables = 0
   Objectives = 2 
   Constraints = 1 
   Values of x are in [0,1] x [-2,2]^(n-1)
   The PF in the objective space is:
   f2 := { 1 -f1          if 0   <= f1 <= 0.5
           -0.5 f1 + 3/4  if 0.5  < f1 <= 0.75
           1 - f1 + 0.125 if 0.75 < f1 <= 1 }. 
*/
void EMO_Benchmark_cf5(EMO_MOP *mop, double *f, double *g, double *x) {
  double sum1, sum2, yj;
  int j;

  sum1 = sum2 = 0.0;

  for(j = 2; j <= mop->nvar; j++) {
    if(j % 2 == 1) {
      yj = x[j-1] - 0.8 * x[0] * cos(6.0 * PI * x[0] + j * PI/(double) mop->nvar);
      sum1 += 2.0 * yj * yj - cos(4.0 * PI * yj) + 1.0;
    } 
    else {
      yj = x[j-1] - 0.8 * x[0] * sin(6.0 * PI * x[0] + j * PI/(double) mop->nvar);

      if(j == 2)
        sum2 += (yj < (1.5 - 0.75 * sqrt(2.0))) ? fabs(yj) : (0.125 + (yj - 1.0) * (yj - 1.0));
      else
        sum2 += 2.0 * yj * yj - cos(4.0 * PI * yj) + 1.0;
    }
  }

  f[0] = x[0] + sum1;
  f[1] = 1.0 - x[0] + sum2;
  g[0] = x[1] - 0.8 * x[0] * sin(6.0 * x[0] * PI + 2.0 * PI/(double) mop->nvar) - 0.5 * x[0] + 0.25;
}

/* Test problem CF6

   Real variables = 10 (scalable)
   Bin variables = 0
   Objectives = 2 
   Constraints = 2 
   Values of x are in [0,1] x [-2,2]^(n-1)
   The PF is:
   f2 := { (1 - f1)^2        if 0   <= f1 <= 0.5
           0.5 (1 - f1)      if 0.5  < f1 <= 0.75
           0.25 sqrt(1 - f1) if 0.75 < f1 <= 1. }
*/
void EMO_Benchmark_cf6(EMO_MOP *mop, double *f, double *g, double *x) {
  double sum1, sum2, yj;
  int j;

  sum1 = sum2 = 0.0;

  for(j = 2; j <= mop->nvar; j++) {
    if(j % 2 == 1) {
      yj = x[j-1] - 0.8 * x[0] * cos(6.0 * PI * x[0] + j * PI/(double) mop->nvar);
      sum1  += yj * yj;
    } 
    else {
      yj = x[j-1] - 0.8 * x[0] * sin(6.0 * PI * x[0] + j * PI/(double) mop->nvar);
      sum2  += yj * yj;
    }
  }

  f[0] = x[0] + sum1;
  f[1] = (1.0 - x[0]) * (1.0 - x[0]) + sum2;
  g[0] = x[1] - 0.8 * x[0] * sin(6.0 * x[0] * PI + 2.0 * PI/(double) mop->nvar) - MYSIGN((x[0] - 0.5) * (1.0 - x[0])) * sqrt(fabs((x[0] - 0.5) * (1.0 - x[0])));
  g[1] = x[3] - 0.8 * x[0] * sin(6.0 * x[0] * PI + 4.0 * PI/(double) mop->nvar) - MYSIGN(0.25 * sqrt(1.0 - x[0]) - 0.5 * (1.0 - x[0])) * sqrt(fabs(0.25 * sqrt(1.0 - x[0]) - 0.5 * (1.0 - x[0])));
}

/* Test problem CF7

   Real variables = 10 (scalable)
   Bin variables = 0
   Objectives = 2
   Constraints = 2
   Values of x are in [0,1] x [-2,2]^(n-1)
   The PF is:
   f2 := { (1 - f1)^2        if 0   <= f1 <= 0.5
           0.5 (1 - f1)      if 0.5  < f1 <= 0.75
           0.25 sqrt(1 - f1) if 0.75 < f1 <= 1. }

*/
void EMO_Benchmark_cf7(EMO_MOP *mop, double *f, double *g, double *x) {
  double sum1, sum2, yj;
  int j;

  sum1 = sum2 = 0.0;

  for(j = 2; j <= mop->nvar; j++) {
    if(j % 2 == 1) {
      yj = x[j-1] - cos(6.0 * PI * x[0] + j * PI/(double) mop->nvar);
      sum1 += 2.0 * yj * yj - cos(4.0 * PI * yj) + 1.0;
    } 
    else {
      yj = x[j-1] - sin(6.0 * PI * x[0] + j * PI / (double) mop->nvar);
      if(j == 2 || j == 4)
        sum2 += yj * yj;
      else
        sum2 += 2.0 * yj * yj - cos(4.0 * PI * yj) + 1.0;
    }
  }

  f[0] = x[0] + sum1;
  f[1] = (1.0 - x[0]) * (1.0 - x[0]) + sum2;
  g[0] = x[1] - sin(6.0 * x[0] * PI + 2.0 * PI/(double) mop->nvar) - MYSIGN((x[0] - 0.5) * (1.0 - x[0])) * sqrt(fabs((x[0] - 0.5) * (1.0 - x[0])));
  g[1] = x[3] - sin(6.0 * x[0] * PI + 4.0 * PI/(double) mop->nvar) - MYSIGN(0.25 * sqrt(1.0 - x[0]) - 0.5 * (1.0 - x[0])) * sqrt(fabs(0.25 * sqrt(1.0 - x[0]) - 0.5 * (1.0 - x[0])));
}

/* Test problem CF8

   Real variables = 10 (scalable)
   Bin variables = 0
   Objectives = 3 
   Constraints = 1 
   Its PF will have 2N + 1 disconnected parts:
    f1 = (i/(2N) (1 - f3^2)^(1/2)
    f2 = (1 - f1^2 - f3^2)^(1/2)
    0 <= f3 <= 1
*/
void EMO_Benchmark_cf8(EMO_MOP *mop, double *f, double *g, double *x) {
  double sum1, sum2, sum3, yj, N, a;
  int j, count1, count2, count3;

  N = 2.0; a = 4.0;
  sum1 = sum2 = sum3 = 0.0;
  count1 = count2 = count3 = 0;

  for(j = 3; j <= mop->nvar; j++) {
    yj = x[j-1] - 2.0 * x[1] * sin(2.0 * PI * x[0] + j * PI/(double) mop->nvar);

    if(j % 3 == 1) {
      sum1  += yj * yj;
      count1++;
    } 
    else if(j % 3 == 2) {
      sum2 += yj * yj;
      count2++;
    }
    else {
      sum3 += yj * yj;
      count3++;
    }
  }

  f[0] = cos(0.5 * PI * x[0]) * cos(0.5 * PI * x[1]) + 2.0 * sum1 / (double) count1;
  f[1] = cos(0.5 * PI * x[0]) * sin(0.5 * PI * x[1]) + 2.0 * sum2 / (double) count2;
  f[2] = sin(0.5 * PI * x[0]) + 2.0 * sum3 / (double) count3;
  g[0] = (f[0] * f[0] + f[1] * f[1])/(1.0 - f[2] * f[2]) - a * fabs(sin(N * PI * ((f[0] * f[0] - f[1] * f[1]) / (1.0 - f[2] * f[2]) + 1.0))) - 1.0;
}

/* Test problem CF9

   Real variables = 10 (scalable)
   Bin variables = 0
   Objectives = 3 
   Constraints = 1 
   Values of x are in [0,1]^2 x [-2,2]^(n-2)

   Its PF consists of:
    + a curve  f1 = 0
               0 <= f2 <= 1
               f3 = (1 - f2^2)^(1/2)

    + N disconnected nonlinear 2-D surfaces, the i-th one is:
               0 <= f3 <= 1
               {(2i - 1)/(2N) (1-f3^2)}^(1/2) <= f1 <= {2i/(2N) (1 - f3^2)}^(1/2)
               f2 = [1 - f1^2 - f2^2 ]^(1/2).
*/
void EMO_Benchmark_cf9(EMO_MOP *mop, double *f, double *g, double *x) {
  double sum1, sum2, sum3, yj, N, a;
  int j, count1, count2, count3;

  N = 2.0; a = 3.0;

  sum1 = sum2 = sum3 = 0.0;
  count1 = count2 = count3 = 0;

  for(j = 3; j <= mop->nvar; j++) {
    yj = x[j-1] - 2.0 * x[1] * sin(2.0 * PI * x[0] + j * PI/(double) mop->nvar);

    if(j % 3 == 1) {
      sum1  += yj * yj;
      count1++;
    } 
    else if(j % 3 == 2) {
      sum2  += yj * yj;
      count2++;
    }
    else {
      sum3  += yj * yj;
      count3++;
    }
  }

  f[0] = cos(0.5 * PI * x[0]) * cos(0.5 * PI * x[1]) + 2.0 * sum1 / (double) count1;
  f[1] = cos(0.5 * PI * x[0]) * sin(0.5 * PI * x[1]) + 2.0 * sum2 / (double) count2;
  f[2] = sin(0.5 * PI * x[0]) + 2.0 * sum3 / (double) count3;
  g[0] = (f[0] * f[0] + f[1] * f[1]) / (1.0 - f[2] * f[2]) - a * sin(N * PI * ((f[0] * f[0] - f[1] * f[1]) / (1.0 - f[2] * f[2]) + 1.0)) - 1.0;
}

/* Test problem CF10

   Real variables = 10 (scalable)
   Bin variables = 0
   Objectives = 3 
   Constraints = 1 
   Values of x are in [0,1]^2 x [-2,2]^(n-2)
   Its PF consists of:
    + a curve
              f1 = 0
              0 <= f2 <= 1
              f3 = (1-f2^2)^(1/2)
    + N disconnected nonlinear 2-D surfaces, the i-th one is:
              0 <= f3 <= 1
              {(2i - 1)/(2N) (1-f3^2)}^(1/2) <= f1 <= {(2i)/(2N)(1-f3^2)}^(1/2)
              f2 = [1 - f1^2 - f2^2]^(1/2).
*/
void EMO_Benchmark_cf10(EMO_MOP *mop, double *f, double *g, double *x) {
  double sum1, sum2, sum3, yj, hj, N, a;
  int j, count1, count2, count3;

  N = 2.0; a = 1.0;

  sum1 = sum2 = sum3 = 0.0;
  count1 = count2 = count3 = 0;

  for(j = 3; j <= mop->nvar; j++) {
    yj = x[j-1] - 2.0 * x[1] * sin(2.0 * PI *x[0] + j * PI/(double) mop->nvar);
    hj = 4.0 * yj * yj - cos(8.0 * PI * yj) + 1.0;

    if(j % 3 == 1) {
      sum1 += hj;
      count1++;
    } 
    else if(j % 3 == 2) {
      sum2 += hj;
      count2++;
    }
    else {
      sum3 += hj;
      count3++;
    }
  }

  f[0] = cos(0.5 * PI * x[0]) * cos(0.5 * PI * x[1]) + 2.0 * sum1 / (double) count1;
  f[1] = cos(0.5 * PI * x[0]) * sin(0.5 * PI * x[1]) + 2.0 * sum2 / (double) count2;
  f[2] = sin(0.5 * PI * x[0]) + 2.0 * sum3 / (double) count3;
  g[0] = (f[0] * f[0] + f[1] * f[1]) / (1.0 - f[2] * f[2]) - a * sin(N * PI * ((f[0] * f[0] - f[1] * f[1]) / (1.0 - f[2] * f[2]) + 1.0)) - 1.0;
}

/* uf1 uf2 uf5 uf6 uf7 cf2 */
void EMO_Benchmark_ufa_range(EMO_MOP *mop) {
  int i;

  mop->xmin[0] = 0.0;
  mop->xmax[0] = 1.0;

  for(i = mop->nvar - 1; i > 0; i--) {
    mop->xmin[i] = -1.0;
    mop->xmax[i] =  1.0;
  }
}

/* uf8, uf9, uf10 cf9 cf10 */
void EMO_Benchmark_ufb_range(EMO_MOP *mop) {
  int i;

  mop->xmin[0] = mop->xmin[1] = 0.0;
  mop->xmax[0] = mop->xmax[1] = 1.0;

  for(i = mop->nvar - 1; i > 1; i--) {
    mop->xmin[i] = -2.0;
    mop->xmax[i] = 2.0;
  }
}

/* cf1 */
void EMO_Benchmark_ufc_range(EMO_MOP *mop) {
  int i;

  for(i = mop->nvar - 1; i > -1; i--) {
    mop->xmin[i] = 0.0;
    mop->xmax[i] = 1.0;
  }
}

/* cf3 cf4 cf5 cf6 cf7 */
void EMO_Benchmark_ufd_range(EMO_MOP *mop) {
  int i;

  mop->xmin[0] = 0.0;
  mop->xmax[0] = 1.0;

  for(i = mop->nvar - 1; i > 0; i--) {
    mop->xmin[i] = -2.0;
    mop->xmax[i] = 2.0;
  }
}

/* cf8 */
void EMO_Benchmark_ufe_range(EMO_MOP *mop) {
  int i;

  mop->xmin[0] = mop->xmin[1] = 0.0;
  mop->xmax[0] = mop->xmax[1] = 1.0;

  for(i = mop->nvar - 1; i > 1; i--) {
    mop->xmin[i] = -4.0;
    mop->xmax[i] = 4.0;
  }
}

#undef MYSIGN

