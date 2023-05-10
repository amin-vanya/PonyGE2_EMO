
#include <math.h>
#include <stdlib.h>
#include "random2.h"

// Random number generation of MOEA/D
/*------Constants for rnd_uni()--------------------------------------------*/

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


//the random generator in [0,1)
double rnd_uni(EMO_Rand *rnd) {
  long j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (rnd->idum <= 0) {
    if (-(rnd->idum) < 1) rnd->idum=1;
    else rnd->idum = -(rnd->idum);
    idum2=(rnd->idum);

    for (j=NTAB+7;j>=0;j--) {
      k=(rnd->idum)/IQ1;
      rnd->idum=IA1*(rnd->idum-k*IQ1)-k*IR1;

      if (rnd->idum < 0) rnd->idum += IM1;
      if (j < NTAB) iv[j] = rnd->idum;
    }

    iy=iv[0];
  }

  k=(rnd->idum)/IQ1;
  rnd->idum=IA1*(rnd->idum-k*IQ1)-k*IR1;

  if (rnd->idum < 0) rnd->idum += IM1;

  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;

  if (idum2 < 0) idum2 += IM2;

  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = rnd->idum;

  if (iy < 1) iy += IMM1;

  if ((temp=AM*iy) > RNMX) return RNMX;

  else return temp;
}

double EMO_Rand_prob3(EMO_Rand *rnd) {
  //the random generator in [0,1)
  return rnd_uni(rnd);
}

/* RHG: generates a random number on [0,a]-interval */
int EMO_Rand_int1(EMO_Rand *rnd, int a, int b) {
  return (int) (rnd_uni(rnd) * b);
}

int EMO_Rand_flip(EMO_Rand *rnd, double prob) {
  if(rnd_uni(rnd) <= prob)
    return 1;
  return 0;
}

void EMO_Rand_alloc(EMO_Rand *rnd, EMO_Debug *dbg, unsigned long long seed) {

  rnd->idum = -(long) ((237 + 111)%1235);
}

int EMO_Rand_next_seed(EMO_Rand *rnd, int skip) {
return 1;
}

void EMO_Rand_alloc_from_file(EMO_Rand *rnd, EMO_Debug *dbg, const char *file, int skip) {

  rnd->idum = -(long) ((237 + 111)%1235);
}

void EMO_Rand_free(EMO_Rand *rnd){ 
}

double EMO_Rand_real1(EMO_Rand *rnd, double a, double b) {
  return a + rnd_uni(rnd) * (b -a);
}

void EMO_Rand_shuffle(EMO_Rand *rnd, int *idx, int size) {
  int i, r, swap;

  for(i = size - 1; i > 0; i--) {
    r = EMO_Rand_int1(rnd, 0, i);
    swap = idx[i];
    idx[i] = idx[r];
    idx[r] = swap;
  }
}

/*
int main() {
  EMO_Rand rnd;
  int i;

  EMO_Rand_alloc(&rnd, NULL, 1);

  for(i = 0; i < 100; i++) {

    printf("%f ", EMO_Rand_prob3(&rnd)); 
    printf("%d ", EMO_Rand_flip(&rnd, 0.5)); 
    printf("%d ", EMO_Rand_int1(&rnd, 2, 50));
    printf("%f\n", EMO_Rand_real1(&rnd, 2, 50)); 
  }
  EMO_Rand_free(&rnd);
  return 0;
}
*/

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


