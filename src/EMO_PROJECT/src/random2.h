#ifndef __RANDOM_H_
#define __RANDOM_H_

#include "debug.h"

typedef struct {
  long idum;
} EMO_Rand;

double EMO_Rand_prob3(EMO_Rand *rnd);
int EMO_Rand_int1(EMO_Rand *rnd, int a, int b);
int EMO_Rand_flip(EMO_Rand *rnd, double prob);

void EMO_Rand_alloc(EMO_Rand *rnd, EMO_Debug *dbg, unsigned long long seed);
int EMO_Rand_next_seed(EMO_Rand *rnd, int skip);
void EMO_Rand_alloc_from_file(EMO_Rand *rnd, EMO_Debug *dbg, const char *file, int skip);
void EMO_Rand_free(EMO_Rand *rnd);
double EMO_Rand_real1(EMO_Rand *rnd, double a, double b);
void EMO_Rand_shuffle(EMO_Rand *rnd, int *idx, int size);
#endif

