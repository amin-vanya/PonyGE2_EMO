#ifndef _PMOEA_H
#define _PMOEA_H

#include "migration.h"
#include "common.h"
#include "island.h"
#include "moea.h"
#include "list.h"

typedef struct {
  EMO_Island comm;      /* MPI communicator */
  EMO_Migration policy;
  EMO_MOEA moea;
  EMO_List src;  /* Source and destination nodes,   */
  EMO_List dest; /* defined by the logical topology */
  EMO_List lst1;   /* Temporary lists */
  EMO_List lst2;
  int interrupt; /* Enable interruption of the program from external nodes */
  int mpolicy;   /* Migration policy */
  int rpolicy;   /* Replacement policy */
  int sync;      /* Enable synchronous migration */
  int nmig;    /* Number of immigrant individuals */
  int epoch;     /* Frequency of migration events   */
  int nproc;
  int myrank;
  int prune_psize;
} EMO_PMOEA;

void EMO_PMOEA_alloc(EMO_PMOEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, const char *str);
void EMO_PMOEA_free(EMO_PMOEA *alg);
void EMO_PMOEA_run(EMO_PMOEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_PMOEA_write(EMO_PMOEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, int run);

#endif

