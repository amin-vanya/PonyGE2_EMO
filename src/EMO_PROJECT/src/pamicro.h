#ifndef _microSMSEMOA_H
#define _microSMSEMOA_H

#include "migration.h"
#include "smsemoa.h"
#include "common.h"
#include "island.h"
#include "niche.h"
#include "vpath.h"
#include "list.h"

typedef struct {
  EMO_Population archive;
  EMO_Island comm;      /* MPI communicator */
  EMO_Migration mpolicy;
  EMO_Migration rpolicy;
  EMO_SMSEMOA moea;
  EMO_VPath vpath;
  EMO_Prune p;
  EMO_List src;  /* Source and destination nodes,   */
  EMO_List dest; /* defined by the logical topology */
  EMO_List lst1;   /* Temporary lists */
  EMO_List lst2;
  int nmpolicy;   /* Migration policy */
  int nrpolicy;   /* Replacement policy */
  int nmig;    /* Number of immigrant individuals */
  int epoch;     /* Frequency of migration events   */
  int nproc;
  int myrank;
  int xrs;
} EMO_PAMICRO;

void EMO_PAMICRO_alloc(EMO_PAMICRO *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_PAMICRO_free(EMO_PAMICRO *alg);
void EMO_PAMICRO_run(EMO_PAMICRO *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_PAMICRO_write(EMO_PAMICRO *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, int run);

#endif

