
#include "common.h"
#include "param.h"
#include "smsemoa.h"
#include "nsga2.h"
#include "nsga3.h"
#include "mombi2.h"
#include "mombi3.h"
#include "ibea.h"
#include "moead.h"
#include "spea2.h"
#include "hype.h"
#include "movap.h"

#ifndef _MOEA_H
#define _MOEA_H

/* Generic MOEA */
typedef void (*EMO_MOEA_falloc)(void *, EMO_Param *, EMO_Population *, EMO_MOP *);
typedef void (*EMO_MOEA_ffree)(void *);
typedef void (*EMO_MOEA_frun)(void *, EMO_Param *, EMO_Population *, EMO_MOP *);

typedef struct {
  EMO_MOEA_falloc alloc;
  EMO_MOEA_ffree free;
  EMO_MOEA_frun run;

  union {
    EMO_SMSEMOA smsemoa;
    EMO_NSGA2 nsga2;
    EMO_NSGA3 nsga3;
    EMO_MOMBI2 mombi2;
    EMO_MOMBI3 mombi3;
    EMO_IBEA ibea;
    EMO_MOEAD moead;
    EMO_SPEA2 spea2;
    EMO_HYPE hype;
    EMO_MOVAP movap;
  } alg;

  void *palg;
} EMO_MOEA;


void EMO_MOEA_alloc(EMO_MOEA *moea, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, const char *str);
void EMO_MOEA_free(EMO_MOEA *moea);
void EMO_MOEA_run(EMO_MOEA *moea, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

#endif

