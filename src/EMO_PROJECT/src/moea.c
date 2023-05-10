
#include <string.h>
#include "moea.h"

const char *EMO_MOEA_list[] = { "SMS_EMOA",
                                "NSGA2",
                                "NSGA3",
                                "MOMBI2",
                                "MOMBI3",
                                "HV_IBEA",
                                "EPS_IBEA",
                                "R2_IBEA",
                                "MOEAD",
                                "SPEA2",
                                "HYPE",
                                "MOVAP",
                                 NULL
                              };


void EMO_MOEA_alloc(EMO_MOEA *moea, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, const char *str) {

  const EMO_MOEA_falloc a[] = { (EMO_MOEA_falloc) EMO_SMSEMOA_alloc, 
                                (EMO_MOEA_falloc) EMO_NSGA2_alloc,
                                (EMO_MOEA_falloc) EMO_NSGA3_alloc,
                                (EMO_MOEA_falloc) EMO_MOMBI2_alloc,
                                (EMO_MOEA_falloc) EMO_MOMBI3_alloc,
                                (EMO_MOEA_falloc) EMO_IBEA_alloc,
                                (EMO_MOEA_falloc) EMO_IBEA_alloc,
                                (EMO_MOEA_falloc) EMO_IBEA_alloc,
                                (EMO_MOEA_falloc) EMO_MOEAD_alloc,
                                (EMO_MOEA_falloc) EMO_SPEA2_alloc,
                                (EMO_MOEA_falloc) EMO_HYPE_alloc,
                                (EMO_MOEA_falloc) EMO_MOVAP_alloc
                              };

  const EMO_MOEA_ffree f[]  = { (EMO_MOEA_ffree) EMO_SMSEMOA_free, 
                                (EMO_MOEA_ffree) EMO_NSGA2_free,
                                (EMO_MOEA_ffree) EMO_NSGA3_free,
                                (EMO_MOEA_ffree) EMO_MOMBI2_free,
                                (EMO_MOEA_ffree) EMO_MOMBI3_free,
                                (EMO_MOEA_ffree) EMO_IBEA_free,
                                (EMO_MOEA_ffree) EMO_IBEA_free,
                                (EMO_MOEA_ffree) EMO_IBEA_free,
                                (EMO_MOEA_ffree) EMO_MOEAD_free,
                                (EMO_MOEA_ffree) EMO_SPEA2_free,
                                (EMO_MOEA_ffree) EMO_HYPE_free,
                                (EMO_MOEA_ffree) EMO_MOVAP_free
                              };

  const EMO_MOEA_frun r[]   = { (EMO_MOEA_frun) EMO_SMSEMOA_run, 
                                (EMO_MOEA_frun) EMO_NSGA2_run,
                                (EMO_MOEA_frun) EMO_NSGA3_run,
                                (EMO_MOEA_frun) EMO_MOMBI2_run,
                                (EMO_MOEA_frun) EMO_MOMBI3_run,
                                (EMO_MOEA_frun) EMO_IBEA_run,
                                (EMO_MOEA_frun) EMO_IBEA_run,
                                (EMO_MOEA_frun) EMO_IBEA_run,
                                (EMO_MOEA_frun) EMO_MOEAD_run,
                                (EMO_MOEA_frun) EMO_SPEA2_run,
                                (EMO_MOEA_frun) EMO_HYPE_run,
                                (EMO_MOEA_frun) EMO_MOVAP_run
                              };


  char *aux;
  int i;

  aux = EMO_toupper(str);
  i = EMO_Dictionary_find(EMO_MOEA_list, aux);

  if(i == -1) {
    printf("Error, unknown MOEA %s in EMO_MOEA_alloc.\n", aux);
    free(aux);
    exit(1);
  }

  moea->alloc = a[i];
  moea->free = f[i];
  moea->run = r[i];

  switch(i) {
    case 0:  moea->palg = &moea->alg.smsemoa; 
             break;
    case 1:  moea->palg = &moea->alg.nsga2; 
             break;
    case 2:  moea->palg = &moea->alg.nsga3; 
             break;
    case 3:  moea->palg = &moea->alg.mombi2; 
             break;
    case 4:  moea->palg = &moea->alg.mombi3; 
             break;
    case 5:
    case 6:
    case 7:  moea->palg = &moea->alg.ibea; 
             break;
    case 8:  moea->palg = &moea->alg.moead; 
             break;
    case 9:  moea->palg = &moea->alg.spea2; 
             break;
    case 10: moea->palg = &moea->alg.hype; 
             break;
    case 11: moea->palg = &moea->alg.movap; 
             break;
 
    default: printf("Error, unknown MOEA (2) %s.\n", aux);
             free(aux);
             exit(1);
  }

  moea->alloc(moea->palg, param, pop, mop);
  free(aux);
}


void EMO_MOEA_free(EMO_MOEA *moea) {
  moea->free(moea->palg);
}

void EMO_MOEA_run(EMO_MOEA *moea, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  moea->run(moea->palg, param, pop, mop);
}

