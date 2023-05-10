/**************************************************************
 * refpoint.c    Determination of the reference points        *
 *                                                            *
 * Computer Science Department, CINVESTAV-IPN                 *
 *                                                            *
 * Professor:   Dr. Carlos A. Coello Coello                   *
 * Author:      Raquel Hernandez Gomez                        *
 *                                                            *
 * Jun 2017                                                   *
 *************************************************************/

#ifndef _REFPOINT_H
#define _REFPOINT_H

#include "utility.h"

typedef struct {
  EMO_Utility utl;
  double *ideal;
  double *nadir;
  int *xtrm;
  double *axis0;
  double *axis1;
  int nobj;
} EMO_Refpoint;

void EMO_Refpoint_alloc(EMO_Refpoint *ref, EMO_Param *param, int nobj);
void EMO_Refpoint_free(EMO_Refpoint *ref);
void EMO_Refpoint_update_ideal(EMO_Refpoint *ref, double *data, int *filter, int size);
void EMO_Refpoint_update_nadir(EMO_Refpoint *ref, double *data, int *filter, int size);

#endif

