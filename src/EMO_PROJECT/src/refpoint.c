
#include <stdlib.h>
#include <float.h>

#include "refpoint.h"
#include "dominance.h"

void EMO_Refpoint_alloc(EMO_Refpoint *ref, EMO_Param *param, int nobj) {
  int i, j;

  EMO_Utility_alloc(&ref->utl, param, nobj, "all");

  ref->ideal = (double *) malloc(sizeof(double) * nobj);
  ref->nadir = (double *) malloc(sizeof(double) * nobj);

  ref->xtrm  = (int *) malloc(sizeof(int) * nobj);
  ref->axis0 = (double *) malloc(sizeof(double) * nobj * nobj);
  ref->axis1 = (double *) malloc(sizeof(double) * nobj * nobj);

  for(i = 0; i < nobj; i++) {
    ref->ideal[i] = DBL_MAX;
    ref->nadir[i] = -DBL_MAX;

    for(j = 0; j < nobj; j++) {
      ref->axis0[i * nobj + j] = (i == j)? 0.0 : 1.0;
      ref->axis1[i * nobj + j] = (i == j)? 1.0 : 0.0;
    }
  }

  ref->nobj = nobj;
}

void EMO_Refpoint_free(EMO_Refpoint *ref) {
  EMO_Utility_free(&ref->utl);
  free(ref->ideal);
  free(ref->nadir);
  free(ref->xtrm);
  free(ref->axis0);
  free(ref->axis1);
}

void EMO_Refpoint_update_ideal(EMO_Refpoint *ref, double *data, int *filter, int size) {
  int i, j, k;

  if(filter == NULL) {
    for(i = 0; i < size; i++) {
      for(j = 0; j < ref->nobj; j++) {
       k = i * ref->nobj + j;

       if(data[k] < ref->ideal[j])
         ref->ideal[j] = data[k];
      }
    }
  }
  else {
    for(i = 0; i < size; i++) {
      if(filter[i] == 0)
        continue;
 
      for(j = 0; j < ref->nobj; j++) {
        k = i * ref->nobj + j;

        if(data[k] < ref->ideal[j])
          ref->ideal[j] = data[k];
      }
    }
  }
}

// shift objectives to the origin (zdt3)
/* It updates the nadir point using the scalarizing functions wpo2 and aasf */
void EMO_Refpoint_update_nadir(EMO_Refpoint *ref, double *data, int *filter, int size) {
  int r, i, j, jmin = -1, h = 0;
  double vmin, v;

  static const EMO_UtilityFunction wdic[] = { EMO_Utility_wpo2,
                                              EMO_Utility_refsf,
                                              NULL
                                            };

  if(filter == NULL) {
    printf("Error, filter cannot be NULL in refpoint.c:EMO_Refpoint_update_nadir\n");
    exit(1);
  } 

  for(i = 0; i < ref->nobj; i++) {
    ref->nadir[i] = -DBL_MAX;
    ref->xtrm[i] = -1;
  }

  while(wdic[h] != NULL) {

    for(i = 0; i < ref->nobj; i++) {
      vmin = DBL_MAX;

      for(j = 0; j < size; j++) {
        if(filter[j] == 1) {

          v = wdic[h](&ref->utl, ref->axis1 + i * ref->nobj, data + j * ref->nobj);

          if(v < vmin) {
            jmin = j;
            vmin = v;
          }
          else if(v == vmin) {
            r = EMO_Dominance_strict(data + j * ref->nobj, data + jmin * ref->nobj, ref->nobj);

            if(r == 1) {
              jmin = j;
              vmin = v;
            }
          }
        }
      }

      if(filter[jmin] == 1 && data[jmin * ref->nobj + i] > ref->nadir[i]) {
        ref->nadir[i] = data[jmin * ref->nobj + i];

        if(ref->xtrm[i] > -1)  // restores values
          filter[ref->xtrm[i]] = 1;

        ref->xtrm[i] = jmin;
        filter[jmin] = 0;  // Avoids to select the individual again
      }
    }


    if(ref->nobj > 2) {
    // inverted POF  (e.g. MIRROR_GAMMA_2.0)
    for(i = 0; i < ref->nobj; i++) {
      vmin = DBL_MAX;

      for(j = 0; j < size; j++) {
        if(filter[j] == 1) {

          v = wdic[h](&ref->utl, ref->axis0 + i * ref->nobj, data + j * ref->nobj);

          if(v < vmin) {
            jmin = j;
            vmin = v;
          }
          else if(v == vmin) {
            r = EMO_Dominance_strict(data + j * ref->nobj, data + jmin * ref->nobj, ref->nobj);

            if(r == 1) {
              jmin = j;
              vmin = v;
            }
          }
        }
      }

      for(j = 0; j < ref->nobj; j++) {
        if(i != j) {

          if(filter[jmin] == 1 && data[jmin * ref->nobj + j] > ref->nadir[j]) {

            if(ref->xtrm[j] > -1)  // restores values
              filter[ref->xtrm[j]] = 1;

            ref->xtrm[j] = jmin;
            filter[jmin] = 0;  // Avoids to select the individual again
            break;
          }
 
        }
      }
    }
    }

    h++;
  }

  /* Checks that the nadir point encloses the extreme points */
  for(i = 0; i < ref->nobj; i++) {
    for(j = 0; j < ref->nobj; j++) {
      h = ref->xtrm[i] * ref->nobj + j;

      if(data[h] > ref->nadir[j])
        ref->nadir[j] = data[h];
    }

    if(ref->xtrm[i] > -1)  // restores values
      filter[ref->xtrm[i]] = 1;
  }
}

