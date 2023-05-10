
#include "dominance.h"
#include "random.h"
#include "list.h"

#ifndef _MIGRATION_H
#define _MIGRATION_H

#define EMO_MIGRATION_RANDOM            0
#define EMO_MIGRATION_ELITIST_RANDOM    1
#define EMO_MIGRATION_ELITIST_RANKING   2
#define EMO_MIGRATION_FRONT             3
#define EMO_MIGRATION_FRONT_RANDOM      4
#define EMO_MIGRATION_FRONT_RANKING     5

#define EMO_REPLACEMENT_RANDOM          0
#define EMO_REPLACEMENT_ELITIST_RANDOM  1
#define EMO_REPLACEMENT_ELITIST_RANKING 2
#define EMO_ELITIST                     3


typedef struct {
  EMO_NDSort nd;
  EMO_Rand *rnd;
  EMO_List best;
  EMO_List worst;
  int *filter;
  int *seq;  /* Sequence for EMO_Migration_random */
  int *tmp;  /* Arbitrary sequence */
  int max_size;
  int nobj;
} EMO_Migration;

void EMO_Migration_alloc(EMO_Migration *m, EMO_Rand *rnd, int max_size, int nobj);
void EMO_Migration_free(EMO_Migration *m);

void EMO_Migration_get_type(int *type, const char *str);
void EMO_Replacement_get_type(int *type, const char *str);

void EMO_Migration_random(EMO_Migration *m, EMO_List *l, int size, int nmig);
void EMO_Migration_elitist_random(EMO_Migration *m, EMO_List *l, double *data, int size, int nmig);
void EMO_Migration_elitist_ranking(EMO_Migration *m, EMO_List *l, double *data, int size, int nmig);
void EMO_Migration_front(EMO_Migration *m, EMO_List *l, double *data, int size);
void EMO_Migration_front_random(EMO_Migration *m, EMO_List *l, double *data, int size, int nmig);
void EMO_Migration_front_ranking(EMO_Migration *m, EMO_List *l, double *data, int size, int nmig);
void EMO_Replacement_random(EMO_Migration *m, EMO_List *migrate, EMO_List *replace, double *data, int size, int nmig);
void EMO_Replacement_elitist_random(EMO_Migration *m, EMO_List *migrate, EMO_List *replace, double *data, int size, int nmig);
void EMO_Replacement_elitist_ranking(EMO_Migration *m, EMO_List *migrate, EMO_List *replace, double *data, int size, int nmig);
void EMO_Replacement_elitist(EMO_Migration *m, EMO_List *migrate, EMO_List *replace, double *data, int size, int nmig);

#endif

