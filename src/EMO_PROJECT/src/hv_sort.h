
#ifndef _HV_SORT_H
#define _HV_SORT_H

//26 07 16:
#include "hv_iwfg.h"

void IWFG_quicksort (void *const pbase, size_t total_elems, int start, EMO_IWFG *hv, size_t size, int(*cmp)(const void *, const void *, int, EMO_IWFG *));

#endif

