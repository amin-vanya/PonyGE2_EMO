
#ifndef _PLOT_H
#define _PLOT_H

#include <stdio.h>

#define GNUPLOT_COMMAND "gnuplot -persist 1>gnuplot.out 2>&1"
#define TERM "x11" /* wxt */

typedef struct {
  FILE **gp;
  int dim;
  long freq;
  char *term;
  char *title;
  char *pftrue;
} EMO_Plot;

void EMO_Plot_alloc(EMO_Plot *p, int dim, const char *term, const char *title, const char *pftrue, long freq);
void EMO_Plot_free(EMO_Plot *p);
void EMO_Plot_run(EMO_Plot *p, double *data, int size, long feval, int end);

#endif

