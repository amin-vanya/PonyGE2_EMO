
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>

#include "plot.h"
#include "io.h"

#define LINECOLOR "3"
#define POINTSIZE "2"

void EMO_Plot_alloc(EMO_Plot *p, int dim, const char *term, const char *title, const char *pftrue, long freq) {
  int i, t;

  p->freq = freq;
  p->dim = 0;

  if(freq >= 0)  {
    signal(SIGPIPE, SIG_IGN);  // Ayuda a que no se salga programa
                               // cuando gnuplot no esta instalado
    p->dim = dim;

    p->term = (char *) malloc(sizeof(char) * 20);
    p->title = (char *) malloc(sizeof(char) * 200);
    p->pftrue = (char *) malloc(sizeof(char) * 200);

    strcpy(p->term, term);
    strcpy(p->title, title);
    strcpy(p->pftrue, pftrue);

    t = 2;  /* Number of plots */
    p->gp = (FILE **) malloc(sizeof(FILE *) * t);

    for(i = 0; i < t; i++) {
      p->gp[i] = popen(GNUPLOT_COMMAND, "w");

      if(p->gp[i] == NULL) {
        printf("Error in EMO_Plot_alloc: %s\n", GNUPLOT_COMMAND);
        exit(1);
      }
    }
  }
  else if (freq < 0) {
    printf("Error, invalid frequency value in plot.c: %ld\n", freq);
    exit(1);
  }
} 

void EMO_Plot_free(EMO_Plot *p) {
  int i, t;

  if(p == NULL) return;

  free(p->term);
  free(p->title);
  free(p->pftrue);

  t = 2;

  for(i = 0; i < t; i++)
    pclose(p->gp[i]);

  free(p->gp);
}

/* Plots the current population */
void EMO_Plot_run(EMO_Plot *p, double *data, int size, long feval, int end) {

  if(p == NULL || (p->freq == 0 && p->dim == 0) || (p->freq == 0 && end == 0)) return;

  if(end == 0 && (feval % p->freq) != 0)  /* Skips */
    return;

  EMO_File_write(data, NULL, size, p->dim, "plot.out", "%f ", 0);
 
  /* Parallel coordinates */
  if(1) {
    fprintf(p->gp[0], "set term %s noenhanced\nset title '%s %dD - Evaluation #%ld'\nunset key\nset xlabel 'Objective'\nset ylabel 'Value'\n", p->term, p->title, p->dim, feval);
    fprintf(p->gp[0], "plot [1:%d] \"<awk '$1 !~ /#/ { for(i = 1; i <= NF; i++) print i,$i; print z}' ./plot.out\" u 1:2:xtic(1) w linesp linecolor %s\n", p->dim, LINECOLOR);
    fflush(p->gp[0]);
  }

  if(p->dim == 2) {
    fprintf(p->gp[1], "set term %s noenhanced\nset title '%s %dD - Evaluation #%ld'\nunset key\nset xlabel 'f1'\nset ylabel 'f2'\nplot ", p->term, p->title, p->dim, feval);

      if(strlen(p->pftrue) != 0) {/* Plots Pareto front true */
        if(p->pftrue[0] == '\'') /* from file */
          fprintf(p->gp[1], "%s w points pointtype 2 linecolor 5, ", p->pftrue);
        else                  /* from expresion */
          fprintf(p->gp[1], "%s with lines linestyle 1 linecolor 3, ", p->pftrue);
      }

      fprintf(p->gp[1], "'plot.out' w points pointtype 6 pointsize %s linecolor %s\n", POINTSIZE, LINECOLOR); // linecolor 2
      fflush(p->gp[1]);
  }
  else if(p->dim == 3) {
      fprintf(p->gp[1], "set term %s noenhanced\nset title '%s %dD - Evaluation #%ld'\nunset key\nset xlabel 'f1'\nset ylabel 'f2'\nset zlabel 'f3'\nset view 30,30\nsplot ", p->term, p->title, p->dim, feval);

      if(strlen(p->pftrue) != 0) { /* Plots Pareto front true */
        if(p->pftrue[0] == '\'') { /* from file */
          fprintf(p->gp[1], "%s w points pointtype 2 linecolor 5, ", p->pftrue);
        }
        else {                /* from parametric expresion */
          fprintf(p->gp[1], "set parametric\nsplot %s w l linecolor 5, ", p->pftrue);
        }
      } 
      fprintf(p->gp[1], "'plot.out' w points pointtype 6 pointsize %s linecolor %s\n", POINTSIZE, LINECOLOR);
      fflush(p->gp[1]);
  }
}

