#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "emo.h"

/**** Definition of a new Multi-objective Optimization Problem (MOP) ***

    minimize   {f_0(x), f_1(x), ..., f_{mop->nobj-1}(x)}

    subject to g_0(x) >= 0
               g_1(x) >= 0
               ...
               g_{mop->ncon - 1}(x) >= 0

               x \in [x_min, x_max]

 f: vector of objective functions
 g: vector of inequality constraints
 x: vector of decision variables

 Note: the following function prototype should be defined for a 
       MOP with inequality constraints:

 void myMOP_eval(EMO_MOP *mop, double *f, double *g, double *x)

 A non feasible solution is that one which for any i g_i(x) < 0.

*/
void myMOP_eval(EMO_MOP *mop, double *f, double *x) {
  f[0] = x[0] * x[0];
  f[1] = pow(x[0] - 1.0, 2.0);
}

void myMOP_alloc(EMO_MOP *mop) {
  int i;

  /* MOP's name */
  mop->name = (char *) malloc(sizeof(char) * 200);
  strcpy(mop->name, "WING_DESIGN");

  /* Number of decision variables */
  mop->nvar = 2;

  /* Number of objectives */
  mop->nobj = 2;

  /* Box constraints */
  mop->xmin = (double *) calloc(sizeof(double), mop->nvar);
  mop->xmax = (double *) calloc(sizeof(double), mop->nvar);

  for(i = 0; i < mop->nvar; i++) {
    mop->xmin[i] = 0.0;
    mop->xmax[i] = 1.0;
  }

  /* Required */
  mop->npos = 0;
  mop->feval = 0;
  mop->coding = EMO_REAL;

  /* MOP without inequality constraints */
  mop->ncon = 0;
  mop->f = myMOP_eval;

  /* MOP with inequality constraints */
  //mop->ncon = #;
  //mop->fc = myMOP_eval;
}

/*********************** End of the definition **********************/


int main(int argc, char **argv) {
  EMO_Population pop;
  EMO_Param param;
  EMO_MOEA moea;
  EMO_MOP mop;
  int r, nrun;

  if(argc != 5) {
    printf("\nSyntax: %s MOEA parameter_file {MOP, default} runs\n\n", argv[0]);
    EMO_Dictionary_print(stdout, EMO_MOEA_list, "MOEA");
    EMO_Dictionary_print(stdout, EMO_Benchmark_list, "\nMOP");
    EMO_Dictionary_print(stdout, EMO_Benchmark_listc, "\nConstrained MOP");
    printf("\nWhen the default option is selected, the MOEA solves the problem defined in %s.c:myMOP_eval\n\n", argv[0]);
    return 1;
  }

  nrun = atoi(argv[4]);

  if(strcmp(argv[3], "default") == 0)
    myMOP_alloc(&mop);
  
  EMO_Param_alloc_from_file(&param, &mop, argv[1], argv[2], argv[3]);
  EMO_MOEA_alloc(&moea, &param, &pop, &mop, argv[1]);

  //EMO_Population_init_from_file(&pop, &mop, "test", 0);

  for(r = 1; r <= nrun; r++) {
    EMO_Stop_start(param.stop);

    EMO_Population_init(&pop, param.rand, &mop);

    EMO_MOEA_run(&moea, &param, &pop, &mop);

    EMO_Population_write(&pop, NULL, &mop, param.prefix, r, 0);
    EMO_Param_save(&param, &pop, &mop, r);

    EMO_Rand_next_seed(param.rand, 0);
  }

  EMO_MOEA_free(&moea);
  EMO_Param_free(&param);  /* free MOP */
  EMO_Population_free(&pop);

  return 0;
}

