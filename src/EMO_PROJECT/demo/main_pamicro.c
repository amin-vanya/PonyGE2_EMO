#include <stdlib.h>
#include <string.h>

#include "emo_mpi.h"

int main(int argc, char **argv) {
  EMO_Population pop;
  EMO_Param param;
  EMO_PAMICRO moea;
  EMO_MOP mop;
  int r, nrun;

  if(argc != 4) {
    printf("Syntax: mpiexec -np N %s parameter_file test_problem runs\n", argv[0]);
    exit(1);
  }

  nrun = atoi(argv[3]);

  MPI_Init(&argc, &argv);
  EMO_Param_alloc_from_file(&param, &mop, "PAMICRO", argv[1], argv[2]);
  EMO_PAMICRO_alloc(&moea, &param, &pop, &mop);

  //EMO_Population_init_from_file(&pop, &mop, "test", 0);

  for(r = 1; r <= nrun; r++) {
    printf("run %d\n", r);
    EMO_Stop_start(param.stop);

    EMO_Population_init(&pop, param.rand, &mop);
    EMO_PAMICRO_run(&moea, &param, &pop, &mop);

    EMO_PAMICRO_write(&moea, &param, &pop, &mop, r);

    EMO_Param_save(&param, &pop, &mop, r);

    EMO_Rand_next_seed(param.rand, 0);

    EMO_Island_reset(&moea.comm);  // barrier
  }

  EMO_PAMICRO_free(&moea);
  EMO_Param_free(&param);
  EMO_Population_free(&pop);
  MPI_Finalize();

  return 0;
}

