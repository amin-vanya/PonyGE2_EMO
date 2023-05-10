
#include <stdio.h>

#include "emo_mpi.h"

int main(int argc, char **argv) {
  EMO_Task task;

  if(argc != 2) {
    printf("Syntax: mpirun -np num_procs %s command_file\n", argv[0]);
    printf("Syntax: mpirun --hostfile file %s command_file\n", argv[0]);
    return 0;
  }

  EMO_Task_alloc(&task, &argc, &argv);
  EMO_Task_run(&task);
  EMO_Task_free(&task);

  return 0;
}

