#include <stdlib.h>
#include <stdio.h>

#include "emo.h"

int main(int argc, char **argv) {
  int row, col, i, j, elem;
  EMO_NDSort nds;
  double *data;

  if(argc != 2) {
    printf("Syntax: %s file\n", argv[0]);
    exit(1);
  }

  elem = row = col = 0;

  data = EMO_File_read(NULL, &row, &col, argv[1], 0);

  EMO_NDSort_alloc(&nds, row);
  EMO_NDSort_run(&nds, data, col, NULL, NULL, row);

  printf("#front objectives\n");
  for(i = 0; i < nds.nfront; i++) {
    for(j = 0; j < nds.front[i].size; j++) {
      EMO_List_get(&(nds.front[i]), &elem, j);
      EMO_vprint(stdout, data + elem*col, col, "%d ", i);
    }
  }

  EMO_NDSort_free(&nds);
  free(data);
  return 0;
}

