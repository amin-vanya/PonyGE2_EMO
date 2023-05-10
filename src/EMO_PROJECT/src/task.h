#ifndef _TASK_H
#define _TASK_H

#include "debug.h"
#include "list.h"

typedef struct {
  char *finput;
  char *foutput;
  char *sbuf, *rbuf;
  int numprocs;
  int myrank;
  int status;
  EMO_Debug dbg;
  EMO_List list;
} EMO_Task;

void EMO_Task_alloc(EMO_Task *task, int *argc, char ***argv);
void EMO_Task_free(EMO_Task *task);
void EMO_Task_run(EMO_Task *task);

#endif

