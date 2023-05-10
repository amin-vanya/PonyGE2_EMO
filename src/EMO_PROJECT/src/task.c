
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "mpi.h"
#include "task.h"

#define MAX_CHAR 2000
#define MASTER 0
#define NO_SENT 40
#define SENT 30

int more_task(EMO_Task *task) {
  static long pos = 0L;
  FILE *fp;
  int l;

  /* Checks if the last command was not sent */
  if(task->status != SENT)
    return 1;

  if((fp = fopen(task->finput, "r")) == NULL) {
    EMO_Debug_printf(&task->dbg, "Error to open file %s", task->finput);
    return 0;
  }

  if(fseek(fp, pos, SEEK_SET) == -1)
    return 0;  /* End of file */
 
  if(fgets(task->sbuf, MAX_CHAR, fp) == NULL) {
    EMO_Debug_printf(&task->dbg, "End of file %s", task->finput);
    return 0;
  }
  l = strlen(task->sbuf);
  task->sbuf[l-1] = '\0';

  pos = ftell(fp);
  fclose(fp);
  task->status = NO_SENT;
  return 1;
}

void report_task(EMO_Task *task) {
  FILE *fp;

  if((fp = fopen(task->foutput, "a")) == NULL) {
    EMO_Debug_printf(&task->dbg, "Error to open file %s", task->foutput);
    return;
  }

  if(fprintf(fp, "%s\n", task->rbuf) == 0) {
    EMO_Debug_printf(&task->dbg, "Error while writing file %s", task->foutput);
    return;
  }

  fclose(fp);
}

void master(EMO_Task *task) {
  int result, tag, dest, source, count, flag, i;
  MPI_Request rreq, sreq;
  MPI_Status status;

  source = MPI_ANY_SOURCE;
  count = MAX_CHAR;
  tag = flag = 0;

  MPI_Irecv(task->rbuf, count, MPI_CHAR, source, tag, MPI_COMM_WORLD, &rreq);

  while(more_task(task)) {

    if(task->list.size > 0) {
      EMO_List_dequeue(&task->list, &dest);
      result = MPI_Isend(task->sbuf, count, MPI_CHAR, dest, tag, MPI_COMM_WORLD, &sreq);
      result += MPI_Wait(&sreq, MPI_STATUS_IGNORE);
      EMO_Debug_printf(&task->dbg, "Master|send to slave %d|%s|result %d", dest, task->sbuf, result);
      task->status = SENT;
    }
    else {
      EMO_Debug_printf(&task->dbg, "Master|there are no available slaves, waiting for one");
      result = MPI_Wait(&rreq, &status);
      flag = 1;
    }

    if(flag == 0) 
      result = MPI_Test(&rreq, &flag, &status);

    if(result == MPI_SUCCESS && flag != 0) {
      EMO_Debug_printf(&task->dbg, "Master|recv from %d|%s", status.MPI_SOURCE, task->rbuf);
      report_task(task);
      EMO_List_queue(&task->list, status.MPI_SOURCE);
      MPI_Irecv(task->rbuf, count, MPI_CHAR, source, tag, MPI_COMM_WORLD, &rreq);
      flag = 0;
    }
  }
 
  sprintf(task->sbuf, "END");

  for(i = 1; i < task->numprocs; i++) {
    result = MPI_Isend(task->sbuf, count, MPI_CHAR, i, tag, MPI_COMM_WORLD, &sreq);
    result += MPI_Wait(&sreq, MPI_STATUS_IGNORE);
    EMO_Debug_printf(&task->dbg, "Master|send to %d|%s|result %d", i, task->sbuf, result);
  }

  /* Waits for the slave responses */
  while(task->list.size != task->numprocs - 1) {
    result = MPI_Wait(&rreq, &status);

    if(result == MPI_SUCCESS) {
      EMO_Debug_printf(&task->dbg, "Master|recv from %d|%s", status.MPI_SOURCE, task->rbuf);
      report_task(task);
      EMO_List_queue(&task->list, status.MPI_SOURCE);

      if(task->list.size != task->numprocs - 1)
        MPI_Irecv(task->rbuf, count, MPI_CHAR, source, tag, MPI_COMM_WORLD, &rreq);
    }
  }
}

void slave(EMO_Task *task) {
  int result, tag, dest, source, count;
  MPI_Request rreq, sreq;

  dest = source = MASTER;
  count = MAX_CHAR;
  tag = 0;

  while(1) {
    result = MPI_Irecv(task->rbuf, count, MPI_CHAR, source, tag, MPI_COMM_WORLD, &rreq);
    result += MPI_Wait(&rreq, MPI_STATUS_IGNORE);
    EMO_Debug_printf(&task->dbg, "Slave|recv from master|%s|result %d", task->rbuf, result);

    if(strncmp(task->rbuf, "END", 3) == 0) {
      EMO_Debug_printf(&task->dbg, "Slave|halt");
      break;
    }

    EMO_Debug_printf(&task->dbg, "Slave|executing %s", task->rbuf);
    result = system(task->rbuf);
    sprintf(task->sbuf, "DONE:result %d:%s", result, task->rbuf);

    result = MPI_Isend(task->sbuf, count, MPI_CHAR, dest, tag, MPI_COMM_WORLD, &sreq);
    result += MPI_Wait(&sreq, MPI_STATUS_IGNORE);
    EMO_Debug_printf(&task->dbg, "Slave|send to master|%s|result %d", task->sbuf, result);
  }
}

void EMO_Task_alloc(EMO_Task *task, int *argc, char ***argv) {
  int version, subversion, i;
  char str[100];

  MPI_Init(argc, argv);
  MPI_Comm_size(MPI_COMM_WORLD, &task->numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &task->myrank);

  sprintf(str, "Proc_%d.log", task->myrank);
  EMO_Debug_alloc(&task->dbg, EMO_DEBUG_ON, task->myrank, str);

  if(MPI_Get_version(&version, &subversion) != MPI_SUCCESS) {
    EMO_Debug_printf(&task->dbg, "Error MPI 3");
    exit(1);
  }

  EMO_Debug_printf(&task->dbg, "Rank: %d, MPI version: %d, subversion: %d", task->myrank, version, subversion);

  if((task->sbuf = (char *) malloc(sizeof(char) * MAX_CHAR)) == NULL) {
    printf("Error when allocating memory in task.c\n");
    exit(1);
  }

  if((task->rbuf = (char *) malloc(sizeof(char) * MAX_CHAR)) == NULL) {
    printf("Error when allocating memory in task.c\n");
    exit(1);
  }

  if(task->myrank == MASTER) {
    if((task->finput = (char *) malloc(sizeof(char) * MAX_CHAR)) == NULL) {
      printf("Error when allocating memory in task.c\n");
      exit(1);
    }

    if((task->foutput = (char *) malloc(sizeof(char) * MAX_CHAR)) == NULL) {
      printf("Error when allocating memory in task.c\n");
      exit(1);
    }

    strcpy(task->finput, (*argv)[1]);
    sprintf(task->foutput, "%s.done", task->finput);

    EMO_List_alloc(&task->list, task->numprocs - 1);

    for(i = 1; i < task->numprocs; i++)
      EMO_List_queue(&task->list, i);

    task->status = SENT;
  }
}

void EMO_Task_free(EMO_Task *task) {

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();

  free(task->sbuf);
  free(task->rbuf);

  if(task->myrank == MASTER) {
    EMO_List_free(&task->list);
    free(task->finput);
    free(task->foutput);
  }

  EMO_Debug_free(&task->dbg);
}

void EMO_Task_run(EMO_Task *task) {
  if(task->myrank == MASTER)
    master(task);
  else
    slave(task);
}

#undef MAX_CHAR
#undef MASTER
#undef NO_SENT
#undef SENT

