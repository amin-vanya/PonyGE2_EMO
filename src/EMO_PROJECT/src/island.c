
#include <stdlib.h>
#include <string.h>

#include "common.h"
#include "island.h"

#define _MAXCHAR 2000

void EMO_Island_alloc(EMO_Island *comm, EMO_Debug *dbg, EMO_MOP *mop, EMO_List *src, EMO_List *dest, int nemig, int nimmig) {
  int n;

  comm->src = src;
  comm->dest = dest;
  comm->nemig = nemig;
  comm->nimmig = nimmig; 

  // msg de interrupcion de origen + migrantes (variables, objetivos y restricciones)
  n = (mop->nvar + mop->nobj + mop->ncon) * sizeof(double);
  comm->nbuffer_recv = sizeof(int) + comm->nimmig * n;
  comm->nbuffer_send = sizeof(int) + comm->nemig  * n;

  comm->enable_dest = (int *) malloc(sizeof(int) * comm->dest->size);
  comm->tag = 0;

  comm->buffer_send = NULL;
  comm->buffer_recv = NULL;
  comm->request_recv = NULL;
  comm->status_recv = NULL;
  comm->request_send = NULL;
  comm->status_send = NULL;

  EMO_Island_reset(comm);
}

void EMO_Island_free(EMO_Island *comm) {
  int i;

  MPI_Barrier(MPI_COMM_WORLD);

  for(i = 0; i < comm->src->size; i++) {
    MPI_Request_free(&comm->request_recv[i]);
  }

  for(i = 0; i < comm->dest->size; i++) {
    MPI_Request_free(&comm->request_send[i]);
  }

  free(comm->buffer_recv);
  free(comm->buffer_send);
  free(comm->request_recv);
  free(comm->request_send);
  free(comm->status_recv);
  free(comm->status_send);
 
  free(comm->enable_dest);
}

void EMO_Island_reset(EMO_Island *comm) {
  int i, r, elem;
  static int res = 0;

  comm->incoming = comm->outcoming = 0;

  printf("Reset %d\n", res++);
  MPI_Barrier(MPI_COMM_WORLD);

  comm->tag++;

  for(i = 0; i < comm->dest->size; i++)
    comm->enable_dest[i] = 1;

  if(comm->request_recv != NULL) {
    for(i = 0; i < comm->src->size; i++) {
      MPI_Request_free(&comm->request_recv[i]);
      comm->request_recv[i] = MPI_REQUEST_NULL;
    }
    free(comm->request_recv);
  }

  if(comm->request_send != NULL) {
    for(i = 0; i < comm->dest->size; i++) {
      MPI_Request_free(&comm->request_send[i]);
      comm->request_send[i] = MPI_REQUEST_NULL;
    }
    free(comm->request_send);
  }

  if(comm->buffer_send != NULL) 
    free(comm->buffer_send);
  if(comm->buffer_recv != NULL) 
    free(comm->buffer_recv);
  if(comm->status_recv != NULL) 
    free(comm->status_recv);
  if(comm->status_send != NULL) 
    free(comm->status_send);

  if((comm->buffer_send = (char *) malloc(sizeof(char) * comm->nbuffer_send)) == NULL) {
    printf("Error, not enough memory in EMO_Island_reset 1.\n");
    exit(1);
  }

  if((comm->buffer_recv = (char *) malloc(sizeof(char) * comm->nbuffer_recv)) == NULL) {
    printf("Error, not enough memory in EMO_Island_reset 2.\n");
    exit(1);
  }

  /* Recibe */
  if((comm->request_recv = (MPI_Request *) malloc(sizeof(MPI_Request) * comm->src->size)) == NULL) {
    printf("Error, not enough memory in EMO_Island_reset 3.\n");
    exit(1);
  }

  if((comm->status_recv = (MPI_Status *) malloc(sizeof(MPI_Status) * comm->src->size)) == NULL) {
    printf("Error, not enough memory in EMO_Island_reset 4.\n");
    exit(1);
  }

  for(i = 0; i < comm->src->size; i++) {
    comm->request_recv[i] = MPI_REQUEST_NULL;
    EMO_List_get(comm->src, &elem, i);
    r = MPI_Recv_init(comm->buffer_recv, comm->nbuffer_recv, MPI_PACKED, elem, comm->tag, MPI_COMM_WORLD, &comm->request_recv[i]);

    if(r == MPI_SUCCESS)
      MPI_Start(&comm->request_recv[i]); /* Lanza solicitud */
  }

  /* Envia */
  if((comm->request_send = (MPI_Request *) malloc(sizeof(MPI_Request) * comm->dest->size)) == NULL) {
    printf("Error, not enough memory in EMO_Island_reset 5.\n");
    exit(1);
  }

  if((comm->status_send = (MPI_Status *) malloc(sizeof(MPI_Status) * comm->dest->size)) == NULL) {
    printf("Error, not enough memory in EMO_Island_reset 6.\n");
    exit(1);
  } 

  for(i = 0; i < comm->dest->size; i++) {
    comm->request_send[i] = MPI_REQUEST_NULL;
    EMO_List_get(comm->dest, &elem, i);
    MPI_Send_init(comm->buffer_send, comm->nbuffer_send, MPI_PACKED, elem, comm->tag, MPI_COMM_WORLD, &comm->request_send[i]);
  }

  printf("after Reset %d\n", res);
}

void EMO_Island_send(EMO_Island *comm, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, EMO_List *emigrant, int stop) {
  int i, r, idx, position, flag = 0, count;

  /**************************** Envia datos a nodos destino activados ****************************/
  position = 0;

  // envia -1 para interrupcion o numero de individuos a migrar
  count = (stop == 1)? -1 : emigrant->size;
  MPI_Pack(&count, 1, MPI_INT, comm->buffer_send, comm->nbuffer_send, &position, MPI_COMM_WORLD);

  EMO_Debug_printf(param->dbg, "send:interrupt %d, count %d", stop, count);

  while(emigrant->size > 0) {
    EMO_List_dequeue(emigrant, &idx);
    MPI_Pack(pop->var + idx * mop->nvar, mop->nvar, MPI_DOUBLE, comm->buffer_send, comm->nbuffer_send, &position, MPI_COMM_WORLD);
    MPI_Pack(pop->obj + idx * mop->nobj, mop->nobj, MPI_DOUBLE, comm->buffer_send, comm->nbuffer_send, &position, MPI_COMM_WORLD);

    if(mop->ncon != 0)
      MPI_Pack(pop->con + idx * mop->ncon, mop->ncon, MPI_DOUBLE, comm->buffer_send, comm->nbuffer_send, &position, MPI_COMM_WORLD);

    comm->outcoming++;
    EMO_Debug_printv(param->dbg, pop->obj + idx * mop->nobj, mop->nobj, "send:individual obj %d", idx);
  }

  for(i = 0; i < comm->dest->size; i++) {
    if(!comm->enable_dest[i]) continue;

    r = MPI_Test(&comm->request_send[i], &flag, &comm->status_send[i]);
    EMO_Debug_printf(param->dbg, "send:result %d, flag %d", r, flag);

    if(r == MPI_SUCCESS && flag != 0) {
      r = MPI_Start(&comm->request_send[i]);
      EMO_Debug_printf(param->dbg, "send:start %d", r);
    }
  }
}

int EMO_Island_receive(EMO_Island *comm, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, int *stop) {
  int i, j, k, r, num, position, flag, count;  //interrupt, 

  *stop = flag = 0;
  k = pop->mu;

  /**************************** Recibe datos de origen ****************************/
  if(comm->request_recv != NULL) {

    EMO_Debug_printf(param->dbg, "receive:0");

    for(i = 0; i < comm->src->size; i++) {
    
      EMO_Debug_printf(param->dbg, "receive:1");
      while(k + comm->nimmig < pop->size) {
        /* Check if there are available messages of interruption and individuals from the internal buffer */
        //EMO_Debug_printf(param->dbg, "before mpi_test data recv %d, %p", i, &comm->request_recv[i]);
        r = MPI_Test(&comm->request_recv[i], &flag, &comm->status_recv[i]);
        EMO_Debug_printf(param->dbg, "receive:result %d, flag %d", r, flag);
 
        if(r == MPI_SUCCESS && flag != 0) {
          r = MPI_Get_count(&comm->status_recv[i], MPI_PACKED, &num);
          EMO_Debug_printf(param->dbg, "receive:MPI_Get_count %d from %d, result %d", num, comm->status_recv[i].MPI_SOURCE, r);

          if(num == 0) {
            EMO_Debug_printf(param->dbg, "receive:no more info from %d", comm->status_recv[i].MPI_SOURCE);
            break;
          }

          position = 0;
          r = MPI_Unpack(comm->buffer_recv, comm->nbuffer_recv, &position, &count, 1, MPI_INT, MPI_COMM_WORLD);
          EMO_Debug_printf(param->dbg, "receive:interrupt/count %d", count);

          if(r == MPI_SUCCESS && count == -1) {
            EMO_Debug_printf(param->dbg, "receive:program interruption from %d", comm->status_recv[i].MPI_SOURCE);
            *stop = *stop || 1;
            break;
          }

          for(j = 0; j < count; j++) {  // RHG, recibe 0's cuando se reciben menos de nmmig individuos
            MPI_Unpack(comm->buffer_recv, comm->nbuffer_recv, &position, pop->var + k * mop->nvar, mop->nvar, MPI_DOUBLE, MPI_COMM_WORLD);
            MPI_Unpack(comm->buffer_recv, comm->nbuffer_recv, &position, pop->obj + k * mop->nobj, mop->nobj, MPI_DOUBLE, MPI_COMM_WORLD);

            if(mop->ncon != 0)
              MPI_Unpack(comm->buffer_recv, comm->nbuffer_recv, &position, pop->con + k * mop->ncon, mop->ncon, MPI_DOUBLE, MPI_COMM_WORLD);

            comm->incoming++;
            EMO_Debug_printv(param->dbg, pop->obj + k * mop->nobj, mop->nobj, "receive:individual obj %d", j);
            k++;
          }

          if(count >= 0)
            MPI_Start(&comm->request_recv[i]);
        }
        else {
          break;
        }
      } /* end-while */
    }
  }
  else
    EMO_Debug_printf(param->dbg, "receive:null");

  return k - pop->mu;
}

