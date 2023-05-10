#ifndef _ISLAND_H
#define _ISLAND_H

#include "common.h" 
#include "param.h"
#include "list.h"
#include "mpi.h"

typedef struct {
  /* Se crea uno por cada nodo de conexion */
  int nbuffer_recv;           /* Tamaño del buffer para recepcion */
  char *buffer_recv;          /* Buffer de recepcion de mensajes */
  MPI_Request *request_recv;  /* Solicitud para atender mensajes entrantes */
  MPI_Status *status_recv;    /* Estatus de mensajes recibidos */

  int nbuffer_send;           /* Tamaño del buffer para envio */
  char *buffer_send;          /* Buffer para el envio de mensajes */
  MPI_Request *request_send;  /* Solicitud para enviar mensajes */
  MPI_Status *status_send;    /* Estatus de mensajes enviados */

  int *enable_dest;
  int tag;

  /* Counters of incoming and outgoing individuals */
  int incoming;
  int outcoming;

  EMO_List *src;  /* Source and destination nodes,   */
  EMO_List *dest; /* defined by the logical topology */
  int nimmig;     /* Number of immigrant individuals */
  int nemig;      /* Number of emigrant individuals  */
} EMO_Island;


void EMO_Island_alloc(EMO_Island *comm, EMO_Debug *dbg, EMO_MOP *mop, EMO_List *src, EMO_List *dest, int nemig, int nimmig);
void EMO_Island_free(EMO_Island *comm);
void EMO_Island_reset(EMO_Island *comm);
void EMO_Island_send(EMO_Island *comm, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, EMO_List *emigrant, int stop);
int EMO_Island_receive(EMO_Island *comm, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, int *stop);

#endif

