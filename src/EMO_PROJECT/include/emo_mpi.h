
#ifndef _EMO_MPI_H
#define _EMO_MPI_H

#include "mpi.h"
#include "emo.h"

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

typedef struct {
  EMO_Island comm;      /* MPI communicator */
  EMO_Migration policy;
  EMO_MOEA moea;
  EMO_List src;  /* Source and destination nodes,   */
  EMO_List dest; /* defined by the logical topology */
  EMO_List lst1;   /* Temporary lists */
  EMO_List lst2;
  int interrupt; /* Enable interruption of the program from external nodes */
  int mpolicy;   /* Migration policy */
  int rpolicy;   /* Replacement policy */
  int sync;      /* Enable synchronous migration */
  int nmig;    /* Number of immigrant individuals */
  int epoch;     /* Frequency of migration events   */
  int nproc;
  int myrank;
} EMO_PMOEA;

void EMO_PMOEA_alloc(EMO_PMOEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, const char *str);
void EMO_PMOEA_free(EMO_PMOEA *alg);
void EMO_PMOEA_run(EMO_PMOEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_PMOEA_write(EMO_PMOEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, int run);

typedef struct {
  EMO_Population archive;
  EMO_Island comm;      /* MPI communicator */
  EMO_Migration mpolicy;
  EMO_Migration rpolicy;
  EMO_SMSEMOA moea;
  EMO_VPath vpath;
  EMO_Prune p;
  EMO_List src;  /* Source and destination nodes,   */
  EMO_List dest; /* defined by the logical topology */
  EMO_List lst1;   /* Temporary lists */
  EMO_List lst2;
  int nmpolicy;   /* Migration policy */
  int nrpolicy;   /* Replacement policy */
  int nmig;    /* Number of immigrant individuals */
  int epoch;     /* Frequency of migration events   */
  int nproc;
  int myrank;
  int xrs;
} EMO_PAMICRO;


void EMO_PAMICRO_alloc(EMO_PAMICRO *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_PAMICRO_free(EMO_PAMICRO *alg);
void EMO_PAMICRO_run(EMO_PAMICRO *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_PAMICRO_write(EMO_PAMICRO *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, int run);

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

