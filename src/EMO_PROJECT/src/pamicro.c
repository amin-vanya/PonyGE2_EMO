
#include <stdlib.h>

#include "topology.h"
#include "pamicro.h"
#include "numeric.h"
#include "island.h"

#define _MAXCHAR 2000

void EMO_PAMICRO_load_param(EMO_PAMICRO *alg, EMO_Param *param) {
  int v, eval;
  char *str;

  str = (char *) malloc(sizeof(char) * _MAXCHAR);

  if(!EMO_Param_get_int(param, &alg->epoch, "epoch")) {
    printf("Error, epoch is not defined in the configuration file.\n");
    exit(1);
  }

  EMO_Stop_set_epoch(param->stop, alg->epoch);

  if(!EMO_Param_get_int(param, &alg->nmig, "migrant")) {
    printf("Error, migrant is not defined in the configuration file.\n");
    exit(1);
  }

  // Define logical topology 
  EMO_Topology_ring(&alg->src, &alg->dest, alg->myrank, alg->nproc, EMO_TOPOLOGY_UNI);

  EMO_List_print(&alg->src, param->dbg->fp,  "source nodes");
  EMO_List_print(&alg->dest, param->dbg->fp, "destiny nodes");


  if(!EMO_Param_get_char(param, str, "mpolicy")) {
    printf("Error, mpolicy is not defined in the configuration file.\n");
    exit(1);
  }

  EMO_Migration_get_type(&alg->nmpolicy, str);

  if(!EMO_Param_get_char(param, str, "rpolicy")) {
    printf("Error, rpolicy is not defined in the configuration file.\n");
    exit(1);
  }

  EMO_Replacement_get_type(&alg->nrpolicy, str);

  if(!EMO_Param_get_int(param, &alg->xrs, "movap_xrs")) {
    printf("Error, movap_xrs is not defined in the configuration file.\n");
    exit(1);
  }

  // Divide el numero de evaluaciones de la funcion objetivo entre el numero
  // de procesadores disponibles
  if(EMO_Param_get_int(param, &v, "feval")) {
    eval = (int) v / alg->nproc;

    if(alg->myrank < v % alg->nproc)
      eval++;

    EMO_Stop_set_feval(param->stop, eval);
  }

  free(str);
}

void EMO_PAMICRO_alloc(EMO_PAMICRO *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int version, subversion, n;

  MPI_Comm_rank(MPI_COMM_WORLD, &alg->myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &alg->nproc);

  if(MPI_Get_version(&version, &subversion) != MPI_SUCCESS) {
    EMO_Debug_error(param->dbg, "MPI error in EMO_Island_alloc");
    exit(1);
  }

  EMO_Debug_printf(param->dbg, "MPI version: %d, subversion: %d", version, subversion);

  EMO_PAMICRO_load_param(alg, param);
  EMO_SMSEMOA_alloc(&alg->moea, param, pop, mop);

  if(alg->nmig > pop->lambda) {
    printf("Error, the number of migrants should be smaller than the population size (%d vs %d)\n", alg->nmig, pop->lambda);
    EMO_SMSEMOA_free(&alg->moea);
    exit(1);    
  }

  EMO_Island_alloc(&alg->comm, param->dbg, mop, &alg->src, &alg->dest, alg->nmig, alg->nmig);

  n = pop->mu * alg->nproc;

  EMO_Migration_alloc(&alg->mpolicy, param->rand, 2*n, mop->nobj);
  EMO_Migration_alloc(&alg->rpolicy, param->rand, pop->size, mop->nobj);

  EMO_Population_alloc(&alg->archive, mop, n, n);

  EMO_VPath_alloc(&alg->vpath, alg->xrs, n, n * alg->nproc, mop->nobj);
  EMO_Prune_alloc(&alg->p, n * alg->nproc, mop->nobj);

  EMO_List_alloc(&alg->lst1, pop->mu);
  EMO_List_alloc(&alg->lst2, pop->mu);
}

void EMO_PAMICRO_free(EMO_PAMICRO *alg) {
  EMO_Population_free(&alg->archive);
  EMO_Island_free(&alg->comm);
  EMO_Migration_free(&alg->mpolicy);
  EMO_Migration_free(&alg->rpolicy);
  EMO_SMSEMOA_free(&alg->moea);
  EMO_VPath_free(&alg->vpath);
  EMO_Prune_free(&alg->p);
  EMO_List_free(&alg->lst1);
  EMO_List_free(&alg->lst2);
  EMO_List_free(&alg->src);
  EMO_List_free(&alg->dest);
}

void EMO_PAMICRO_run(EMO_PAMICRO *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int nmig, stop;
  int aux, aux2;

  /* Update current epoch according to the function evaluations performed during the
   * initialization of the population */
  EMO_Stop_update_epoch(param->stop);

  EMO_Debug_printf(param->dbg, "microSMSEMOA:run");

  alg->archive.asize = 0;  /* Initialize external archive */

  while(!EMO_Stop_end(param->stop)) {

    /* Run algorithm */
    EMO_SMSEMOA_run(&alg->moea, param, pop, mop);

    /* Receive migrants at position pop->mu */
    nmig = EMO_Island_receive(&alg->comm, param, pop, mop, &stop);

    /* Attach current population plus migrants to the archive population */
    EMO_Population_copy2(&alg->archive, pop, mop, 0, pop->mu + nmig);

    memset(alg->p.filter, 0, sizeof(int) * alg->archive.size);

    aux = EMO_VPath_prune(&alg->vpath, &alg->p, alg->archive.obj, alg->archive.asize, alg->archive.mu);
    aux2 = alg->archive.mu;
    alg->archive.mu = aux;

    EMO_Population_survive(&alg->archive, NULL, NULL, mop, &alg->lst1, &alg->lst2, alg->p.filter);

    alg->archive.mu = aux2;
    alg->archive.asize = aux;

    if(alg->archive.asize >= alg->nmig) {
      /* Select individuals to migrate */
      switch(alg->nmpolicy) {
        case EMO_MIGRATION_RANDOM:
          EMO_Migration_random(&alg->mpolicy, &alg->lst1, alg->archive.asize, alg->nmig);
          break;
        case EMO_MIGRATION_ELITIST_RANDOM:
          EMO_Migration_elitist_random(&alg->mpolicy, &alg->lst1, alg->archive.obj, alg->archive.asize, alg->nmig);
          break;
        case  EMO_MIGRATION_ELITIST_RANKING:
          EMO_Migration_elitist_ranking(&alg->mpolicy, &alg->lst1, alg->archive.obj, alg->archive.asize, alg->nmig);
          break;
        case EMO_MIGRATION_FRONT:
          EMO_Migration_front(&alg->mpolicy, &alg->lst1, alg->archive.obj, alg->archive.asize);
          break;
        case EMO_MIGRATION_FRONT_RANDOM:
          EMO_Migration_front_random(&alg->mpolicy, &alg->lst1, alg->archive.obj, alg->archive.asize, alg->nmig);
          break;
        case EMO_MIGRATION_FRONT_RANKING:
          EMO_Migration_front_ranking(&alg->mpolicy, &alg->lst1, alg->archive.obj, alg->archive.asize, alg->nmig);
          break;
      }

      /* Send migrants */
      EMO_Island_send(&alg->comm, param, &alg->archive, mop, &alg->lst1, param->stop->interrupt);
    }

    /* Select individuals to replace and migrants to be incorporated into population (lst1 = migrate, lst2 = replace) */
    switch(alg->nrpolicy) {
      case EMO_REPLACEMENT_RANDOM:
        EMO_Replacement_random(&alg->rpolicy, &alg->lst1, &alg->lst2, pop->obj, pop->mu, nmig);
        break;
      case EMO_REPLACEMENT_ELITIST_RANDOM:
        EMO_Replacement_elitist_random(&alg->rpolicy, &alg->lst1, &alg->lst2, pop->obj, pop->mu, nmig);
        break;
      case EMO_REPLACEMENT_ELITIST_RANKING:
        EMO_Replacement_elitist_ranking(&alg->rpolicy, &alg->lst1, &alg->lst2, pop->obj, pop->mu, nmig);
        break;
      case EMO_ELITIST:
        EMO_Replacement_elitist(&alg->rpolicy, &alg->lst1, &alg->lst2, pop->obj, pop->mu, nmig);
    }

    /* Order population swaping lists (lst1 = missing, lst2 = available) */
    EMO_Population_survive(pop, NULL, NULL, mop, &alg->lst1, &alg->lst2, NULL);
  }

  EMO_Debug_printf(param->dbg, "Total immigrants: %d", alg->comm.incoming);
  EMO_Debug_printf(param->dbg, "Total emigrants: %d", alg->comm.outcoming);
  EMO_Debug_printf(param->dbg, "==========================================");
}


void EMO_PAMICRO_write(EMO_PAMICRO *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, int run) {
  int i, j, mu, aux, aux2, size;
  EMO_Population all;
  char *str, *file;

  mu = alg->nproc * alg->archive.mu;

  str = (char *) malloc(sizeof(char) * _MAXCHAR);
  file = (char *) malloc(sizeof(char) * _MAXCHAR);

  if(!EMO_Param_get_char(param, str, "output")) {
    printf("Error, output is not defined in the configuration file.\n");
    exit(1);
  }

  aux = alg->archive.mu;
  alg->archive.mu = alg->archive.asize;
  sprintf(file, "%s/%s_archive_%d_%s_%.2dD_R", str, param->algorithm, alg->myrank, mop->name, mop->nobj);
  EMO_Population_write(&alg->archive, NULL, mop, file, run, 0);
  alg->archive.mu = aux;

  if(param->dbg != NULL && param->dbg->level != EMO_DEBUG_OFF) {
    sprintf(file, "%s/%s_%d_%s_%.2dD_R", str, param->algorithm, alg->myrank, mop->name, mop->nobj);
    EMO_Population_write(pop, NULL, mop, file, run, 0);
  }

  // Solo lo hace el proceso 0
  // Se asume que todos los algoritmos tienen el mismo tamaño de poblacion
  // Tambien, se asume que existe un file system compartido

  MPI_Barrier(MPI_COMM_WORLD);  // necesita que todos los archivos esten listos

  if(alg->myrank == 0) {
 
    EMO_Population_alloc(&all, mop, mu, 0);

     size = 0;

     // 0 se puede omitir 
     for(i = 0; i < alg->nproc; i++) {
      sprintf(file, "%s/%s_archive_%d_%s_%.2dD_R%.2d", str, param->algorithm, i, mop->name, mop->nobj, run);

      // Construye poblacion unica a partir de diferentes archivos
      j = EMO_Population_init_from_file(&all, mop, file, size);
      size += j;
    }

    aux = EMO_VPath_prune(&alg->vpath, &alg->p, all.obj, size, alg->archive.mu);
    aux2 = all.mu;
    all.mu = size;

    sprintf(file, "%s/%s_%s_%.2dD_R", str, param->algorithm, mop->name, mop->nobj);
    EMO_Population_write(&all, alg->p.filter, mop, file, run, 0);
    all.mu = aux2;

    EMO_Population_free(&all);
  }

  free(str);
  free(file);
}

#undef _MAXCHAR

