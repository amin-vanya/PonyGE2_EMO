
#include <stdlib.h>

#include "topology.h"
#include "island.h"
#include "pmoea.h"
#include "niche.h"

#define _MAXCHAR 2000

void EMO_PMOEA_load_param(EMO_PMOEA *alg, EMO_Param *param) {
  int topo, type, i1, i2, eval, psize;
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

  if(!EMO_Param_get_char(param, str, "topology")) {
    printf("Error, topology is not defined in the configuration file.\n");
    exit(1);
  }

  EMO_Topology_get_type(&topo, str);

  if(topo != EMO_TOPOLOGY_FULL) {
    if(!EMO_Param_get_char(param, str, "comm")) {
      printf("Error, comm is not defined in the configuration file.\n");
      exit(1);
    }

    EMO_Topology_get_flow(&type, str);
  }

  if(topo == EMO_TOPOLOGY_TREE) {
    if(!EMO_Param_get_int(param, &i1, "degree")) {
      printf("Error, degree is not defined in the configuration file.\n");
      exit(1);
    }
  }

  if(topo == EMO_TOPOLOGY_TORUS || topo == EMO_TOPOLOGY_TORUSD) {
    if(!EMO_Param_get_int(param, &i1, "torus_row")) {
      printf("Error, torus_row is not defined in the configuration file.\n");
      exit(1);
    }

    if(!EMO_Param_get_int(param, &i2, "torus_col")) {
      printf("Error, torus_col is not defined in the configuration file.\n");
      exit(1);
    }
  }

  if(topo == EMO_TOPOLOGY_MESH) {
    if(!EMO_Param_get_char(param, str, "mesh_file")) {
      printf("Error, mesh_file is not defined in the configuration file.\n");
      exit(1);
    }
  }

  // Define the logical topology 
  switch(topo) {
    case EMO_TOPOLOGY_LINE:
      EMO_Topology_line(&alg->src, &alg->dest, alg->myrank, alg->nproc, type);
      break;
    case EMO_TOPOLOGY_RING: 
      EMO_Topology_ring(&alg->src, &alg->dest, alg->myrank, alg->nproc, type);
      break;
    case EMO_TOPOLOGY_STAR: 
      EMO_Topology_star(&alg->src, &alg->dest, alg->myrank, alg->nproc, type);
      break;
    case EMO_TOPOLOGY_TREE: 
      EMO_Topology_tree(&alg->src, &alg->dest, alg->myrank, alg->nproc, type, i1);
      break;
    case EMO_TOPOLOGY_FULL: 
      EMO_Topology_full(&alg->src, &alg->dest, alg->myrank, alg->nproc);
      break;
    case EMO_TOPOLOGY_TORUS: 
      EMO_Topology_torus(&alg->src, &alg->dest, alg->myrank, i1, i2, type);
      break;
    case EMO_TOPOLOGY_TORUSD: 
      EMO_Topology_torus_diagonal(&alg->src, &alg->dest, alg->myrank, i1, i2, type);
      break;
    case EMO_TOPOLOGY_MESH: 
      EMO_Topology_mesh(&alg->src, &alg->dest, alg->myrank, alg->nproc, str);
  }

  EMO_List_print(&alg->src, param->dbg->fp,  "source nodes");
  EMO_List_print(&alg->dest, param->dbg->fp, "destiny nodes");

  if(!EMO_Param_get_int(param, &alg->sync, "sync")) {
    printf("Error, sync is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_int(param, &alg->interrupt, "interrupt")) {
    printf("Error, interrupt is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_char(param, str, "mpolicy")) {
    printf("Error, mpolicy is not defined in the configuration file.\n");
    exit(1);
  }

  EMO_Migration_get_type(&alg->mpolicy, str);

  if(!EMO_Param_get_char(param, str, "rpolicy")) {
    printf("Error, rpolicy is not defined in the configuration file.\n");
    exit(1);
  }

  EMO_Replacement_get_type(&alg->rpolicy, str);

  // Divide el numero de evaluaciones de la funcion objetivo entre el numero
  // de procesadores
  if(!EMO_Param_get_int(param, &i1, "feval")) {
    printf("Error, feval is not defined in the configuration file.\n");
    exit(1);
  } 

  eval = (int) i1 / alg->nproc;
  eval = (alg->myrank < i1 % alg->nproc)? eval + 1 : eval;
  EMO_Stop_set_feval(param->stop, eval);

  // Divide el tamaño de la poblacion entre el numero de procesadores
  if(!EMO_Param_get_int(param, &psize, "psize")) {
    printf("Error, psize is not defined in the configuration file.\n");
    exit(1);
  }

  param->mu = psize / alg->nproc;
  param->mu += (alg->myrank < (psize % alg->nproc ))? 1 : 0;

  sprintf(str, "%d", param->mu);
  EMO_Param_set(param, "psize", str);


  if(!EMO_Param_get_int(param, &alg->prune_psize, "prune_psize")) {
    alg->prune_psize = psize;  // coleccion de individuos de todas las islas
    printf("Warning, prune_psize is not defined in the configuration file, setting the default value %d.\n", psize);
  }

  if(alg->prune_psize < 0 || alg->prune_psize > psize) {
    printf("Error, invalid value for prune_psize in the configuration file (%d vs %d).\n", alg->prune_psize, psize);
    exit(1);
  }

  // ojo psize no se lee en hype
  free(str);
}

void EMO_PMOEA_alloc(EMO_PMOEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, const char *str) {
  int version, subversion;

  MPI_Comm_rank(MPI_COMM_WORLD, &alg->myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &alg->nproc);

  if(MPI_Get_version(&version, &subversion) != MPI_SUCCESS) {
    EMO_Debug_error(param->dbg, "MPI error in EMO_Island_alloc");
    exit(1);
  }

  EMO_Debug_printf(param->dbg, "MPI version: %d, subversion: %d", version, subversion);

  EMO_PMOEA_load_param(alg, param);

  EMO_MOEA_alloc(&alg->moea, param, pop, mop, str); 

  if(alg->nmig > pop->lambda) {
    printf("Error, the number of migrants should be smaller than the population size (%d vs %d)\n", alg->nmig, pop->lambda);
    EMO_MOEA_free(&alg->moea);
    exit(1);    
  }

  EMO_Island_alloc(&alg->comm, param->dbg, mop, &alg->src, &alg->dest, alg->nmig, alg->nmig);
  EMO_Migration_alloc(&alg->policy, param->rand, pop->size, mop->nobj); // pop->mu

  EMO_List_alloc(&alg->lst1, pop->mu);
  EMO_List_alloc(&alg->lst2, pop->mu);
}

void EMO_PMOEA_free(EMO_PMOEA *alg) {
  EMO_Island_free(&alg->comm);
  EMO_Migration_free(&alg->policy);
  EMO_MOEA_free(&alg->moea);
  EMO_List_free(&alg->lst1);
  EMO_List_free(&alg->lst2);
  EMO_List_free(&alg->src);
  EMO_List_free(&alg->dest);
}

void EMO_PMOEA_run(EMO_PMOEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop) {
  int nmig, stop;
  //int i;

  /* Update current epoch according to the function evaluations performed during the
   * initialization of the population */
  EMO_Stop_update_epoch(param->stop);

  while(!EMO_Stop_end(param->stop)) {
    EMO_Debug_printf(param->dbg, "MOEA:run");

    /* Run algorithm */
    EMO_MOEA_run(&alg->moea, param, pop, mop);

    /* Receive migrants */
    nmig = EMO_Island_receive(&alg->comm, param, pop, mop, &stop);

    if(alg->interrupt) param->stop->interrupt = stop;

    /* Select emigrant individuals */
    switch(alg->mpolicy) {
      case EMO_MIGRATION_RANDOM:
        EMO_Migration_random(&alg->policy, &alg->lst1, pop->mu, alg->nmig);
        break;
      case EMO_MIGRATION_ELITIST_RANDOM:
        EMO_Migration_elitist_random(&alg->policy, &alg->lst1,  pop->obj, pop->mu, alg->nmig);
        break;
      case  EMO_MIGRATION_ELITIST_RANKING:
        EMO_Migration_elitist_ranking(&alg->policy, &alg->lst1, pop->obj, pop->mu, alg->nmig);
        break;
      case EMO_MIGRATION_FRONT:
        EMO_Migration_front(&alg->policy, &alg->lst1, pop->obj, pop->mu);
        break;
      case EMO_MIGRATION_FRONT_RANDOM:
        EMO_Migration_front_random(&alg->policy, &alg->lst1, pop->obj, pop->mu, alg->nmig);
        break;
      case EMO_MIGRATION_FRONT_RANKING:
        EMO_Migration_front_ranking(&alg->policy, &alg->lst1, pop->obj, pop->mu, alg->nmig);
        break;
    }

    /* Send migrants */
    EMO_Island_send(&alg->comm, param, pop, mop, &alg->lst1, param->stop->interrupt);

    /* Select individuals to replace and migrants to be incorporated into population (lst1 = migrate, lst2 = replace) */

    switch(alg->rpolicy) {
      case EMO_REPLACEMENT_RANDOM:
        EMO_Replacement_random(&alg->policy, &alg->lst1, &alg->lst2, pop->obj, pop->mu, nmig);
        break;
      case EMO_REPLACEMENT_ELITIST_RANDOM:
        EMO_Replacement_elitist_random(&alg->policy, &alg->lst1, &alg->lst2, pop->obj, pop->mu, nmig);
        break;
      case EMO_REPLACEMENT_ELITIST_RANKING:
        EMO_Replacement_elitist_ranking(&alg->policy, &alg->lst1, &alg->lst2, pop->obj, pop->mu, nmig);
        break;
      case EMO_ELITIST:
        EMO_Replacement_elitist(&alg->policy, &alg->lst1, &alg->lst2, pop->obj, pop->mu, nmig);
    }

    /* Order population swaping lists (lst1 = missing, lst2 = available) */
    EMO_Population_survive(pop, NULL, NULL, mop, &alg->lst1, &alg->lst2, NULL);

    if(alg->sync)
      MPI_Barrier(MPI_COMM_WORLD); 
  }

  EMO_Debug_printf(param->dbg, "Total immigrants: %d", alg->comm.incoming);
  EMO_Debug_printf(param->dbg, "Total emigrants: %d", alg->comm.outcoming);
  EMO_Debug_printf(param->dbg, "==========================================");
}

void EMO_PMOEA_write(EMO_PMOEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, int run) {
  int i, j, mu, *filter, size, aux;
  char file[_MAXCHAR], *str;
  EMO_Population all;
  EMO_KNN knn;

  str = (char *) malloc(sizeof(char) * _MAXCHAR);

  if(!EMO_Param_get_char(param, str, "output")) {
    printf("Error, output is not defined in the configuration file.\n");
    exit(1);
  }

  EMO_Population_write(pop, NULL, mop, param->prefix, run, 0);

  // Solo lo hace el proceso 0
  // Se asume que todos los algoritmos tienen el mismo tamaño de poblacion
  // Tambien, se asume que existe un file system compartido

  MPI_Barrier(MPI_COMM_WORLD);  // necesita que todos los archivos esten listos

  if(alg->myrank == 0) {
    //printf("algoritmo: %s\n", param->algorithm);
    //printf("subalgoritmo: %s\n", param->subalgorithm);

    mu = alg->nproc * pop->mu;

    filter = (int *) malloc(sizeof(int) * mu);
    EMO_Population_alloc(&all, mop, mu, 0);

    EMO_KNN_alloc(&knn, mu, mop->nobj);

    size = 0;

    for(i = 0; i < alg->nproc; i++) {
      sprintf(file, "%s/%s_%d_%s_%.2dD_R%.2d", str, param->algorithm, i, mop->name, mop->nobj, run);

      // Construye poblacion unica a partir de diferentes archivos
      j = EMO_Population_init_from_file(&all, mop, file, size);
      size += j;

      if(j != pop->mu)
        printf("Warning, %d individuals were read instead of %d\n", j, pop->mu);
    }

    mu = EMO_Dominance_ndset(filter, all.obj, NULL, size, mop->nobj, EMO_Dominance_strict);
    aux = all.mu; 
    all.mu = size;

    if(param->dbg != NULL && param->dbg->level != EMO_DEBUG_OFF) {
      sprintf(file, "%s/%s_ND_%s_%.2dD_R", str, param->algorithm, mop->name, mop->nobj);
      EMO_Population_write(&all, filter, mop, file, run, 0);
    }

    /* truncate procedure */
    EMO_KNN_prune(&knn, filter, alg->prune_psize, all.obj, all.mu);
    sprintf(file, "%s/%s_%s_%.2dD_R", str, param->algorithm, mop->name, mop->nobj);
    EMO_Population_write(&all, filter, mop, file, run, 0);

    all.mu = aux;

    free(filter);
    EMO_Population_free(&all);
    EMO_KNN_free(&knn);
  }

  free(str);
}

#undef _MAXCHAR

