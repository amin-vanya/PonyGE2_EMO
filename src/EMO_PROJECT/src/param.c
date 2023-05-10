
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#ifdef EMO_MPI
  #include "mpi.h"
#endif

#include "benchmark.h"
#include "param.h"
#include "io.h"

#define MAX_CHAR 1000

void EMO_Param_alloc(EMO_Param *param, int max_param) {

  param->parser = (EMO_Parser *) malloc(sizeof(EMO_Parser));
  EMO_Parser_alloc(param->parser, NULL, max_param);
 
  param->mop = NULL;
  param->dbg = NULL;
  param->rand = NULL;
  param->stop = NULL;
  param->plot = NULL;
  param->fp = NULL;

  param->prefix = NULL;
  param->flog = NULL;
  param->fsum = NULL;
  param->algorithm = NULL;
  param->subalgorithm = NULL;
}

void EMO_Param_init(EMO_Param *param, EMO_MOP *mop, char *alg, char *problem) {
  char *str1, *str2, *str3, *algorithm;
  int v = 0, flag, proc;
  double fitness = 0;

  param->mop = mop;
  param->dbg = (EMO_Debug *) malloc(sizeof(EMO_Debug));
  param->rand = (EMO_Rand *) malloc(sizeof(EMO_Rand));
  param->stop = (EMO_Stop *) malloc(sizeof(EMO_Stop));

  param->parser->dbg = param->dbg;

  param->prefix = (char *) malloc(sizeof(char) * MAX_CHAR);
  param->flog = (char *) malloc(sizeof(char) * MAX_CHAR);
  param->fsum = (char *) malloc(sizeof(char) * MAX_CHAR);
  param->algorithm = (char *) malloc(sizeof(char) * MAX_CHAR);
  param->subalgorithm = (char *) malloc(sizeof(char) * MAX_CHAR);

  str1 = (char *) malloc(sizeof(char) * MAX_CHAR);
  str2 = (char *) malloc(sizeof(char) * MAX_CHAR);
  str3 = (char *) malloc(sizeof(char) * MAX_CHAR);

  strcpy(param->algorithm, alg);

  algorithm = EMO_toupper(alg);

  #ifdef EMO_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    sprintf(param->algorithm, "p%s", algorithm);
  #else
    proc = 0;
    strcpy(param->algorithm, algorithm);
  #endif

  sprintf(param->flog, "%s_%s_%d_%d.log", param->algorithm, problem, proc, getpid()); // random number
  EMO_Debug_alloc(param->dbg, EMO_DEBUG_ON, proc, param->flog);

  strcpy(param->subalgorithm, algorithm); // required by IBEA

  if(strcmp(problem, "DEFAULT") != 0 && strcmp(problem, "default") != 0)
    EMO_Benchmark_alloc(mop, param, problem);

  if(!EMO_Param_get_char(param, str1, "output")) {
    printf("Error, output is not defined in the configuration file.\n");
    exit(1);
  }

  #ifdef EMO_MPI
    EMO_Debug_printf(param->dbg, "EMO_MPI is defined");
    sprintf(param->flog, "%s/%s_%d_%s_%.2dD.log", str1, param->algorithm, proc, mop->name, mop->nobj);
    sprintf(param->prefix, "%s/%s_%d_%s_%.2dD_R", str1, param->algorithm, proc, mop->name, mop->nobj);
    sprintf(param->fsum, "%s/%s_%d_%s_%.2dD.sum", str1, param->algorithm, proc, mop->name, mop->nobj);
  #else
    EMO_Debug_printf(param->dbg, "EMO_MPI is not defined");
    sprintf(param->flog, "%s/%s_%s_%.2dD.log", str1, param->algorithm, mop->name, mop->nobj);
    sprintf(param->prefix, "%s/%s_%s_%.2dD_R", str1, param->algorithm, mop->name, mop->nobj);
    sprintf(param->fsum, "%s/%s_%s_%.2dD.sum", str1, param->algorithm, mop->name, mop->nobj);
  #endif

  EMO_Debug_rename(param->dbg, param->flog);

  if(!EMO_Param_get_char(param, str1, "seed")) {
    printf("Error, seed is not defined in the configuration file.\n");
    exit(1);
  }

  EMO_Rand_alloc_from_file(param->rand, param->dbg, str1, proc);
  EMO_Stop_alloc(param->stop, mop, param->dbg);
  param->fp = fopen(param->fsum, "w");
  fprintf(param->fp, "#Run|Feval|Fitness|Time (s)\n");

  flag = 0;

  if(EMO_Param_get_int(param, &v, "feval")) {
    EMO_Stop_set_feval(param->stop, v);
    flag = 1;
  }

  if(EMO_Param_get_int(param, &v, "time")) { // 10 * 60
    EMO_Stop_set_time(param->stop, v);
    flag = 1;
  }

  if(EMO_Param_get_double(param, &fitness, "fitness")) {
    EMO_Stop_set_fitness(param->stop, fitness);
    flag = 1;
  }

  if(flag == 0) {
    printf("Error, missing stopping condition in configuration file\n");
    exit(1);
  }

  param->mu = 0;

  if(!EMO_Param_get_int(param, &v, "plot_freq")) {
    printf("Error, plot_freq is not defined in the configuration file.\n");
    exit(1);
  }

  param->plot = NULL;

  if(v != -1) {
    if(!EMO_Param_get_char(param, str1, "plot_pftrue")) {
      printf("Error, plot_pftrue is not defined in the configuration file.\n");
      exit(1);
    }

    if(!EMO_Param_get_char(param, str2, "plot_term")) {
      printf("Error, plot_term is not defined in the configuration file.\n");
      exit(1);
    }

    sprintf(str3, "%s, %s", param->algorithm, mop->name);

    #ifndef EMO_MPI
      if(v >= 0) {
        param->plot = (EMO_Plot *) malloc(sizeof(EMO_Plot));
        EMO_Plot_alloc(param->plot, mop->nobj, str2, str3, str1, v);
      }
    #endif
  }

  if(!EMO_Param_get_int(param, &v, "debug")) {
    printf("Error, debug is not defined in the configuration file.\n");
    exit(1);
  }

  param->dbg->level = (v == EMO_DEBUG_OFF)? v : EMO_DEBUG_ON;

  printf("DEBUG level %d\n", param->dbg->level);

  free(str1);
  free(str2);
  free(str3);
  free(algorithm);
}

void EMO_Param_alloc_from_file(EMO_Param *param, EMO_MOP *mop, char *alg, char *file, char *problem) {
  char *str1, *str2, *str3, *algorithm;
  int v = 0, flag, proc;
  double fitness = 0;

  param->mop = mop;
  param->parser = (EMO_Parser *) malloc(sizeof(EMO_Parser));
  param->dbg = (EMO_Debug *) malloc(sizeof(EMO_Debug));
  param->rand = (EMO_Rand *) malloc(sizeof(EMO_Rand));
  param->stop = (EMO_Stop *) malloc(sizeof(EMO_Stop));

  param->prefix = (char *) malloc(sizeof(char) * MAX_CHAR);
  param->flog = (char *) malloc(sizeof(char) * MAX_CHAR);
  param->fsum = (char *) malloc(sizeof(char) * MAX_CHAR);
  param->algorithm = (char *) malloc(sizeof(char) * MAX_CHAR);
  param->subalgorithm = (char *) malloc(sizeof(char) * MAX_CHAR);

  str1 = (char *) malloc(sizeof(char) * MAX_CHAR);
  str2 = (char *) malloc(sizeof(char) * MAX_CHAR);
  str3 = (char *) malloc(sizeof(char) * MAX_CHAR);

  strcpy(param->algorithm, alg);

  algorithm = EMO_toupper(alg);

  #ifdef EMO_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &proc);
    sprintf(param->algorithm, "p%s", algorithm);
  #else
    proc = 0;
    strcpy(param->algorithm, algorithm);
  #endif

  sprintf(param->flog, "%s_%s_%d.log", param->algorithm, problem, proc); // random number

  printf("log tmp %s\n", param->flog);
  EMO_Debug_alloc(param->dbg, EMO_DEBUG_ON, proc, param->flog);

  strcpy(param->subalgorithm, algorithm); // required by IBEA

  EMO_Parser_alloc_from_file(param->parser, param->dbg, file);

  if(strcmp(problem, "DEFAULT") != 0 && strcmp(problem, "default") != 0)
    EMO_Benchmark_alloc(mop, param, problem);

  if(!EMO_Param_get_char(param, str1, "output")) {
    printf("Error, output is not defined in the configuration file.\n");
    exit(1);
  }

  #ifdef EMO_MPI
    EMO_Debug_printf(param->dbg, "EMO_MPI is defined");
    sprintf(param->flog, "%s/%s_%d_%s_%.2dD.log", str1, param->algorithm, proc, mop->name, mop->nobj);
    sprintf(param->prefix, "%s/%s_%d_%s_%.2dD_R", str1, param->algorithm, proc, mop->name, mop->nobj);
    sprintf(param->fsum, "%s/%s_%d_%s_%.2dD.sum", str1, param->algorithm, proc, mop->name, mop->nobj);
  #else
    EMO_Debug_printf(param->dbg, "EMO_MPI is not defined");
    sprintf(param->flog, "%s/%s_%s_%.2dD.log", str1, param->algorithm, mop->name, mop->nobj);
    sprintf(param->prefix, "%s/%s_%s_%.2dD_R", str1, param->algorithm, mop->name, mop->nobj);
    sprintf(param->fsum, "%s/%s_%s_%.2dD.sum", str1, param->algorithm, mop->name, mop->nobj);
  #endif

  EMO_Debug_rename(param->dbg, param->flog);

  if(!EMO_Param_get_char(param, str1, "seed")) {
    printf("Error, seed is not defined in the configuration file.\n");
    exit(1);
  }

  EMO_Rand_alloc_from_file(param->rand, param->dbg, str1, proc);
  EMO_Stop_alloc(param->stop, mop, param->dbg);
  param->fp = fopen(param->fsum, "w");
  fprintf(param->fp, "#Run|Feval|Fitness|Time (s)\n");

  flag = 0;

  if(EMO_Param_get_int(param, &v, "feval")) {
    EMO_Stop_set_feval(param->stop, v);
    flag = 1;
  }

  if(EMO_Param_get_int(param, &v, "time")) { // 10 * 60
    EMO_Stop_set_time(param->stop, v);
    flag = 1;
  }

  if(EMO_Param_get_double(param, &fitness, "fitness")) {
    EMO_Stop_set_fitness(param->stop, fitness);
    flag = 1;
  }

  if(flag == 0) {
    printf("Error, missing stopping condition in configuration file\n");
    exit(1);
  }

  param->mu = 0;

  if(!EMO_Param_get_int(param, &v, "plot_freq")) {
    printf("Error, plot_freq is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_char(param, str1, "plot_pftrue")) {
    printf("Error, plot_pftrue is not defined in the configuration file.\n");
    exit(1);
  }

  if(!EMO_Param_get_char(param, str2, "plot_term")) {
    printf("Error, plot_term is not defined in the configuration file.\n");
    exit(1);
  }

  sprintf(str3, "%s, %s", param->algorithm, mop->name);

  #ifdef EMO_MPI
    param->plot = NULL;
  #else
    if(v >= 0) {
      param->plot = (EMO_Plot *) malloc(sizeof(EMO_Plot));
      EMO_Plot_alloc(param->plot, mop->nobj, str2, str3, str1, v);
    }
    else
      param->plot = NULL;
  #endif

  if(!EMO_Param_get_int(param, &v, "debug")) {
    printf("Error, debug is not defined in the configuration file.\n");
    exit(1);
  }

  param->dbg->level = (v == EMO_DEBUG_OFF)? v : EMO_DEBUG_ON;

  free(str1);
  free(str2);
  free(str3);
  free(algorithm);
}

void EMO_Param_free(EMO_Param *param) {

  EMO_Parser_free(param->parser);
  free(param->parser);

  if(param->prefix != NULL) free(param->prefix);
  if(param->flog != NULL) free(param->flog);
  if(param->fsum != NULL) free(param->fsum);
  if(param->algorithm != NULL) free(param->algorithm);
  if(param->subalgorithm != NULL) free(param->subalgorithm);
  if(param->fp != NULL) fclose(param->fp);

  if(param->dbg != NULL) {
    EMO_Debug_free(param->dbg);
    free(param->dbg);
  }

  if(param->dbg != NULL) {
    EMO_Stop_free(param->stop);
    free(param->stop);
  }

  if(param->rand != NULL) {
    EMO_Rand_free(param->rand);
    free(param->rand);
  }

  #ifndef EMO_MPI
    if(param->plot != NULL) {
      EMO_Plot_free(param->plot);
      free(param->plot);
    }
  #endif

  if(param->mop != NULL)
    EMO_Benchmark_free(param->mop);
}

void EMO_Param_save(EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, int run) {
  int d;

  EMO_Plot_run(param->plot, pop->obj, pop->mu, mop->feval, 1);

  d = param->dbg->level;
  param->dbg->level = EMO_DEBUG_ON;

  param->stop->fitness = 5;

  fprintf(param->fp, "%d %ld %e %ld.%04ld\n", run, mop->feval, param->stop->fitness, param->stop->total_time.tv_sec, param->stop->total_time.tv_usec);
  fflush(param->fp);

  EMO_Debug_printf(param->dbg, "Total execution time: %lds.%04ld", param->stop->total_time.tv_sec, param->stop->total_time.tv_usec);
  EMO_Debug_printf(param->dbg, "Number of objective function evaluations: %d", mop->feval);
  EMO_Debug_printf(param->dbg, "Performance indicator of the Pareto front: %f", param->stop->fitness);

  param->dbg->level = d;
}

void EMO_Param_set(EMO_Param *param, const char *name, const char *value) {
  EMO_Parser_set(param->parser, name, value);
}

/* Obtiene un parametro de acuerdo al tipo de dato definido */
int EMO_Param_get_int(EMO_Param *param, int *v, const char *s) {
  int i;

  if((i = EMO_Parser_find(param->parser, s)) == -1)
    return 0;

  *v = atoi(param->parser->value[i]);

  EMO_Debug_printf(param->dbg, "Parser:get_int:%s=%d", s, *v);
  return 1;
}

int EMO_Param_get_double(EMO_Param *param, double *v, const char *s) {
  int i;

  if((i = EMO_Parser_find(param->parser, s)) == -1)
    return 0;

  *v = atof(param->parser->value[i]);

  EMO_Debug_printf(param->dbg, "Parser:get_double:%s=%f", s, *v);
  return 1;
}

int EMO_Param_get_char(EMO_Param *param, char *v, const char *s) {
  int i;

  if(v == NULL || (i = EMO_Parser_find(param->parser, s)) == -1)
    return 0;

  strcpy(v, param->parser->value[i]);

  EMO_Debug_printf(param->dbg, "Parser:get_char:%s=%s", s, v);
  return 1;
}

int EMO_Param_get_vector_double(EMO_Param *param, double *v, int *size, const char *s) {
  char *tmp, *p;
  int n, m, i;

  if(v == NULL || (i = EMO_Parser_find(param->parser, s)) == -1)
    return 0;

  p = param->parser->value[i];
  n = EMO_Parser_word_count(p, ',');

  if(n > *size) {
    *size = n;
    return 0;
  }

  tmp = (char *) malloc(sizeof(char) * (strlen(p)+ 1));

  for(i = 0; i < n; i++) {
    m = EMO_Parser_get_token(tmp, &p, ',');

    if(m > 0) {
      v[i] = atof(tmp);
    }
    else {
      free(tmp);
      return 0;
    }
  }

  EMO_Debug_printv(param->dbg, v, n, "Parser:get_vector_double:%s", s);
  free(tmp);
  *size = n;

  return 1;
}

#undef MAX_CHAR

