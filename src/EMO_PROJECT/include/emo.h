
#ifndef _EMO_H
#define _EMO_H

#include <sys/time.h>
#include <stdio.h>

typedef struct _EMO_Node {
  int index;
  int item;
  struct _EMO_Node *ant;
  struct _EMO_Node *sig;
  int inblock;  /* specifies whether the node is taken from memblock or not */
} _EMO_Node;

typedef struct EMO_List {
  int size;
  int max_size;
  _EMO_Node *ini;
  _EMO_Node *fin;
  _EMO_Node **memblock;
  int available;
} EMO_List;

void EMO_List_alloc(EMO_List *l, int size);
void EMO_List_free(EMO_List *l);
void EMO_List_clear(EMO_List *l);
void EMO_List_queue(EMO_List *l, int item);
void EMO_List_dequeue(EMO_List *l, int *item);
void EMO_List_add(EMO_List *l, int item, int pos);
int  EMO_List_remove(EMO_List *l, int item);
int  EMO_List_retrieve(EMO_List *l, int *item, int pos);
int  EMO_List_get(EMO_List *l, int *item, int pos);
int  EMO_List_seek(EMO_List *l, int item);
int  EMO_List_count(EMO_List *l, int item);
void EMO_List_append(EMO_List *dest, EMO_List *src);
int  EMO_List_move(EMO_List *from, EMO_List *to, int item);
int EMO_List_move_all(EMO_List *from, EMO_List *to);
void EMO_List_print(EMO_List *l, FILE *fp, const char *s);

#define EMO_DEBUG_OFF 0
#define EMO_DEBUG_ON  1

typedef struct {
  unsigned int rank;
  unsigned int level;
  FILE *fp;
  unsigned int line;
  char *name, buf[200];
} EMO_Debug;

void EMO_Debug_alloc(EMO_Debug *dbg, unsigned int level, unsigned int rank, const char *file);
void EMO_Debug_free(EMO_Debug *dbg);
void EMO_Debug_rename(EMO_Debug *dbg, const char *file);
void EMO_Debug_printf(EMO_Debug *dbg, const char *format, ...);
void EMO_Debug_printv(EMO_Debug *dbg, double *x, int n, const char *fmt, ...);
void EMO_Debug_error(EMO_Debug *dbg, const char *fmt, ...);

#define NN 312

typedef struct {
  unsigned long long mt[NN]; /* The array for the state vector */
  int mti; // = NN + 1;      /* mti==NN+1 means mt[NN] is not initialized */
  unsigned long long seed;   /* Current seed */
  FILE *fp;                  /* Pointer to the seed file */
  char *token;               /* Temporary array for parsing the seed file */
  EMO_Debug *dbg;
} EMO_Rand;


void EMO_Rand_alloc(EMO_Rand *rnd, EMO_Debug *dbg, unsigned long long seed);
int EMO_Rand_next_seed(EMO_Rand *rnd, int skip);
void EMO_Rand_alloc_from_file(EMO_Rand *rnd, EMO_Debug *dbg, const char *file, int skip);
void EMO_Rand_free(EMO_Rand *rnd);
double EMO_Rand_prob1(EMO_Rand *rnd);
double EMO_Rand_prob2(EMO_Rand *rnd);
double EMO_Rand_prob3(EMO_Rand *rnd);
double EMO_Rand_real1(EMO_Rand *rnd, double a, double b);
double EMO_Rand_real2(EMO_Rand *rnd, double a, double b);
double EMO_Rand_real3(EMO_Rand *rnd, double a, double b);
int EMO_Rand_int1(EMO_Rand *rnd, int a, int b);
int EMO_Rand_int2(EMO_Rand *rnd, int a, int b);
int EMO_Rand_int3(EMO_Rand *rnd, int a, int b);  // me salio b
int EMO_Rand_flip(EMO_Rand *rnd, double prob);
void EMO_Rand_shuffle(EMO_Rand *rnd, int *idx, int size);
double EMO_Rand_box_muller(EMO_Rand *rnd, double m, double s);

/*typedef struct {
  long idum;
} EMO_Rand;

double EMO_Rand_prob3(EMO_Rand *rnd);
int EMO_Rand_int1(EMO_Rand *rnd, int a, int b);
int EMO_Rand_flip(EMO_Rand *rnd, double prob);
void EMO_Rand_alloc(EMO_Rand *rnd, EMO_Debug *dbg, unsigned long long seed);
int EMO_Rand_next_seed(EMO_Rand *rnd, int skip);
void EMO_Rand_alloc_from_file(EMO_Rand *rnd, EMO_Debug *dbg, const char *file, int skip);
void EMO_Rand_free(EMO_Rand *rnd);
double EMO_Rand_real1(EMO_Rand *rnd, double a, double b);
void EMO_Rand_shuffle(EMO_Rand *rnd, int *idx, int size);*/


double EMO_vnorm(double *x, double p, int n);
void EMO_vmul(double *y, double *x, double alpha, int n);
void EMO_vsum(double *y, double *x, double alpha, int n);
double EMO_vdot(double *x, double *y, int n);
void EMO_vproj(double *r, double *x, double *y, int n);
void EMO_vadd(double *r, double *x, double *y, int n);
void EMO_vdiff(double *r, double *x, double *y, int n);
void EMO_vorth(double *r, double *p, double *v, int n);
double EMO_vdist(double *x, double *y, int n);
void EMO_vzero(double *x, int n);
int EMO_isvzero(double *x, int n);
void EMO_vaxes(double *x, int n);
void EMO_vcopy(double *dest, double *src, int n);
void EMO_vprint(FILE *fp, double *x, int n, const char *fmt, ...);
void EMO_vprinti(FILE *fp, int *x, int n, const char *fmt, ...);

void EMO_mprint(FILE *fp, double *A, int n, int m, const char *fmt, ...);
void EMO_matmul(double *a, double *b, double *c, int n, int m, int l);
int EMO_minverse(double *B, double *A, int n, int nthread, int metodo, int tridiag);


#define PI 3.1415926535897932384626433832795
#define max(A,B) ((A)>(B) ? (A) : (B))  //ver funcion fmin:  double fmin(double x, double y); en math
#define min(A,B) ((A)<(B) ? (A) : (B))

double EMO_dmin(int *idx, double *arr, int *filter, int size);
double EMO_dmax(int *idx, double *arr, int *filter, int size);
int EMO_min(int *idx, int *arr, int *filter, int size);
int EMO_max(int *idx, int *arr, int *filter, int size);
void EMO_maxBound(double *max, double *data, int *filter, int row, int col);
void EMO_maxBound2(double *max, int *idx, double *data, int *filter, int row, int col);
void EMO_minBound(double *min, double *data, int *filter, int row, int col);
void EMO_maxminBound(double *max, double *min, double *data, int *filter, int row, int col);
void EMO_findminBound(int *min, double *data, int *filter, int row, int col);
void EMO_findmaxminBound(int *max, int *min, double *data, int *filter, int row, int col);
void EMO_shift(double *shift, double *data, int *filter, int size, double *z, int dim);
void EMO_normalize(double *norm, double *data, int *filter, int size, double *min, double *max, int dim);
void EMO_normalize2(double *norm, double *data, int *filter, int size, double *max, int dim);
void EMO_normalize3(double *norm, double *data, int *filter, int size, double *min, double *max, int dim);
void EMO_lsq(double *a, double *tmp1, double *tmp2, double *y, int n);
int EMO_find(int *data, int size, int elem);

double EMO_mean(double *data, int *filter, int size);
double EMO_quartile(double *q, double *data, int *filter, double **sort, int size);
double EMO_median(int *idx, double *data, int *filter, double **sort, int size);
double EMO_var(double *data, int *filter, double mean, int size);
double EMO_std(double *data, int *filter, double mean, int size);

#define EMO_BINARY 0
#define EMO_REAL 1

typedef struct EMO_MOP {
  int nvar;            /* number of decision variables  */
  int nobj;            /* number of objective functions */
  int ncon;            /* number of constraints */

  int npos;            /* number of position related parameters (WFG test problems) */
  double gamma;        /* Lame super-spheres */

  double *xmin, *xmax; /* valid ranges of decision variables */

  void (*f)(struct EMO_MOP *mop, double *, double *);  /* pointer to the evaluation function */
  void (*fc)(struct EMO_MOP *mop, double *, double *g, double *x); /* pointer to the evaluation function with constraints */

  char *name;          /* name of the problem */
  int coding;          /* representation: real or binary */
  int L;               /* chromosome length (binary representation) */
  unsigned long feval; /* number of objective function evaluations */

                       /* Auxiliary variables */
  double *t;           /* theta for DTLZ5, DTLZ6 and transition vector for WFG test problems */
  double *x;           /* underlying parameters */
  double *h;           /* shape */
  double *S;           /* constant */
  double *A;           /* constant */
  double *y;           /* temporary vector */
  double *w;           /* weight vector */
  double (*g)(double );  /* g function in test problems based on Lame Superspheres */
} EMO_MOP;


typedef struct {
  double *var;  // RHG, representacion binaria, union
  double *obj;
  double *con;
  int *vio;    /* violation of constraints */
  double *cv;  /* constraint violation value */
  int mu;      /* parent population size */
  int lambda;  /* offspring population size */
  int size;    /* total population size */
  double *vdummy, *odummy, *cdummy; // Internal use
  int asize;  /* current archive size (only when the population is secondary) */
  char *format;   /* writing format */
} EMO_Population;

void EMO_Population_alloc(EMO_Population *pop, EMO_MOP *mop, int mu, int lambda);
void EMO_Population_init(EMO_Population *pop, EMO_Rand *rand, EMO_MOP *mop);
int EMO_Population_init_from_file(EMO_Population *pop, EMO_MOP *mop, const char *prefix, int start);
void EMO_Population_write(EMO_Population *pop, int *filter, EMO_MOP *mop, const char *prefix, int run, unsigned long itera);
void EMO_Population_free(EMO_Population *pop);
void EMO_Population_copy(EMO_Population *pop, int *fit1, double *fit2, EMO_MOP *mop, int i, int j);
void EMO_Population_copy2(EMO_Population *dest, EMO_Population *src, EMO_MOP *mop, int start, int size);
void EMO_Population_swap(EMO_Population *pop, int *fit1, double *fit2, EMO_MOP *mop, int i, int j);
void EMO_Population_survive(EMO_Population *pop, int *fit1, double *fit2, EMO_MOP *mop, EMO_List *missing, EMO_List *available, int *filter);
void EMO_Population_evaluate(EMO_Population *pop, EMO_MOP *mop, int start, int size);
void EMO_Population_constrain(EMO_Population *pop, EMO_MOP *mop, int start, int size);

// Funciones de utileria, manejo de cadenas
char *EMO_toupper(const char *src);
char *EMO_tolower(const char *src);
int EMO_Dictionary_find(const char *dicc[], char *pattern);
void EMO_Dictionary_print(FILE *fp, const char *dicc[], const char *s);

typedef struct {
  int type;               // define the type of stop, by evals, time, fitness, epoch or a combination
  unsigned long *feval, max_feval;             // evaluations
  unsigned long epoch, max_epoch;              // epoch (for pMOEAs)
  struct timeval t1, t2, total_time, max_time; // execution time
  double fitness, max_fitness;                 // indicator value: hv, igd, delta_p
  int interrupt;          // program interruption
  int msg;                // print summary
  EMO_Debug *dbg;         // associated log
} EMO_Stop;

void EMO_Stop_alloc(EMO_Stop *stop, EMO_MOP *mop, EMO_Debug *dbg);
void EMO_Stop_free(EMO_Stop *stop);
void EMO_Stop_start(EMO_Stop *stop);
void EMO_Stop_set_feval(EMO_Stop *stop, unsigned long max_feval);
void EMO_Stop_set_time(EMO_Stop *stop, time_t seconds);
void EMO_Stop_set_fitness(EMO_Stop *stop, double max_fitness);
void EMO_Stop_set_epoch(EMO_Stop *stop, int epoch);
void EMO_Stop_update_epoch(EMO_Stop *stop);
void EMO_Stop_now(EMO_Stop *stop);
int EMO_Stop_end(EMO_Stop *stop);


typedef struct {
  char **name;        /* Parameter names */
  char **value;       /* Parameter value */
  char *tmp1, *tmp2;  /* Temporary variables */
  int size;           /* Number of parameters */
  int current;        /* Current parameter (only for set function) */
  EMO_Debug *dbg;
} EMO_Parser;

int  EMO_Parser_word_count(const char *s, const char delim);
int  EMO_Parser_get_token(char *dest, char **src, const char delim);
void EMO_Parser_alloc(EMO_Parser *parser, EMO_Debug *dbg, int size);
void EMO_Parser_alloc_from_file(EMO_Parser *parser, EMO_Debug *dbg, char *file);
void EMO_Parser_free(EMO_Parser *parser);
void EMO_Parser_set(EMO_Parser *parser, const char *name, const char *value);
int  EMO_Parser_find(EMO_Parser *parser, const char *s);
void EMO_Parser_print(EMO_Parser *parser);


typedef struct {
  FILE **gp;
  int dim;
  long freq;
  char *term;
  char *title;
  char *pftrue;
} EMO_Plot;

void EMO_Plot_alloc(EMO_Plot *p, int dim, const char *term, const char *title, const char *pftrue, long freq);
void EMO_Plot_free(EMO_Plot *p);
void EMO_Plot_run(EMO_Plot *p, double *data, int size, long feval, int end);

typedef struct {
  EMO_Parser *parser;
  EMO_Debug *dbg;
  EMO_Rand *rand;
  EMO_Stop *stop;
  EMO_Plot *plot;
  EMO_MOP *mop;  // for freeing resources

  char *algorithm;
  char *subalgorithm;
  FILE *fp; // summary
  int mu;  // population size
  double Pc;
  double Pm;
  double Nc;
  double Nm;

  char *prefix;
  char *flog;
  char *fsum;
} EMO_Param;

void EMO_Param_alloc(EMO_Param *param, int max_param);
void EMO_Param_alloc_from_file(EMO_Param *param, EMO_MOP *mop, char *algorithm, char *file, char *problem);
void EMO_Param_free(EMO_Param *param);
void EMO_Param_init(EMO_Param *param, EMO_MOP *mop, char *alg, char *problem);
void EMO_Param_save(EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, int run);
void EMO_Param_set(EMO_Param *param, const char *name, const char *value);
int EMO_Param_get_int(EMO_Param *param, int *v, const char *s);
int EMO_Param_get_double(EMO_Param *param, double *v, const char *s);
int EMO_Param_get_char(EMO_Param *param, char *v, const char *s);
int EMO_Param_get_vector_double(EMO_Param *param, double *v, int *size, const char *s);


void EMO_crossSBX(double *child1, double *child2, double *parent1, double *parent2, EMO_Rand *rand, EMO_MOP *mop, double Pc, double eta_c);
void EMO_mutatePolynom(double *child, EMO_Rand *rand, EMO_MOP *mop, double Pm, double eta_m);
  
// Representacion binaria
void EMO_one_point_crossover(double *child1, double *child2, double *parent1, double *parent2, EMO_Rand *rand, EMO_MOP *mop, double Pc);
void EMO_one_bit_mutation(double *child, EMO_Rand *rand, EMO_MOP *mop);
   
void tournament_selection_dominance(int *ind, EMO_Rand *rand, int *seedp, double *obj, int n, int nobj);
void tournament_selection_rank(int *ind, EMO_Rand *rand, int *seedp, int *rank, int n);
void tournament_selection_rank_crowding(int *ind, EMO_Rand *rand, int *seedp, int *rank, double *cd, int n);
void tournament_selection_fitness(int *ind, EMO_Rand *rand, int *seedp, double *fit, int n);
void tournament_selection_r2(int *ind, EMO_Rand *rand, int *seedp, double *r2, int n, int m);
int tournament_selection_fitness2(double *fit, EMO_Rand *rand, int size);

typedef struct EMO_Utility {
  char *name;
  int nobj;
  int inverted;            /* inverted scalarizing function */
  int tlpbi_H;
  int qpbi_H;
  double wzero;            /* Small value of wi, used to avoid division by zero */
  double wcp_p;
  double wcp2_p;
  double wpo_p;
  double wpo2_p;
  double ewc_p;
  double ewc2_p;
  double wn_p;
  double wn2_p;
  double vads_p;
  double vads2_p;
  double ache_alpha;
  double ache2_alpha;
  double mche_alpha;
  double mche2_alpha;
  double aasf_alpha;
  double tlpbi_alpha;
  double qpbi_alpha;
  double gsf_alpha;
  double gsf2_alpha;
  double nsf_alpha;
  double nsf2_alpha;
  double cs_alpha;
  double cs2_alpha;
  double pbi_theta;
  double qpbi_theta;
  double tlpbi_theta1;
  double tlpbi_theta2;
  double gsf_beta;
  double gsf2_beta;
  double didass_beta;
  double didass2_beta;
  double tlpbi_dstar;
  double qpbi_dstar;
  double *pbi_v;      /* temporary vectors */
  double *tlpbi_v;
  double *qpbi_v;
  double *didass_v;
  double *didass2_v;
  double (*uf)(struct EMO_Utility *u, double *w, double *x);
} EMO_Utility;

typedef double (*EMO_UtilityFunction)(EMO_Utility *u, double *w, double *x);

extern const char *EMO_Utility_list[];
extern const EMO_UtilityFunction fdic[];

void EMO_Utility_alloc(EMO_Utility *u, EMO_Param *param, int nobj, const char *str);
void EMO_Utility_free(EMO_Utility *u);

double EMO_Utility_wcp(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wcp2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wpo(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wpo2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ewc(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ewc2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wpr(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wpr2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ws(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ws2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wn(EMO_Utility *u, double *w, double *x);
double EMO_Utility_wn2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ls(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ls2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_che(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ache(EMO_Utility *u, double *w, double *x);
double EMO_Utility_ache2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_mche(EMO_Utility *u, double *w, double *x);
double EMO_Utility_mche2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_asf(EMO_Utility *u, double *w, double *x);
double EMO_Utility_aasf(EMO_Utility *u, double *w, double *x);
double EMO_Utility_pbi(EMO_Utility *u, double *w, double *x);
void EMO_Utility_update_dstar(EMO_Utility *u, double *zmin, double *zmax);
double EMO_Utility_tlpbi(EMO_Utility *u, double *w, double *x);
double EMO_Utility_qpbi(EMO_Utility *u, double *w, double *x);
double EMO_Utility_gsf(EMO_Utility *u, double *w, double *x);
double EMO_Utility_gsf2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_nsf(EMO_Utility *u, double *w, double *x);
double EMO_Utility_nsf2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_cs(EMO_Utility *u, double *w, double *x);
double EMO_Utility_cs2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_vads(EMO_Utility *u, double *w, double *x);
double EMO_Utility_vads2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_didass(EMO_Utility *u, double *w, double *x);
double EMO_Utility_didass2(EMO_Utility *u, double *w, double *x);
double EMO_Utility_refsf(EMO_Utility *u, double *w, double *x);

typedef struct {
  EMO_Utility utl;
  double *ideal;
  double *nadir;
  int *xtrm;
  double *axis0;
  double *axis1;
  int nobj;
} EMO_Refpoint;

void EMO_Refpoint_alloc(EMO_Refpoint *ref, EMO_Param *param, int nobj);
void EMO_Refpoint_free(EMO_Refpoint *ref);
void EMO_Refpoint_update_ideal(EMO_Refpoint *ref, double *data, int *filter, int size);
void EMO_Refpoint_update_nadir(EMO_Refpoint *ref, double *data, int *filter, int size);

double *EMO_File_read(double *data, int *row, int *col, const char *str, int start);
double *EMO_File_read_skip(double *data, int *row, int *col, const char *str, int start, int skip);
void EMO_File_write(double *data, int *filter, int row, int col, const char *str, const char *fmt, unsigned long itera);


double EMO_hv2(double *data, int n, const double *ref, int d);
double EMO_hv2_contribution(double *deltahv, double *data, int *enable, int n, const double *ref, int d);


typedef double OBJECTIVE;

typedef struct
{
	OBJECTIVE *objectives;
        int idx;  // RHG, se usa en hv_iwfg.h
} POINT;

typedef struct
{
	int nPoints;
	int n;    // RHG, se usa en hv_iwfg.h
 	POINT *points;
} FRONT;

typedef struct
{
	double width;
	FRONT front;
	int index;
} SLICE;


typedef struct {
  int col;      // the number of objectives 
  int row;
  FRONT *fs;    // memory management stuff 
  FRONT f;
  int fr;       // current depth 
  int safe;     // the number of points that don't need sorting 
} EMO_HV;



void EMO_HV_alloc(EMO_HV *hv, int max_row, int col);
void EMO_HV_free(EMO_HV *hv);
double EMO_HV_run(EMO_HV *hv, double *data, int *enable, int row, const double *ref);
double EMO_HV_run2(EMO_HV *hv, double *data, int *enable, int row, const double *ref);
double EMO_HV_contribution(EMO_HV *hv, double *deltahv, double *data, int *enable, int row, const double *ref, int col);


typedef struct
{
	int	n;		// The number of objectives 
	//POINT	ref;		// The refernce point	
	//POINT	dirs;		// Records the directions of the objectives 

	FRONT	*fs;		// Memory manage stuff 
        FRONT f;
	int	fr;		// Current depth 
	int	maxm;		// Maximum number of points 
	int	maxn;		// Maximum number of objectives 
	int	safe;		// The number of points that don't need sorting

	double* partial;	// Partial exclusive hipervolumes
	int*	heap;		// Heap-based priority queue
	int	heapsize;	// Number of points in queue
	SLICE	**stacks;	// Set of slices per point per slicing depth
	int	*stacksize;	// Current slicing depth per point

	int*	gorder;		// Objective order used by comparison funtions
	int**	torder;		// Order of objectives per point
	int**	tcompare;	
	FRONT*	fsorted;	// Front sorted in each objective
        int maxStackSize;       // Maximum stack size
} EMO_IWFG;
 
void EMO_IWFG_alloc(EMO_IWFG *hv, int maxm, int maxn);
void EMO_IWFG_free(EMO_IWFG *hv);
int EMO_IWFG_run(EMO_IWFG *hv, double *data, int *enable, int row, const double *ref, double *hv_worst);
void IWFG_quicksort (void *const pbase, size_t total_elems, int start, EMO_IWFG *hv, size_t size, int(*cmp)(const void *, const void *, int, EMO_IWFG *));

double EMO_Indicator_gd(double *f, int *filter, int n, double *pf, int m, int dim, double p);
double EMO_Indicator_igd(double *f, int *filter, int n, double *pf, int m, int dim, double p);
double EMO_Indicator_igd2(double *f, int *filter, int n, double *pf, int m, int dim);
double EMO_Indicator_igd_plus(double *f, int *filter, int n, double *pf, int m, int dim);
double EMO_Indicator_deltap(double *f, int *filter, int n, double *pf, int m, int dim, double p);
double EMO_Indicator_r2(double *data, int *filter, int size, double *W, int wsize, EMO_Utility *utl);
void EMO_Indicator_r2_ranking(double *rank, double **sort, double *norm, double *tmp, double *data, int size, double *W, int wsize, EMO_Utility *utl);
double EMO_Indicator_epsilon_multiplicative(double *a, int *filter, int na, double *b, int nb, int dim);
double EMO_Indicator_epsilon_additive(double *a, int *filter, int na, double *b, int nb, int dim);
double EMO_Indicator_maximin(double *fit, double *f, int *filter, int n, int dim);
double EMO_Indicator_onvg(double *f, int *filter, int n, int dim);
double EMO_Indicator_onvgr(double *a, int *filter, int na, double *b, int nb, int dim);
double EMO_Indicator_c(double *a, int *filter, int na, double *b, int nb, int dim);
double EMO_Indicator_spacing(double *f, int *filter, int n, int dim);
double EMO_Indicator_senergy(double *fit, double *f, int *filter, int n, int dim);
double EMO_Indicator_senergy_update(double *fit, double *f, int *filter, int n, int dim, int idx, int action);
double EMO_Indicator_solow_polasky(double *fit, double *f, int *filter, int n, int dim);

double EMO_WFG_correct_to_01(double a);
double EMO_WFG_sh_linear_1(double *x, int M);
double EMO_WFG_sh_linear_m(double *x, int M, int m);
double EMO_WFG_sh_linear_M(double *x);
void EMO_WFG_sh_linear(double *h, double *x, int M);
double EMO_WFG_sh_convex_1(double *x, int M);
double EMO_WFG_sh_convex_m(double *x, int M, int m);
double EMO_WFG_sh_convex_M(double *x);
void EMO_WFG_sh_convex(double *h, double *x, int M);
double EMO_WFG_sh_concave_1(double *x, int M);
double EMO_WFG_sh_concave_m(double *x, int M, int m);
double EMO_WFG_sh_concave_M(double *x);
void EMO_WFG_sh_concave(double *h, double *x, int M);
double EMO_WFG_sh_mixed_M(double *x, double A, double alpha);
double EMO_WFG_sh_disc_M(double *x, double A, double alpha, double beta);
double EMO_WFG_tr_b_poly(double y, double alpha);
double EMO_WFG_tr_b_flat(double y, double A, double B, double C);
double EMO_WFG_tr_b_param(double y, double u, double A, double B, double C);
double EMO_WFG_tr_s_linear(double y, double A);
double EMO_WFG_tr_s_decept(double y, double A, double B, double C);
double EMO_WFG_tr_s_multi(double y, double A, double B, double C);
double EMO_WFG_tr_r_sum(double *y, double *w, int n);
double EMO_WFG_tr_r_nonsep(double *y, double A, int n);
void EMO_WFG_wfg1_t1(EMO_MOP *mop, double *z);
void EMO_WFG_wfg1_t2(EMO_MOP *mop, double *z);
void EMO_WFG_wfg1_t3(EMO_MOP *mop, double *z);
void EMO_WFG_wfg1_t4(EMO_MOP *mop, double *z);
void EMO_WFG_wfg2_t2(EMO_MOP *mop, double *z);
void EMO_WFG_wfg2_t3(EMO_MOP *mop, double *z);
void EMO_WFG_wfg4_t1(EMO_MOP *mop, double *z);
void EMO_WFG_wfg4_t2(EMO_MOP *mop, double *z);
void EMO_WFG_wfg5_t1(EMO_MOP *mop, double *z);
void EMO_WFG_wfg6_t2(EMO_MOP *mop, double *z);
void EMO_WFG_wfg7_t1(EMO_MOP *mop, double *z);
void EMO_WFG_wfg8_t1(EMO_MOP *mop, double *z);
void EMO_WFG_wfg9_t1(EMO_MOP *mop, double *z);
void EMO_WFG_wfg9_t2(EMO_MOP *mop, double *z);
void EMO_WFG_calc_x(double *x, double *t, double *A, int M);
void EMO_WFG_calc_f(double *f, double *x, double *S, double *h, int M);
void EMO_WFG_wfg1(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg2(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg3(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg4(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg5(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg6(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg7(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg8(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg9(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg1_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg2_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg3_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg4_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg5_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg6_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg7_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg8_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_wfg9_minus(EMO_MOP *mop, double *f, double *z);
void EMO_WFG_range(EMO_MOP *mop);
int EMO_WFG_alloc(EMO_MOP *mop);
void EMO_WFG_free(EMO_MOP *mop);

void EMO_SZF_f1(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f2(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f3(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f4(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f5(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f6(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f7(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f8(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f9(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f10(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f11(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f12(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f13(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f14(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f15(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f16(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f17(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f18(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f19(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_f20(EMO_MOP *mop, double *F, double *y);
void EMO_SZF_range(EMO_MOP *mop);

extern const char *EMO_Benchmark_list[];
extern const char *EMO_Benchmark_listc[];

void EMO_Benchmark_alloc(EMO_MOP *mop, EMO_Param *param, const char *problem);
void EMO_Benchmark_free(EMO_MOP *mop);

typedef int (*EMO_Dominance)(double *, double *, int);

int EMO_Dominance_weak(double *x, double *y, int n);
int EMO_Dominance_strict(double *x, double *y, int n);
int EMO_Dominance_strong(double *x, double *y, int n);
int EMO_Dominance_alpha(double *x, double *y, double *tmp, double *alpha, int n);
int EMO_Dominance_favor(double *x, double *y, int n);

int EMO_Dominance_incomparable(double *x, double *y, int n, EMO_Dominance r);
int EMO_Dominance_indifferent(double *x, double *y, int n);
int EMO_Dominance_constraint(double *x, double *y, int n, double *gx, double *gy, int k, EMO_Dominance r);
int EMO_Dominance_constraint2(double *x, double *y, int n, double *gx, double *gy, int k, EMO_Dominance r);
int EMO_Dominance_constraint3(double *x, double *y, int n, double cvx, double cvy, EMO_Dominance r);
int EMO_Dominance_feasible(double *g, int k);
int EMO_Dominance_ndset(int *nd, double *data, int *filter, int row, int col, EMO_Dominance r);

typedef struct {
  int size;         // numero de soluciones a evaluar
  int *n;           // auxiliar, numero de veces que una solucion i ha sido dominada
  int *rank;        // jerarquia de los individuos
  int nfront;       // numero de frentes
  EMO_List *S;      // auxiliar, lista de soluciones que domina el i-esimo individuo
  EMO_List *front;  // Lista de frentes
} EMO_NDSort;

void EMO_NDSort_alloc(EMO_NDSort *nd, int max_size);
void EMO_NDSort_free(EMO_NDSort *nd);
void EMO_NDSort_run(EMO_NDSort *nd, double *obj, int nobj, double *cv, int *filter, int size);
void EMO_NDSort_run2(EMO_NDSort *nd, double *obj, int nobj, double *con, int ncon, int *filter, int size);

void EMO_Dominance_ranking_sum(double *rank, double **sort, double *data, int size, int nobj);
void EMO_Dominance_ranking_min(double *rank, double **sort, double *data, int size, int nobj);

typedef struct {
  EMO_List lst;
  EMO_NDSort nd;
  EMO_Refpoint ref;
  double *min;
  double *ideal;
  double *nadir;
  double **sort;    /* temporary array for sorting population */
  double *norm; 
  double *vnorm;    /* norm value for each individual */
  int *filter;
  int *tmp;
  int ssize;        /* temporal population size (for memory free) */
} EMO_Prune;

void EMO_Prune_alloc(EMO_Prune *p, int size, int nobj);
void EMO_Prune_free(EMO_Prune *p);


void EMO_nicheCount(double *niche, double *data, int row, int col, double alpha, double sigma_share);
void EMO_crowdingDistance(double *cd, double **sort, double *data, int row, double *max, double *min, int col);
void EMO_crowdingDistance2(double *cd, double **sort, double *data, EMO_List *front, double *max, double *min, int col);
void EMO_wdist(double *wd, double *data, int size, double *W, int wsize, double *tmp, int dim);

// Metodo de truncado de soluciones basado en SPEA2
typedef struct {
  EMO_List *lnn;          /* array of lists for storing neighbors */
  EMO_List copy;
  double **sort;          /* temporary array for sorting */
  double **dist;
  int *min, *max;         /* individuals that are part of the extreme points */
  int dim, ssize;
} EMO_KNN;

void EMO_KNN_alloc(EMO_KNN *knn, int size, int dim);
void EMO_KNN_free(EMO_KNN *knn);
void EMO_KNN_prune(EMO_KNN *knn, int *filter, int max_size, double *data, int size);

void EMO_knn(EMO_List *l, double **sort, double *W, int size, int dim, int k, int itself);
void EMO_knn2(EMO_List *l, EMO_List *copy, double **dist, double **sort, double *W, int size, int dim);
void EMO_kNN2(EMO_List *lst, double **sort, double *data, int row, int col);

void EMO_Topology_get_type(int *topo, const char *str);
void EMO_Topology_get_flow(int *flow, const char *str);
void EMO_Topology_line(EMO_List *in, EMO_List *out, int rank, int size, int type);
void EMO_Topology_ring(EMO_List *in, EMO_List *out, int rank, int size, int type);
void EMO_Topology_star(EMO_List *in, EMO_List *out, int rank, int size, int type);
void EMO_Topology_tree(EMO_List *in, EMO_List *out, int rank, int size, int type, int degree);
void EMO_Topology_full(EMO_List *in, EMO_List *out, int rank, int size);
void EMO_Topology_torus(EMO_List *in, EMO_List *out, int rank, int row, int col, int type);
void EMO_Topology_torus_diagonal(EMO_List *in, EMO_List *out, int rank, int row, int col, int type);
void EMO_Topology_mesh(EMO_List *in, EMO_List *out, int rank, int size, const char *file);

#define EMO_MIGRATION_RANDOM            0
#define EMO_MIGRATION_ELITIST_RANDOM    1
#define EMO_MIGRATION_ELITIST_RANKING   2
#define EMO_MIGRATION_FRONT             3
#define EMO_MIGRATION_FRONT_RANDOM      4
#define EMO_MIGRATION_FRONT_RANKING     5

#define EMO_REPLACEMENT_RANDOM          0
#define EMO_REPLACEMENT_ELITIST_RANDOM  1
#define EMO_REPLACEMENT_ELITIST_RANKING 2
#define EMO_ELITIST                     3

typedef struct {
  EMO_NDSort nd;
  EMO_Rand *rnd;
  EMO_List best;
  EMO_List worst;
  int *filter;
  int *seq;  /* Sequence for EMO_Migration_random */
  int *tmp;  /* Arbitrary sequence */
  int max_size;
  int nobj;
} EMO_Migration;

void EMO_Migration_alloc(EMO_Migration *m, EMO_Rand *rnd, int max_size, int nobj);
void EMO_Migration_free(EMO_Migration *m);

void EMO_Migration_get_type(int *type, const char *str);
void EMO_Replacement_get_type(int *type, const char *str);

void EMO_Migration_random(EMO_Migration *m, EMO_List *l, int size, int nmig);
void EMO_Migration_elitist_random(EMO_Migration *m, EMO_List *l, double *data, int size, int nmig);
void EMO_Migration_elitist_ranking(EMO_Migration *m, EMO_List *l, double *data, int size, int nmig);
void EMO_Migration_front(EMO_Migration *m, EMO_List *l, double *data, int size);
void EMO_Migration_front_random(EMO_Migration *m, EMO_List *l, double *data, int size, int nmig);
void EMO_Migration_front_ranking(EMO_Migration *m, EMO_List *l, double *data, int size, int nmig);
void EMO_Replacement_random(EMO_Migration *m, EMO_List *migrate, EMO_List *replace, double *data, int size, int nmig);
void EMO_Replacement_elitist_random(EMO_Migration *m, EMO_List *migrate, EMO_List *replace, double *data, int size, int nmig);
void EMO_Replacement_elitist_ranking(EMO_Migration *m, EMO_List *migrate, EMO_List *replace, double *data, int size, int nmig);
void EMO_Replacement_elitist(EMO_Migration *m, EMO_List *migrate, EMO_List *replace, double *data, int size, int nmig);

int EMO_compare_asc(const void **a, const void **b);
int EMO_compare_desc(const void **a, const void **b);
int EMO_compare_rows(const void *a, const void *b, int start, int n);
void EMO_quicksort (void *const pbase, size_t total_elems, int start, int dimension, size_t size, int(*cmp)(const void *, const void *, int, int));
int EMO_binsert(const void *A, const void *key, int total_elems, int start, int dimension, size_t size, int(*cmp)(const void *, const void *, int, int));

typedef struct {
  EMO_NDSort nd;
  EMO_HV hv;
  EMO_IWFG iwfg;
  double *max, *min;
  int *filter;
  double *chv;

  /* mpirun */
  double **sort;       // temporary array for sorting population
  int ssize;           // temporal population size (for memory free)
  EMO_List lst1, lst2; // temporary lists
  int iwfg_flag;       // incremental IWFG algorithm
  double sfactor;      // reference point for the hypervolume indicator
} EMO_SMSEMOA;

void EMO_SMSEMOA_alloc(EMO_SMSEMOA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_SMSEMOA_free(EMO_SMSEMOA *alg);
void EMO_SMSEMOA_run(EMO_SMSEMOA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_SMSEMOA_prun(EMO_SMSEMOA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

typedef struct {
  EMO_NDSort nd;
  EMO_List lst1, lst2; // temporary lists
  double *norm;        // normalized objective functions
  double *a;           // Intercepts of the hyper-plane
  double *one;         // Vector of ones
  double *axis;        // coordinate axes
  int wsize;           // number of weight vectors
  double *W;           // weight vectors
  double *waxis;       // weight vectors that are closest to the axes
  int *niche;          // niche count
  double *xtrm;        // matrix of extreme point
  double *inv;         // Inverse matrix
  double *min;         // ideal point
  int *pi;             // closest reference point for each individual
  double *dist;        // distance to the reference point
  int *filter;         // selected individuals
  char *wfile;         // file of weight vectors
  double *cv;          // constraint value of individuals
} EMO_NSGA3;

void EMO_NSGA3_alloc(EMO_NSGA3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_NSGA3_free(EMO_NSGA3 *alg);
void EMO_NSGA3_run(EMO_NSGA3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);


typedef struct {
  EMO_NDSort nd;
  EMO_List lst1, lst2;    /* temporary lists */
  int *filter;            /* selected individuals */
  double *cd;             /* crowding distance */
  double *min;            /* ideal point */
  double *max;            /* nadir point */
  int *parent1, *parent2; /* candidates for reproduction */
  int *seedp;             /* enumeration of parent population */ 
  double **sort;          /* temporary array for sorting */
  int ssize;
} EMO_NSGA2;

void EMO_NSGA2_alloc(EMO_NSGA2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_NSGA2_free(EMO_NSGA2 *alg);
void EMO_NSGA2_run(EMO_NSGA2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);


typedef struct {
  EMO_List lst1, lst2; // temporary lists
  double *norm;        // normalized objective functions
  double *rank;        // rank of individuals
  double *l2;          // l2 norm
  double **sort;       // temporary array for sorting population
  int ssize;           // temporal population size (for memory free)
  double *tmp;         // temporary array, borrar RHG
  int wsize;           // number of weight vectors
  double *W;           // weight vectors
  double *min;         // ideal point
  double *max;         // nadir point
  double *ideal;       // ideal point
  double *new_min;
  double *new_max;
  double *hist;
  int *update;
  char *utl_name;
  EMO_Utility utl;
  int *filter;         // selected individuals
  int max_hist;        /* parameters of the algorithm */
  double epsilon;
  double alpha;
  char *wfile;
  int dm;
  int (*fcomp)(const void *, const void *);
  void (*fnorm)(double *, double *, int *, int, double *, double *, int);
} EMO_MOMBI2;

void EMO_MOMBI2_alloc(EMO_MOMBI2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_MOMBI2_free(EMO_MOMBI2 *alg);
void EMO_MOMBI2_run(EMO_MOMBI2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_MOMBI2_prun(EMO_MOMBI2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

typedef struct {
  EMO_List lst1, lst2; // temporary lists
  EMO_Refpoint ref;    // Reference points
  double *norm;        // normalized objective functions
  double *rank;        // rank of individuals
  double *rank2;       // rank of individuals
  double *l2;          // l2 norm
  double **sort;       // temporary array for sorting population
  int ssize;           // temporal population size (for memory free)
  double *tmp;         // temporary array, borrar RHG
  int wsize;           // number of weight vectors
  double *W;           // weight vectors
  char **H;            // hyper-heuristic
  char **T;            // target directions
  double *min;         // ideal point
  double *max;         // nadir point
  double *ideal;       // ideal point
  double *new_min;
  double *new_max;
  double *hist;
  int *update;
  EMO_Utility utl;
  int *filter0;        // selected individuals
  int *filter;         // selected individuals
  int *filter2;        // selected individuals
  int max_hist;        /* parameters of the algorithm */
  double epsilon;
  double alpha;
  char *wfile;
  int dm;
  int (*fcomp)(const void *, const void *);
  void (*fnorm)(double *, double *, int *, int, double *, double *, int);
} EMO_MOMBI3;


void EMO_MOMBI3_alloc(EMO_MOMBI3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_MOMBI3_free(EMO_MOMBI3 *alg);
void EMO_MOMBI3_run(EMO_MOMBI3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_MOMBI3_prun(EMO_MOMBI3 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

typedef struct EMO_IBEA {
  int indicator;
  double (*indf)(struct EMO_IBEA *, double *, double *, int); /* indicator-pointer function */
  double *indv;           /* matrix of indicator values */
  double *fit;            /* fitness */
  int *filter;            /* selected individuals */
  int *seedp;             /* enumeration of parent population */
  int *parent1, *parent2; /* candidates for reproduction */
  EMO_List lst1, lst2;    /* temporary lists */

  double kappa;           /* scaling factor */
  double rho;             /* scaling factor for the reference point (hv) */

  double *min;            /* ideal point (hv) */
  double *max;            /* nadir point (hv) */
  double *diff;           /* difference between max and min (eps, hv) */
  double *zref;           /* reference point (hv, r2) */
  double max_vol;         /* maximum volume (hv) */
  double *norm;           /* normalized objective functions (r2) */
  int wsize;              /* number of weight vectors (r2) */
  double *W;              /* weight vectors (r2) */
  EMO_Utility utl;        /* utility function (r2) */

  /* mpirun */
  double **sort;          /* temporary array for sorting population */
  int ssize;
  char *wfile;
} EMO_IBEA;

void EMO_IBEA_alloc(EMO_IBEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_IBEA_free(EMO_IBEA *alg);
void EMO_IBEA_run(EMO_IBEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_IBEA_prun(EMO_IBEA *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

typedef struct {
  EMO_NDSort nd;
  int *filter;        /* selected individuals */
  //double *fit;
  double *min;            /* ideal point */
  double *max;            /* nadir point */
  int *parent1, *parent2; /* candidates for reproduction */
  int *seedp;             /* enumeration of parent population */
  EMO_List lst1, lst2;    /* temporary lists */
 
  int samples;
  double bound;
  //int alpha;  /* number of individuals in initial population */ 

  double *tmp;            //arreglo temporal para copiar obj del ultimo frente
  double *rho;            //double rho[ pop->mu ];   >>hypeIndicator
  int *hitstat;           //int hitstat[ popsize ];  >>hypeSampling
  double *sample;         //double sample[ dim ];    >>hypeSampling
  double *val;            // double val[ pop->size];  >>Parameter of hypeIndicator()
  double *boundsVec;       //double boundsVec[ nobj ];  >>hypeExact
  int *indices;           //int indices[ popsize ];   >>hypeExact
  //double extrusion;       //double extrusion;         >>hypeExactRecursive
  int *pvec;              //int pvec[pnts];           >>hypeExactRecursive
  double *p;              //double p[pnts*dim];       >>hypeExactRecursive
  int *beg; 	          //int beg[MAX_LEVELS];      >>rearrangeIndicesByColumn
  int *end;               //int end[MAX_LEVELS];      >>rearrangeIndicesByColumn
  double *ref;             //double ref[rows];         >>rearrangeIndicesByColumn
} EMO_HYPE;

void EMO_HYPE_alloc(EMO_HYPE *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_HYPE_free(EMO_HYPE *alg);
void EMO_HYPE_run(EMO_HYPE *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

typedef struct {
  int niche;         /* neighborhood size */
  int wsize;         /* number of weight vectors = population size */
  double *W;         /* weight vectors */
  double *min;       /* ideal point */
  double *max;       /* nadir point */
  double *diff1;     /* temporary vector for subtraction */
  double *diff2;     /* temporary vector for subtraction */
  EMO_Utility utl;   /* Utility function: tchebycheff, pbi, etc. */
  EMO_List lp;       /* list for storing population members */
  EMO_List lt1;      /* temporary list */
  EMO_List lt2;      /* temporary list */
  EMO_List *lnn;     /* array of lists for storing neighbors */
  double **sort;     /* temporary array for sorting */
  int ssize;         /* temporary variable for the population size */
  char *wfile;       /* file containing the weight vectors */
  char *utl_name;    /* scalarizing or utility function */
  int norm;          /* flag that indicates if objectives should be normalized */
  int ver;           /* MOEA/D version */
  double delta;         /* probability that parent solutions are selected from the neighborhood */
  int nr;            /* maximal number of solutions replaced by each child solution */
  void (*fdiff)(double *, double *, double *, double *, int);
} EMO_MOEAD;


void EMO_MOEAD_alloc(EMO_MOEAD *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_MOEAD_free(EMO_MOEAD *alg);
void EMO_MOEAD_run(EMO_MOEAD *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

typedef struct {
  int k;                  /**/
  int niche;
  int wsize;              /* number of weight vectors = population size */
  double *W;              /* weight vectors */
  int *parent1, *parent2; /* candidates for reproduction */
  double *fit;            /* fitness */ 
  double *S;
  //double *R;
  //double *D;
  double *F;
  EMO_List *lnn;          /* array of lists for storing neighbors */
  EMO_List copy;
  double **sort;          /* temporary array for sorting */
  double **dist;
  int *filter;            /* selected individuals */
  int *seedp;             /* enumeration of parent population */
  EMO_List lst1, lst2;    /* temporary lists */ 
  int ssize;
} EMO_SPEA2;

void EMO_SPEA2_alloc(EMO_SPEA2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_SPEA2_free(EMO_SPEA2 *alg);
void EMO_SPEA2_run(EMO_SPEA2 *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

typedef struct {
  int *M;  /* matrix of value path */
  double *c;
  int x, y;
  int theta;
  int dim;
} EMO_VPath;

void EMO_VPath_alloc(EMO_VPath *vpath, int resolution, int mu, int max_mu, int dim);
void EMO_VPath_free(EMO_VPath *vpath);
void EMO_VPath_write(EMO_VPath *vpath, const char *str);
int EMO_VPath_run(EMO_VPath *vpath, double *data, int *filter, int size, double z);
void EMO_VPath_update(EMO_VPath *vpath, double *v, double z, int val);
void EMO_VPath_contribution(EMO_VPath *vpath, double *data, int *filter, int size, double z, EMO_List *lst);
int EMO_VPath_prune(EMO_VPath *vpath, EMO_Prune *p, double *data, int size, int new_size);

typedef struct {
  EMO_NDSort nd;
  EMO_VPath vpath;
  EMO_List lst;
  int xrs;          /* x-axis resolution of vpath */
  int *filter;      /* auxiliar for enabling/disabling individuals */
  double *tmp;      /* temporary array for storing an individual */
  double *min;      /* Reference points */
  double *max0;
  double *ideal;
  double *nadir;
  double **sort;    /* temporary array for sorting population */
  double *norm;     /* normalized objective functions */
  double *vnorm;    /* norm value for each individual */
  int ssize;        /* temporal population size (for memory free) */
} EMO_MOVAP;

void EMO_MOVAP_alloc(EMO_MOVAP *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);
void EMO_MOVAP_free(EMO_MOVAP *alg);
void EMO_MOVAP_run(EMO_MOVAP *alg, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

/* Al final, Generic MOEA */
typedef void (*EMO_MOEA_falloc)(void *, EMO_Param *, EMO_Population *, EMO_MOP *);
typedef void (*EMO_MOEA_ffree)(void *);
typedef void (*EMO_MOEA_frun)(void *, EMO_Param *, EMO_Population *, EMO_MOP *);

typedef struct {
  EMO_MOEA_falloc alloc;
  EMO_MOEA_ffree free;
  EMO_MOEA_frun run;

  union {
    EMO_SMSEMOA smsemoa;
    EMO_NSGA2 nsga2;
    EMO_NSGA3 nsga3;
    EMO_MOMBI2 mombi2;
    EMO_MOMBI3 mombi3;
    EMO_IBEA ibea;
    EMO_MOEAD moead;
    EMO_SPEA2 spea2;
    EMO_HYPE hype;
    EMO_MOVAP movap;
  } alg;

  void *palg;
} EMO_MOEA;

extern const char *EMO_MOEA_list[];

void EMO_MOEA_alloc(EMO_MOEA *moea, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop, const char *str);
void EMO_MOEA_free(EMO_MOEA *moea);
void EMO_MOEA_run(EMO_MOEA *moea, EMO_Param *param, EMO_Population *pop, EMO_MOP *mop);

#endif

