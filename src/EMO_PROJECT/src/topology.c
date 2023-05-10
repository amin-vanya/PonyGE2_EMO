/* Defines a network topology */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "topology.h"
#include "common.h"
#include "io.h"

void EMO_Topology_get_type(int *topo, const char *str) {
  const char *name [] = { "line",
                          "ring",
                          "star",
                          "tree",
                          "full",
                          "torus",
                          "torusd",
                          "mesh",
                          NULL };

   const int id []  = {  EMO_TOPOLOGY_LINE,
                         EMO_TOPOLOGY_RING,
                         EMO_TOPOLOGY_STAR,
                         EMO_TOPOLOGY_TREE,
                         EMO_TOPOLOGY_FULL,
                         EMO_TOPOLOGY_TORUS,
                         EMO_TOPOLOGY_TORUSD,
                         EMO_TOPOLOGY_MESH,
                       };

  char *topology;
  int i;

  /* Find the function's name in dictionary */
  topology = EMO_tolower(str);
  i = EMO_Dictionary_find(name, topology);

  if(i == -1) {
    printf("Error, unknown topology %s in EMO_Topology_get_type.\n", str);
    free(topology);
    exit(1);
  }

  *topo = id[i];
  free(topology);
}

void EMO_Topology_get_flow(int *flow, const char *str) {
  const char *name [] = { "uni",
                          "bidi",
                          "scatter",
                          "gather",
                          NULL };

  const int id []  = { EMO_TOPOLOGY_UNI,
                       EMO_TOPOLOGY_BIDI,
                       EMO_TOPOLOGY_SCATTER,
                       EMO_TOPOLOGY_GATHER
                      };

  char *comm;
  int i;

  /* Find the function's name in dictionary */
  comm = EMO_tolower(str);
  i = EMO_Dictionary_find(name, comm);

  if(i == -1) {
    printf("Error, unknown communication %s in EMO_Topology_get_flow.\n", str);
    free(comm);
    exit(1);
  }

  *flow = id[i];
  free(comm);
}

/* Line topology 
 *
 *  0 ---> 1 ---> 2 ---> ... ---> size - 1  (unidirectional)
 *  0 <--> 1 <--> 2 <--> ... <--> size - 1  (bidirectional)
 */
void EMO_Topology_line(EMO_List *in, EMO_List *out, int rank, int size, int type) {
  int n;

  if(type == EMO_TOPOLOGY_UNI)
    n = (rank % size - 1)? 0 : 1;
  else
    n = (rank == 0 || rank == size - 1)? 1 : 2;

  EMO_List_alloc(in, n); 
  EMO_List_alloc(out, n); 

  if(rank > 0) {
    EMO_List_queue(in, rank - 1);

    if(type != EMO_TOPOLOGY_UNI)
      EMO_List_queue(out, rank - 1);
  }

  if(rank < size - 1) {
    EMO_List_queue(out, rank + 1);

    if(type != EMO_TOPOLOGY_UNI)
      EMO_List_queue(in, rank + 1);
  }
}

/* Ring topology
 *
 *    0 --->  1  ---> 2
 *    ^               |   (unidirectional)
 *    |               |
 *    |               v 
 *      <--- ... <--- 3
 *  size - 1 
 *
 *
 *    0 <-->  1  <--> 2
 *    ^               ^ 
 *    |               |   (bidirectional)
 *    v               v 
 *      <--> ... <--> 3
 *  size - 1 
 */
void EMO_Topology_ring(EMO_List *in, EMO_List *out, int rank, int size, int type) {
  int n;

  n = (type == EMO_TOPOLOGY_UNI)? 1 : 2;

  EMO_List_alloc(in, n); 
  EMO_List_alloc(out, n); 

  EMO_List_queue(in, (size + rank - 1) % size);
  EMO_List_queue(out, (rank + 1) % size);

  if(type != EMO_TOPOLOGY_UNI) {
    EMO_List_queue(out, (size + rank - 1) % size);
    EMO_List_queue(in, (rank + 1) % size);
  }
}

/* Star topology
 *
 *              1
 *              ^ 
 *              |
 * size - 1 <-- 0 --> 2  (scatter, unidirectional)
 *              |
 *              v 
 *             ...
 *
 *
 *              1
 *              |
 *              v
 * size - 1 --> 0 <-- 2  (gather)
 *              ^
 *              | 
 *             ...
 *
 *
 *               1
 *               ^ 
 *               |
 *               v
 * size - 1 <--> 0 <--> 2  (bidirectional)
 *               ^
 *               |
 *               v 
 *              ...
 */
void EMO_Topology_star(EMO_List *in, EMO_List *out, int rank, int size, int type) {
  int nin, nout, i;

  switch(type) {
    case EMO_TOPOLOGY_SCATTER:
         nin  = (rank == 0)? 0 : 1;
         nout = (rank == 0)? size - 1 : 0;
         break;

    case EMO_TOPOLOGY_GATHER:
         nin  = (rank == 0)? size  - 1: 0;
         nout = (rank == 0)? 0 : 1;
         break;

    case EMO_TOPOLOGY_BIDI:
         nin = nout = (rank == 0)? size  - 1: 1;
         break;

    default:
         printf("Error, unknowkn topology %d\n", type);
         exit(1);
  }

  EMO_List_alloc(in, nin); 
  EMO_List_alloc(out, nout); 

  if(rank == 0) {
    for(i = 1; i < size; i++) {
      if(type == EMO_TOPOLOGY_GATHER || type == EMO_TOPOLOGY_BIDI)
        EMO_List_queue(in, i);

      if(type == EMO_TOPOLOGY_SCATTER || type == EMO_TOPOLOGY_BIDI)
        EMO_List_queue(out, i);
    }
  }
  else {
    if(type == EMO_TOPOLOGY_SCATTER || type == EMO_TOPOLOGY_BIDI)
      EMO_List_queue(in, 0);

    if(type == EMO_TOPOLOGY_GATHER || type == EMO_TOPOLOGY_BIDI)
      EMO_List_queue(out, 0);
  }
}

/* Tree topology
 *
 * degree: number of children of a node.
 *
 *  size = 15
 *  degree = 2
 *                 0
 *          1              2
 *       3     4       5       6
 *      7 8   9 10   11 12   13 14
 * 
 *     0: 1,2    4: 9,10
 *     1: 3,4    2: 5,6
 *     2: 5,6    5: 11,12
 *     3: 7,8    6: 13,14
 */
void EMO_Topology_tree(EMO_List *in, EMO_List *out, int rank, int size, int type, int degree) {
  int cmin, cmax, nchild, parent, nin, nout, i;

  /* children interval */
  cmin = degree * rank + 1;
  cmax = degree * (rank + 1);

  if(cmax > size - 1) /* incomplete children */
   cmax = size - 1;

  /* parent node (root has no parent) */
  parent = (rank == 0)? -1 : (int) (rank / degree);

  if(rank % degree == 0)
    parent -= 1;

  nchild = (cmin > size - 1)? 0 : cmax - cmin + 1;

  switch(type) {
    case EMO_TOPOLOGY_SCATTER:
         nin  = (rank == 0)? 0 : 1;
         nout = nchild;
         break;

    case EMO_TOPOLOGY_GATHER:
         nin  = nchild;
         nout = (rank == 0)? 0 : 1;
         break;

    case EMO_TOPOLOGY_BIDI:
         nin = nout = (rank == 0)? nchild: nchild + 1;
         break;

    default:
         printf("Error, unknowkn topology %d\n", type);
         exit(1);
  }

  EMO_List_alloc(in, nin); 
  EMO_List_alloc(out, nout); 

  if(rank == 0) {
    for(i = cmin; i <= cmax; i++) {
      if(type == EMO_TOPOLOGY_GATHER || type == EMO_TOPOLOGY_BIDI)
        EMO_List_queue(in, i);

      if(type == EMO_TOPOLOGY_SCATTER || type == EMO_TOPOLOGY_BIDI)
        EMO_List_queue(out, i);
    }
  }
  else {
    if(type == EMO_TOPOLOGY_SCATTER || type == EMO_TOPOLOGY_BIDI) {
      EMO_List_queue(in, parent);

      for(i = cmin; i <= cmax; i++)
        EMO_List_queue(out, i);
     }

    if(type == EMO_TOPOLOGY_GATHER || type == EMO_TOPOLOGY_BIDI) {
      for(i = cmin; i <= cmax; i++)
        EMO_List_queue(in, i);

      EMO_List_queue(out, parent);
    }
  }
}

/* Fully connected mesh, by default it uses bidirectional communication
 *
 *  0 - 1
 *  | \ |
 *  2 - 3
 */
void EMO_Topology_full(EMO_List *in, EMO_List *out, int rank, int size) {
  int i;

  EMO_List_alloc(in, size - 1); 
  EMO_List_alloc(out, size - 1); 

  for(i = 0; i < size; i++) {
    if(rank != i) {
      EMO_List_queue(in, i);
      EMO_List_queue(out, i);
    }
  }
}

/* Torus topology
 *
 *  3 - 0 - 1 - 2 - 3 - 0 ...
 *  |   |   |   |   |
 *  7 - 4 - 5 - 6 - 7 - 4 ...
 *  |   |   |   |   |
 * 11 - 8 - 9 -10 -11 - 8 ...
 *  |   |   |   |   |
 *  3 - 0 - 1 - 2 - 3
 *  |   |   |   |   |
 *  ...
 *
 *  Neighbors
 *
 *      d
 *      |
 *   b- w -a    (unidirectional: in = {b, d}, out = {a, c})
 *      |
 *      c
 */ 
void EMO_Topology_torus(EMO_List *in, EMO_List *out, int rank, int row, int col, int type) {
  int x, y, a, b, c, d;

  if(row <= 0 || col <= 0) {
    printf("Error, invalid dimensions row %d, col %d\n", row, col);
    exit(1);
  }

  /* Get x and y dimensions */
  x = rank % col;
  y = (int) rank / col;

  EMO_List_alloc(in, 4);
  EMO_List_alloc(out, 4);

  /* rank at position (x,y) -> col * y + x */
  a = col * y + (x + 1) % col;
  b = col * y + (col + x - 1) % col;
  c = col * ((y + 1) % row) + x;
  d = col * ((row + y - 1) % row) + x;

  if(!EMO_List_seek(in, b)) EMO_List_queue(in, b);
  if(!EMO_List_seek(in, d)) EMO_List_queue(in, d);
  if(!EMO_List_seek(out,a)) EMO_List_queue(out,a);
  if(!EMO_List_seek(out,c)) EMO_List_queue(out,c);

  if(type != EMO_TOPOLOGY_UNI) {
    if(!EMO_List_seek(in, a)) EMO_List_queue(in, a);
    if(!EMO_List_seek(in, c)) EMO_List_queue(in, c);
    if(!EMO_List_seek(out,b)) EMO_List_queue(out,b);
    if(!EMO_List_seek(out,d)) EMO_List_queue(out,d);
  }
}

/* Diagonal torus topology
 *
 *  3 - 0 - 1 - 2 - 3 - 0 ...
 *  | X | X | X | X | X
 *  7 - 4 - 5 - 6 - 7 - 4 ...
 *  | X | X | X | X | X
 * 11 - 8 - 9 -10 -11 - 8 ...
 *  | X | X | X | X | X
 *  3 - 0 - 1 - 2 - 3  ...
 *  | X | X | X | X |
 *  ...
 *
 *  Neighbors
 *
 *   e  d  f
 *    \ | /
 *   b- w -a    (unidirectional: in = {b, d, e, f}, out = {a, c, g, h})
 *    / | \
 *   g  c  h
 */ 
void EMO_Topology_torus_diagonal(EMO_List *in, EMO_List *out, int rank, int row, int col, int type) {
  int x, y, e, f, g, h;

  EMO_Topology_torus(in, out, rank, row, col, type);

  /* Get x and y dimensions */
  x = rank % col;
  y = (int) rank / col;
 
  /* rank at position (x,y) -> col * y + x */
  /* Neighbors along diagonals */ 
  e = col * ((row + y - 1) % row) + (col + x - 1) % col;
  f = col * ((row + y - 1) % row) + (x + 1) % col;
  g = col * ((y + 1) % row) + (col + x - 1) % col;
  h = col * ((y + 1) % row) + (x + 1) % col;

  if(!EMO_List_seek(in, e)) EMO_List_queue(in, e);
  if(!EMO_List_seek(in, f)) EMO_List_queue(in, f);
  if(!EMO_List_seek(out,g)) EMO_List_queue(out,g);
  if(!EMO_List_seek(out,h)) EMO_List_queue(out,h);

  if(type != EMO_TOPOLOGY_UNI) {
    if(!EMO_List_seek(in, g)) EMO_List_queue(in, g);
    if(!EMO_List_seek(in, h)) EMO_List_queue(in, h);
    if(!EMO_List_seek(out,e)) EMO_List_queue(out,e);
    if(!EMO_List_seek(out,f)) EMO_List_queue(out,f);
  }
}

/* Mesh topology
 *
 * The file must have a square matrix of zeros and ones that represents a graph.
 * Given a rank, its row gives the output nodes and its column gives the input nodes.
 * For simplicity the diagonals are omitted.
 *
 * Example:  file.dat:
 *    # 4 4
 *    0 1 0 1
 *    0 0 1 0
 *    0 0 0 1
 *    0 1 0 0
 *
 *    0 --> 1
 *    |   7 |
 *    v  /  v
 *    3 <-- 2 
 */
void EMO_Topology_mesh(EMO_List *in, EMO_List *out, int rank, int size, const char *file) {
  double *data, v;
  int row, col, j;

  row = col = 0;
  data = EMO_File_read(NULL, &row, &col, file, 0);

  if(row != col) {
    printf("Error, dimensions of file %s must be the same (row %d, col %d) in EMO_Topology_mesh\n", file, row, col);
    exit(1);
  }

  if(row != size) {
    printf("Error, dimensions of file %s must agree with size %d in EMO_Topology_mesh\n", file, size);
    exit(1);
  }

  EMO_List_alloc(in, size);
  EMO_List_alloc(out, size);

  for(j = 0; j < col; j++) {
    v = data[rank * col + j];

    if(v != 0 && v != 1) {
      printf("Error, elements of file %s must be {0,1} at (%d,%d) = %.0f in EMO_Topology_mesh\n", file, rank, j, v);
      exit(1);
    }

    if(rank != j && v) EMO_List_queue(out, j);
  }

  for(j = 0; j < row; j++) {
    v = data[j * col + rank];

    if(v != 0 && v != 1) {
      printf("Error, elements of file %s must be {0,1} at (%d,%d) = %.0f in EMO_Topology_mesh\n", file, j, rank, v);
      exit(1);
    }

    if(rank != j && v) EMO_List_queue(in, j);
  }
}

