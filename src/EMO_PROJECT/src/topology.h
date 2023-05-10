
#ifndef _TOPOLOGY_H
#define _TOPOLOGY_H

#include "list.h"

#define EMO_TOPOLOGY_UNI     0 /* Unidirectional communication  */
#define EMO_TOPOLOGY_BIDI    1 /* Bidirectional communication  */

#define EMO_TOPOLOGY_SCATTER 0  /* Star and tree topologies */
#define EMO_TOPOLOGY_GATHER  2

#define EMO_TOPOLOGY_LINE   0
#define EMO_TOPOLOGY_RING   1
#define EMO_TOPOLOGY_STAR   2
#define EMO_TOPOLOGY_TREE   3
#define EMO_TOPOLOGY_FULL   4
#define EMO_TOPOLOGY_TORUS  5
#define EMO_TOPOLOGY_TORUSD 6
#define EMO_TOPOLOGY_MESH   7

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

#endif

