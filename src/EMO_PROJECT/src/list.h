/**************************************************************
 * list.h       Implementation of doublly linked list         *
 *                                                            *
 * Computer Science Department, CINVESTAV-IPN                 *
 *                                                            *
 * Professor:   Dr. Carlos A. Coello Coello                   *
 * Author:      Raquel Hernandez Gomez                        *
 *                                                            *
 * September 2013                                             *
 *************************************************************/

#ifndef _LIST_H
#define _LIST_H

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

/* Non-public functions */
_EMO_Node * _EMO_Node_alloc(int item);
_EMO_Node * _EMO_Node_get(EMO_List *l, int item); 
void _EMO_Node_free(EMO_List *l, _EMO_Node *n);
//void        _EMO_Node_free(_EMO_Node *n); 
_EMO_Node * _EMO_Node_find(EMO_List *l, int item);
_EMO_Node * _EMO_Node_retrieve(EMO_List *l, int pos);

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

#endif

