/**************************************************************
 * list.c       Implementation of doubly linked list.         *
 *                                                            *
 * Computer Science Department, CINVESTAV-IPN                 *
 *                                                            *
 * Professor:   Dr. Carlos A. Coello Coello                   *
 * Author:      Raquel Hernandez Gomez                        *
 *                                                            *
 * September 2014                                             *
 *************************************************************/

#include <stdlib.h>
#include "list.h"

/* Create a new node */
_EMO_Node *_EMO_Node_alloc(int item) {
  _EMO_Node *n;

  if((n = (_EMO_Node *) malloc(sizeof(_EMO_Node))) == NULL) {
    printf("Error, not enough memory in list.c:_EMO_Node_alloc.\n");
    exit(1);
  }
  n->index = -1;
  n->item = item;
  n->ant = n->sig = NULL;
  n->inblock = -1;
  return n;
}

/* Get a node from the block of memory.
   If there is not available memory, it creates a new one */
_EMO_Node *_EMO_Node_get(EMO_List *l, int item) {
  _EMO_Node *n;
  int i;

  if(l->available == -1)
    return _EMO_Node_alloc(item);

  i = l->available;
  n = l->memblock[i];
  n->item = item;
 
  while(1) {
    i = (i + 1) % l->max_size;

    if(i == l->available) {
      l->available = -1;
      break;
    }

    if(l->memblock[i]->index == -1) {
      l->available = i;
      break;
    }
  }
  return n;
}

/* Free a node */
void _EMO_Node_free(EMO_List *l, _EMO_Node *n) {

  l->size--; 

  if(n->inblock > -1) {
    n->index = -1;
    n->item = 0;
    n->ant = n->sig = NULL;

    // Updates 'available'
    if(l->available == -1)
      l->available = n->inblock;
  }
  else
    free(n);
}

/* Look for the node which contains the given item */ 
_EMO_Node *_EMO_Node_find(EMO_List *l, int item) {
  _EMO_Node *n;

  n = l->ini;
    
  while(n != NULL) {
    if(n->item == item)
      break;
    n = n->sig;
  }
  return n;
}

/* Look for the node at the given position */ 
_EMO_Node *_EMO_Node_retrieve(EMO_List *l, int pos) {
  _EMO_Node *n;

  if(pos < 0 || pos >= l->size)
    return NULL;

  if(l->size - pos - 1 >= pos) {
    n = l->ini;
    
    while(n != NULL) {
      if(n->index == pos)
        break;
      n = n->sig;
    }
  }
  else {
    n = l->fin;
  
    while(n != NULL) {
      if(n->index == pos)
        break;
      n = n->ant;
    }
  }
  return n;
}

/* Initialize the structure of a list and reserves
   the memory from which the nodes will be allocated. */
void EMO_List_alloc(EMO_List *l, int size) {
  _EMO_Node *n;
  int i;

  l->ini = l->fin = NULL;
  l->size = 0;
  l->max_size = size;

  if(size == 0) 
    l->available = -1;
  else {
    l->memblock = (_EMO_Node **) malloc(sizeof(_EMO_Node *) * size);

    if(l->memblock == NULL) {
      printf("Error, not enough memory in list.c:initList.\n");
      exit(1);
    }
 
    for(i = 0; i < size; i++) {
      n = _EMO_Node_alloc(0);
      n->inblock = i;
      l->memblock[i] = n;
    }

    l->available = 0;
  }
}

/* Free the list */
void EMO_List_free(EMO_List *l) {
  int i;

  EMO_List_clear(l);

  if(l->max_size > 0) {
    for(i = 0; i < l->max_size; i++) {
      if(l->memblock[i]->index == -1) {
        free(l->memblock[i]);
      }
    }
    free(l->memblock);
  }
}

/* Clean the list */
void EMO_List_clear(EMO_List *l) {
  int i, n;

  n = l->size;

  for(i = 0; i < n; i++)  
    EMO_List_dequeue(l, NULL);

  l->ini = l->fin = NULL;
}

/* Add an item at the end of the list */
void EMO_List_queue(EMO_List *l, int item) {
  _EMO_Node *n = _EMO_Node_get(l, item);

  if (l->size == 0) {
    l->ini = l->fin = n;
  }
  else {
    n->ant = l->fin;
    l->fin->sig = n;
    l->fin = n;
  }
  n->index = l->size++; 
}

/* Delete the first node from the list */
void EMO_List_dequeue(EMO_List *l, int *item) {
  _EMO_Node *tmp;

  if(l->size == 0) 
    return;

  if(l->ini->sig != NULL)
    l->ini->sig->ant = NULL;

  tmp = l->ini;

  if(item != NULL) *item = tmp->item;

  l->ini = l->ini->sig; 
  _EMO_Node_free(l, tmp);
}

/* Insert an item at the given position */
void EMO_List_add(EMO_List *l, int item, int pos) {
  _EMO_Node *n, *ap;

  if(pos < 0 || pos > l->size) {
    printf("Error, invalid position in list.c: EMO_List_add\n");
    exit(1);
  }

  if(pos == l->size || pos == l->size - 1) {
    EMO_List_queue(l, item); 
    return;
  }
 
  n = _EMO_Node_get(l, item);
 
  if (l->size == 0 && pos == 0) {
    l->ini = l->fin = n;
    n->index = l->size;
  }
  else { 
    ap = _EMO_Node_retrieve(l, pos);
 
    n->ant = ap->ant;
    n->sig = ap;
 
    if(ap->ant != NULL)
      ap->ant->sig = n;
 
    ap->ant = n;

    n->index = ap->index;

    while(ap != NULL) {
      ap->index++;
      ap = ap->sig;
    }

    if(pos == 0)
      l->ini = n;

    if(pos == l->size)
      l->fin = n;
  }
  l->size++;
}

/* Delete the first occurrence that matches with the given item from a list */
int EMO_List_remove(EMO_List *l, int item) {
  _EMO_Node *ap, *tmp;

  if(l->size == 0)
    return 0;

  ap = _EMO_Node_find(l, item);

  if(ap == NULL)
    return 0;

  if(ap->ant != NULL)
    ap->ant->sig = ap->sig;

  if(ap->sig != NULL)
    ap->sig->ant = ap->ant;

  if(ap == l->ini)
    l->ini = ap->sig;

  if(ap == l->fin)
    l->fin = ap->ant;

  tmp = ap;
  ap = ap->sig;

  while(ap != NULL) {
    ap->index--;
    ap = ap->sig;
  }

  _EMO_Node_free(l, tmp);
  return 1;
}

/* Delete the node from the given position and stores its value in item */
int EMO_List_retrieve(EMO_List *l, int *item, int pos) {
  _EMO_Node *ap, *tmp;

  if(l->size == 0)
    return 0;

  ap = _EMO_Node_retrieve(l, pos);

  if(ap == NULL)
    return 0;

  if(item != NULL)
    *item = ap->item;

  if(ap->ant != NULL)
    ap->ant->sig = ap->sig;

  if(ap->sig != NULL)
    ap->sig->ant = ap->ant;

  if(ap == l->ini)
    l->ini = ap->sig;

  if(ap == l->fin)
    l->fin = ap->ant;

  tmp = ap;
  ap = ap->sig;

  while(ap != NULL) {
    ap->index--;
    ap = ap->sig;
  }

  _EMO_Node_free(l, tmp);
  return 1;
}

/* Get the item at the given position */
int EMO_List_get(EMO_List *l, int *item, int pos) {
  _EMO_Node *n = _EMO_Node_retrieve(l, pos);

  if(n != NULL && item != NULL) {
    *item = n->item; 
    return 1;
  }

  printf("Warning, item at %d position not found in list.c: EMO_List_get\n", pos);
  return 0;
}

/* Look for the node which contains a given item */
int EMO_List_seek(EMO_List *l, int item) {
  _EMO_Node *n;

  n = l->ini;

  while(n != NULL) {
    if(n->item == item)
      return 1;
    n = n->sig;
  }

  return 0;
}

/* Count how many times item is in the list */
int EMO_List_count(EMO_List *l, int item) {
  int c = 0;
  _EMO_Node *n;

  n = l->ini;

  while(n != NULL) {
    if(n->item == item)
      c++;
    n = n->sig;
  }
  return c;
}

// Add the elements of src list to the dest list
void EMO_List_append(EMO_List *dest, EMO_List *src) {
  int i, n, elem = 0;

  n = src->size;

  for(i = 0; i < n; i++) {
    EMO_List_get(src, &elem, i);
    EMO_List_queue(dest, elem);
  }
}

/* Removes an item from a list and adds this item to another list */
int EMO_List_move(EMO_List *from, EMO_List *to, int item) {
  if(EMO_List_remove(from, item) == 0)
    return 0;

  EMO_List_queue(to, item);
  return 1;
}

/* Removes all elements from a list and adds them to another list */
int EMO_List_move_all(EMO_List *from, EMO_List *to) {
  int item;

  while(from->size > 0) {
    if(EMO_List_retrieve(from, &item, 0) == 0)
      return 0;
    EMO_List_queue(to, item);
  }
  return 1;
}


/* Print the list */
void EMO_List_print(EMO_List *l, FILE *fp, const char *s) {
  _EMO_Node *n;

  n = l->ini;

  if(fp == NULL)
    fp = stdout;

  if(s != NULL && s[0] != '\n')
    fprintf(fp, "%s: ", s);

  while(n != NULL) {
    fprintf(fp, "[%d]=%d ", n->index, n->item);
    n = n->sig;
  }

  fprintf(fp, " size=%d\n", l->size);
}

