
/**************************************************************
 * parameter.h  Function prototypes for processing            *
 *              configuration files.                          *
 *                                                            *
 * Computer Science Department, CINVESTAV-IPN                 *
 *                                                            *
 * Professor:   Dr. Carlos A. Coello Coello                   *
 * Author:      Raquel Hernandez Gomez                        *
 *                                                            *
 * March 2013                                                 *
 *************************************************************/

#ifndef _PARSER_H
#define _PARSER_H

#include <stdio.h>
#include "debug.h"

typedef struct {
  char **name;        /* Parameter names */
  char **value;       /* Parameter value */
  char *tmp1, *tmp2;  /* Temporary variables */
  int size;           /* Number of parameters */
  int current;        /* Current parameter (only for set function) */
  EMO_Debug *dbg;
} EMO_Parser;

/* Non-public functions */
int  _EMO_Parser_count(FILE *fp);
int  _EMO_hash(const char *str, int n);
int  _EMO_rhash(int h, int j, int n);
void _EMO_Parse_run(EMO_Parser *parser, FILE *fp);

int  EMO_Parser_word_count(const char *s, const char delim);
int  EMO_Parser_get_token(char *dest, char **src, const char delim);
void EMO_Parser_alloc(EMO_Parser *parser, EMO_Debug *dbg, int size);
void EMO_Parser_alloc_from_file(EMO_Parser *parser, EMO_Debug *dbg, char *file);
void EMO_Parser_free(EMO_Parser *parser);
void EMO_Parser_set(EMO_Parser *parser, const char *name, const char *value);
int  EMO_Parser_find(EMO_Parser *parser, const char *s);
void EMO_Parser_print(EMO_Parser *parser);

#endif

