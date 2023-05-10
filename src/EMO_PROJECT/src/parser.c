/**************************************************************
 * parameter.c  Process the parameters of a configuration     *
 *              file and makes them available for a program.  * 
 *                                                            *
 * Computer Science Department, CINVESTAV-IPN                 *
 *                                                            *
 * Professor:   Dr. Carlos A. Coello Coello                   *
 * Author:      Raquel Hernandez Gomez                        *
 *                                                            *
 * March 2013                                                 *
 *************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "parser.h"
#include "vector.h"

#define EMO_CSIZE  1000  /* Longitud de cadenas */
#define EMO_COMMENT 1    /* Estados de lectura posibles en archivo de configuracion */
#define EMO_PARAM   2
#define EMO_VALUE   3
#define EMO_NEWLN   4

/* Contabiliza parametros de un archivo */
int _EMO_Parser_count(FILE *fp) {
  int c, flag, i = 0;

  flag = EMO_NEWLN;
 
  while((c = fgetc(fp)) != EOF) {
    if(c == '\n') {
      flag = EMO_NEWLN;
      continue;
    }

    if(c == ' ' || c == '\t' || c == '\r') /* Ignora espacios */
      continue;

    if(flag == EMO_NEWLN) {
      flag = 0;

      if(c != '#') 
        i++;
    }
  }

  rewind(fp);
  return i;
}

int _EMO_hash(const char *str, int n) {
  long h = 0;
  int i = 0;

  // key is the sum of ascii codes of str
  while(str[i] != '\0') {
    h += str[i];
    i++;
  }

  return (int) (h % n);
}

int _EMO_rhash(int h, int j, int n) {
  return (h + (int) sqrt((double) j)) % n;
}

/* Procesa archivo de configuracion */
void _EMO_Parser_run(EMO_Parser *parser, FILE *fp) {
  int j, c, flag;
  char *name, *value;

  name = parser->tmp1;
  value = parser->tmp2;

  flag = EMO_PARAM;
  j = 0;
 
  while((c = fgetc(fp)) != EOF) {

    if(c == '\n') {  /* Nueva linea */

      if(flag == EMO_VALUE) {  /* Termina de leer el valor de un parametro */
        value[j] = '\0';
        EMO_Parser_set(parser, name, value);
      }

      flag = EMO_PARAM;
      j = 0;
      continue;
    }

    if(c == ' ' || c == '\t' || c == '\r') /* Ignora espacios */
      continue;

    if(flag == EMO_COMMENT)  /* Ignora comentarios */
     continue;

    if(c == '#') {  /* Comentarios inician con simbolo # */

      if(flag == EMO_VALUE) {  //param = value #comentario
        value[j] = '\0';
        EMO_Parser_set(parser, name, value);
      }

      flag = EMO_COMMENT;
      continue;
    }

    if(c == '=') {  /* Finaliza lectura del nombre de parametro */
      name[j] = '\0';
      flag = EMO_VALUE;
      j = 0;
      continue;
    }

    /* Establece caracter por caracter nombres de parametros y valores */
    if(flag == EMO_PARAM) 
      name[j++] = c;

    if(flag == EMO_VALUE)
      value[j++] = c;
  }
}

int EMO_Parser_word_count(const char *s, const char delim) {
  int c, i;

  c = i = 0;

  if(s[0] == '\0') 
    return 0;

  while(s[i] != '\0') {
    if(s[i] == delim)
      c++;

    i++;
  }

  return c + 1;
}

int EMO_Parser_get_token(char *dest, char **src, const char delim) {
  char *p;
  int n;

  while(**src == delim)
    (*src)++;

  p = *src;

  while(*p != '\0' && *p != delim) 
    p++;

  n = p - *src;

  if(n > 0)
    strncpy(dest, *src, n);

  dest[n] = '\0';
  *src = p;

  return n;
}

void EMO_Parser_alloc(EMO_Parser *parser, EMO_Debug *dbg, int size) {
  int i;

  parser->dbg = dbg;
  parser->current = 0;
  parser->size = 3 * size;  // size of the hash table

  if(parser->size % 2 == 0) parser->size++;

  if((parser->name = (char **) malloc(sizeof(char *) * parser->size)) == NULL) {
    printf("Error, not enough memory in EMO_Parser 1.\n");
    exit(1);
  }

  if((parser->value = (char **) malloc(sizeof(char *) * parser->size)) == NULL) {
    printf("Error, not enough memory in EMO_Parser 2.\n");
    exit(1);
  }

  if((parser->tmp1 = (char *) malloc(sizeof(char) * EMO_CSIZE)) == NULL) {
    printf("Error, not enough memory in EMO_Parser 3.\n");
    exit(1);
  }

  if((parser->tmp2 = (char *) malloc(sizeof(char) * EMO_CSIZE)) == NULL) {
    printf("Error, not enough memory in EMO_Parser 4.\n");
    exit(1);
  }

  for(i = parser->size-1; i > -1; i--) {
    if((parser->name[i] = (char *) malloc(sizeof(char) * EMO_CSIZE)) == NULL) {
      printf("Error, not enough memory in EMO_Parser 5.\n");
      exit(1);
    }

    if((parser->value[i] = (char *) malloc(sizeof(char) * EMO_CSIZE)) == NULL) {
      printf("Error, not enough memory in EMO_Parser 6.\n");
      exit(1);
    }

    strcpy(parser->name[i], "");
    strcpy(parser->value[i], "");
  }
}

/* Reserva memoria 
 * file: Nombre del archivo de configuracion */
void EMO_Parser_alloc_from_file(EMO_Parser *parser, EMO_Debug *dbg, char *file) {
  FILE * fp;

  parser->dbg = dbg;
  EMO_Debug_printf(parser->dbg, "Reading parameter file %s", file);

  if((fp = fopen(file, "r")) == NULL) {
    printf("Error to open the configuration file %s.\n", file);
    exit(1);
  }

  EMO_Parser_alloc(parser, dbg, _EMO_Parser_count(fp));
  _EMO_Parser_run(parser, fp);
  fclose(fp);
}

/* Libera memoria */
void EMO_Parser_free(EMO_Parser *parser) {
  int i;

  for(i = parser->size-1; i > -1; i--) {
    free(parser->name[i]);
    free(parser->value[i]);
  }
  free(parser->name);
  free(parser->value);
  free(parser->tmp1);
  free(parser->tmp2);
}

void EMO_Parser_set(EMO_Parser *parser, const char *name, const char *value) {
  int i, j;

  if(parser->current >= parser->size) {
    printf("Error, not enough space in Parser (%d vs %d).\n", parser->current, parser->size);
    exit(1);
  }

  i = _EMO_hash(name, parser->size);
  j = 0;

  while(strlen(parser->name[i]) != 0) {
    if(strcmp(parser->name[i], name) == 0) {
      EMO_Debug_printf(parser->dbg, "Warning, parameter %s was updated (%s vs %s)", name, parser->value[i], value);
      break;
    }
    i = _EMO_rhash(i, j++, parser->size);
  }
  
  strcpy(parser->name[i], name);
  strcpy(parser->value[i], value);
  parser->current++;
}

/* Identifica la posicion que ocupa un parametro en arreglo */
int EMO_Parser_find(EMO_Parser *parser, const char *s) {
  int i, j;

  i = _EMO_hash(s, parser->size);
  j = 0;

  while(strlen(parser->name[i]) != 0) {
    if(strcmp(parser->name[i], s) == 0)
      return i;

    i = _EMO_rhash(i, j++, parser->size);
  }

  return -1; 
}

void EMO_Parser_print(EMO_Parser *parser) {
  int i;

  for(i = 0; i < parser->size; i++)
    printf("EMO_Parser [%d]: %s = %s\n", i, parser->name[i], parser->value[i]);
}

#undef EMO_CSIZE
#undef EMO_COMMENT
#undef EMO_PARAM
#undef EMO_VALUE
#undef EMO_NEWLN

