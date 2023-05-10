/**************************************************************
 * io.h   Toolkit functions for input/output.                 *
 *                                                            *
 * Computer Science Department, CINVESTAV-IPN                 *
 *                                                            *
 * Professor:   Dr. Carlos A. Coello Coello                   *
 * Author:      Raquel Hernandez Gomez                        *
 *                                                            *
 * September 2014                                             *
 *************************************************************/

#ifndef _IO_H
#define _IO_H

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "io.h"

/* Lee un archivo, si row y col son cero las establece
   de acuerdo al archivo, si col es diferente de cero
   verifica dimensiones con el archivo, si row es
   diferente de cero y es menor que el numero de registros
   del archivo, trunca el archivo (lee cierta cantidad de registros).
   Si no alcanzan marca error
   Si no se reservo memoria escribe a partir de la posicion start */
double *EMO_File_read(double *data, int *row, int *col, const char *str, int start) {
  double v;
  FILE *fp;
  int i, j;

  if((fp = fopen(str, "r")) == NULL) {
    printf("Error to open file %s\n", str);
    exit(1);
  }

  if(!fscanf(fp, "# %d %d\n", &i, &j)) {
    printf("Error, missing head %s.\n", str); 
    exit(1);
  }

  if(*row != 0 && (start + i) > *row) {
    printf("Error, not enough space %s %d vs %d.\n", str, (start+i), *row);
    exit(1);
  }
 
  if(*row != 0 && *row < i)
    printf("Warning, truncated file %s: %d rows were read, instead of %d.\n", str, *row, i);
 
  if(*col != 0 && *col != j) {
    printf("Error, mismatch dimension %s col (%d vs %d).\n", str, *col, j);
    exit(1);
  }

  if(i == 0 || j == 0) {

    printf("Error in io.c:read, zero rows or columns in file %s (row %d, col %d) \n", str, i, j);
    exit(1);
  }

  *row = i;
  *col = j;

  if(data == NULL) {
    if((data = (double *) malloc(sizeof(double)*(*row)*(*col))) == NULL) {
      printf("Error, not enough memory in EMO_File_read.\n");
      exit(1);
    }

    if(start != 0) {
      printf("EMO_File_read %s: storing at the beginning of the array.\n", str);
      start = 0;
    }
  }

  for(i = 0; i < *row; i++) {
    for(j = 0; j < *col; j++) {
      if(!fscanf(fp, "%lf", &v)) {
        printf("Error to read file %s in row %d, col %d\n", str, i, j);
        exit(1);
      }

      data[(i + start) * (*col) + j] = v;
    }
  }

  fclose(fp);

  return data; 
}

/* Lee un archivo desde start saltando skip lineas */
double *EMO_File_read_skip(double *data, int *row, int *col, const char *str, int start, int skip) {
  int i, j, size = 0, flag, cont = 0;
  double v;
  FILE *fp;

  if(skip < 0) {
    printf("Error, invalid value for skip %d in io.c:EMO_File_read_skip.\n", skip);
    exit(1);
  }

  skip++;

  printf("start %d, skip %d\n", start, skip);

  if((fp = fopen(str, "r")) == NULL) {
    printf("Error to open file %s\n", str);
    exit(1);
  }

  if(!fscanf(fp, "# %d %d\n", &i, &j)) {
    printf("Error, missing head %s.\n", str); 
    exit(1);
  }

  if(*row != 0 && i > *row) {
    printf("Error, not enough space %s (%d vs %d).\n", str, i, *row);
    exit(1);
  }
 
  if(*row != 0 && *row < i)
    printf("Warning, truncated file %s: %d rows were read, instead of %d.\n", str, *row, i);
 
  if(*col != 0 && *col != j) {
    printf("Error, mismatch dimension %s col (%d vs %d).\n", str, *col, j);
    exit(1);
  }

  if(i == 0 || j == 0) {

    printf("Error in io.c:read, zero rows or columns in file %s (row %d, col %d) \n", str, i, j);
    exit(1);
  }

  *row = i;
  *col = j;

  if(data == NULL) {
    size = (*row) / skip;
    size += (start < (*row) % skip)? 1 : 0;

    printf("start %d, skip %d, size %d\n", start, skip, size);

    if((data = (double *) malloc(sizeof(double)*size*(*col))) == NULL) {
      printf("Error, not enough memory in EMO_File_read.\n");
      exit(1);
    }
  }

  for(i = 0; i < *row; i++) {
    flag = 0;

    for(j = 0; j < *col; j++) {
      if(!fscanf(fp, "%lf", &v)) {
        printf("Error to read file %s in row %d, col %d\n", str, i, j);
        exit(1);
      }

      if(i >= start) {
        if(((i-start) % skip) == 0) {

          if(cont >= size) {
            printf("Error io.c, start %d, skip %d, size %d, file %s\n", start, skip, size, str); 
            exit(1);
          }

          data[cont * (*col) + j] = v;
          flag = 1; 
        }
      }

    }

    if(flag)
      cont++;
  }

  *row = cont;

  fclose(fp);

  return data; 
}

void EMO_File_write(double *data, int *filter, int row, int col, const char *str, const char *fmt, unsigned long itera) {
  static const char *mode[] = {"w", "a"};
  int i, j, count;
  FILE *fp;


  if(str == NULL || strcmp(str, "stdout") == 0) {
    fp = stdout;
    fflush(stdout);
  }
  else {
    i = (itera > 0) ? 1 : 0;

    if((fp = fopen(str, mode[i])) == NULL) {
      printf("Error to open file %s\n", str);
      exit(1);
    }
  }

  if(filter == NULL) {

    if(row == 0 || col == 0) {
      printf("Error in io.c:write, zero rows or columns in file %s (row %d, col %d) \n", str, row, col);
      exit(1);
    }

    if(fprintf(fp, "# %d %d\n", row, col) == 0) {
      printf("Error to write file %s\n", str);
      exit(1);
    }

   if(itera > 0 && fprintf(fp, "# iteration %ld\n", itera) == 0) {
     printf("Error to write file %s\n", str);
     exit(1);
   }
 
    for(i = 0; i < row; i++) {
      for(j = 0; j < col; j++) {
        if(fprintf(fp, fmt, data[i*col+j]) == 0) { // fmt = %.20lf for example
          printf("Error to write file %s\n", str);
          exit(1);
        }
      }
      fprintf(fp, "\n");
    }
  }
  else {
    count = 0;

    for(i = 0; i < row; i++)
      if(filter[i]) count++; 

    if(count == 0 || col == 0) {
      printf("Error in io.c:write, zero rows or columns in file %s (row %d, col %d) \n", str, count, col);
      exit(1);
    }

    if(fprintf(fp, "# %d %d\n", count, col) == 0) {
      printf("Error to write file %s\n", str);
      exit(1);
    }

    for(i = 0; i < row; i++) {
      if(filter[i]) { 
        for(j = 0; j < col; j++) {
          if(fprintf(fp, fmt, data[i*col+j]) == 0) {
            printf("Error to write file %s\n", str);
            exit(1);
          }
        }
        fprintf(fp, "\n");
      }
    }
  }

  if(strcmp(str, "stdout") != 0) 
    fclose(fp); 
}

#endif
