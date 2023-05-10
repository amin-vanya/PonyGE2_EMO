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

double *EMO_File_read(double *data, int *row, int *col, const char *str, int start);
double *EMO_File_read_skip(double *data, int *row, int *col, const char *str, int start, int skip);
void EMO_File_write(double *data, int *filter, int row, int col, const char *str, const char *fmt, unsigned long itera);

