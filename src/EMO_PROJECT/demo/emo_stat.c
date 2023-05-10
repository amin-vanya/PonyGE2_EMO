
// Inspired by  http://des-testbed.net/node/166
// http://des-testbed.net/node/166
// https://en.wikipedia.org/wiki/Quartile
// https://en.wikipedia.org/wiki/Interquartile_range
// https://en.wikipedia.org/wiki/Quantile

// http://stackoverflow.com/questions/18707029/how-to-combine-two-box-whisker-plots-into-one-using-gnuplot
// http://stackoverflow.com/questions/15404628/how-can-i-generate-box-and-whisker-plots-with-variable-box-width-in-gnuplot

// fractales
// http://chartgnuplot.sourceforge.net/advanced.html
// http://benice-equation.blogspot.mx/2013_01_01_archive.html
// http://ayapin-film.sakura.ne.jp/Gnuplot/Tips/fractal.html


#include <stdlib.h>
#include "emo.h"

int main(int argc, char **argv) {
  double vmin, vmax, vmedian, vmean, vvar, vstd, iqr;
  double *data, q[3], **sort;
  int i, j, imin, imax, imedian, size = 0, col = 1;

  if(argc < 2) {
    printf("\nSyntax: %s data_file(s)\n", argv[0]);
    return 1;
  }
 
  data = EMO_File_read(NULL, &size, &col, argv[1], 0);

  // temporary array for calculating median
  if((sort = (double **) malloc(sizeof(double *) * size)) == NULL) {
    printf("Error, not enough memory in %s.\n", argv[0]);
    return 1;
  }

  for(i = 0; i < size; i++) {
    if((sort[i] = (double *) malloc(sizeof(double) * 2)) == NULL) {
      printf("Error, not enough memory in %s (2)\n", argv[0]);
      return 1;
    }
  }

  printf("# file|min|max|median|mean|var|std|q1|q2|q3|iqr|imin|imax|imedian\n");

  for(i = 1; i < argc; i++) {
    vmin    = EMO_dmin(&imin, data, NULL, size);
    vmax    = EMO_dmax(&imax, data, NULL, size);
    vmedian = EMO_median(&imedian, data, NULL, sort, size);
    vmean   = EMO_mean(data, NULL, size);
    vvar    = EMO_var(data, NULL, vmean, size);
    vstd    = EMO_std(data, NULL, vmean, size);
    iqr = EMO_quartile(q, data, NULL, sort, size);

    printf("%s|%e|%e|%e|%e|%e|%e|%e|%e|%e|%e|%d|%d|%d\n", argv[i],vmin, vmax, vmedian, vmean, vvar, vstd, q[0], q[1], q[2], iqr, imin+1, imax+1, imedian+1);

    if(i+1 < argc) {
      free(data);
      j = 0;
      data = EMO_File_read(NULL, &j, &col, argv[i+1], 0);

      if(j != size) {
        printf("Warning, different number of rows in the file %s (%d vs %d)\n", argv[i+1], j, size);
      }
    }
  }

  free(data);

  for(i = 0; i < size; i++)
    free(sort[i]);

  free(sort);

  return 0;
}

