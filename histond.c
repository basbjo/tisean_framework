/*Author: Rainer Hegger. Last modified: May 20, 2014 */
/*Changes by Bjoern Bastian:
    2014/07/07: adapted version for n-dimensional histograms
    2014/07/07: option -F for relative frequencies
    2014/08/28: print out intervals and binwidths as comment line
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "routines/tsa.h"

#ifndef _MATH_H
#include <math.h>
#endif

#define WID_STR "Creates a n-d-histogram of a multivariate time series [2014/08/28: n-d-histogram]"

unsigned long length=ULONG_MAX;
unsigned long exclude=0;
unsigned int dim=2;
char *columns=NULL,dimset=0;
unsigned int base=16;
unsigned int verbosity=0xff;
unsigned int stout=1;
char density=1;
char *outfile=NULL;
char *infile=NULL;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted as a"
          " possible datafile.\nIf no datafile is given stdin is read. "
          " Just - also means stdin\n");
  fprintf(stderr,"\t-l length of file [default whole file]\n");
  fprintf(stderr,"\t-x # of lines to ignore [default %ld]\n",exclude);
  fprintf(stderr,"\t-m # of components to be read [default %u]\n",dim);
  fprintf(stderr,"\t-c columns to read [default 1,2]\n");
  fprintf(stderr,"\t-b # of intervals per dim [default %u]\n",base);
  fprintf(stderr,"\t-F output relative frequencies not densities"
	  " [default not set]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.his ;"
          " If no -o is given: stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h show these options\n");
  exit(0);
}

void scan_options(int n,char **argv)
{
  char *out;

  if ((out=check_option(argv,n,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(argv,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(argv,n,'m','u')) != NULL) {
    sscanf(out,"%u",&dim);
    dimset=1;
  }
  if ((out=check_option(argv,n,'c','s')) != NULL)
    columns=out;
  if ((out=check_option(argv,n,'b','u')) != NULL)
    sscanf(out,"%u",&base);
  if ((out=check_option(argv,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(argv,n,'F','n')) != NULL)
    density=0;
  if ((out=check_option(argv,n,'o','o')) != NULL) {
    stout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

int main(int argc,char **argv)
{
  char stdi=0;
  double base_1,norm2;
  double min[dim],interval[dim];
  double **series;
  unsigned long i,j;
  unsigned int n,bi[dim];
  unsigned long *box;
  FILE *fout=NULL;

  if (scan_help(argc,argv))
    show_options(argv[0]);

  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,NULL,verbosity);
  if (infile == NULL)
    stdi=1;

  if (!stout && (outfile == NULL)) {
    if (!stdi) {
      check_alloc(outfile=calloc(strlen(infile)+5,(size_t)1));
      sprintf(outfile,"%s.his",infile);
    }
    else {
      check_alloc(outfile=calloc((size_t)10,(size_t)1));
      sprintf(outfile,"stdin.his");
    }
  }

  if (columns == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&dim,"",dimset,
                                      verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&dim,columns,
                                      dimset,verbosity);

  for (n=0;n<dim;n++) {
    min[n]=interval[n]=series[n][0];
    for (i=1;i<length;i++) {
      if (series[n][i] < min[n]) min[n]=series[n][i];
      else if (series[n][i] > interval[n]) interval[n]=series[n][i];
    }
    interval[n] -= min[n];

    for (i=0;i<length;i++) {
      series[n][i]=(series[n][i]-min[n]);
    }
  }

  check_alloc(box=(unsigned long*)malloc(sizeof(unsigned long*)*pow(base,dim)));
  for (i=0;i<pow(base,dim);i++)
  {
    box[i]=1;
  }
  base_1=(double)base;
  norm2=(double)(length+pow(base,dim));
  if (density) {
    for (n=0;n<dim;n++) {
      norm2*=(double)interval[n]/base_1;
    }
  }

  for (i=0;i<length;i++) {
    j=0;
    for (n=0;n<dim;n++) {
      bi[n]=(unsigned int)(series[n][i]*base_1/interval[n]);
      bi[n]=(bi[n]>=base)? base-1:bi[n];
      j+=bi[n]*pow(base,dim-1-n);
    }
    box[j]++;
  }

  if (!stout)
    test_outfile(outfile);

  /* print out intervals */
  fout=fopen(outfile,"w");
  if (stout) {
    fprintf(stdout,"#intervals: ");
    for (n=0;n<dim;n++) {
      fprintf(stdout,"%e ",interval[n]);
    }
    fprintf(stdout,"\n");
  }
  else {
    fprintf(fout,"#intervals: ");
    for (n=0;n<dim;n++) {
      fprintf(fout,"%e ",interval[n]);
    }
    fprintf(fout,"\n");
  }

  for (n=0;n<dim;n++) {
    interval[n] /= base_1;
  }
  /* print out binwidths */
  fout=fopen(outfile,"w");
  if (stout) {
    fprintf(stdout,"#binwidths: ");
    for (n=0;n<dim;n++) {
      fprintf(stdout,"%e ",interval[n]);
    }
    fprintf(stdout,"\n");
  }
  else {
    fprintf(fout,"#intervals: ");
    for (n=0;n<dim;n++) {
      fprintf(fout,"%e ",interval[n]);
    }
    fprintf(fout,"\n");
  }

  /* print out histogram */
  for (i=0;i<pow(base,dim);i++) {
    if (stout) {
      for (n=0;n<dim;n++) {
        fprintf(stdout,"%e ",
                ((double)(((i)/(unsigned int)
                    pow(base,dim-1-n))%base)+0.5)*interval[n]+min[n]);
      }
      fprintf(stdout,"%e\n",(double)box[i]/norm2);
    }
    else {
      for (n=0;n<dim;n++) {
        fprintf(fout,"%e ",
                ((double)(((i)/(unsigned int)
                    pow(base,dim-1-n))%base)+0.5)*interval[n]+min[n]);
      }
      fprintf(fout,"%e\n",(double)box[i]/norm2);
    }
    if (!((i+1)%base)) {
      if (stout)
        fprintf(stdout,"\n");
      else
        fprintf(fout,"\n");
    }
  }
  if (!stout)
    fclose(fout);

  /*Freeing all allocated arrays*/
  if (outfile != NULL) free(outfile);
  if (infile != NULL) free(infile);
  if (columns != NULL) free(columns);
  free(box);
  for (n=0;n<dim;n++) {
    free(series[n]);
  }
  free(series);

  return 0;
}
