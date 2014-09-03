/*Derived from source code by Rainer Hegger. */
/*Author: Bjoern Bastian.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "routines/tsa.h"

#ifndef _MATH_H
#include <math.h>
#endif

#define WID_STR "Calculates column-wise negentropies from 1d-histogram"

unsigned long length=ULONG_MAX;
unsigned int dim=1;
unsigned long minmaxlength=3;
unsigned long exclude=0;
unsigned int base=300;
unsigned int verbosity=0xff;
char *columns=NULL,dimset=0;
unsigned int stout=1;
char *outfile=NULL;
char *infile=NULL;
char *minmaxfile=NULL;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted as a"
          " possible datafile.\nIf no datafile is given stdin is read. "
          " Just - also means stdin\n");
  fprintf(stderr,"\t-l # of lines to use [default is whole file]\n");
  fprintf(stderr,"\t-x # of lines to ignore [default %ld]\n",exclude);
  fprintf(stderr,"\t-m # of components to be read [default %u]\n",dim);
  fprintf(stderr,"\t-c column selection [default 1,...,# of components]\n");
  fprintf(stderr,"\t-b # of intervals per dim [default %u]\n",base);
  fprintf(stderr,"\t-r reference file for binning range [optional]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.nen ;"
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
  if ((out=check_option(argv,n,'r','o')) != NULL) {
    if (strlen(out) > 0)
      minmaxfile=out;
  }
  if ((out=check_option(argv,n,'o','o')) != NULL) {
    stout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

int main(int argc,char **argv)
{
  char stdi=0;
  char *cols=NULL;
  unsigned long i,j,k;
  unsigned long *offset,*negoffset,*range;
  double entropygauss;
  double pi=3.14159265358979;
  double e=2.71828182845905;
  double x,norm,*size;
  double *min,*interval,*refmin,*refinterval;
  double **series,**minmax;
  double *average,*std,*entropy;
  unsigned long **box;
  FILE *fout=NULL,*test=NULL;

  if (scan_help(argc,argv))
    show_options(argv[0]);

  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  /*Get reference range for option '-r'*/
  check_alloc(refmin=(double*)malloc(sizeof(double)*dim));
  check_alloc(refinterval=(double*)malloc(sizeof(double)*dim));
  if (minmaxfile != NULL) {
    test=fopen(minmaxfile,"r");
    if (test == NULL) {
      fprintf(stderr,"File %s not found!\n",minmaxfile);
      exit(HISTOGRAM__MINMAX_MISSING_OR_WRONG_FORMAT);
    }
    if (verbosity&VER_INPUT) {
      fprintf(stderr,"Get reference ranges from file %s\n",minmaxfile);
    }

    if (columns == NULL) {
      minmax=(double**)get_multi_series(minmaxfile,&minmaxlength,0,&dim,"",
                                          dimset,verbosity);
    }
    else {
      check_alloc(cols=calloc(strlen(columns),(size_t)1));
      strcpy(cols,columns);
      minmax=(double**)get_multi_series(minmaxfile,&minmaxlength,0,&dim,cols,
                                        dimset,verbosity);
    }

    if(minmaxlength!=2) {
      fprintf(stderr,"Wrong format in file '%s'. Needs exactly two lines"
          " with minima and maxima for each column.\n",minmaxfile);
      exit(HISTOGRAM__MINMAX_MISSING_OR_WRONG_FORMAT);
    }
    for (i=0;i<dim;i++) {
      refmin[i]=minmax[i][0];
      refinterval[i]=minmax[i][1]-refmin[i];
    }
  }

  /*Read data*/
  infile=search_datafile(argc,argv,NULL,verbosity);
  if (infile == NULL)
    stdi=1;

  if (!stout && (outfile == NULL)) {
    if (!stdi) {
      check_alloc(outfile=calloc(strlen(infile)+5,(size_t)1));
      sprintf(outfile,"%s.nen",infile);
    }
    else {
      check_alloc(outfile=calloc((size_t)10,(size_t)1));
      sprintf(outfile,"stdin.nen");
    }
  }

  if (columns == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&dim,"",dimset,
                    verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&dim,columns,
                    dimset,verbosity);

  /*Get data minimum, interval and variance*/
  check_alloc(min=(double*)malloc(sizeof(double)*dim));
  check_alloc(interval=(double*)malloc(sizeof(double)*dim));
  check_alloc(average=(double*)malloc(sizeof(double)*dim));
  check_alloc(std=(double*)malloc(sizeof(double)*dim));
  for (i=0;i<dim;i++) {
    min[i]=interval[i]=series[i][0];
    for (j=1;j<length;j++) {
      if (series[i][j] < min[i]) min[i]=series[i][j];
      else if (series[i][j] > interval[i]) interval[i]=series[i][j];
    }
    interval[i] -= min[i];
  /*variance calculates standard deviation!*/
    variance(series[i],length,&average[i],&std[i]);
  /*use the unbiased estimator*/
  std[i]*=pow((double)length/(double)(length-1),0.5);
  }

  /*Settings*/
  check_alloc(size=(double*)malloc(sizeof(double)*dim));
  check_alloc(offset=(long*)malloc(sizeof(long)*dim));
  check_alloc(negoffset=(long*)malloc(sizeof(long)*dim));
  check_alloc(range=(long*)malloc(sizeof(long)*dim));
  if (minmaxfile != NULL) {
    for (i=0;i<dim;i++) {
      size[i]=refinterval[i]/base;
      if (refmin[i] > min[i]) {
        offset[i]=(long)((refmin[i]-min[i])/size[i]);
        negoffset[i]=0;
      }
      else {
        offset[i]=0;
        negoffset[i]=(long)((min[i]-refmin[i])/size[i]);
      }
      range[i]=(long)((min[i]+interval[i]-refmin[i])/size[i])+offset[i];
    }
  }
  else {
    for (i=0;i<dim;i++) {
      refmin[i]=min[i];
      refinterval[i]=interval[i];
      size[i]=interval[i]/base;
      offset[i]=0;
      negoffset[i]=0;
      range[i]=base;
    }
  }

  /*Binning*/
  check_alloc(box=(unsigned long**)malloc(sizeof(unsigned long*)*dim));
  for (i=0;i<dim;i++) {
    if (range[i] > 0) {
      check_alloc(box[i]=(long*)malloc(sizeof(long)*range[i]));
      for (j=negoffset[i];j<range[i];j++)
        box[i][j]=0;
      for (j=0;j<length;j++) {
        k=(long)((series[i][j]-refmin[i])*base/refinterval[i]+offset[i]);
        if (k >= range[i]) {
          k=range[i]-1;
        }
        box[i][k]++;
      }
    }
  }

  /*Entropy*/
  check_alloc(entropy=(double*)malloc(sizeof(double)*dim));
  norm=1.0/(double)length;
  for (i=0;i<dim;i++) {
    entropy[i]=0.0;
    for (j=negoffset[i];j<range[i];j++) {
      if (box[i][j]>0) {
        x=norm*(double)box[i][j];
        entropy[i]-=log(x/size[i])*x;
      }
    }
  }

  if (!stout) {
    test_outfile(outfile);
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    fprintf(fout,"#component negentropy entropy_of_gaussian entropy_of_system"
              " (%d bins)\n",base);
    for (i=0;i<dim;i++) {
      /*exact entropy of a Gaussian with same variance*/
      entropygauss=0.5*log(2.0*pi*e*pow(std[i],2));
      fprintf(fout,"%ld %e %e %e\n",i+1,entropygauss-entropy[i],
        entropygauss,entropy[i]);
    }
    fclose(fout);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    fprintf(stdout,"#component negentropy entropy_of_gaussian entropy_of_system"
              " (%d bins)\n",base);
    for (i=0;i<dim;i++) {
      /*exact entropy of a Gaussian with same variance*/
      entropygauss=0.5*log(2.0*pi*e*pow(std[i],2));
      fprintf(stdout,"%ld %e %e %e\n",i+1,entropygauss-entropy[i],
        entropygauss,entropy[i]);
      fflush(stdout);
    }
  }

  /*Freeing all allocated arrays*/
  if (outfile != NULL) free(outfile);
  if (infile != NULL) free(infile);
  if (columns != NULL) free(columns);
  for (i=0;i<dim;i++) {
      free(box[i]);
  }
  free(box);
  for (i=0;i<dim;i++) {
    free(series[i]);
  }
  free(series);
  free(entropy);
  if (minmaxfile != NULL) {
    for (i=0;i<dim;i++) {
      free(minmax[i]);
    }
    free(minmax);
    free(minmaxfile);
  }

  return 0;
}
