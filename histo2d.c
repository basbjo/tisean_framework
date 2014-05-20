/*Author: Rainer Hegger. Last modified: May 20, 2014 */
/*Changes by Bjoern Bastian:
    2014/05/20: option -r to set reference binning range
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "routines/tsa.h"

#ifndef _MATH_H
#include <math.h>
#endif

#define WID_STR "Creates a 2d-histogram of an bivariate time series [2014/05/20: option -r added]"

unsigned long length=ULONG_MAX;
unsigned long minmaxlength=3;
unsigned long exclude=0;
char *column=NULL;
unsigned int base=16;
unsigned int verbosity=0xff;
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
  fprintf(stderr,"\t-l length of file [default whole file]\n");
  fprintf(stderr,"\t-x # of lines to ignore [default %ld]\n",exclude);
  fprintf(stderr,"\t-c columns to read [default 1,2]\n");
  fprintf(stderr,"\t-b # of intervals per dim [default %u]\n",base);
  fprintf(stderr,"\t-r reference file for binning range [optional]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.dat ;"
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
  if ((out=check_option(argv,n,'c','s')) != NULL)
    column=out;
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
  unsigned int dim=2;
  char stdi=0;
  char *col=NULL;
  double base_1,sx,sy,logmax,logout,norm1,norm2;
  double min[2],interval[2],refmin[2],refinterval[2];
  double **series,**minmax;
  unsigned long i,j,lmax;
  unsigned int bi,bj;
  unsigned long **box,*box1d;
  FILE *fout=NULL,*test=NULL;

  if (scan_help(argc,argv))
    show_options(argv[0]);

  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  /*Get reference range for option '-r'*/
  if (minmaxfile != NULL) {
    test=fopen(minmaxfile,"r");
    if (test == NULL) {
      fprintf(stderr,"File %s not found!\n",minmaxfile);
      exit(HISTOGRAM__MINMAX_MISSING_OR_WRONG_FORMAT);
    }
    if (verbosity&VER_INPUT) {
      fprintf(stderr,"Get reference range from file %s\n",minmaxfile);
    }

    if (column == NULL) {
      minmax=(double**)get_multi_series(minmaxfile,&minmaxlength,0,&dim,"",1,
                                        verbosity);
    }
    else {
      check_alloc(col=calloc(strlen(column),(size_t)1));
      strcpy(col,column);
      minmax=(double**)get_multi_series(minmaxfile,&minmaxlength,0,&dim,col,
                                        1,verbosity);
    }

    if(minmaxlength!=2) {
      fprintf(stderr,"Wrong format in file '%s'. Needs exactly two lines"
          " with minima and maxima for each column.\n",minmaxfile);
      exit(HISTOGRAM__MINMAX_MISSING_OR_WRONG_FORMAT);
    }
    refmin[0]=minmax[0][0];
    refinterval[0]=minmax[0][1]-refmin[0];
    refmin[1]=minmax[1][0];
    refinterval[1]=minmax[1][1]-refmin[1];
  }

  /*Read data*/
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

  if (column == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&dim,"",1,
                                      verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&dim,column,
                                      1,verbosity);

  /*Get intervals and shift data to zero minimum*/
  min[0]=interval[0]=series[0][0];
  min[1]=interval[1]=series[1][0];
  for (i=1;i<length;i++) {
    if (series[0][i] < min[0]) min[0]=series[0][i];
    else if (series[0][i] > interval[0]) interval[0]=series[0][i];
    if (series[1][i] < min[1]) min[1]=series[1][i];
    else if (series[1][i] > interval[1]) interval[1]=series[1][i];
  }
  interval[0] -= min[0];
  interval[1] -= min[1];

  for (i=0;i<length;i++) {
    series[0][i]=(series[0][i]-min[0]);
    series[1][i]=(series[1][i]-min[1]);
  }

  /*Binning*/
  check_alloc(box1d=(unsigned long*)malloc(sizeof(unsigned long)*base));
  for (i=0;i<base;i++)
    box1d[i]=1;

  check_alloc(box=(unsigned long**)malloc(sizeof(unsigned long*)*base));
  for (i=0;i<base;i++) {
    check_alloc(box[i]=(unsigned long*)malloc(sizeof(unsigned long)*base));
    for (j=0;j<base;j++)
      box[i][j]=1;
  }
  base_1=(double)base;
  sx=interval[0]/base_1;
  sy=interval[1]/base_1;
  norm1=(double)(length+base)*sx;
  norm2=(double)(length+base*base)*sx*sy;

  for (i=0;i<length;i++) {
    bi=(unsigned int)(series[0][i]*base_1/interval[0]);
    bj=(unsigned int)(series[1][i]*base_1/interval[1]);
    bi=(bi>=base)? base-1:bi;
    bj=(bj>=base)? base-1:bj;
    box[bi][bj]++;
    box1d[bi]++;
  }

  lmax=0;
  for (i=0;i<base;i++)
    for (j=0;j<base;j++)
      if (box[i][j] > 0)
        lmax=(box[i][j]>lmax)? box[i][j]:lmax;
  logmax=log((double)lmax/norm2);

  if (!stout)
    test_outfile(outfile);

  fout=fopen(outfile,"w");

  interval[0] /= base_1;
  interval[1] /= base_1;
  for (i=0;i<base;i++) {
    for (j=0;j<base;j++) {
      logout=log((double)box[i][j]/norm2)-logmax;
      if (stout) {
        fprintf(stdout,"%e %e %e %e %e\n",((double)(i)+0.5)*interval[0]+min[0],
                ((double)(j)+0.5)*interval[1]+min[1],
                (double)box[i][j]/norm2,
                (double)box[i][j]/(double)box1d[i]/norm2*norm1,-logout);
      }
      else {
        fprintf(fout,"%e %e %e %e %e\n",((double)(i)+0.5)*interval[0]+min[0],
                ((double)(j)+0.5)*interval[1]+min[1],
                (double)box[i][j]/norm2,
                (double)box[i][j]/(double)box1d[i]/norm2*norm1,-logout);
      }
    }
    if (stout)
      fprintf(stdout,"\n");
    else
      fprintf(fout,"\n");
  }
  if (!stout)
    fclose(fout);

  /*Freeing all allocated arrays*/
  if (outfile != NULL) free(outfile);
  if (infile != NULL) free(infile);
  if (minmaxfile != NULL) free(minmaxfile);
  if (column != NULL) free(column);
  free(box1d);
  for (i=0;i<base;i++)
    free(box[i]);
  free(box);
  free(series[0]);
  free(series[1]);
  free(series);
  if (minmaxfile != NULL) {
    free(minmax[0]);
    free(minmax[1]);
    free(minmax);
  }

  return 0;
}
