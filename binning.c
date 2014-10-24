/*
 *   This file is part of TISEAN
 *
 *   Copyright (c) 1998-2007 Rainer Hegger, Holger Kantz, Thomas Schreiber
 *
 *   TISEAN is free software; you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation; either version 2 of the License, or
 *   (at your option) any later version.
 *
 *   TISEAN is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with TISEAN; if not, write to the Free Software
 *   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
/*Author: Rainer Hegger. Last modified May 16, 2014*/
/*Changes by Bjoern Bastian:
    2014/09/22: fork for binning instead of histogram creation
    2014/10/22: option -R to set reference binning range and output range
    2014/10/22: option -s to set reference binning range by argument
    2014/10/22: option -S to set reference binning range and output range
*/

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "routines/tsa.h"

#define WID_STR "Averages second column with respect to binning of first column"

unsigned long length=ULONG_MAX;
unsigned long minmaxlength=3;
unsigned long base=50;
unsigned long exclude=0;
unsigned int verbosity=0xff;
char *columns=NULL;
char my_stdout=1,gotsize=0,cropoutput=0;
char *outfile=NULL;
char *infile=NULL;
char *minmaxfile=NULL,*minmaxstring=NULL;

void show_options(char *progname)
{
  what_i_do(progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted as a"
	  " possible datafile.\nIf no datafile is given stdin is read. "
	  " Just - also means stdin.\n");
  fprintf(stderr,"A minmax file contains column wise minima and maxima on the"
          " first and second row.\n");
  fprintf(stderr,"\t-l length of file [default whole file]\n");
  fprintf(stderr,"\t-x # of lines to ignore [default %ld]\n",exclude);
  fprintf(stderr,"\t-c column selection [default 1,2]\n");
  fprintf(stderr,"\t-b # of intervals [default %ld]\n",base);
  fprintf(stderr,"\t-r minmax file to set reference range with # of intervals [optional]\n");
  fprintf(stderr,"\t-R minmax file to set reference range and resctrict output [optional]\n");
  fprintf(stderr,"\t-s num,num to set reference range with # of intervals [optional]\n");
  fprintf(stderr,"\t-S num,num to set reference range and resctrict output [optional]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.bins ;"
	  " If no -o is given: stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h show these options\n");
  exit(0);
}

void scan_options(int n,char **str)
{
  char *out;

  if ((out=check_option(str,n,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(str,n,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(str,n,'c','s')) != NULL)
    columns=out;
  if ((out=check_option(str,n,'b','u')) != NULL)
    sscanf(out,"%lu",&base);
  if ((out=check_option(str,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(str,n,'r','o')) != NULL) {
    if (strlen(out) > 0)
      minmaxfile=out;
  }
  if ((out=check_option(str,n,'R','o')) != NULL) {
    if (strlen(out) > 0) {
      minmaxfile=out;
      cropoutput=1;
    }
  }
  if ((out=check_option(str,n,'s','o')) != NULL) {
    if (strlen(out) > 0)
      minmaxstring=out;
  }
  if ((out=check_option(str,n,'S','o')) != NULL) {
    if (strlen(out) > 0) {
      minmaxstring=out;
      cropoutput=1;
    }
  }
  if ((out=check_option(str,n,'o','o')) != NULL) {
    my_stdout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

int main(int argc,char **argv)
{
  char stdi=0;
  int refcolumn;
  unsigned long i,j,N;
  unsigned int dummy=2;
  unsigned long offset,negoffset,range,fullrange;
  double x,size;
  double min,interval,refmin,refinterval;
  double **series,*minmax=NULL;
  long *box;
  double *sum,*sumsq;
  FILE *fout,*test;

  if (scan_help(argc,argv))
    show_options(argv[0]);

  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,0L,verbosity);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,1));
      strcpy(outfile,infile);
      strcat(outfile,".bins");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,1));
      strcpy(outfile,"stdin.bins");
    }
  }
  if (!my_stdout)
    test_outfile(outfile);

  /*Get reference range for options '-r' and '-R'*/
  if (minmaxfile != NULL) {
    test=fopen(minmaxfile,"r");
    if (test == NULL) {
      fprintf(stderr,"File %s not found!\n",minmaxfile);
      exit(HISTOGRAM__MINMAX_MISSING_OR_WRONG_FORMAT);
    }
    if (columns == NULL) {
      refcolumn=1;
    }
    else {
      sscanf(columns,"%d",&refcolumn);
    }
    if (verbosity&VER_INPUT) {
      fprintf(stderr,"Get reference range from %s, reading column %u\n",
          minmaxfile,refcolumn);
    }

    minmax=(double*)get_series(minmaxfile,&minmaxlength,0,refcolumn,verbosity);
    if (minmaxlength != 2) {
      fprintf(stderr,"Wrong format in file '%s'. Needs exactly two lines"
          " with minima and maxima for each column.\n",minmaxfile);
      exit(HISTOGRAM__MINMAX_MISSING_OR_WRONG_FORMAT);
    }
    refmin=minmax[0];
    refinterval=minmax[1]-refmin;
  }

  /*Get reference range for options '-s' and '-S'*/
  if (minmaxstring != NULL) {
    if (minmax == NULL) {
      check_alloc(minmax=(double*)malloc(sizeof(double)*2));
    }
    sscanf(minmaxstring,"%lf,%lf",&minmax[0],&minmax[1]);
    refmin=minmax[0];
    refinterval=minmax[1]-refmin;
  }

  /*Read data*/
  if (columns == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&dummy,"",
        (char)1,verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&dummy,columns,
        (char)1,verbosity);

  /*Get data minimum and interval*/
  min=interval=series[0][0];
  for (i=1;i<length;i++) {
    if (series[0][i] < min) min=series[0][i];
    else if (series[0][i] > interval) interval=series[0][i];
  }
  interval -= min;

  /*Settings*/
  if (minmaxfile != NULL || minmaxstring != NULL) {
    size=refinterval/base;
    if (refmin > min) {
      offset=(long)((refmin-min)/size);
      if (!cropoutput) {
        negoffset=0;
      }
      else {
        negoffset=(long)((refmin-min)/size);
      }
    }
    else {
      offset=0;
      negoffset=(long)((min-refmin)/size);
    }
    range=(long)((min+interval-refmin)/size)+offset;
    fullrange=range;
    if (cropoutput && ((min+interval) > (refmin+refinterval))) {
      range-=((long)((min+interval-refmin-refinterval)/size));
    }
  }
  else {
    refmin=min;
    refinterval=interval;
    size=interval/base;
    offset=0;
    negoffset=0;
    range=base;
    fullrange=range;
  }

  /*Binning*/
  if (range > 0) {
    check_alloc(box=(long*)malloc(sizeof(long)*fullrange));
    check_alloc(sum=(double*)malloc(sizeof(double)*fullrange));
    check_alloc(sumsq=(double*)malloc(sizeof(double)*fullrange));
    for (i=negoffset;i<range;i++) {
      box[i]=0;
      sum[i]=0.0;
      sumsq[i]=0.0;
    }
    for (i=0;i<length;i++) {
      j=(long)((series[0][i]-refmin)*base/refinterval+offset);
      if ((!cropoutput) || ((long)(min+interval-refmin-refinterval) == 0)) {
        if (j >= range) {
          j=range-1;
        }
      }
      box[j]++;
      sum[j]+=series[1][i];
      sumsq[j]+=pow(series[1][i],2);
    }
  }

  if (!my_stdout) {
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    fprintf(fout,"#binning range: [%e:%e]\n",min,min+interval);
    if (minmaxfile != NULL || minmaxstring != NULL)
      fprintf(fout,"#reference interval: [%e:%e]\n",refmin,refmin+refinterval);
    fprintf(fout,"#bin_center mean stddev_of_mean n_bin_entries\n");
    for (i=negoffset;i<range;i++) {
      x=(double)(i*size-offset*size);
      N=box[i];
      if(N>0) {
        fprintf(fout,"%e %e",(x+size/2.0)+refmin,sum[i]/N);
        if(N>1) {
          fprintf(fout," %e %ld\n",pow((sumsq[i]-pow(sum[i],2)/N)/(N-1)/N,0.5),N);
        }
        else {
          fprintf(fout," nan %ld\n",N);
        }
      }
      else {
        fprintf(fout,"%e %s %s %ld\n",(x+size/2.0)+refmin,"nan","nan",N);
      }
    }
    fclose(fout);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    fprintf(stdout,"#binning range: [%e:%e]\n",min,min+interval);
    if (minmaxfile != NULL || minmaxstring != NULL)
      fprintf(stdout,"#reference interval: [%e:%e]\n",refmin,refmin+refinterval);
    fprintf(stdout,"#bin_center mean stddev_of_mean n_bin_entries\n");
    for (i=negoffset;i<range;i++) {
      x=(double)(i*size-offset*size);
      N=box[i];
      if(N>0) {
        fprintf(stdout,"%e %e",(x+size/2.0)+refmin,sum[i]/N);
        if(N>1) {
          fprintf(stdout," %e %ld\n",pow((sumsq[i]-pow(sum[i],2)/N)/(N-1)/N,0.5),N);
        }
        else {
          fprintf(stdout," nan %ld\n",N);
        }
      }
      else {
        fprintf(stdout,"%e %s %s %ld\n",(x+size/2.0)+refmin,"nan","nan",N);
      }
      fflush(stdout);
    }
  }
  return 0;
}
