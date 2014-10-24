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
    2014/05/16: option -r to set reference binning range by minmax file
    2014/10/21: option -R to set reference binning range and output range
    2014/10/21: option -s to set reference binning range by argument
    2014/10/21: option -S to set reference binning range and output range
*/

#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "routines/tsa.h"

#define WID_STR "Creates a histogram of a onedimensional dataset [2014/05/16: option -r added]"

unsigned long length=ULONG_MAX;
unsigned long minmaxlength=3;
unsigned long base=50;
unsigned long exclude=0;
unsigned int column=1;
unsigned int verbosity=0xff;
char my_stdout=1,gotsize=0,density=0,cropoutput=0;
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
  fprintf(stderr,"\t-c column to read [default %d]\n",column);
  fprintf(stderr,"\t-b # of intervals [default %ld]\n",base);
  fprintf(stderr,"\t-D output densities not relative frequencies"
	  " [default not set]\n");
  fprintf(stderr,"\t-r minmax file to set reference range with # of intervals [optional]\n");
  fprintf(stderr,"\t-R minmax file to set reference range and resctrict output [optional]\n");
  fprintf(stderr,"\t-s num,num to set reference range with # of intervals [optional]\n");
  fprintf(stderr,"\t-S num,num to set reference range and resctrict output [optional]\n");
  fprintf(stderr,"\t-o output file [default 'datafile'.his ;"
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
  if ((out=check_option(str,n,'c','u')) != NULL)
    sscanf(out,"%u",&column);
  if ((out=check_option(str,n,'b','u')) != NULL)
    sscanf(out,"%lu",&base);
  if ((out=check_option(str,n,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(str,n,'D','n')) != NULL)
    density=1;
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
  unsigned long i,j;
  unsigned long offset,negoffset,range,fullrange;
  double x,norm,size;
  double min,interval,refmin,refinterval;
  double *series,*minmax=NULL;
  double average,var;
  long *box;
  FILE *fout,*test;

  if (scan_help(argc,argv))
    show_options(argv[0]);
  
  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    what_i_do(argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,&column,verbosity);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,1));
      strcpy(outfile,infile);
      strcat(outfile,".his");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,1));
      strcpy(outfile,"stdin.his");
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
    if (verbosity&VER_INPUT) {
      fprintf(stderr,"Get reference range from %s, reading column %u\n",
          minmaxfile,column);
    }

    minmax=(double*)get_series(minmaxfile,&minmaxlength,0,column,verbosity);
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
  series=(double*)get_series(infile,&length,exclude,column,verbosity);
  variance(series,length,&average,&var);

  /*Get data minimum and interval*/
  min=interval=series[0];
  for (i=1;i<length;i++) {
    if (series[i] < min) min=series[i];
    else if (series[i] > interval) interval=series[i];
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
    for (i=negoffset;i<range;i++)
      box[i]=0;
    for (i=0;i<length;i++) {
      j=(long)((series[i]-refmin)*base/refinterval+offset);
      if ((!cropoutput) || ((long)(min+interval-refmin-refinterval) == 0)) {
        if (j >= range) {
          j=range-1;
        }
      }
      box[j]++;
    }
  }

  if (!density)
    norm=1.0/(double)length;
  else
    norm=1.0/(double)length/size;

  if (!my_stdout) {
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    fprintf(fout,"#interval of data:   [%e:%e]\n",min,min+interval);
    if (minmaxfile != NULL || minmaxstring != NULL)
      fprintf(fout,"#reference interval: [%e:%e]\n",refmin,refmin+refinterval);
    fprintf(fout,"#average= %e\n",average);
    fprintf(fout,"#standard deviation= %e\n",var);
    for (i=negoffset;i<range;i++) {
      x=(double)(i*size-offset*size);
      fprintf(fout,"%e %e\n",(x+size/2.0)+refmin,(double)box[i]*norm);
    }
    fclose(fout);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    fprintf(stdout,"#interval of data:   [%e:%e]\n",min,min+interval);
    if (minmaxfile != NULL || minmaxstring != NULL)
      fprintf(stdout,"#reference interval: [%e:%e]\n",refmin,refmin+refinterval);
    fprintf(stdout,"#average= %e\n",average);
    fprintf(stdout,"#standard deviation= %e\n",var);
    for (i=negoffset;i<range;i++) {
      x=(double)(i*size-offset*size);
      fprintf(stdout,"%e %e\n",(x+size/2.0)+refmin,(double)box[i]*norm);
      fflush(stdout);
    }
  }
  return 0;
}
