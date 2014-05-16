/*
 *   This file is derived from corr.c of TISEAN,
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
/*Author: Rainer Hegger. Last modified: Sep 3, 1999
 * Adapted by Bjoern Bastian. Last modified: May 16, 2014 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "routines/tsa.h"

#define WID_STR "Enumerates and formats a data set"

unsigned long length=ULONG_MAX,exclude=0;
unsigned int column=1;
unsigned int verbosity=0xff;
char *outfile=NULL,stout=1;
double *array;
char *infile=NULL;

void show_options(char *progname)
{
  fprintf(stderr, "\n%s: %s\n\n",progname,WID_STR);
  fprintf(stderr," Usage: %s [options]\n",progname);
  fprintf(stderr," Options:\n");
  fprintf(stderr,"Everything not being a valid option will be interpreted"
          " as a possible"
          " datafile.\nIf no datafile is given stdin is read. Just - also"
          " means stdin.\n");
  fprintf(stderr,"\t-l # of lines to use [default is whole file]\n");
  fprintf(stderr,"\t-x # of lines to ignore [default 0]\n");
  fprintf(stderr,"\t-c column to be read [default 1]\n");
  fprintf(stderr,"\t-o output file name [default 'datafile'.same; no -o"
  " means stdout]\n");
  fprintf(stderr,"\t-V verbosity level [default 1]\n\t\t"
          "0='only panic messages'\n\t\t"
          "1='+ input/output messages'\n");
  fprintf(stderr,"\t-h show these options\n");
  fprintf(stderr,"\n");
  exit(0);
}

void scan_options(int argc,char **argv)
{
  char *out;

  if ((out=check_option(argv,argc,'l','u')) != NULL)
    sscanf(out,"%lu",&length);
  if ((out=check_option(argv,argc,'x','u')) != NULL)
    sscanf(out,"%lu",&exclude);
  if ((out=check_option(argv,argc,'c','u')) != NULL)
    sscanf(out,"%u",&column);
  if ((out=check_option(argv,argc,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(argv,argc,'o','o')) != NULL) {
    stout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

double same(long i)
{
  double c=0.0;

  c = array[i];

  return c;
}

int main(int argc,char** argv)
{
  char stdi=0;
  long i;
  FILE *fout=NULL;

  /* get options */
  if (scan_help(argc,argv))
    show_options(argv[0]);

  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    fprintf(stderr, "\n%s: %s\n\n",argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,&column,verbosity);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,(size_t)1));
      strcpy(outfile,infile);
      strcat(outfile,".same");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.same");
    }
  }
  if (!stout)
    test_outfile(outfile);

  /* get data */
  array=(double*)get_series(infile,&length,exclude,column,verbosity);

  /* processing */

  /* write results */
  if (!stout) {
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    for (i=0;i<length;i++) {
      fprintf(fout,"%ld %e\n",i,same(i));
      fflush(fout);
    }
    fclose(fout);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    for (i=0;i<length;i++) {
      fprintf(stdout,"%ld %e\n",i,same(i));
      fflush(stdout);
    }
  }

  if (outfile != NULL)
    free(outfile);
  if (infile != NULL)
    free(infile);
  free(array);

  return 0;
}
