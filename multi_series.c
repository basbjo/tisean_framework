/*
 *   This file is derived from rescale.c of TISEAN,
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
/*Author: Rainer Hegger. Last modified: Nov 23, 2000
 * Adapted by Bjoern Bastian. Last modified: May 16, 2014 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include "routines/tsa.h"

#define WID_STR "Rescales the data"

unsigned long length=ULONG_MAX,exclude=0;
unsigned int dim=1;
unsigned int verbosity=0xff;
char *columns=NULL,dimset=0;
char *outfile=NULL,stout=1;
double **series;
double xmin=0.0,xmax=1.0;
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
  fprintf(stderr,"\t-m # of components to be read [default %u]\n",dim);
  fprintf(stderr,"\t-c columns to be read [default 1,...,# of components]\n");
  fprintf(stderr,"\t-z minimum of the new series [default 0.0]\n");
  fprintf(stderr,"\t-Z maximum of the new series [default 1.0]\n");
  fprintf(stderr,"\t-o output file name [default 'datafile'.res; no -o"
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
  if ((out=check_option(argv,argc,'m','u')) != NULL) {
    sscanf(out,"%u",&dim);
    dimset=1;
  }
  if ((out=check_option(argv,argc,'c','s')) != NULL)
    columns=out;
  if ((out=check_option(argv,argc,'V','u')) != NULL)
    sscanf(out,"%u",&verbosity);
  if ((out=check_option(argv,argc,'z','f')) != NULL)
    sscanf(out,"%lf",&xmin);
  if ((out=check_option(argv,argc,'Z','f')) != NULL)
    sscanf(out,"%lf",&xmax);
  if ((out=check_option(argv,argc,'o','o')) != NULL) {
    stout=0;
    if (strlen(out) > 0)
      outfile=out;
  }
}

int main(int argc,char** argv)
{
  char stdi=0;
  long i,n;
  FILE *fout=NULL;
  double min,max;
  double av,varianz;

  if (scan_help(argc,argv))
    show_options(argv[0]);

  scan_options(argc,argv);
#ifndef OMIT_WHAT_I_DO
  if (verbosity&VER_INPUT)
    fprintf(stderr, "\n%s: %s\n\n",argv[0],WID_STR);
#endif

  infile=search_datafile(argc,argv,NULL,verbosity);
  if (infile == NULL)
    stdi=1;

  if (outfile == NULL) {
    if (!stdi) {
      check_alloc(outfile=(char*)calloc(strlen(infile)+5,(size_t)1));
      strcpy(outfile,infile);
      strcat(outfile,".res");
    }
    else {
      check_alloc(outfile=(char*)calloc((size_t)10,(size_t)1));
      strcpy(outfile,"stdin.res");
    }
  }
  if (!stout)
    test_outfile(outfile);

  if (xmin >= xmax) {
    fprintf(stderr,"Choosing the minimum larger or equal the maximum\n"
        "makes no sense. Exiting!\n");
    exit(RESCALE__WRONG_INTERVAL);
  }

  if (columns == NULL)
    series=(double**)get_multi_series(infile,&length,exclude,&dim,"",dimset,
                    verbosity);
  else
    series=(double**)get_multi_series(infile,&length,exclude,&dim,columns,
                    dimset,verbosity);

  for (n=0;n<dim;n++) {
    variance(series[n],length,&av,&varianz);
  }

  if (!stout) {
    fout=fopen(outfile,"w");
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Opened %s for writing\n",outfile);
    for (i=0;i<length;i++) {
      fprintf(fout,"%e",series[0][i]);
      for (n=1;n<dim;n++)
    fprintf(fout," %e",series[n][i]);
      fprintf(fout,"\n");
    }
    fclose(fout);
  }
  else {
    if (verbosity&VER_INPUT)
      fprintf(stderr,"Writing to stdout\n");
    for (i=0;i<length;i++) {
      fprintf(stdout,"%e",series[0][i]);
      for (n=1;n<dim;n++)
    fprintf(stdout," %e",series[n][i]);
      fprintf(stdout,"\n");
    }
  }

  return 0;
}
