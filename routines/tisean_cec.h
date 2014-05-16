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
/*Author: Rainer Hegger Last modified: May 26, 2000*/
/*Changes: Bjoern Bastian Last modified: May 14, 2014 */

/* These definitions give the exit codes for the C part of the Tisean package.
   Typically the name is build up of, first, the name of the routine creating
   the exception, secondly, sort of an description of the exception.
   */

#ifndef _TISEAN_CEC_H
#define _TISEAN_CEC_H

/* These are the codes for the routines subtree */
#define RESCALE_DATA_ZERO_INTERVAL 11
#define CHECK_ALLOC_NOT_ENOUGH_MEMORY 12
#define CHECK_OPTION_NOT_UNSIGNED 13
#define CHECK_OPTION_NOT_INTEGER 14
#define CHECK_OPTION_NOT_FLOAT 15
#define CHECK_OPTION_NOT_TWO 16
#define CHECK_OPTION_C_NO_VALUE 17
#define TEST_OUTFILE_NO_WRITE_ACCESS 18
#define GET_SERIES_NO_LINES 20
#define GET_MULTI_SERIES_WRONG_TYPE_OF_C 21
#define GET_MULTI_SERIES_NO_LINES 22
#define VARIANCE_VAR_EQ_ZERO 23
#define CHECK_OPTION_NOT_THREE 25

/* These are the codes for the main routines */

#endif
