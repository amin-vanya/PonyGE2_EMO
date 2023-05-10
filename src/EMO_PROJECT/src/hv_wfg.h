/* hv_wfg.h  wfg algorithm for calculating hypervolume indicator (see \cite{While12})
 *
 * Lyndon While, Lucas Bradstreet, Luigi Barone 
 * Email lyndon.while@uwa.edu.au for any queries regarding usage/bugs/improvements. 
 *
 * This code includes a high performance implementation of the WFG algorithm, 
 * used to calculate the hypervolume indicator for a set of non-dominated points. 
 * 
 * COPYRIGHT: 
 * This program is free software (software libre); you can redistribute
 *  it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, you can obtain a copy of the GNU
 *  General Public License at:
 *                  http://www.gnu.org/copyleft/gpl.html
 *  or by writing to:
 *            Free Software Foundation, Inc., 59 Temple Place,
 *                  Suite 330, Boston, MA 02111-1307 USA
 * 
 * ----------------------------------------------------------------------
 * Updates
 * Original source code taken in april 2, 2016:
 * http://www.wfg.csse.uwa.edu.au/hypervolume/#code
 * version: metric WFG implementation (22 November 2015)
 * Raquel Hernandez Gomez
 *
 */

#ifndef _HV_WFG_H_
#define _HV_WFG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef double OBJECTIVE;

typedef struct
{
	OBJECTIVE *objectives;
        int idx;  // RHG, se usa en hv_iwfg.h
} POINT;

typedef struct
{
	int nPoints;
	int n;    // RHG, se usa en hv_iwfg.h
 	POINT *points;
} FRONT;

typedef struct {
  int col;      // the number of objectives 
  int row;
  FRONT *fs;    // memory management stuff 
  FRONT f;
  int fr;       // current depth 
  int safe;     // the number of points that don't need sorting 
} EMO_HV;


/*typedef struct
{
	int nFronts;
	FRONT *fronts;
} FILECONTENTS;

FILECONTENTS *readFile(char[]);

extern void printContents(FILECONTENTS *);*/

void EMO_HV_alloc(EMO_HV *hv, int max_row, int col);
void EMO_HV_free(EMO_HV *hv);
double EMO_HV_run(EMO_HV *hv, double *data, int *enable, int row, const double *ref);
double EMO_HV_run2(EMO_HV *hv, double *data, int *enable, int row, const double *ref);
double EMO_HV_contribution(EMO_HV *hv, double *deltahv, double *data, int *enable, int row, const double *ref, int col);

#endif

