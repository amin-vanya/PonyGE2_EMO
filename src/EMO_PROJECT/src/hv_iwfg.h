/* hv_iwfg.h  incremental IWFG algorithm (see \cite{Cox16})
 *
 * Lyndon While, Lucas Bradstreet, Wesley Cox 
 * Email lyndon.while@uwa.edu.au for any queries regarding usage/bugs/improvements. 
 *
 * This code includes a high performance implementation of the IWFG algorithm, 
 * used to identify the smallest hypervolume-contributor of a set of non-dominated points. 
 *
 * COPYRIGHT: 
 * This software is Copyright (C) 2015 Lyndon While, Lucas Bradstreet, Wesley Cox. 
 * This program is free software (software libre). You can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software Foundation; 
 * either version 2 of the License, or (at your option) any later version. 
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * See the GNU General Public License for more details. 
 *
 * Taken from version: IWFG_1.01 (24 November 2015)
 *
 * http://www.wfg.csse.uwa.edu.au/hypervolume/#code
 *
 * Updates:
 * Fixed two bugs:  binarySearch returns -1
 *                  different allocs for maxStackSize and not enough size
 *
 * 20 May 2017, Raquel Hernandez Gomez
 */

#ifndef _HV_IWFG_H_
#define _HV_IWFG_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "hv_wfg.h"

//typedef double OBJECTIVE;

//typedef struct
//{
//	OBJECTIVE *objectives;
//} POINT;

//typedef struct
//{
//	int	nPoints;
//	int	n;
//	POINT	*points;
//} FRONT;

typedef struct
{
	double width;
	FRONT front;
	int index;
} SLICE;

typedef struct {
	int	n;		// The number of objectives 
	// RHG POINT	ref;		// The refernce point	
	// RHG POINT	dirs;		// Records the directions of the objectives 

	FRONT	*fs;		// Memory manage stuff 
        FRONT f;
	int	fr;		// Current depth 
	int	maxm;		// Maximum number of points 
	int	maxn;		// Maximum number of objectives 
	int	safe;		// The number of points that don't need sorting

	double* partial;	// Partial exclusive hipervolumes
	int*	heap;		// Heap-based priority queue
	int	heapsize;	// Number of points in queue
	SLICE	**stacks;	// Set of slices per point per slicing depth
	int	*stacksize;	// Current slicing depth per point

	int*	gorder;		// Objective order used by comparison funtions
	int**	torder;		// Order of objectives per point
	int**	tcompare;	
	FRONT*	fsorted;	// Front sorted in each objective
        int maxStackSize;       // Maximum stack size
} EMO_IWFG;
 
void EMO_IWFG_alloc(EMO_IWFG *hv, int maxm, int maxn);
void EMO_IWFG_free(EMO_IWFG *hv);
int EMO_IWFG_run(EMO_IWFG *hv, double *data, int *enable, int row, const double *ref, double *hv_worst);

#endif
