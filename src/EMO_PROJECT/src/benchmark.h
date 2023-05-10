/**************************************************************
 * benchmark.h   Definition of test functions.                *
 *                                                            *
 * Computer Science Department, CINVESTAV-IPN                 *
 *                                                            *
 * Professor:   Dr. Carlos A. Coello Coello                   *
 * Author:      Raquel Hernandez Gomez                        *
 *                                                            *
 * March 2013                                                 *
 *************************************************************/
/*  Reference (chapter 4):

    Carlos A. Coello Coello, Gary B. Lamont, and David A. Van Veldhuizen
    Evolutionary Algorithms for Solving Multi-Objective Problems,
    Second edition, Springer, New York
    September, 2007,
    ISBN 978-0-387-33254-3
 */

#ifndef _BENCHMARK_
#define _BENCHMARK_

#include "common.h"
#include "param.h"

void EMO_Benchmark_alloc(EMO_MOP *mop, EMO_Param *param, const char *problem);
void EMO_Benchmark_free(EMO_MOP *mop);

#endif



