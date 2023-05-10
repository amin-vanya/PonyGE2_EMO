# gcc -g main_dtlz.c -lm -o main_dtlz.o
###############################################################
#  Makefile     Compilation file.                             #
#                                                             #
#  Computer Science Department, CINVESTAV-IPN                 #
#                                                             #
#  Professor:   Dr. Carlos A. Coello Coello                   #
#  Author:      Raquel Hernandez Gomez                        #
#                                                             #
#  September 2014                                             #
###############################################################/


OBJS_NOMPI	= emo_moea emo_indicator emo_stat emo_ndset emo_ndsort 
OBJS_MPI	= emo_pmoea main_pamicro.o emo_task

LPATH	= ../lib
IPATH	= ../include
CC	= gcc
MPICC	= mpicc
#CFLAGS	= -Wall -g -L$(LPATH) -I$(IPATH)#-pedantic
CFLAGS	= -Wall -O3 -L$(LPATH) -I$(IPATH)
MPIDEF	= -DEMO_MPI

all:			$(OBJS_NOMPI)

mpi:			$(OBJS_MPI)

emo_moea:		emo_moea.c
			$(CC) $(CFLAGS) emo_moea.c -o emo_moea -lemo -lm

emo_task:		emo_task.c
			$(MPICC) $(CFLAGS) emo_task.c -o emo_task -lemo_mpi -lm

emo_indicator:		emo_indicator.c
			$(CC) $(CFLAGS) emo_indicator.c -o emo_indicator -lemo -lm

emo_stat:		emo_stat.c
			$(CC) $(CFLAGS) emo_stat.c -o emo_stat -lemo -lm

main_pamemo.o:		main_pamemo.c
			$(MPICC) $(CFLAGS) main_pamemo.c -o main_pamemo.o -lemo_mpi -lm

emo_pmoea:		emo_pmoea.c
			$(MPICC) $(CFLAGS) emo_pmoea.c -o emo_pmoea -lemo_mpi -lm

main_pamicro.o:		main_pamicro.c
			$(MPICC) $(CFLAGS) main_pamicro.c -o main_pamicro.o -lemo_mpi -lm

emo_ndset:		emo_ndset.c
			$(CC) $(CFLAGS) emo_ndset.c -o emo_ndset -lemo -lm

emo_ndsort:		emo_ndsort.c
			$(CC) $(CFLAGS) emo_ndsort.c -o emo_ndsort -lemo -lm


clean:		
		rm -f $(OBJS_NOMPI)  $(OBJS_MPI)

