
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

#OBJS_BASE	= io.o parser.o list.o avl.o random.o stat.o numeric.o sort.o vector.o matrix.o common.o evop.o benchmark.o lame.o wfg.o cec09.o knapsack.o dominance.o niche.o utility.o refpoint.o hv.o hv_wfg.o hv_iwfg.o hv_sort.o indicator.o plot.o topology.o migration.o nsga2.o spea2.o moea.o hype.o mombi2.o mombi3.o vpath.o
OBJS_BASE	= io.o parser.o list.o avl.o random.o stat.o numeric.o sort.o vector.o matrix.o common.o evop.o benchmark.o lame.o wfg.o cec09.o dominance.o niche.o utility.o refpoint.o hv.o hv_wfg.o hv_iwfg.o hv_sort.o indicator.o plot.o topology.o migration.o nsga2.o spea2.o moea.o hype.o mombi2.o mombi3.o vpath.o
OBJS_NOMPI	= stop.o debug.o param.o ibea.o smsemoa.o moead.o nsga3.o movap.o
OBJS_MPI	= task.o debug_mpi.o stop_mpi.o param_mpi.o ibea_mpi.o smsemoa_mpi.o moead_mpi.o nsga3_mpi.o movap_mpi.o island.o pmoea.o pamicro.o

LPATH	= ../lib
CC	= gcc
MPICC	= mpicc
AR	= ar
PROG	= libemo.a
PPROG	= libemo_mpi.a
#CFLAGS	= -Wall -g #-pedantic #-DEMO_MPI
CFLAGS	= -Wall -O3
MPIDEF	= -DEMO_MPI

all:			$(OBJS_BASE) $(OBJS_NOMPI)
			$(AR) -cvr $(PROG) $(OBJS_BASE) $(OBJS_NOMPI)
			mv $(PROG) $(LPATH)

mpi:			$(OBJS_BASE) $(OBJS_MPI)
			$(AR) -cvr $(PPROG) $(OBJS_BASE) $(OBJS_MPI)
			mv $(PPROG) $(LPATH)

io.o:			io.h io.c
			$(CC) $(CFLAGS) -c io.c

parser.o:		parser.h parser.c
			$(CC) $(CFLAGS) -c parser.c

param.o:		param.h param.c
			$(CC) $(CFLAGS) -c param.c

param_mpi.o:		param.h param.c
			$(MPICC) $(CFLAGS) $(MPIDEF) -c param.c -o param_mpi.o

list.o:			list.h list.c
			$(CC) $(CFLAGS) -c list.c

avl.o:			avl.h avl.c
			$(CC) $(CFLAGS) -c avl.c

#random.o:		random2.h random2.c
#			$(CC) $(CFLAGS) -c random2.c -o random.o

random.o:		random.h random.c
			$(CC) $(CFLAGS) -c random.c

stat.o:			stat.h  stat.c
			$(CC) $(CFLAGS) -c stat.c

numeric.o:		numeric.h  numeric.c
			$(CC) $(CFLAGS) -c numeric.c

sort.o:			sort.h  sort.c
			$(CC) $(CFLAGS) -c sort.c

vector.o:		vector.h  vector.c
			$(CC) $(CFLAGS) -c vector.c

matrix.o:		matrix.h  matrix.c
			$(CC) $(CFLAGS) -c matrix.c

common.o:		common.h common.c
			$(CC) $(CFLAGS) -c common.c

moea.o:			moea.h moea.c
			$(CC) $(CFLAGS) -c moea.c

evop.o:			evop.h evop.c
			$(CC) $(CFLAGS) -c evop.c

benchmark.o:		benchmark.h benchmark.c
			$(CC) $(CFLAGS) -c benchmark.c

lame.o:			lame.h lame.c
			$(CC) $(CFLAGS) -c lame.c

wfg.o:			wfg.h wfg.c
			$(CC) $(CFLAGS) -c wfg.c

cec09.o:		cec09.h cec09.c
			$(CC) $(CFLAGS) -c cec09.c

knapsack.o:		knapsack.h knapsack.c
			$(CC) $(CFLAGS) -c knapsack.c

dominance.o:		dominance.h dominance.c list.h sort.h
			$(CC) $(CFLAGS) -c dominance.c

niche.o:		niche.h niche.c vector.h sort.h
			$(CC) $(CFLAGS) -c niche.c

refpoint.o:		refpoint.h refpoint.c
			$(CC) $(CFLAGS) -c refpoint.c

utility.o:		utility.h utility.c
			$(CC) $(CFLAGS) -c utility.c

hv.o:			hv.h hv.c avl.h
			$(CC) $(CFLAGS) -c hv.c

hv_wfg.o:		hv_wfg.h hv_wfg.c
			$(CC) $(CFLAGS) -c hv_wfg.c

hv_iwfg.o:		hv_iwfg.h hv_iwfg.c
			$(CC) $(CFLAGS) -c hv_iwfg.c

hv_sort.o:		hv_sort.h hv_sort.c
			$(CC) $(CFLAGS) -c hv_sort.c

indicator.o:		indicator.h indicator.c
			$(CC) $(CFLAGS) -c indicator.c

vpath.o:		vpath.h vpath.c
			$(CC) $(CFLAGS) -c vpath.c

vpath_mpi.o:		vpath.h vpath.c
			$(MPICC) $(CFLAGS) $(MPIDEF) -c vpath.c -o vpath_mpi.o

plot.o:			plot.h plot.c
			$(CC) $(CFLAGS) -c plot.c

topology.o:		topology.h topology.c
			$(CC) $(CFLAGS) -c topology.c

migration.o:		migration.h migration.c
			$(CC) $(CFLAGS) -c migration.c

debug.o:		debug.h debug.c
			$(CC) $(CFLAGS) -c debug.c

debug_mpi.o:		debug.h debug.c
			$(MPICC) $(CFLAGS) $(MPIDEF) -c debug.c -o debug_mpi.o

task.o:			task.h task.c
			$(MPICC) $(CFLAGS) $(MPIDEF) -c task.c -o task.o

stop.o:			stop.h stop.c
			$(CC) $(CFLAGS) -c stop.c

stop_mpi.o:		stop.h stop.c
			$(MPICC) $(CFLAGS) $(MPIDEF) -c stop.c -o stop_mpi.o

smsemoa.o:		smsemoa.h smsemoa.c
			$(CC) $(CFLAGS) -c smsemoa.c

smsemoa_mpi.o:		smsemoa.h smsemoa.c
			$(MPICC) $(CFLAGS) $(MPIDEF) -c smsemoa.c -o smsemoa_mpi.o

movap.o:		movap.h movap.c
			$(CC) $(CFLAGS) -c movap.c

movap_mpi.o:		movap.h movap.c
			$(MPICC) $(CFLAGS) $(MPIDEF) -c movap.c -o movap_mpi.o

hype.o:			hype.h hype.c
			$(CC) $(CFLAGS) -c hype.c

nsga2.o:		nsga2.h nsga2.c
			$(CC) $(CFLAGS) -c nsga2.c

mombi2.o:		mombi2.h mombi2.c
			$(CC) $(CFLAGS) -c mombi2.c

mombi3.o:		mombi3.h mombi3.c
			$(CC) $(CFLAGS) -c mombi3.c

ibea.o:			ibea.h ibea.c
			$(CC) $(CFLAGS) -c ibea.c

ibea_mpi.o:		ibea.h ibea.c
			$(MPICC) $(CFLAGS) $(MPIDEF) -c ibea.c -o ibea_mpi.o

moead.o:		moead.h moead.c
			$(CC) $(CFLAGS) -c moead.c

moead_mpi.o:		moead.h moead.c
			$(MPICC) $(CFLAGS) $(MPIDEF) -c moead.c -o moead_mpi.o

nsga3.o:		nsga3.h nsga3.c
			$(CC) $(CFLAGS) -c nsga3.c

nsga3_mpi.o:		nsga3.h nsga3.c
			$(MPICC) $(CFLAGS) $(MPIDEF) -c nsga3.c -o nsga3_mpi.o

spea2.o:		spea2.h spea2.c
			$(CC) $(CFLAGS) -c spea2.c

island.o:		island.h island.c
			$(MPICC) $(CFLAGS) -c island.c

pmoea.o:		pmoea.h pmoea.c
			$(MPICC) $(CFLAGS) -c pmoea.c

pamicro.o:		pamicro.h pamicro.c
			$(MPICC) $(CFLAGS) -c pamicro.c

clean:
	rm -f $(LPATH)/$(PROG)
	rm -f $(OBJS_BASE) $(OBJS_NOMPI) $(OBJS_MPI)
