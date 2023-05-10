
EMO Project


Requierements:

  gcc 4.1.2 
  MPI 3.0, MPICH 2.2 (optional)
  GNU Make 3.81

For compilation:

  make clean
  make

If MPI is installed
  make mpi 

For running the programs:

  cd demo
  ./program

A help message will appear about the syntaxis.

Examples:

+ Execution of MOEAs:
  ./emo_moea SMS_EMOA input/Param_02D.cfg zdt1 1
  ./emo_moea NSGA2 input/Param_02D.cfg zdt1 1
  ./emo_moea NSGA3 input/Param_02D.cfg zdt1 1
  ./emo_moea MOMBI2 input/Param_02D.cfg zdt1 1
  ./emo_moea HV_IBEA input/Param_02D.cfg zdt1 1
  ./emo_moea EPS_IBEA input/Param_02D.cfg zdt1 1
  ./emo_moea R2_IBEA input/Param_02D.cfg zdt1 1
  ./emo_moea MOEAD input/Param_02D.cfg zdt1 1
  ./emo_moea SPEA2 input/Param_02D.cfg zdt1 1
  ./emo_moea HYPE input/Param_02D.cfg zdt1 1
  ./emo_moea MOVAP input/Param_02D.cfg zdt1 1

The output files are generated in demo/output.

+ Execution of parallel tasks:
  mpirun -np 3 ./emo_task cmd.txt

where cmd.txt contains:
  ./emo_moea NSGA2 input/Param_02D.cfg zdt1 1
  ./emo_moea NSGA3 input/Param_02D.cfg zdt1 1
  ./emo_moea MOMBI2 input/Param_02D.cfg zdt1 1

+ Filtering nondominated solutions:
  ./emo_ndset output/NSGA2_ZDT1_02D_R01.pof

It creates a filed named output/NSGA2_ZDT1_02D_R01.pof.nd

+ Performance indicators (HV, GD, IGD, IGD+, \Delta_p, R2, SP):
  ./emo_indicator HV output/NSGA2_ZDT1_02D_R01.pof.nd 1 1.2 1.2
  Running WFG program
  1 output/NSGA2_ZDT1_02D_R01.pof.nd 1.100411

+ Parallel execution of a MOEA using the island model:
  mpiexec -np 2 ./emo_pmoea NSGA2 input/Param_02D.cfg zdt1 1

The output files are generated in demo/output.



Support, comments and suggestions:

rhernandez@computacion.cs.cinvestav.mx

