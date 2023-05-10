#python ponyge.py --parameters mombi_sf.txt
from fitness.base_ff_classes.base_ff import base_ff
#from subprocess import call
import subprocess
# to check if best functions file exists
from os import path

from algorithm.parameters import params

class scalarizingFunction_multipleDims(base_ff):
    """
    Fitness function for generating new scalarizing functions used in MOMBI-II
    to solve multi-objective optimization problems.

    Example execution: python ponyge.py --parameters ge_sf.txt --sf_problems DTLZ1,WFG1 --random_seed 1
    """

    # Fitness is directly based on the average hypervolume obtained, hence it should be maximized
    maximise = True

    def __init__(self):
        # Initialize base fitness function class.
        super().__init__()

        self.MAX_FITNESS_FEW = -1;
        self.MAX_FITNESS_EXT = -1;

        self.NMULTIPLEDIMS = 3
        self.NUM_DIMENSIONS = [2, 3, 5]
        self.TXT_DIMENSIONS = ["02D", "03D", "05D"]
        self.TOTAL_N_FEW = ["20", "15", "10"]
        self.TOTAL_N_EXT = ["30", "30", "30"]
        self.PARAM_FILES_FEW = ["input/Param_02D_50k.cfg", "input/Param_03D_25k.cfg", "input/Param_05D_25k.cfg"]
        self.PARAM_FILES_EXT = ["input/Param_02D_100k.cfg", "input/Param_03D_100k.cfg", "input/Param_05D_100k.cfg"]
        self.PROBLEMS = []
        self.AVG_HV = []
        self.REFERENCE_POINTS = []
        self.AVG_HVS = [0.0,0.0,0.0]

        self.TEST_TXT_DIM = self.TXT_DIMENSIONS[0]
        self.TEST_NUM_DIM = self.NUM_DIMENSIONS[0]
        self.TEST_PARAM_FILE = "input/Param_02D_25k.cfg"
        self.TEST_TOTAL_N = "5"

        self.EMO_DIR_FULL = "EMO_PROJECT/"
        self.EMO_LINE_1 = 804
        self.EMO_LINE_2 = 810
        self.THRESHOLD = 0.5

        reference_points = {
            "DTLZ1_02D":[1.6,1.6],"DTLZ2_02D":[2,2],"DTLZ3_02D":[2.1,2.1],"DTLZ4_02D":[2,2],"DTLZ5_02D":[2,2],"DTLZ6_02D":[2.5,2.5],"DTLZ7_02D":[1.9,5],"WFG1_02D":[3.7,4.1],"WFG2_02D":[2.4,5.1],"WFG3_02D":[3.1,5.1],"WFG4_02D":[3.1,5.1],"WFG5_02D":[3.1,5.1],"WFG6_02D":[3.1,5.1],"WFG7_02D":[3.1,5.3],"WFG8_02D":[3.3,5.2],"WFG9_02D":[3.1,5.1],"DTLZ1_03D":[1.5,1.5,1.5],"DTLZ2_03D":[2,2,2],"DTLZ3_03D":[2.1,2.1,2.1],"DTLZ4_03D":[2,2,2],"DTLZ5_03D":[1.8,1.8,2],"DTLZ6_03D":[2,2,2.2],"DTLZ7_03D":[1.9,1.9,7],"WFG1_03D":[3.5,5.6,6.5],"WFG2_03D":[2.8,4.6,6.9],"WFG3_03D":[4,3.1,7.1],"WFG4_03D":[3.1,5.1,7.1],"WFG5_03D":[3.1,5.1,7.1],"WFG6_03D":[3.1,5.1,7.1],"WFG7_03D":[3.1,5.1,7.1],"WFG8_03D":[3.3,5.1,7.1],"WFG9_03D":[3.2,5.1,7.1],"DTLZ1_04D":[1.5,1.5,1.5,1.5],"DTLZ2_04D":[2,2,2,2],"DTLZ3_04D":[2.1,2.1,2.1,2.1],"DTLZ4_04D":[2,2,2,2],"DTLZ5_04D":[1.8,1.8,4.5,2],"DTLZ6_04D":[1.8,1.8,11.5,2.2],"DTLZ7_04D":[1.9,1.9,1.9,9],"WFG1_04D":[3.5,3.6,5,5],"WFG2_04D":[2.6,4.3,6,8.8],"WFG3_04D":[3.6,4.7,5.2,9],"WFG4_04D":[3.1,5.1,7.1,9.1],"WFG5_04D":[3.1,5.1,7.1,9.1],"WFG6_04D":[3.1,5.1,7.1,9.1],"WFG7_04D":[3.1,5.1,7.1,9.1],"WFG8_04D":[3.6,5.2,7.1,9.1],"WFG9_04D":[3.2,5.2,7.1,9.1],"DTLZ1_05D":[1.5,1.5,1.5,1.5,1.5],"DTLZ2_05D":[2,2,2,2,2],"DTLZ3_05D":[2.1,2.1,2.1,2.1,2.1],"DTLZ4_05D":[2,2,2,2,2],"DTLZ5_05D":[1.8,1.8,4.4,4.5,2],"DTLZ6_05D":[4,3.5,8.5,9.7,2.8],"DTLZ7_05D":[1.9,1.9,2,2,11.8],"WFG1_05D":[3.6,2.2,2.7,2.9,5.6],"WFG2_05D":[2.5,4,5.5,7.2,9.3],"WFG3_05D":[3.9,5.3,7.6,8.6,11.1],"WFG4_05D":[3.1,5.1,7.1,9.1,11.1],"WFG5_05D":[3.1,5.1,7.1,9.1,11.1],"WFG6_05D":[3.1,5.1,7.1,9.1,11.1],"WFG7_05D":[3.1,5.1,7.1,9.1,11.1],"WFG8_05D":[4,5.5,7.1,9.1,11.1],"WFG9_05D":[3.2,5.2,7.2,9.2,11.2],"DTLZ1_06D":[1.5,1.5,1.5,1.5,1.5,1.5],"DTLZ2_06D":[2,2.1,2.1,2,2.1,2.1],"DTLZ3_06D":[2.1,2.1,2.1,2.1,2.1,2.1],"DTLZ4_06D":[2,2.1,2,2.1,2.1,2.1],"DTLZ5_06D":[1.8,1.8,1.9,4.4,4.5,2.1],"DTLZ6_06D":[8.3,8.4,8.4,8.5,9.2,3.2],"DTLZ7_06D":[1.9,2,2,2,1.9,13.2],"WFG1_06D":[3.5,2.1,2.5,2.5,3,6.5],"WFG2_06D":[2.4,3.8,5.5,6.9,8.6,13.1],"WFG3_06D":[3.9,5.4,7.6,9.6,10.8,12.9],"WFG4_06D":[3.1,5.1,7.1,9.1,11.1,13.1],"WFG5_06D":[3.1,5.1,7.1,9.1,11.1,13.1],"WFG6_06D":[3.1,5.1,7.1,9.1,11.1,13.1],"WFG7_06D":[3.1,5.1,7.1,9.1,11.1,13.1],"WFG8_06D":[4,5.9,7.3,9.1,11.1,13.1],"WFG9_06D":[3.2,5.2,7.2,9.2,11.2,13.2],"DTLZ1_07D":[1.5,1.5,1.5,1.5,1.5,1.5,1.5],"DTLZ2_07D":[2,2,2,2,2,2,2],"DTLZ3_07D":[2.1,2.1,2.1,2.1,2.1,2.1,2.1],"DTLZ4_07D":[2,2,2,2,2,2,2],"DTLZ5_07D":[1.8,1.8,1.8,1.8,3.7,4.5,2.1],"DTLZ6_07D":[7.3,3.2,3.9,8,4,8.9,3.1],"DTLZ7_07D":[2,1.9,1.9,2,1.9,2,15.3],"WFG1_07D":[3.4,2.2,2.5,2.7,2.7,3.6,7.3],"WFG2_07D":[2.3,3.5,5.4,6.2,8.1,9.8,14.9],"WFG3_07D":[3.8,5.7,7.7,9.6,11.7,12.7,15],"WFG4_07D":[3.1,5.1,7.1,9.1,11.1,13.1,15.1],"WFG5_07D":[3.1,5.1,7.1,9.1,11.1,13.1,15.1],"WFG6_07D":[3.1,5.1,7.1,9.1,11.1,13.1,15.1],"WFG7_07D":[3.1,5.1,7.1,9.1,11.1,13.1,15.1],"WFG8_07D":[4,6,7.9,9.1,11.1,13.1,15.1],"WFG9_07D":[3.2,5.2,7.2,9.2,11.2,13.2,15.2]
        }

        test_avg_hypervolumes = {
            "DTLZ1_02D":2.4326586,"DTLZ2_02D":3.209867,"DTLZ3_02D":2.11129528,"DTLZ4_02D":3.2100918,"DTLZ5_02D":3.209867,"DTLZ6_02D":4.6621132,"DTLZ7_02D":4.146444,"WFG1_02D":4.723919,"WFG2_02D":7.3807536,"WFG3_02D":10.74322,"WFG4_02D":8.0314748,"WFG5_02D":8.0440914,"WFG6_02D":8.297435,"WFG7_02D":7.3033126,"WFG8_02D":7.7019586,"WFG9_02D":8.6084454
        }
        few_avg_hypervolumes = {
            "DTLZ1_02D":2.43359825,"DTLZ2_02D":3.2101245,"DTLZ3_02D":3.6097688,"DTLZ4_02D":3.21012195,"DTLZ5_02D":3.2101245,"DTLZ6_02D":5.25184545,"DTLZ7_02D":4.14773,"WFG1_02D":4.9933269,"WFG2_02D":7.4126962,"WFG3_02D":10.9427225,"WFG4_02D":8.480721,"WFG5_02D":8.39357855,"WFG6_02D":8.59799495,"WFG7_02D":7.5376156,"WFG8_02D":8.2369505,"WFG9_02D":8.83736125,"DTLZ1_03D":3.3380598,"DTLZ2_03D":7.393965,"DTLZ3_03D":8.13361173333333,"DTLZ4_03D":7.40960746666667,"DTLZ5_03D":4.73071513333333,"DTLZ6_03D":4.653516,"DTLZ7_03D":11.8851986666667,"WFG1_03D":57.348802,"WFG2_03D":79.3938706666667,"WFG3_03D":57.7790786666667,"WFG4_03D":81.5478913333333,"WFG5_03D":78.0861633333334,"WFG6_03D":78.9473766666667,"WFG7_03D":82.6069726666666,"WFG8_03D":81.3396066666667,"WFG9_03D":79.1981666666667,"DTLZ1_05D":7.5900398,"DTLZ2_05D":31.661133,"DTLZ3_05D":2.35336827,"DTLZ4_05D":31.682748,"DTLZ5_05D":106.89187,"DTLZ6_05D":2617.2881,"DTLZ7_05D":82.715266,"WFG1_05D":58.050499,"WFG2_05D":3401.5117,"WFG3_05D":9912.6035,"WFG4_05D":9132.4202,"WFG5_05D":9042.7801,"WFG6_05D":9174.8116,"WFG7_05D":9642.8794,"WFG8_05D":10856.824,"WFG9_05D":9018.4258
        }
        ext_avg_hypervolumes = {
            "DTLZ1_02D":2.43370673333333,"DTLZ2_02D":3.21012963333333,"DTLZ3_02D":3.61879233333333,"DTLZ4_02D":3.21012953333333,"DTLZ5_02D":3.21012963333333,"DTLZ6_02D":5.2603026,"DTLZ7_02D":4.14769643333333,"WFG1_02D":5.30806463333333,"WFG2_02D":7.44443453333333,"WFG3_02D":11.31074,"WFG4_02D":8.89683666666667,"WFG5_02D":8.6680239,"WFG6_02D":8.82691336666667,"WFG7_02D":7.8360852,"WFG8_02D":8.8897825,"WFG9_02D":8.944584,"DTLZ1_03D":3.32222066666667,"DTLZ2_03D":7.35335696666667,"DTLZ3_03D":8.63676886666666,"DTLZ4_03D":7.37043966666667,"DTLZ5_03D":4.7312296,"DTLZ6_03D":6.52167766666667,"DTLZ7_03D":11.878287,"WFG1_03D":71.6920733333333,"WFG2_03D":80.6450653333333,"WFG3_03D":59.0454706666667,"WFG4_03D":83.1890413333334,"WFG5_03D":79.5154213333333,"WFG6_03D":80.542609,"WFG7_03D":83.5000383333333,"WFG8_03D":83.4252366666667,"WFG9_03D":79.928651,"DTLZ1_05D":7.5651614,"DTLZ2_05D":31.5983703333333,"DTLZ3_05D":40.4840766666667,"DTLZ4_05D":31.685935,"DTLZ5_05D":104.979836666667,"DTLZ6_05D":2997.45546666667,"DTLZ7_05D":86.182004,"WFG1_05D":81.055235,"WFG2_05D":3492.88923333333,"WFG3_05D":10353.8539,"WFG4_05D":9901.8265,"WFG5_05D":9437.53446666667,"WFG6_05D":9506.95406666667,"WFG7_05D":9989.46146666667,"WFG8_05D":11786.0653333333,"WFG9_05D":9195.0749
        }

        # Read custom parameter specifying which problems are to be used in the GE search
        try:
            custom_param = params['SF_PROBLEMS']
            self.SF_PROBLEMS = custom_param.split(",")
            print("Problem(s) to be used:", self.SF_PROBLEMS)

            n_problems = len(self.SF_PROBLEMS)
            # Adjust auxiliar arrays for multiple problems
            if(n_problems>1):
                for i in range(n_problems-1):
                    for j in range(self.NMULTIPLEDIMS):
                        self.NUM_DIMENSIONS.append(self.NUM_DIMENSIONS[j])
                        self.TXT_DIMENSIONS.append(self.TXT_DIMENSIONS[j])
                        self.TOTAL_N_FEW.append(self.TOTAL_N_FEW[j])
                        self.TOTAL_N_EXT.append(self.TOTAL_N_EXT[j])
                        self.PARAM_FILES_FEW.append(self.PARAM_FILES_FEW[j])
                        self.PARAM_FILES_EXT.append(self.PARAM_FILES_EXT[j])
                        self.AVG_HVS.append(self.AVG_HVS[j])
            try:
                for i in range(n_problems):
                    for j in range(self.NMULTIPLEDIMS):
                        self.PROBLEMS.append( self.SF_PROBLEMS[i] )

                        aux_key = self.SF_PROBLEMS[i]+"_"+self.TXT_DIMENSIONS[j]
                        self.AVG_HV.append(ext_avg_hypervolumes[aux_key])
                        self.REFERENCE_POINTS.append(reference_points[aux_key])

                # Define test problem parameters
                aux_key = self.SF_PROBLEMS[0]+"_"+self.TXT_DIMENSIONS[0]
                self.TEST_PROBLEM = self.SF_PROBLEMS[0]
                self.TEST_AVG_HV = test_avg_hypervolumes[aux_key]
                self.TEST_REFPOINT = reference_points[aux_key]
            except (KeyError):
                print("Error: Unknown problem specified.")
                exit(1)

        except (KeyError):
            print("No problem was specified. Default DTLZ4 will be used.")
            self.SF_PROBLEMS = ['DTLZ4']
            self.PROBLEMS = self.SF_PROBLEMS * self.NMULTIPLEDIMS
            try:
                for i in range(self.NMULTIPLEDIMS):
                    aux_key = self.SF_PROBLEMS[0]+"_"+self.TXT_DIMENSIONS[i]
                    self.AVG_HV.append(ext_avg_hypervolumes[aux_key])
                    self.REFERENCE_POINTS.append(reference_points[aux_key])
                # Define test problem parameters
                aux_key = self.SF_PROBLEMS[0]+"_"+self.TXT_DIMENSIONS[0]
                self.TEST_PROBLEM = self.SF_PROBLEMS[0]
                self.TEST_AVG_HV = test_avg_hypervolumes[aux_key]
                self.TEST_REFPOINT = reference_points[aux_key]
            except (KeyError):
                print("Error: Unknown problem specified.")
                exit(1)

        # Create best functions log file
        self.FINAL_LOG_FILENAME = path.join(params['FILE_PATH'], "best_functions_log.txt")
        with open(self.FINAL_LOG_FILENAME, "w") as log_file:
            header_str = "function\tfitness\t"

            for p in self.SF_PROBLEMS:
                for i in range(self.NMULTIPLEDIMS):
                    header_str += p+'_'+self.TXT_DIMENSIONS[i]+'(HV)\t'
            header_str += '\n'

            log_file.write(header_str)

    # splits phenotype in two parts
    # it also checks if x and y variables are present in the first part
    # when both variables are present, flag = 0
    # when either variable is missing, flag = 1
    def process_scalarizing_function(self, phenotype):
        flag = 0
        sf_p1 = ''
        sf_p2 = ''

        # Process individual phenotype to split it in two parts
        utility_func = phenotype

        # Assuming first four characters are 'max{', store the substring inside operator max{}
        aux_index = utility_func.find('}')
        sf_p1 = utility_func[4:aux_index]
        if('x[i]' in sf_p1 and 'y' in sf_p1):
            print('\tpart1: '+sf_p1)
            sf_p1 = '\tv = ' + sf_p1 + ';\n'

            # Store the substring after operator max{}
            sf_p2 = utility_func[aux_index+1:len(utility_func)]
            print('\tpart2: '+sf_p2)
            sf_p2 = '\t\tvmax = vmax' + sf_p2 + ';\n'
        else:
            flag = 1

        return flag, sf_p1, sf_p2

    def recompile_emo_project(self, sf_p1, sf_p2, emo_dir, emo_line1, emo_line2):
        # Replace utility function in utility.c
        # Open utility file
        aux_file = open(emo_dir+"src/utility.c",'r')
        aux_lines = aux_file.readlines()
        aux_file.close()
        # Replace utility function line (inside loop)
        aux_lines[emo_line1] = sf_p1
        # Replace utility function line (outside loop)
        aux_lines[emo_line2] = sf_p2
        # Write file again
        aux_file = open(emo_dir+"src/utility.c",'w')
        aux_file.writelines(aux_lines)
        aux_file.close()

        # Re-compile emo project
        process1 = subprocess.Popen(["make", "clean"], cwd=emo_dir, stdout=subprocess.DEVNULL)
        process1.wait()
        process2 = subprocess.Popen("make", cwd=emo_dir, stdout=subprocess.DEVNULL)
        process2.wait()

    def perform_moea_executions(self, parameters, problem, n, emo_dir):
        flag = 0
        try:
            # Execute moea using the cmd line:
            #./emo_moea MOEA PARAMETERS PROBLEM N_RUNS
            # For instance:
            # ./emo_moea MOMBI2 input/Param_02D.cfg DTLZ1 15
            tmp=subprocess.check_output(["./emo_moea", "MOMBI2", parameters, problem, n], cwd=emo_dir+"/demo")
        except subprocess.CalledProcessError as e:
            flag = 1
            print("signal error: "+str(e.output))

        return flag

    def get_avg_hv(self, problem, dimensions, n, refpoint, emo_dir):
        avg_hv_val = 0

        # Execute hypervolume calculation using the cmd line:
        #./emo_indicator HV output/TEST_PROBLEM N REFPOINT REFPOINT WFG
        # For instance:
        #./emo_indicator HV output/MOMBI2_DTLZ1 10 1 1 WFG
        cmd_instruction = ["./emo_indicator", "HV", "output/MOMBI2_"+problem+"_"+dimensions, n]
        # append reference point to instruction
        for item in refpoint:
            cmd_instruction.append(str(item))
        cmd_instruction.append("WFG")
        subprocess.call(cmd_instruction,stdout=subprocess.DEVNULL, cwd=emo_dir+"/demo")

        # Average resulting hypervolumes
        hv_filename = emo_dir+"demo/output/MOMBI2_"+problem+"_"+dimensions+".hv"
        file = open(hv_filename,'r')
        lines = file.readlines()
        file.close()

        # Retrieve the number of lines from the second token from the first line
        file_size = int(lines[0].split()[1])

        for i in range(file_size):
            line = lines[i+1].split()
            current_hv = float(line[0])
            #print("#"+str(i)+" "+str(current_hv))
            avg_hv_val += current_hv

        avg_hv_val /= file_size
        return avg_hv_val

    def get_avg_se(self, problem, dimensions, numeric_dimension, n, emo_dir):
        avg_se_val = 0

        # Execute hypervolume calculation using the cmd line:
        #./main_se.o output/MOEA_PROBLEM_DIM N_RUNS S_PARAM
        # For instance:
        #./main_se.o output/MOMBI2_DTLZ1_02D 15 1
        cmd_instruction = ["./main_se.o", "output/MOMBI2_"+problem+"_"+dimensions, n, str(numeric_dimension-1)]
        # append reference point to instruction
        subprocess.call(cmd_instruction,stdout=subprocess.DEVNULL, cwd=emo_dir+"/demo")

        # Average resulting hypervolumes
        se_filename = emo_dir+"demo/output/MOMBI2_"+problem+"_"+dimensions+".se"
        file = open(se_filename,'r')
        lines = file.readlines()
        file.close()

        # Retrieve the number of lines from the second token from the first line
        file_size = int(lines[0].split()[1])

        for i in range(file_size):
            line = lines[i+1].split()
            current_se = float(line[0])
            #print("#"+str(i)+" "+str(current_se))
            avg_se_val += current_se

        avg_se_val /= file_size
        return avg_se_val

    def evaluate(self, ind, **kwargs):
        """
        :param ind: An individual to be evaluated.
        :param kwargs: Optional extra arguments.
        :return: The fitness of the evaluated individual.
        """

        print("f(x,y):\t"+ind.phenotype)
        error_flag = 0

        ################################################################################
        # Obtain final scalarizing function from phenotype
        ################################################################################
        flag, sf_p1, sf_p2 = self.process_scalarizing_function(ind.phenotype)

        if(flag):
            fitness = 0
            print("\tError: " + ind.phenotype + " does not contain x and y variables in max operator")
            print("\tFinal fitness: "+str(fitness)+"\n")
            return fitness

        ################################################################################
        # Recompile EMO_PROJECT using new scalarizing function
        ################################################################################
        self.recompile_emo_project(sf_p1, sf_p2, self.EMO_DIR_FULL, self.EMO_LINE_1, self.EMO_LINE_2)

        ################################################################################
        # Execute modified EMO_PROJECT in test problem, with a small amount of executions and few function evaluations
        ################################################################################
        flag = self.perform_moea_executions(self.TEST_PARAM_FILE, self.TEST_PROBLEM, self.TEST_TOTAL_N, self.EMO_DIR_FULL)

        if(flag):
            fitness = 0
            print("\n\t\tError executing EMO_PROJECT")
            print("\t\tFinal fitness: "+str(fitness)+"\n")
            return fitness

        ################################################################################
        # Obtain hypervolume values of the resulting pareto fronts
        ################################################################################
        avg_hv = self.get_avg_hv(self.TEST_PROBLEM, self.TEST_TXT_DIM, self.TEST_TOTAL_N, self.TEST_REFPOINT, self.EMO_DIR_FULL)

        #Scale current average hypervolume
        #current_fitness = avg_hv/self.TEST_AVG_HV
        current_fitness = (avg_hv-0.5*self.TEST_AVG_HV) / self.TEST_AVG_HV

        # Check if function is worth more evaluations
        if(avg_hv < self.TEST_AVG_HV * self.THRESHOLD):
        	fitness = current_fitness
        	print("Final_fitness:\t"+str(fitness)+" (only "+str(self.TEST_TOTAL_N)+" evals done)")
        else:
            # Reset fitness as more evaluations are about to be performed
            fitness = 0

            for index in range(len(self.PROBLEMS)):

                ################################################################################
                # Execute modified EMO_PROJECT in every problem,
                # with a decent amount of executions and few function evaluations
                ################################################################################
                flag = self.perform_moea_executions(self.PARAM_FILES_FEW[index], self.PROBLEMS[index], self.TOTAL_N_FEW[index], self.EMO_DIR_FULL)

                if(flag):
                    fitness = 0
                    print("\n\t\tError executing EMO_PROJECT")
                    print("\t\tFinal fitness: "+str(fitness)+"\n")
                    return fitness

                ################################################################################
                # Obtain hypervolume values of the resulting pareto fronts
                ################################################################################
                avg_hv = self.get_avg_hv(self.PROBLEMS[index], self.TXT_DIMENSIONS[index], self.TOTAL_N_FEW[index], self.REFERENCE_POINTS[index], self.EMO_DIR_FULL)
                adjusted_fitness = (avg_hv-0.5*self.AVG_HV[index]) / self.AVG_HV[index]

                print("\t("+self.PROBLEMS[index]+"_"+self.TXT_DIMENSIONS[index]+")\tavg_HV:\t"+str(avg_hv)+"\tadj_fitness: \t"+str(adjusted_fitness))
                #fitness += avg_hv / self.AVG_HV[index]
                fitness += adjusted_fitness

            print("\tFinal fitness:\t"+str(fitness)+"\n")

        # Check if a better fitness has been found (comparing against the whole run) and perform more function evaluations
        fitness_ext = 0
        if(fitness > self.MAX_FITNESS_FEW):
            print("\tCurrent fitness ("+str(fitness)+") is greater than current max ("+str(self.MAX_FITNESS_FEW)+")")
            print("\tPerforming more evaluations...")

            for index in range(len(self.PROBLEMS)):

                ################################################################################
                # Execute modified EMO_PROJECT in every problem,
                # with a large amount of executions and more function evaluations
                ################################################################################
                flag = self.perform_moea_executions(self.PARAM_FILES_EXT[index], self.PROBLEMS[index], self.TOTAL_N_EXT[index], self.EMO_DIR_FULL)

                if(flag):
                    fitness = 0
                    print("\n\t\tError executing EMO_PROJECT")
                    print("\t\tFinal fitness: "+str(fitness)+"\n")
                    return fitness

                ################################################################################
                # Obtain hypervolume values of the resulting pareto fronts
                ################################################################################
                avg_hv = self.get_avg_hv(self.PROBLEMS[index], self.TXT_DIMENSIONS[index], self.TOTAL_N_EXT[index], self.REFERENCE_POINTS[index], self.EMO_DIR_FULL)
                adjusted_fitness = (avg_hv-0.5*self.AVG_HV[index]) / self.AVG_HV[index]

                print("\t("+self.PROBLEMS[index]+"_"+self.TXT_DIMENSIONS[index]+")\tavg_HV:\t"+str(avg_hv)+"\tadj_fitness: \t"+str(adjusted_fitness))
                #fitness_ext += avg_hv / self.AVG_HV[index]
                fitness_ext += adjusted_fitness

                self.AVG_HVS[index] = avg_hv

            print("\tFinal extended fitness:\t"+str(fitness_ext)+"\n")

        if(fitness_ext > self.MAX_FITNESS_EXT):
            self.MAX_FITNESS_FEW = fitness
            self.MAX_FITNESS_EXT = fitness_ext
            print("\tUpdating MAX_FITNESS_EXT: "+str(self.MAX_FITNESS_EXT)+" MAX_FITNESS_FEW: "+str(self.MAX_FITNESS_FEW)+" \n")
            fitness = fitness_ext

            # Register best value found
            with open(self.FINAL_LOG_FILENAME, "a") as log_file:
                aux_str = ""
                for index in range(len(self.PROBLEMS)):
                    aux_str += str(self.AVG_HVS[index])+"\t"
                log_file.write(ind.phenotype+"\t"+str(fitness_ext)+"\t"
                    +aux_str+"\n")

        return fitness
