from fitness.base_ff_classes.base_ff import base_ff
import os
from os import path
import numpy as np
import scipy.stats as ss
from algorithm.parameters import params

class HVC_8D(base_ff):
    """
    Basic fitness function first_linelate for writing new fitness functions. This
    basic first_linelate inherits from the base fitness function class, which
    contains various checks and balances.

    Note that all fitness functions must be implemented as a class.

    Note that the class name must be the same as the file name.

    Important points to note about base fitness function class from which
    this first_linelate inherits:

      - Default Fitness values (can be referenced as "self.default_fitness")
        are set to NaN in the base class. While this can be over-written,
        PonyGE2 works best when it can filter solutions by NaN values.

      - The standard fitness objective of the base fitness function class is
        to minimise fitness. If the objective is to maximise fitness,
        this can be over-written by setting the flag "maximise = True".

    """

    # The base fitness function class is set up to minimise fitness.
    # However, if you wish to maximise fitness values, you only need to
    # change the "maximise" attribute here to True rather than False.
    # Note that if fitness is being minimised, it is not necessary to
    # re-define/overwrite the maximise attribute here, as it already exists
    # in the base fitness function class.
    maximise = False

    def __init__(self):

        # Initialise base fitness function class.
        super().__init__()
        self.population = None
        self.hvr = None
        self.worst_inds = None
        self.sums = None
        self.n_files = None
        self.file_sizes = None

        # User defined parameters
        self.max_file_size = 120
        self.dim = 8

        self.get_data()

    def evaluate(self, ind, **kwargs):
        """
        Default fitness execution call for all fitness functions. When
        implementing a new fitness function, this is where code should be added
        to evaluate target phenotypes.

        There is no need to implement a __call__() method for new fitness
        functions which inherit from the base class; the "evaluate()" function
        provided here allows for this. Implementing a __call__() method for new
        fitness functions will over-write the __call__() method in the base
        class, removing much of the functionality and use of the base class.

        :param ind: An individual to be evaluated.
        :param kwargs: Optional extra arguments.
        :return: The fitness of the evaluated individual.
        """

        # Evaluate the fitness of the phenotype
        fitness_sum = 0
        evaluation_values = np.zeros(self.max_file_size)
        for file_index in range(self.n_files):
            pop_size = self.file_sizes[file_index]
            sums = self.sums[file_index]        # should be defined in the grammar as sums[0] and sums[1]

            # Evaluate and store each individual value
            for ind_index in range(pop_size):
                x = self.population[file_index,ind_index] # in the grammar it should be refered as x[0], x[1]...
                x_min = min(x)
                x_max = max(x)
                evaluation_values[ind_index] = eval(ind.phenotype)
            
            # Rank the obtained values
            ranked_values = ss.rankdata(evaluation_values[0:pop_size])
            
            #Determine the fitness based on the difference in the worst individual
            goal_hvr = self.hvr[file_index][0:pop_size]
            worst_index = np.argmax(ranked_values)
            actual_worst_index = self.worst_inds[file_index]
            wi_original_rank = goal_hvr[worst_index]            
            awi_original_rank = goal_hvr[actual_worst_index]
            ind_fitness = abs(wi_original_rank-awi_original_rank)*100 / pop_size
            if ind_fitness<25:
                fitness_sum += ind_fitness
            elif ind_fitness<50:
                fitness_sum += ind_fitness*1.25
            else:
                fitness_sum += ind_fitness*1.5
            
        fitness = fitness_sum/self.n_files

        return fitness


    def get_data(self):
        """ Return the training data for the current experiment.
        """
        datasetTrain = params['DATASET_TRAIN']
        trainRoot = path.join("..","datasets",datasetTrain)

        # List all files within the directory
        content = os.listdir(trainRoot)
        self.n_files = len(content)

        # Initialize size of storage variables
        self.file_sizes = np.zeros(self.n_files,dtype=int)
        self.population = np.zeros((self.n_files, self.max_file_size, self.dim))
        self.sums = np.zeros((self.n_files,self.dim))
        self.hvr = np.zeros((self.n_files, self.max_file_size))
        self.worst_inds = np.ones(self.n_files,dtype=int)

        # Process each file
        for file_index in range(self.n_files): # Recorre los archivos

            filename = path.join(trainRoot,content[file_index])
            file = open(filename)
            lines = file.readlines()
            file.close()

            first_line = lines[0].split()
            n = int(first_line[1]) # Number of points

            self.file_sizes[file_index] = n

            # Store points and other information
            worst_index = -1
            worst_hvr = -1
            for ind_index in range(n):
                data_row = np.array(lines[ind_index+1].split())
                for d in range(self.dim):
                    data = float(data_row[d])
                    self.population[file_index, ind_index, d] = data
                    self.sums[file_index,d] += data
                hvr = float(data_row[d+2])   # this column should contain the hypervolume contribution rank
                if hvr > worst_hvr:
                    worst_index = ind_index
                    worst_hvr = hvr
                self.hvr[file_index, ind_index] = hvr
            self.worst_inds[file_index] = worst_index