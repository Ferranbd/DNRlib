import numpy as np
import networkx as nx
import random
import sys
import os
import logging
import warnings
import pandas as pd
import math

import GC_utils 

import GridCalEngine.api as gce  # For interfacing with the GridCal API
from GridCalEngine.IO.file_handler import FileOpen, FileSave

warnings.filterwarnings('ignore')  # Ignore warnings during execution

class BPSO:
    """Binary Particle Swarm Optimization (BPSO) class for optimizing network configurations.

    Attributes:
        net: An object representing the network to be optimized.
        NumCandidates: Number of candidate solutions to generate.
        NumOffLines: Number of lines that can be taken out of service.
        order: Optional parameter for ordering.
        flagUsed: Flag indicating whether to use certain properties.
        shuffled: Flag indicating whether to shuffle candidates.
        verbose_logging: Logging level for debug messages.
    """

    def __init__(self, grid, NumCandidates=10, TieLines=None, fitness_ratio=1, loss_factor = 0.08, verbose_logging=logging.WARNING):
        """Initializes the BPSO with the given parameters.

        Args:
            net: Network object for optimization.
            NumCandidates: Number of candidate solutions.
            order: Order for configurations, if any.
            flagUsed: Whether to use specific properties.
            shuffled: Whether to shuffle candidates.
            verbose_logging: Logging level for debug messages.
        """
        # PSO parameters
        self.grid = grid
        self.NumCandidates : int = NumCandidates
        self.TieLines : list = TieLines
        self.gbest : list = []
        self.fgbest : int = 0
        self.NumPF=0
        self.totalLoad = abs(sum([self.grid.loads[i].P for i in range(len(self.grid.loads))]))
        self.loss_factor = loss_factor
        self.fitness_ratio = fitness_ratio
#        self.NumOffLines : int = self.net.NumLines - self.net.NumBuses + 1  # Calculate out-of-service lines
        GC_utils.NetworkReconfiguration(self.grid, all=True, selected_configuration=None, value_all=True, value_configuration=False)
        self.NetworkLoops = GC_utils.SearchLoopsLines(self.grid, flagUsed=True, shuffle=False)      #following the paper, it should be shuffle=True, but not works properly
        self.NumOffLines: int = len(self.NetworkLoops)
        self.ConvergenceList : list = []
        # Set logging level for the BPSO class
        self.logger = logging.getLogger('bpso.py')  # Create a logger object for debug messages
        self.logger.setLevel(verbose_logging)  # Set the logging level for debug messages
        # configurate the algorithm
        self.config()
        self.logger.info(f"candidates: {self.NetworkLoops}")
        
    def __FitnessCalculation(self, config=None):
            self.NumPF+=1
    
            if config:
                GC_utils.NetworkReconfiguration(self.grid, all=True, value_all=True, selected_configuration= config, value_configuration=False)
            fitness = GC_utils.GC_FitnessCalculation(self.grid, self.totalLoad, loss_factor=self.loss_factor, split_factor=self.fitness_ratio)
            return fitness    
            
    def __powerflow(self, config=None):
            self.NumPF+=1
            resultPF,old_losses = GC_utils.GC_PowerFlow(self.grid, config=config)
            return resultPF,old_losses

    def RandomCandidateGenerator(self, n=1):
        #the Jakus paper proposed to create the candidates by generating MST with random weights
        RandomCandidates = []
        idx =0 
        while idx < n:
            candidate = [random.choice(sublist) for sublist in self.NetworkLoops] 
            #res, loss = GC_utils.GC_PowerFlow(self.grid, config=candidate)
            #res, loss = self.__powerflow(config= candidate)
            loss = self.__FitnessCalculation(config = candidate)
            radiality,_,_ = GC_utils.CheckRadialConnectedNetwork(self.grid)
            logging.getLogger('bpso.py').debug(f"Candidate generator (random){idx}/{n} : {loss,radiality} -- {candidate}")
            if radiality and loss < np.inf:
                RandomCandidates.append(candidate)
                logging.getLogger('bpso.py').debug(f"Candidate generator:{loss},{radiality} -- {candidate}")

                idx+=1

        return RandomCandidates

    def DefineActionSpace(self):
        """Defines the action space for the optimization based on network configurations."""
        # Generate the list of loops
        self.ActionSpaceSize = self.NumOffLines  # Set action space size based on out-of-service lines

    def config(self):
        """Configures the BPSO instance by defining the action space and initializing parameters."""
        self.DefineActionSpace()  # Define action space based on the current network configuration
        # Create initial list of candidates and log their fitness
        self.fpbest = np.full(self.NumCandidates, np.inf)  # Initialize personal best fitness values
        self.fitness = np.full(self.NumCandidates, np.inf)  # Initialize fitness values for all candidates
        self.pbest = []
    
        self.CreateCandidateList()
        for idx, candidate in enumerate(self.Candidates):
            self.logger.info(f"candidate[{idx}]:{candidate} : {self.fitness[idx]}")

        idx_best = np.argmin(self.fitness)
        self.gbest = self.Candidates[idx_best].copy()  # Initialize global best candidate as the first candidate
        self.fgbest = self.fpbest[idx_best]  # Initialize global best fitness as the one for the first candidate

        self.r1 = np.random.rand(self.NumCandidates, self.ActionSpaceSize)  # Random values for local search
        self.r2 = np.random.rand(self.NumCandidates, self.ActionSpaceSize)  # Random values for global search

    def getConvergeceList(self):
        return self.ConvergenceList


    def SortCandidates(self, Candidates, Loops):
        logging.getLogger('bpso.py').debug(f"SortCandidates")
        #logging.getLogger('bpso.py').debug(f"Loops to match {Loops}")
        SortedCandidates = []
        for candidate in Candidates:
            #print("######## Candidate ", candidate)
            sorted_candidate = []
            persistence = False
            occurrences = [[1 if num in sublist else 0 for sublist in Loops] for num in candidate]
            occurrences_array = np.array(occurrences)
            #logging.getLogger('bpso.py').debug(f"Candidate {candidate} occurences_array:{occurrences_array}")
            while not persistence:
                persistence = True
                column_summary = occurrences_array.sum(axis=0)
                #logging.getLogger('bpso.py').debug(f"occurences_array:{occurrences_array}")
                #logging.getLogger('bpso.py').debug(f"column_summary {column_summary}")
                for idx, element in enumerate(column_summary):
                    if element == 1:
                        persistence = False
                        idx_candidate = np.where(occurrences_array[:, idx] == 1)[0]
                        sorted_candidate.append(candidate[idx_candidate[0]])
                        occurrences_array[:, idx] = 0
                        occurrences_array[idx_candidate[0], :] = 0
                row_summary = occurrences_array.sum(axis=1)
                #logging.getLogger('bpso.py').debug(f"occurences_array:{occurrences_array}")
                #logging.getLogger('bpso.py').debug(f"sorted candidate: {sorted_candidate}")
                #logging.getLogger('bpso.py').debug(f"row_summary {row_summary}")
                for idx, element in enumerate(row_summary):
                    if element == 1:
                        persistence = False
                        idx_candidate = np.where(occurrences_array[idx, :] == 1)[0][0]
                        #logging.getLogger('bpso.py').debug(f"Candidate:{int(candidate[idx])} row:{idx} - {occurrences_array[idx, :]} found:{idx_candidate}")
                        logging.getLogger('bpso.py').debug(f"idx_candidate:{idx_candidate}")
                        sorted_candidate[idx_candidate[0]] = candidate[idx]
                        occurrences_array[idx, :] = 0
                        occurrences_array[:, idx_candidate] = 0
                #print("sorted candidate:", sorted_candidate)
                if persistence:
                    #print("persistence")
                    for idx, element in enumerate(column_summary):
                        if element:
                            persistence = False
                            idx_candidate = np.where(occurrences_array[:, idx] == 1)[0]
                            #print(f"Candidate:{int(candidate[idx_candidate[0]])} column:{idx} - {occurrences_array[:, idx]} found:{idx_candidate}")
                            sorted_candidate[idx] = candidate[idx_candidate[0]]
                            occurrences_array[:, idx] = 0
                            occurrences_array[idx_candidate[0], :] = 0
                            break
                #logging.getLogger('bpso.py').debug(f"FINAL sorted candidate: {sorted_candidate}")
            SortedCandidates.append(list(sorted_candidate))
        return SortedCandidates


    def CalculateFitness(self, configuration):
        """Calculates the fitness of a candidate configuration.

        Args:
            candidate: List of selected configurations to evaluate.

        Returns:
            A tuple containing the fitness value and a boolean indicating if the configuration is valid.
        """
        GC_utils.NetworkReconfiguration(self.grid, all=True, value_all=True,
                                        selected_configuration=list(configuration),
                                        value_configuration=False)
        radialconnected, _, _ = GC_utils.CheckRadialConnectedNetwork(self.grid)

        if radialconnected:
            #_,losses = GC_utils.GC_PowerFlow(self.grid)  # Compute power flow for the network
            #_, losses = self.__powerflow()
            losses = self.__FitnessCalculation()
            fitness_value = losses
            logging.getLogger('network.py').debug(
                    f"The configuration is complaint with the radiality constraint and has a fitness={fitness_value}")
            return (fitness_value, True)
        logging.getLogger('network.py').debug(f"The configuration does not comply with the radiality constraint")
        return (np.inf, False)  # Return infinite fitness for invalid configurations


    def UpdateCandidateFitness(self, candidate_idx):
        """Updates the fitness of a specific candidate and updates personal and global bests if needed.

        Args:
            candidate_idx: Index of the candidate to update.

        Returns:
            A boolean indicating if the candidate configuration is valid.
        """
        fitness_Value, result = self.CalculateFitness(self.Candidates[candidate_idx])
        self.fitness[candidate_idx] = fitness_Value  # Update the fitness of the candidate

        if result:
            # Updating personal best (pbest)
            if self.fitness[candidate_idx] < self.fpbest[candidate_idx]:
                self.pbest[candidate_idx] = self.Candidates[candidate_idx].copy()
                self.fpbest[candidate_idx] = self.fitness[candidate_idx].copy()

            # Updating global best (gbest)
            if self.fitness[candidate_idx] < self.fgbest:
                self.gbest = self.Candidates[candidate_idx].copy()
                self.fgbest = self.fitness[candidate_idx].copy()
                logging.getLogger('bpso.py').debug(f"\t\t new gbest= {self.gbest} with fitness={self.fgbest}")
        else:
            self.logger.debug(
                f"       bad valid candidate  {candidate_idx} -- {self.Candidates[candidate_idx]} -- {self.fitness[candidate_idx]}  ")
        return result

    def CreateCandidateList(self):
        """Creates a list of candidates and evaluates their fitness."""
        self.logger.debug("create candidates")

        self.Candidates = np.zeros((self.NumCandidates, self.ActionSpaceSize))  # Initialize candidate configurations
        if len(self.TieLines) > self.NumCandidates:
            self.TieLines = self.TieLines[:self.NumCandidates]

        # Add, as first candidates, the ones passed as "TieLines" argument to the algorithm
        # Generate candidates until the specified number is reached
        self.Candidates : list = self.RandomCandidateGenerator(n=self.NumCandidates-len(self.TieLines))
        self.Candidates.extend(self.TieLines)

        self.logger.debug(f"TieLines: {self.TieLines}")
        self.logger.debug(f"Candidates: {self.Candidates}")

        #convert and adapt the Candidates list
        #self.Candidates = np.array(self.Candidates)  # Convert list to numpy array for efficient processing
        #self.Candidates = [[int(element) for element in sublist] for sublist in self.Candidates]
        self.Candidates = self.SortCandidates(self.Candidates, self.NetworkLoops)
        self.pbest = self.Candidates.copy()     #initialize the particular best to the initial value of the candidates

        for idx_candidates, candidate in enumerate(self.Candidates):
            self.UpdateCandidateFitness(idx_candidates)
        self.fpbest = self.fitness.copy()

        #for candidate in self.Candidates:
        #    self.logger.debug(candidate)

    def CalculateNewPosition(self, candidate_idx, c1, c2, w):
        """Calculates the new position of a candidate based on its velocity and best known positions.

        Args:
            candidate_idx: Index of the candidate whose position is being updated.
            c1: Cognitive coefficient for personal best influence.
            c2: Social coefficient for global best influence.
            w: Inertia weight for velocity.

        Returns:
            A tuple of the new candidate position and the updated velocity.
        """
        newValue = self.Candidates[candidate_idx].copy()  # Copy the current candidate position
        oldValue = self.Candidates[candidate_idx].copy()  # Copy for reference during updates

        for j in range(self.ActionSpaceSize):
            old_value = oldValue[j]
            gbest_idx = self.NetworkLoops[j].index(self.gbest[j])  # Index of the global best in the action space
            #self.logger.debug(f"\taa {candidate_idx},{j}")
            #self.logger.debug(f"\tbb {self.pbest[candidate_idx][j]}")
            #self.logger.debug(f"\tcc {self.NetworkLoops[j]}")
            pbest_idx = self.NetworkLoops[j].index(self.pbest[candidate_idx][j])  # Index of the personal best
            cand_idx = self.NetworkLoops[j].index(oldValue[j])  # Index of the current candidate

            # Calculate distances for velocity updates
            global_distance = gbest_idx - cand_idx
            local_distance = pbest_idx - cand_idx

            # Calculate corrections based on distances and randomness
            global_correction = c2 * self.r2[candidate_idx][j] * global_distance
            local_correction = c1 * self.r1[candidate_idx][j] * local_distance

            prev_velocity = self.velocity[candidate_idx][j]  # Store the previous velocity
            self.velocity[candidate_idx][j] = w * self.velocity[candidate_idx][j] + local_correction + global_correction

            # Ensure velocity is within specified limits
            self.velocity[candidate_idx][j] = min(max(self.velocity[candidate_idx][j], self.vmin), self.vmax)

            # Randomly alter velocity if it hasn't changed
            if abs(self.velocity[candidate_idx][j]) == abs(prev_velocity):
                self.velocity[candidate_idx][j] = random.random() * self.velocity[candidate_idx][j]

            # Updating particles' coordinate
            scope = len(np.nonzero(self.NetworkLoops[j])[0]) - 1  # Scope of available lines on the loop

            # Sigmoid function to map the velocity to a new position
            sigmoid = 1 / (1 + math.exp(
                -self.beta * self.velocity[candidate_idx][j])) - 0.5  # Sigmoid function centered at 0
            newPosition = cand_idx + int(np.floor(scope * sigmoid))  # Calculate new position based on sigmoid

            newPosition = min(max(newPosition, 0), scope - 1)  # Clamp new position within bounds
            newValue[j] = self.NetworkLoops[j][newPosition]  # Update new value for the candidate
            # Log details of the update process
            self.logger.debug(
                f"\tCandidate[{candidate_idx},{j}] - loc.dist.//global.dist ({self.pbest[candidate_idx][j]}/{pbest_idx}-{oldValue[j]}/{cand_idx})={local_distance} // {self.gbest[j]}/{gbest_idx}-{oldValue[j]}/{cand_idx})={global_distance}  \t\t w//local correction//global correction//wnew//scope//sigmoid//newPosition//old//new value (  {w:.3f} {local_correction:.4f}  {global_correction:.4f} {self.velocity[candidate_idx][j]:.4f}  {scope}  {sigmoid}  {newPosition} {old_value}  {newValue[j]}")

        return newValue, self.velocity[candidate_idx]  # Return the new position and updated velocity

    def UpdateCandidateNewPosition(self, candidate_idx, w, c1, c2):
        """Updates the position of a candidate based on the BPSO algorithm.

        Args:
            candidate_idx: Index of the candidate to update.
            w: Inertia weight for velocity calculation.
            c1: Cognitive coefficient for personal best influence.
            c2: Social coefficient for global best influence.
        """
        self.logger.debug(
            f"\tCandidate[{candidate_idx}]: {self.Candidates[candidate_idx]}  Pbest {self.pbest[candidate_idx]}   Gbest {self.gbest}  -self.velocity={self.velocity[candidate_idx]}")

        # Calculate new position and velocity for the candidate
        old_candidate = self.Candidates[candidate_idx].copy()
        self.Candidates[candidate_idx], self.velocity[candidate_idx] = self.CalculateNewPosition(candidate_idx, c1, c2, w)
        # Update fitness for the new candidate position
        if self.UpdateCandidateFitness(candidate_idx):
            self.logger.debug(f"\t\tnew Candidate[{candidate_idx}]: {self.Candidates[candidate_idx]}")
        else :
            self.logger.debug(f"\t\tnew Candidate[{candidate_idx}]: {self.Candidates[candidate_idx]} is not valid, so we dismiss the changes which creates the invalid solution")
            self.Candidates[candidate_idx] = old_candidate
            self.UpdateCandidateFitness(candidate_idx)  #update gbest and fgbest, pbest, fpbest
            self.logger.debug(f"\t\t\t keep the previuos Candidate[{candidate_idx}]: {self.Candidates[candidate_idx]} ")

    # @jit(target_backend='cuda')
    def Solve(self, wmax=0.9, wmin=0.5, beta=1, maxiter=30, vmax=4, c1=2, c2=2):
        """Solves the optimization problem using the BPSO algorithm.

        Args:
            wmax: Maximum inertia weight.
            wmin: Minimum inertia weight.
            beta: Parameter for the sigmoid function.
            maxiter: Maximum number of iterations for the optimization.
            vmax: Maximum velocity for the particles.

        Returns:
            The best global configuration found during optimization.
        """
        # Main loops
        self.logger.info("Start solving by BPSO")
        self.wmax = wmax
        self.wmin = wmin
        self.beta = beta
        self.vmax = vmax
        self.vmin = -vmax  # Set minimum velocity to the negative of max velocity
        self.maxiter = maxiter
        self.c1 = c1
        self.c2 = c2
        self.velocity = np.random.uniform(self.vmin, self.vmax,
                                          size=(self.NumCandidates, self.NumOffLines))  # Initialize velocities

        self.ConvergenceList = []
        iter = 0

        # Iterate for a maximum number of iterations
        while iter < self.maxiter:
            iter += 1
            self.logger.debug(
                f"Iter {iter} \n \t Candidates {self.Candidates} \n\t iter{iter}  gbest[{self.gbest}]{self.fgbest}   ")

            w = self.wmax - (self.wmax - self.wmin) * iter / self.maxiter  # Calculate current inertia weight
#            self.c1 = 2  # Cognitive coefficient (could be randomized for variability)
#            self.c2 = 2  # Social coefficient (could be randomized for variability)

            for loop in self.NetworkLoops:
                            self.logger.debug(f"loop : {loop}")

            # Update the position for each candidate
            for candidate_idx in range(self.NumCandidates):
                self.UpdateCandidateNewPosition(candidate_idx, w, self.c1, self.c2)

            self.ConvergenceList.append((self.gbest, self.fgbest))

            self.logger.debug(
                f"    iter{iter}  gbest[{self.gbest}]{self.fgbest}    candidate[0]:{self.fitness[1]}-{self.Candidates[1]}")

        return self.gbest  # Return the best global candidate found



if __name__ == '__main__':
    print('Algorithms to find the optimal distribution network configuration')

    logging.basicConfig(
        level=logging.ERROR,  # Set the log level to DEBUG
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Set the log format
        datefmt='%Y-%m-%d %H:%M:%S'  # Set the date format
    )

    case=3
    # Create a network object
    if case==1:
        gridGC = FileOpen("D:\\15_Thesis-code\\DistributionNetwork_libraries\\NetworkExamples\\gridcal\\case33.gridcal").open()
        TieLinesName=['line 32','line 33','line 34','line 35','line 36']
    if case==2:
        gridGC = FileOpen("D:\\15_Thesis-code\\DistributionNetwork_libraries\\NetworkExamples\\gridcal\\case69.gridcal").open()
        TieLinesName = ['line 57','line 10','line 69','line 14','line 19']
    if case==3:
        gridGC = FileOpen("D:\\15_Thesis-code\\DistributionNetwork_libraries\\NetworkExamples\\gridcal\\case118.m").open()
        TieLinesName = ['75_77_1', '69_75_1', '77_80_1', '80_97_1', '94_95_1', '92_94_1', '105_106_1', '100_103_1', '100_104_1', '103_110_1', '103_104_1', '92_102_1', '80_99_1', '80_98_1', '92_93_1', '89_90_1', '85_88_1', '82_83_1', '83_84_1', '68_81_1', '62_66_1', '60_61_1', '64_65_1', '59_60_1', '63_64_1', '54_55_1', '55_56_1', '54_56_1', '49_51_1', '51_52_1', '49_50_1', '49_54_1', '49_69_1', '45_46_1', '46_47_1', '34_37_1', '37_39_1', '40_42_1', '40_41_1', '15_19_1', '15_17_1', '27_32_1', '23_25_1', '17_31_1', '17_113_1', '27_28_1', '23_24_1', '24_72_1', '19_20_1', '4_5_1', '8_30_1', '3_5_1', '12_14_1', '12_16_1', '5_6_1', '1_2_1', '17_18_1', '34_36_1', '47_69_1', '77_78_1', '70_74_1', '69_70_1']
    if case==4:
        import pandapower as pp
        import simbench as sb
        import GC_PandaPowerImporter
        sb_code1 = "1-HVMV-urban-2.203-0-no_sw"
        gridPP = sb.get_simbench_net(sb_code1)
        gridPP.switch.drop([232,234,236,238,240, 242,244,246], inplace=True)
        gridPP.trafo.drop([1,3,4], inplace=True)
        gridPP.line.drop(set([123,226,139,140,151,161,166,170,173,178,180,186,187,188,208,223,225,123,226,227,232,228,229,230,231,227,232,233]), inplace=True)
        gridPP.ext_grid.at[0,'name']="grid_ext"
        gridPP.line['in_service'] = True
        pp.runpp(gridPP)
        gridGC = GC_PandaPowerImporter.PP2GC(gridPP)
        TieLinesName=['1_2_1', '1_24_1', '1_36_1', '1_47_1', '51_52_1', '1_60_1', '1_74_1', '1_85_1', '117_181_1', '171_117_1', '117_125_1', '127_164_1', '121_188_1', '146_147_1', '171_181_1', '116_196_1', '116_154_1']

    TieLinesID=GC_utils.GC_Line_Name2idtag_array(gridGC, TieLinesName)

    _, loss = GC_utils.GC_PowerFlow(gridGC, config=TieLinesID)
    radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
    print(f"Original network:{loss}, radiality:{radiality} is {GC_utils.GC_Line_idtag2name_array(gridGC,TieLinesID)}" )

    for hyper in [1]:
        # Create an MSTgreedy object
        khalil = BPSO(gridGC, TieLines=[], NumCandidates=16, fitness_ratio=1, loss_factor=0.08, verbose_logging=logging.ERROR)

        # Solve the Minimum Spanning Tree problem
        disabled_lines = khalil.Solve(wmax=0.9, wmin=0.5, beta=1, maxiter=20, vmax=4, c1=2, c2=2)

        # Print the list of disabled line indices
        _,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        print(f"{hyper} : The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{khalil.NumPF} ") #is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )