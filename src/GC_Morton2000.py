import numpy as np
import networkx as nx
import random
import sys
import os
import logging
import warnings
import pandas as pd

import GridCalEngine.api as gce  # For interfacing with the GridCal API
from GridCalEngine.IO.file_handler import FileOpen, FileSave
import GC_utils 

warnings.filterwarnings('ignore')  # Ignore warnings during execution


class Morton2000:
    """Implementation of the Morton 2000 algorithm for optimizing network configurations.

    Attributes:
        net: The network object containing distribution network data.
        candidates: A list to store candidate network configurations.
        initial_loops: Initial loop configurations for the network.
        init_config: Initial configuration of the network.
        best_power_loss: Variable to track the best (minimum) power loss during optimization.
    """

    def __init__(self, grid, init_config, fitness_ratio=1, loss_factor = 0.08, verbose_logging=logging.WARNING) -> None:
        """Initializes the Morton2000 class with required parameters.

        Args:
            net: The network object for optimization.
            init_config: The initial configuration of the network.
            verbose_logging: Logging level for debug messages.
        """
        self.grid = grid  # Set the network object
        self.candidates = []  # Initialize an empty list for candidates
        self.init_config = init_config  # Store the initial configuration
        self.best_power_loss = np.inf  # Initialize best power loss to infinity
        logging.getLogger('morton2000.py').setLevel(verbose_logging)  # Set logging level
        self.NumPF=0
        self.totalLoad = abs(sum([self.grid.loads[i].P for i in range(len(self.grid.loads))]))
        self.loss_factor = loss_factor
        self.fitness_ratio = fitness_ratio
        GC_utils.NetworkReconfiguration(self.grid, all=True, value_all=True, selected_configuration=[], value_configuration=False)
        self.initial_loops = GC_utils.SearchLoopsLines(self.grid, flagUsed=False)

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
    
    def _search_trees(self, present_candidate, loops, rec_level, total_candidates):
        """Recursively searches for trees in the network based on given loops.

        Args:
            present_candidate: The current candidate configuration being explored.
            loops: The available loops to explore.
            rec_level: The current recursion level.
            total_candidates: A set of all candidate configurations found so far.

        Returns:
            A set of all candidate configurations found during the search.
        """
        present_candidate = present_candidate.copy()  # Create a copy of the present candidate
        available_loops = loops.copy()  # Create a copy of available loops
        total_candidates = total_candidates.copy()  # Create a copy of total candidates

        # If the number of total candidates exceeds the maximum, return the current candidates
        if len(total_candidates) > self.max_candidates:
            return total_candidates

        any_loop = available_loops.pop(0)  # Remove the first loop to explore
        for idx,any_line in enumerate(any_loop):
            if len(total_candidates) > self.max_candidates:
                return total_candidates
            present_candidate[rec_level] = any_line  # Set the current line in the candidate configuration
            added_candidated = tuple(sorted(tuple(present_candidate.copy())))  # Create a sorted tuple for uniqueness
            total_candidates.add(added_candidated)  # Add the candidate to the total candidates set
            logging.getLogger('morton2000.py').debug(f"[{rec_level}] Total : {len(total_candidates)} - New candidate {added_candidated}")
            if len(available_loops) > 0:  # If there are more loops to explore
                total_candidates = self._search_trees(present_candidate, available_loops, rec_level + 1, total_candidates)
            #print(f"aqui:({idx}): {any_line}")
            self.grid.lines[idx].active = True  # Enable the current line in the network
        return total_candidates  # Return the total candidates found

    def SearchTrees(self):
        """Starts the search for trees in the network based on initial loops.

        This method reconfigures the network, initializes candidates, and starts the recursive search for trees.
        """
        logging.getLogger('morton2000.py').debug("Start SearchTrees")
        self.tot_candidates = set()  # Initialize total candidates as a set
        self.tot_candidates.add(tuple(self.init_config.copy()))  # Add initial configuration as a candidate
        GC_utils.NetworkReconfiguration(self.grid, all=True, value_all=True, selected_configuration=self.init_config,value_configuration=False)
        logging.getLogger('morton2000.py').debug(f"loops: {self.initial_loops}")
        # Start the recursive search for trees
        self.tot_candidates = self._search_trees(self.init_config, self.initial_loops, 0, self.tot_candidates)

    def RemoveDuplicates(self, tot_candidates):
        """Removes duplicate candidates from the list of total candidates.

        This method populates the candidates list with unique configurations.
        """
        logging.getLogger('morton2000.py').debug("remove duplicated candidates")
        clean_candidates=[]
        len_tot = len(tot_candidates)
        for idx, item in enumerate(tot_candidates):
            if item not in clean_candidates:  # Check if the item is already in the unique candidates list
                if len(item) == len(set(item)): #if the network does not have any repeated element (if it has it, it means not radial)
                    clean_candidates.append(item)  # Add unique item to the candidates list
            if idx % 100 == 0:
                logging.getLogger('morton2000.py').debug(f"reviewed candidates:{idx}/{len_tot}")

        logging.getLogger('morton2000.py').debug(f"found {len(tot_candidates)} candidates")  # Print the total number of candidates found
        logging.getLogger('morton2000.py').debug(f"found {len(clean_candidates)} unique candidates")  # Print the number of unique candidates found
    
        return clean_candidates
    
    def ChooseBestTree(self):
        """Chooses the best tree configuration based on the lowest power loss.

        This method tests each candidate configuration and tracks the configuration with the lowest power loss.

        Returns:
            The best tree configuration and the associated power flow results.
        """
        logging.getLogger('morton2000.py').debug("choose best tree")
        best_tree = []  # Initialize the best tree configuration
        pf_result = []  # Initialize a list to store power flow results
        idx = 0  # Index to track progress
        #best_power_loss = np.inf  # Initialize best power loss to infinity
        best_fitness = np.inf  # Initialize best power loss to infinity

        len_candidates = len(self.candidates)
        for candidate_tree in self.candidates:  # Iterate through each candidate configuration
            idx += 1
            if idx % 100 == 0:
                logging.getLogger('morton2000.py').debug(f"tested candidates:{idx} / {len_candidates}")
            # Reconfigure the network with the current candidate configuration
            #_,losses = GC_utils.GC_PowerFlow(self.grid, config=candidate_tree)
            fitness  =self.__FitnessCalculation(config=candidate_tree)
            validTree,_,_ = GC_utils.CheckRadialConnectedNetwork(self.grid)
            logging.getLogger('morton2000.py').debug(f"tested candidate:{candidate_tree} : {fitness} / {validTree}")
            if validTree:  # If the configuration is valid
                fitness  =self.__FitnessCalculation()
                ##_, losses = GC_utils.GC_PowerFlow(self.grid)# Calculate power flow for the current configuration
                pf_result.append(fitness)  # Store the power flow result
                # Check if the current power loss is less than the best found
                if best_fitness > fitness:
                    best_fitness = fitness  # Update best power loss
                    best_tree = candidate_tree  # Update the best tree configuration
            else:
                pf_result.append(np.inf)  # Append infinity if the configuration is invalid
        logging.getLogger('morton2000.py').debug(f"new best: {best_tree} with {best_fitness}")  # Print the best configuration found

        return best_tree, pf_result  # Return the best tree and power flow results

    def Solve(self, max_candidates=np.inf):
        """Main method to solve the network optimization problem.

        Args:
            max_candidates: The maximum number of candidate configurations to consider.

        Returns:
            The best tree configuration and associated results.
        """
        logging.getLogger('morton2000.py').info("Start solving Morton 2000 (brute force)")
        self.max_candidates = max_candidates  # Set the maximum number of candidates
        # The following methods can be uncommented to run the full optimization process
        self.SearchTrees()  # Search for candidate trees
        self.candidates=self.RemoveDuplicates(self.tot_candidates)  # Remove duplicates from the candidates
        print(f"found {len(self.candidates)} possible networks")
        best, fitness = self.ChooseBestTree()
        return best  # Choose the best tree from the candidates



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

    # Create an Morton2000 object
    morton = Morton2000(gridGC, init_config=TieLinesID, fitness_ratio=1, loss_factor = 0.08, verbose_logging=logging.INFO)

    # Solve the Minimum Spanning Tree problem
    disabled_lines = morton.Solve() #max_candidates=100)
    
    # Print the list of disabled line indices
    print(disabled_lines)

    _,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
    radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
    print(f"The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{morton.NumPF}") # is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )
    #the solution for the case33 is ('line 13', 'line 31', 'line 36', 'line 6', 'line 8') 139551.3472210437