import numpy as np
import networkx as nx
import random
import sys
import os
import logging
import warnings
import pandas as pd
import GC_utils

import GridCalEngine.api as gce  # For interfacing with the GridCal API
from GridCalEngine.IO.file_handler import FileOpen, FileSave

# Add the directory containing utility_script.py to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

import GC_utils


warnings.filterwarnings('ignore')  # Ignore warnings during execution

## based on baran 1989
class Baran1989:
    """
    Class for solving the Baran 1989 test feeder network problem.

    Attributes:
        net (DistributionNetwork): The distribution network object.
        initial_config (list, optional): The initial tie line configuration. Defaults to None.
        logger (logging.Logger): Logger object for debug messages.
    """

    def __init__(self, grid=None, TieLines=None, fitness_ratio=1, loss_factor = 0.08, verbose_logging=logging.WARNING):
        """
        Initializes the Baran1989 class.

        Args:
            net (DistributionNetwork): The distribution network object.
            TieLines (list, optional): The initial tie line configuration. Defaults to None.
            verbose_logging (int, optional): The logging level for debug messages. Defaults to logging.WARNING.
        """
        self.grid = grid  # Store the distribution network object
        self.logger = logging.getLogger('baran1989.py')  # Create a logger object for debug messages
        self.logger.setLevel(verbose_logging)  # Set the logging level for debug messages
        self.initial_config = TieLines  # Store the initial tie line configuration
        self.NumPF=0
        self.totalLoad = 1 #abs(sum([self.grid.loads[i].P for i in range(len(self.grid.loads))]))
        self.loss_factor = loss_factor
        self.fitness_ratio = fitness_ratio
        
    def __FitnessCalculation(self, config=None):
            self.NumPF+=1
            if config:
                GC_utils.NetworkReconfiguration(self.grid, all=True, value_all=True, selected_configuration= config, value_configuration=False)
            #fitness = GC_utils.GC_FitnessCalculation(self.grid, self.totalLoad, loss_factor=self.loss_factor, split_factor=self.fitness_ratio)
            fitness = GC_utils.GC_FitnessCalculation(self.grid, 1, loss_factor=1, split_factor=1)
            return fitness
            
    def Solve(self, order=None):
        """
        Solves the Baran 1989 problem by finding the best tie line configuration that minimizes power losses.

        This function iterates through all loops in the initial network and tries opening/closing tie lines
        within each loop. It then checks for radial connectivity and performs power flow analysis to determine
        the power losses. The configuration with the minimum losses is returned.

        Args:
            order (list, optional): An order for iterating through loops (not used in this implementation). Defaults to None.

        Returns:
            list: The best tie line configuration that minimizes power losses.
        """
        self.logger.info("Start solving Baran 1989")

        # Perform network reconfiguration (consider using a more descriptive method name)
        #for line in self.grid.lines:
        #    line.active = True

        # Find initial loops
        GC_utils.NetworkReconfiguration(self.grid, all=True, value_all=True)
        init_loops = GC_utils.SearchLoopsLines(self.grid,flagUsed=True)
        self.logger.debug(f"Initial loops =")
        for loop in init_loops:
            self.logger.debug(f"\t\t{loop}")

        #set initial configuration
        if self.initial_config is None:
            self.initial_config = [sublist[0] for sublist in init_loops if sublist]
        #GC_utils.NetworkReconfiguration(self.grid,all=True,value_all=True,selected_configuration=self.initial_config,value_configuration=False)
        #self.logger.debug(f"initial config: {self.initial_config}")
        #self.logger.debug(f"radiality: {GC_utils.CheckRadialConnectedNetwork(self.grid)}")
        #_, losses = GC_utils.GC_PowerFlow(grid=self.grid)
        #_, losses = self.__powerflow()
        #self.logger.debug(f"losses: {losses}")

        GC_utils.NetworkReconfiguration(self.grid, all=True, selected_configuration=self.initial_config, value_all=True, value_configuration=False)
        # Check for radial connectivity
        radial, _, _ = GC_utils.CheckRadialConnectedNetwork(self.grid)
        #_,init_losses = GC_utils.GC_PowerFlow(self.grid)
        #_, init_losses = self.__powerflow()
        bestTieLines = self.initial_config
        if radial:
            #bestlosses = init_losses
            best_fitness  = self.__FitnessCalculation()
        else:
            #bestlosses = np.inf
            best_fitness = np.inf

        #self.logger.debug(f"initial configuration {bestTieLines}, losses:{bestlosses},radiality:{radial} ")
        self.logger.debug(f"initial configuration {bestTieLines}, fitness:{best_fitness},radiality:{radial} ")

        for idx, loop in enumerate(init_loops):
            self.logger.debug(f"new loop {loop}")
            for tie in loop:
                self.logger.debug(f" new tie: '{tie}' while loop(idx)={idx}")
                orig_tie = bestTieLines[idx]
                newTieLines = [tie if x == orig_tie else x for x in bestTieLines]
                # Perform network reconfiguration with the new tie line configuration
                GC_utils.NetworkReconfiguration(self.grid, all=True, selected_configuration=newTieLines, value_all=True, value_configuration=False)
                # Check for radial connectivity
                radial, _, _ = GC_utils.CheckRadialConnectedNetwork(self.grid)
                # Log debug information
                self.logger.debug(f" temporary newTieLines {newTieLines} : radial:{radial}")
                # If the network is radial, perform power flow analysis and check for improved losses
                if radial:
                    #_,losses = GC_utils.GC_PowerFlow(self.grid)
                    #_, losses = self.__powerflow()
                    fitness  = self.__FitnessCalculation()

                    #self.logger.debug(f"  TieLine {newTieLines} : radial = {radial} losses={self.net.losses.real:.4f} outofservice={OutofServiceLines_end}")

                    # If the new losses are better than the best found so far, update the best tie lines and losses
                    #if (losses < bestlosses):
                    if (fitness < best_fitness):
                        bestTieLines = newTieLines
                    #    bestlosses = losses
                        best_fitness = fitness
                        #self.logger.debug(f" losses are smaller {bestlosses} than previous tieline, so new one is :{bestTieLines}")
                        self.logger.debug(f" fitness is smaller {best_fitness} than previous tieline, so new one is :{bestTieLines}")
                    else:
                        #self.logger.debug(f" losses are bigger {losses}  than the previous besttielines {bestlosses}")
                        self.logger.debug(f" fitness is bigger {fitness}  than the previous besttielines {best_fitness}")

        # Log the best tie lines and losses
        #self.logger.debug(f"Best TieLines ={bestTieLines} with losses={bestlosses}")
        self.logger.debug(f"Best TieLines ={bestTieLines} with fitness={best_fitness}")

        # Return the best tie lines
        return list(bestTieLines)
    
if __name__ == '__main__':
    print('Algorithms to find the optimal distribution network configuration')

    logging.basicConfig(
        level=logging.ERROR,  # Set the log level to DEBUG
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Set the log format
        datefmt='%Y-%m-%d %H:%M:%S'  # Set the date format
    )

    case=1
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
    print("Original network: ",GC_utils.GC_Line_idtag2name_array(gridGC,TieLinesID), loss, radiality )
   
    fitness_ratios = [1]

    for fitness_ratio in fitness_ratios:
        # Create a Baran object
        baran = Baran1989(gridGC, TieLines=TieLinesID, fitness_ratio=fitness_ratio, loss_factor = 0.08, verbose_logging=logging.ERROR)
        disabled_lines = baran.Solve()

        # Print the list of disabled line indices
        res,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        Vmin = res.results.get_bus_df().Vm.min()
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        print(f"The new optimal configuration losses:{loss}, Vmin:{Vmin}, radiality:{radiality}, numPF:{baran.NumPF} ")#is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )






