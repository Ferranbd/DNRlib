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

warnings.filterwarnings('ignore')  # Ignore warnings during execution


class Merlin1975:
    """Implementation of the Merlin 1975 algorithm for optimizing network configurations.

    Attributes:
        net: The network object that contains information about the distribution network.
    """

    def __init__(self, grid=None, verbose_logging=logging.WARNING):
        """Initializes the Merlin1975 class with a network object.

        Args:
            net: The network object for optimization.
            verbose_logging: Logging level for debug messages.
        """
        logging.getLogger('merlin1975.py').setLevel(verbose_logging)  # Set logging level
        self.grid = grid  # Store the network object
        # Set logging level for the Merlin1975 class   
        self.NumPF=0
        
    def __powerflow(self, config=None):
            self.NumPF+=1
            resultPF,old_losses = GC_utils.GC_PowerFlow(self.grid, config=config)
            return resultPF,old_losses
    
    def __FitnessCalculation(self, config=None):
            self.NumPF+=1

            self.loss_factor = 0.08
            self.fitness_ratio = 1
            self.totalLoad = abs(sum([self.grid.loads[i].P for i in range(len(self.grid.loads))]))
            
            if config:
                GC_utils.NetworkReconfiguration(self.grid, all=True, value_all=True, selected_configuration= config, value_configuration=False)
            fitness = GC_utils.GC_FitnessCalculation(self.grid, self.totalLoad, loss_factor=self.loss_factor, split_factor=self.fitness_ratio)
            return fitness
            
            
    def SearchMinFlow(self, banned):
        """Searches for the line with the minimum current flow that is not banned.

        Args:
            banned: List of indices of banned lines.

        Returns:
            The index of the line with minimum current flow if available; otherwise, -1.
        """

        # Iterate through the sorted lines to find a suitable line
        for idx, line in self.sortedlines.iterrows():
            #logging.getLogger('Merlin1975.py').debug(f"SearchMinFlow:, idx:{idx} line_current:{line.If}, line_Enabled:{line.Enabled}, line_controllable:{line.controllable}")
            if idx not in banned and (line['active']):
                return idx  # Return the index of the line if it is not banned and is controllable
        return -1  # Return -1 if no suitable line is found

    def Solve(self):
        """Solves the optimization problem using the Merlin 1975 algorithm.

        Returns:
            The list of lines that are out of service after optimization.
        """
        logging.getLogger('Merlin1975.py').info("Start solving Merlin 1975")
        logging.getLogger('Merlin1975.py').debug("\t\t Reconfiguration all True")
        banned = []  # Initialize a list for banned lines

        for line in self.grid.lines:
            line.active = True
        for trafo in self.grid.transformers2w:
            trafo.active = True

        graph = GC_utils.GC2Graph(self.grid)  # Initialize the Minimum Spanning Tree as the network graph      
        Cycles = GC_utils.EdgeCycles(graph)
        NumCycles = len(Cycles)  # Count the number of cycles


        logging.getLogger('merlin1975.py').debug(
            f"inici : numBus={len(self.grid.buses)}, cycles={NumCycles} ")
        logging.getLogger('merlin1975.py').debug(f"NumCycles {NumCycles} - Cycles:{Cycles}")

        mu = NumCycles  # Set mu to the number of cycles

        while mu:
            # Calculate power flow for the current network configuration
            active = list([line.active for idx, line in enumerate(self.grid.lines) ])
            name = list([line.idtag for idx, line in enumerate(self.grid.lines) ])
            #self.__FitnessCalculation()
            resultsPF, _ = self.__powerflow()
            current = resultsPF.results.If[:len(self.grid.lines)].real

            logging.getLogger('merlin1975.py').debug(f"currents({len(current)}):{current}")
            logging.getLogger('merlin1975.py').debug(f"name({len(name)}):{name}")

            lines = pd.DataFrame({"active":active, "name":name, "current":current})
            self.sortedlines = lines.sort_values(by="current", ascending=True)
            # Search for the line with the minimum flow that is not banned
            name_line_min_flow = self.SearchMinFlow(banned)
            if name_line_min_flow < 0:
                logging.getLogger('merlin1975.py').debug(f"not find any line available, I am leaving the solver")
                break
            logging.getLogger('merlin1975.py').debug(
                f"hem trobat un bus de mínims i not banned, with idx={name_line_min_flow}, the network will be reconfigured")
            # Reconfigure the network with the found line
            self.grid.lines[name_line_min_flow].active = False

            #graph = GC_utils.GC2Graph(self.grid)  # Initialize the Minimum Spanning Tree as the network graph      
            NumCycles = len(GC_utils.SearchLoopsLines(self.grid))
            _, connected, _ = GC_utils.CheckRadialConnectedNetwork(self.grid)  # Check if the network remains connected
            #Cycles = list(nx.simple_cycles(graph))
            #NumCycles = len(Cycles)  # Count the number of cycles
            logging.getLogger('merlin1975.py').debug(f"connected : {connected} cycles={NumCycles}")

            if not connected:
                logging.getLogger('merlin1975.py').debug(f"The line {name_line_min_flow} is banned because do not create a tree [cycles {NumCycles}], the network will be reconfigured to the previous state")
                # If the network becomes disconnected, revert to the previous configuration
                self.grid.lines[name_line_min_flow].active = True
                banned.append(name_line_min_flow)  # Add the line to the banned list
                logging.getLogger('merlin1975.py').debug(f"\t\tbanned={banned}")
            else:
                # If the network remains connected, continue optimizing
                mu -= 1  # Decrement mu as a line has been successfully included
                logging.getLogger('merlin1975.py').debug(
                    f"Line {name_line_min_flow} is included in the solution :  ...")
                logging.getLogger('merlin1975.py').debug(f"\t\tNumCycles {NumCycles} - mu={mu}")
                #logging.getLogger('merlin1975.py').debug(f"\t\tCycles {Cycles}")

            logging.getLogger('merlin1975.py').debug(f"finished: {GC_utils.LinesOutofService(self.grid)}")
        return([line.idtag for idx, line in enumerate(self.grid.lines) if not line.active])



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


    # Create an MSTgreedy object
    merlin = Merlin1975(gridGC, verbose_logging=logging.INFO)

    # Solve the Minimum Spanning Tree problem
    disabled_lines = merlin.Solve()
    # Print the list of disabled line indices
    _,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
    radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
    print(f"The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{merlin.NumPF} ")#is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )