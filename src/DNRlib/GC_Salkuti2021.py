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

## based on Saltuki 2021
class Saltuki2021:
    """Class to implement the algorithm described in Saltuki 2021 for network reconfiguration.

    Attributes:
        net: The network object containing distribution network data.
        initial_config: The initial configuration of tie lines in the network.
    """

    def __init__(self, grid= None, TieLines= None, init_config= None, verbose_logging=logging.WARNING):
        """Initializes the Saltuki2021 class with required parameters.

        Args:
            net: The network object for optimization.
            init_config: The initial configuration of tie lines.
            verbose_logging: Logging level for debug messages.
        """
        self.grid = grid  # Store the network object
        logging.getLogger('saltuki2021.py').setLevel(verbose_logging)  # Set logging level
        self.initial_config = TieLines # Store initial configuration of tie lines
        self.NumPF=0
        
    def __powerflow(self, config=None):
            self.NumPF+=1
            resultPF,old_losses = GC_utils.GC_PowerFlow(self.grid, config=config)
            return resultPF,old_losses


    def _searchLoopLine(self, bus, openedTieLine):
        """Search for a line connected to the given bus that is not already in the tie lines.

        Args:
            bus: The bus object containing bus information.
            TieLines: List of currently used tie lines.
            openedTieLine: The tie line that has been opened.

        Returns:
            The line object connected to the bus that is not in TieLines.
        """
        logging.getLogger('saltuki2021.py').debug(f"to be found bus={bus}, with previous line={openedTieLine}")
        lines=[]
        
        for line in self.grid.lines:
            if line.bus_from==bus and line.idtag != openedTieLine:
                lines.append(line.idtag)
            
        for line in self.grid.lines:
            if line.bus_to==bus and line.idtag != openedTieLine:
                lines.append(line.idtag)
            
        logging.getLogger('saltuki2021.py').debug(f"found {lines}")
        return lines

    def FindBusVoltages(self, target_line, resultPF):
        logging.getLogger('saltuki2021.py').debug(f"FindBusVoltages:{target_line}")
        vsender=0
        vreceiver=0
        for line in self.grid.lines:
            if line.idtag == target_line:
                busFrom = line.bus_from
                busTo = line.bus_to
                vreceiver = resultPF.results.get_bus_df().loc[busFrom.name,'Vm']
                vsender = resultPF.results.get_bus_df().loc[busTo.name,'Vm']
                logging.getLogger('saltuki2021.py').debug(f"line: {line.idtag} busF,busT: {line.bus_from} {line.bus_to} vr,vs= {vreceiver},{vsender}")
                break
        return line, vsender, vreceiver

    def Solve(self, order=None):
        """Solves the network reconfiguration problem based on Saltuki's approach.

        Args:
            order: The order of processing tie lines (not currently used).

        Returns:
            List of tie lines in the new configuration after solving.
        """
        logging.getLogger('saltuki2021.py').info("Start solving Saltuki 2021")
        newTieLines = self.initial_config  # Initialize new tie lines with the initial configuration
        logging.getLogger('saltuki2021.py').debug(f"with tielines:{self.initial_config} - {newTieLines}")

        for tie in newTieLines:  # Iterate over each tie line in the initial configuration
            prevTieLines = newTieLines  # Keep track of the previous configuration
            resultPF,old_losses = self.__powerflow(config=newTieLines)
            gctie, vsender,vreceiver = self.FindBusVoltages(tie, resultPF=resultPF)
            # Decision-making based on voltage levels
            if vsender > vreceiver:     #Add tie switch to the sending end node. Open the branch that is connected to the receiving end node
                possibleNewTieLines = self._searchLoopLine(gctie.bus_from, tie)
            else:                       #Add tie switch to the receiving end node. Open the branch that is connected to the sending end node
                possibleNewTieLines = self._searchLoopLine(gctie.bus_to, tie)

            logging.getLogger('saltuki2021.py').debug(f"   the possible TieLines are  {possibleNewTieLines}")
            for newtieline in possibleNewTieLines:
                newTieLines = [newtieline if x == tie else x for x in newTieLines]  # Update tie lines
                logging.getLogger('saltuki2021.py').debug(f"   newtieline vs>vr ={newtieline} ==> Newtielines :{newTieLines}")
                # Run the power flow analysis and calculate losses
                GC_utils.NetworkReconfiguration(self.grid, all=True, value_all=True, selected_configuration=newTieLines, value_configuration=False)
                radiality, _, _ = GC_utils.CheckRadialConnectedNetwork(self.grid)  # Check if the network remains radial
                # If the new configuration is not valid (non-radial), revert to the previous configuration
                if not radiality:
                    logging.getLogger('saltuki2021.py').debug(f"   new config is not valid, let's go back to the previous TieLines")
                    newTieLines = prevTieLines
                else:
                    # Update voltage values for further calculations
                    resultPF, losses = self.__powerflow(config=newTieLines)
                    gctie, vsender,vreceiver = self.FindBusVoltages(newtieline, resultPF=resultPF)
                    logging.getLogger('saltuki2021.py').debug(f"   valid new configuration : {newTieLines} losses{losses} old_losses{old_losses}")
                    prevTieLines = newTieLines  # Update previous tie lines for the next iteration

                    # Continue to optimize while sending voltage is greater than receiving voltage
                    while (vsender > vreceiver) and (losses > old_losses):
                        # Open the branch connected to the receiving end node
                        possibleNewTieLines = self._searchLoopLine(gctie.bus_from, tie)
                        for newtieline in possibleNewTieLines:
                            logging.getLogger('saltuki2021.py').debug(f"   *newtieline={newtieline}")
                            newTieLines = [newtieline if x == tie else x for x in newTieLines]  # Update tie lines
                            logging.getLogger('saltuki2021.py').debug(f"   *new tie lines :{newTieLines} with losses={losses:.2f}")
                            old_losses = losses  # Update old losses
                           # Run power flow analysis again
                            GC_utils.NetworkReconfiguration(self.grid, all=True, value_all=True, selected_configuration=newTieLines, value_configuration=False)
                            radiality, _, _ = GC_utils.CheckRadialConnectedNetwork(self.grid)  # Check if the network is still radial
                            # If the new configuration is not valid, revert to the previous configuration
                            if not radiality:
                                logging.getLogger('saltuki2021.py').debug(f"   new config is not valid, let's go back to the previous TieLines")
                                newTieLines = prevTieLines

        return list(newTieLines)  # Return the final list of tie lines in the new configuration



    
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
    print(f"Original network:{loss}, radiality:{radiality} ") #is {GC_utils.GC_Line_idtag2name_array(gridGC,TieLinesID)}" )

    # Create an MSTgreedy object
    salkuti = Saltuki2021(gridGC, TieLines=TieLinesID,verbose_logging=logging.ERROR)

    # Solve the Minimum Spanning Tree problem
    disabled_lines = salkuti.Solve()

    # Print the list of disabled line indices
    _,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
    radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
    print(f"The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{salkuti.NumPF} ")#is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )