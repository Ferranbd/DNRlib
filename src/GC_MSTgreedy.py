import numpy as np
import networkx as nx
import random
import sys
import os
import logging
import warnings

import GC_utils 

import GridCalEngine.api as gce  # For interfacing with the GridCal API
from GridCalEngine.IO.file_handler import FileOpen, FileSave

warnings.filterwarnings('ignore')  # Ignore warnings during execution


class MSTgreedy:
    """Class to implement the greedy algorithm for finding the Minimum Spanning Tree (MST).

    Attributes:
        net: The network object containing distribution network data.
    """

    def __init__(self, grid=None, verbose_logging=logging.WARNING):
        """Initializes the MSTgreedy class with required parameters.

        Args:
            net: The network object for optimization.
            verbose_logging: Logging level for debug messages.
        """
        logging.getLogger('mstgreedy.py').setLevel(verbose_logging)  # Set logging level
        self.grid = grid  # Store the network object
        self.NumPF=0
        
    def __powerflow(self, config=None):
            self.NumPF+=1
            resultPF,old_losses = GC_utils.GC_PowerFlow(self.grid, config=config)
            return resultPF,old_losses
    
    def Solve(self, algorithm="kruskal", randomMST=False, one=False, current_power=False, inverted=False):
        """Solves the Minimum Spanning Tree problem using a greedy algorithm.

        Args:
            randomMST: If True, assigns random weights to edges.
            algorithm: The algorithm to use for finding the MST ('kruskal', 'prim', or 'boruvka').
            one: If True, assigns a weight of 50 to all edges.
            current_power: If True, weights are based on current power values.

        Returns:
            List of disabled line indices in the network after reconfiguration.
        """
        logging.getLogger('mstgreedy.py').info("Start solving MSTgreedy")

        for line in self.grid.lines:
            line.active = True
        for trafo in self.grid.transformers2w:
            trafo.active = True

        gridGraph = GC_utils.GC2Graph(self.grid)  # Initialize the Minimum Spanning Tree as the network graph
        
        if randomMST:  # If randomMST is True, assign random weights to edges
            #print("randomMST")
            for _, (_, _, attributes) in enumerate(gridGraph.edges(data=True)):
                attributes['weight'] = random.random() * 100.0  # Random weight between 0 and 100
            MinSpanningTree = nx.minimum_spanning_tree(gridGraph, weight='weight',
                                                       algorithm=algorithm)  # Calculate MST with random weights

        elif one:  # If one is True, assign a weight of 50 to all edges
            #print("one")
            for _, (_, _, attributes) in enumerate(gridGraph.edges(data=True)):
                attributes['weight'] = 50
            MinSpanningTree = nx.minimum_spanning_tree(gridGraph, weight='weight',
                                                       algorithm=algorithm)  # Calculate MST with weight 50

        else:  # If a specific algorithm is provided
            #print("current/power")
            logging.getLogger('mstgreedy.py').debug(f"current_power: {current_power}")
            power_flow, loss = self.__powerflow()

            #options = gce.PowerFlowOptions(gce.SolverType.NR, verbose=False)
            #power_flow = gce.PowerFlowDriver(self.grid, options)
            #power_flow.run()
            maxCurrent = power_flow.results.If.real.max()  # Get the maximum current
            maxLosses = power_flow.results.losses.real.max()  # Get the maximum losses
            for _, (u, v, attributes) in enumerate(gridGraph.edges(data=True)):
                # Log the current values for debugging
                logging.getLogger('mstgreedy.py').debug("line with buses: %s %s", u, v)

                if current_power:  # If current_power is True, set weights based on current
                    current_tmp1 = [power_flow.results.losses.real[idx].real for idx, line in enumerate(self.grid.lines) 
                            if (((line.bus_to.idtag==v) & (line.bus_from.idtag==u)) | (line.bus_to.idtag==u) & (line.bus_from.idtag==v)) ]
                    current_tmp2 = [power_flow.results.losses.real[idx].real for idx, line in enumerate(self.grid.transformers2w) 
                            if (((line.bus_to.idtag==v) & (line.bus_from.idtag==u)) | (line.bus_to.idtag==u) & (line.bus_from.idtag==v)) ]
                    current = list(set(current_tmp1).union(current_tmp2))
                    if len(current)>1:
                        current = [sum(current)]
                    attributes['weight'] = 100 * maxCurrent / current
                else:  # Otherwise, set weights based on losses
                    losses_tmp1 = [power_flow.results.If[idx].real for idx, line in enumerate(self.grid.lines)
                            if (((line.bus_to.idtag==v) & (line.bus_from.idtag==u)) | (line.bus_to.idtag==u) & (line.bus_from.idtag==v)) ]
                    losses_tmp2 = [power_flow.results.If[idx].real for idx, line in enumerate(self.grid.transformers2w)
                            if (((line.bus_to.idtag==v) & (line.bus_from.idtag==u)) | (line.bus_to.idtag==u) & (line.bus_from.idtag==v)) ]
                    losses = list(set(losses_tmp1).union(losses_tmp2))
                    if len(losses)>1:
                        losses = [sum(losses)]
                    attributes['weight'] = 100 * maxLosses / losses
                if inverted:
                    attributes['weight']=1/attributes['weight']
            # Calculate the Minimum Spanning Tree using the defined weights
            MinSpanningTree = nx.minimum_spanning_tree(gridGraph, weight='weight', algorithm=algorithm)

        # Return the list of disabled line indices after reconfiguration
        #print("MinSpanningTree: ", MinSpanningTree)


        # Set all lines to inactive
        self.grid = GC_utils.Graph2GC(MinSpanningTree, self.grid)

        #for line in self.grid.lines:
        #    line.active = False
        # Loop through each edge in the graph (each edge corresponds to a line between two buses).
        #for _, (u, v, attributes) in enumerate(MinSpanningTree.edges(data=True)):
            # Enable the line where Bus1 is connected to Bus2 (i.e., the direction from Bus1 to Bus2).
        #    for idx, line in enumerate(self.grid.lines) :
        #        if (((line.bus_to.name==v) & (line.bus_from.name==u)) | (line.bus_to.name==u) & (line.bus_from.name==v)):
        #            line.active=True
        #            #print("line with buses",line.bus_to.name,line.bus_from.name, "compared to u,v:",u,v)
        #            continue
            

        return GC_utils.LinesOutofService(self.grid)


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

    mstgreedy = MSTgreedy(gridGC)
    disabled_lines = mstgreedy.Solve(randomMST=True)
    _,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
    radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
    print(f"randomMST The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{mstgreedy.NumPF} ") #is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )

    mstgreedy = MSTgreedy(gridGC)
    disabled_lines = mstgreedy.Solve(randomMST=True, algorithm="prim")
    _,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
    radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
    print(f"randomMST/Prim The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{mstgreedy.NumPF} ") #is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )

    mstgreedy = MSTgreedy(gridGC)
    disabled_lines = mstgreedy.Solve(one=True)
    _,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
    radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
    print(f"One The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{mstgreedy.NumPF} ") #is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )

    mstgreedy = MSTgreedy(gridGC)
    disabled_lines = mstgreedy.Solve(one=True, algorithm="prim")
    _,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
    radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
    print(f"One/Prim The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{mstgreedy.NumPF} ") #is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )

    mstgreedy = MSTgreedy(gridGC)
    disabled_lines = mstgreedy.Solve(current_power=False, algorithm="kruskal")
    _,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
    radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
    print(f"Current/Kruskal The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{mstgreedy.NumPF} ") #is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )
    
    mstgreedy = MSTgreedy(gridGC)
    disabled_lines = mstgreedy.Solve(current_power=False, algorithm="prim")
    _,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
    radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
    print(f"Current/Prim The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{mstgreedy.NumPF} ") #is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )
    
    mstgreedy = MSTgreedy(gridGC)
    disabled_lines = mstgreedy.Solve(current_power=True, algorithm="prim")
    _,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
    radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
    print(f"Power/Prim The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{mstgreedy.NumPF}") #is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )