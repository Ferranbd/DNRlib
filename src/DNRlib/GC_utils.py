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
from GridCalEngine.Simulations.PowerFlow.power_flow_worker import multi_island_pf_nc
from GridCalEngine.enumerations import ExternalGridMode

warnings.filterwarnings('ignore')  # Ignore warnings during execution

def GC_FitnessCalculation(grid, totalLoad, loss_factor=2, voltage_factor=1, split_factor=1):
	#Obj_tensió = (MaxAllowedVoltage - minVoltage@iteration) / MaxAllowedVoltage 
	#   if minVoltage == MaxVoltage --> ObjVoltage = 0 : Perfect
	#   if minVoltage == 0 --> ObjVoltage = 1 : wrong
	#Obj_Losses = Losses / 0,08*Loads 
	#   if Losses == 0,08*Loads --> Obj_Losses = 1 : bad
	#   if Losses == 0 --> Obj_Losses =0 : good
    #Where 0.08 is a parameter to be adjusted

    #calculate the power flow without reconfiguring the network
    #calculate the max and min voltage
    #calculate the losses
    res, loss  = GC_PowerFlow(grid)

    AbsVmax = max([grid.buses[i].Vmax for i in range(len(grid.buses))])
    Vmin = res.results.get_bus_df().Vm.min()
    objVoltage = abs((AbsVmax-Vmin))*voltage_factor

    objLoss = loss / totalLoad * loss_factor

    fitness = split_factor * objLoss + (1-split_factor) * objVoltage
    return fitness*100
    #return loss


def GC2Graph(grid, idtag_name=False, weights = None):
    # Create a NetworkX graph
    G = nx.Graph()

    # Add nodes (buses) to the graph
    for bus in grid.buses:
        if idtag_name :
            G.add_node(bus.name, voltage=bus.Vnom)
        else :
            G.add_node(bus.idtag, voltage=bus.Vnom)

    # Add edges (lines) to the graph
    for idx,line in enumerate(grid.lines):
        if line.active:
            if weights == None:
                weight=0
            else:
                print(weights[0], line.name)
                weight=weights[line.name]
            if idtag_name :
                G.add_edge(line.bus_from.name, line.bus_to.name, r=line.R, x=line.X, length=line.length, type="line", weight=weight)
            else :
                G.add_edge(line.bus_from.idtag, line.bus_to.idtag, r=line.R, x=line.X, length=line.length, type="line", weight=weight)
    for trafo in grid.transformers2w:
        if trafo.active:
            if idtag_name :
                G.add_edge(trafo.bus_from.name, trafo.bus_to.name, r=trafo.R, x=trafo.X, type="trafo")
            else :
                G.add_edge(trafo.bus_from.idtag, trafo.bus_to.idtag, r=trafo.R, x=trafo.X, type="trafo")
    return G

def Graph2GC(graph, grid):
    # Set all lines to inactive
    for line in grid.lines:
        line.active = False
    for trafo in grid.transformers2w:
        trafo.active = False
    # Loop through each edge in the graph (each edge corresponds to a line between two buses).
    for _, (u, v, attributes) in enumerate(graph.edges(data=True)):
        # Enable the line where Bus1 is connected to Bus2 (i.e., the direction from Bus1 to Bus2).
        for idx, line in enumerate(grid.lines) :
            if (((line.bus_to.idtag==v) & (line.bus_from.idtag==u)) | (line.bus_to.idtag==u) & (line.bus_from.idtag==v)):
                line.active=True
                #print("line with buses",line.bus_to.name,line.bus_from.name, "compared to u,v:",u,v)
                continue
        for idx, trafo in enumerate(grid.transformers2w) :
            if (((trafo.bus_to.idtag==v) & (trafo.bus_from.idtag==u)) | (trafo.bus_to.idtag==u) & (trafo.bus_from.idtag==v)):
                trafo.active=True
                #print("line with buses",line.bus_to.name,line.bus_from.name, "compared to u,v:",u,v)
                continue
    return grid

def GCNC_GenerateNC(grid):
    nc = gce.compile_numerical_circuit_at(grid)
    opt = gce.PowerFlowOptions(solver_type=gce.SolverType.NR)
    res_nc = multi_island_pf_nc(nc, options=opt)

    return nc, res_nc

def GCNC_CreateBranchDF(nc):
    data = pd.DataFrame({
        'index':nc.passive_branch_data.names,
        'bus_from':nc.passive_branch_data.F,
        'bus_to':nc.passive_branch_data.T,
        'R':nc.passive_branch_data.R,
        'X':nc.passive_branch_data.X,
        'B':nc.passive_branch_data.B,
        'active': nc.passive_branch_data.active,
    }, copy=False)

    data.set_index('index', inplace=True)
    return data

def GCNC_CreateBusDF(nc):
    data = pd.DataFrame({
        'index':nc.bus_data.names,
        'Vnom':nc.bus_data.Vnom,
    }, copy=False)

    data.set_index('index', inplace=True)
    return data

def GC_PowerFlow(grid, config=None):
    if config:
        if len(config)>0:
            NetworkReconfiguration(grid, all=True, value_all=True, selected_configuration=config, value_configuration=False)
    
    options = gce.PowerFlowOptions(gce.SolverType.NR, verbose=False)
    power_flow = gce.PowerFlowDriver(grid, options)
    power_flow.run()
    return power_flow, power_flow.results.losses.sum().real

def LinesOutofService(grid):
    return([line.idtag for idx, line in enumerate(grid.lines) if not line.active] + [line.idtag for idx, line in enumerate(grid.transformers2w) if not line.active])

def EdgeCycles(graph):

#    cycles = nx.cycle_basis(graph)
#    cycles_edges = []
#
#    for cycle in cycles:
#        cycle_edges = []
#        for i in range(len(cycle)):
#            edge = (cycle[i], cycle[(i + 1) % len(cycle)])
#            if graph.has_edge(*edge):
#                cycle_edges.append(edge)
#            else:
#                cycle_edges.append((edge[1], edge[0]))  # Add the edge in the correct direction
#        cycles_edges.append(cycle_edges)
#
#    return cycles_edges

    cycles = nx.cycle_basis(graph)
    cycles_edges = []

    for cycle in cycles:
        cycle_edges = []
        for edge in graph.edges:
            if edge[0] in cycle and edge[1] in cycle:
                cycle_edges.append(edge)
        cycles_edges.append(cycle_edges)
    return cycles_edges


def GC_Line_Name2idtag_single(grid, LineName):
    for line in grid.lines:
        if line.name == LineName:
            return line.idtag
    for trafo in grid.transformers2w:
        if trafo.name == LineName:
            return trafo.idtag

def GC_Line_idtag2name_single(grid, idtag):
    for line in grid.lines:
        if line.idtag == idtag:
            return line.name
    for trafo in grid.transformers2w:
        if trafo.idtag == idtag:
            return trafo.name
def GC_Line_Name2idtag_array(grid, LineNames):
    LineIdtags=[]
    for name in LineNames:
        LineIdtags.append(GC_Line_Name2idtag_single(grid,name))
    return LineIdtags

def GC_Line_idtag2name_array(grid, idtags):
    LineNames=[]
    for idtag in idtags:
        LineNames.append(GC_Line_idtag2name_single(grid,idtag))
    return LineNames



def remove_repeated_elements(list_of_lists):
    seen = set()
    for sublist in list_of_lists:
        i = 0
        while i < len(sublist):
            if sublist[i] in seen:
                sublist.pop(i)
            else:
                seen.add(sublist[i])
                i += 1
    return list_of_lists

def SearchLoopsLines(grid, flagUsed=False, shuffle=False):
    G = GC2Graph(grid)
    loops = EdgeCycles(G)
    solutionA=[]
    used=set()
    if (shuffle):
        random.shuffle(loops)
    for loop in loops:
        lines_in_loop=[]
        for edge in loop:
            for line in grid.lines:
                if line.bus_from.idtag in edge and line.bus_to.idtag in edge:
                    lines_in_loop.append(line.idtag)
            for line in grid.transformers2w:
                if line.bus_from.idtag in edge and line.bus_to.idtag in edge:
                    lines_in_loop.append(line.idtag)
                    #print("searchloops, trafo found:",line.idtag)
        solutionA.append(lines_in_loop)
    if (flagUsed):
        solutionB = remove_repeated_elements(solutionA)
        solutionB = [sublist for sublist in solutionB if sublist]
        return solutionB
    return solutionA

def CheckRadialConnectedNetwork(grid, units=False):
    """Check if the network is both radial and connected.

    This method verifies if the network represented by the graph is a radial network
    and whether it is connected. A radial network is defined as a tree structure
    (i.e., acyclic and connected).

    Returns:
        tuple: A tuple containing three boolean values:
            - RadialConnected (bool): True if the network is both radial and connected, False otherwise.
            - connected (bool): True if the network is connected, False otherwise.
            - radial (bool): True if the network is a tree (radial), False otherwise.
    """
    graph = GC2Graph(grid)

    if units:
        isolated_nodes = list(nx.isolates(graph))
        num_isolated_nodes = len(isolated_nodes)
        cycles = nx.cycle_basis(graph)
        num_cycles = len(cycles)
        RadialConnected=num_isolated_nodes==0 and num_cycles==0
        return RadialConnected, num_isolated_nodes, num_cycles 
    else:
        # Check if the graph is a tree (radial)
        radial = nx.is_tree(graph)
        # Check if the graph is connected
        connected = nx.is_connected(graph)
        # Determine if the network is radial and connected
        RadialConnected = connected and radial
        return RadialConnected, connected, radial 

def find_value_in_arrays(self, arrays, value):
    """
    Finds the index and array containing a specific value.

    Args:
        arrays (list): A list of arrays to search.
        value: The value to search for.

    Returns:
        tuple: A tuple containing the index of the array (if found) and the array itself, or None if not found.
    """
    for i, array in enumerate(arrays):
        if value in array:
            return i, array  # Return the index and array if found
    return None, None  # Return None if not found

# Function to substitute the value in the array
def substitute_value(self, array, old_value, new_value):
    """
    Substitutes a value in an array with a new value.

    Args:
        array (list): The array to modify.
        old_value: The value to be replaced.
        new_value: The new value to replace with.

    Returns:
        list: The modified array with the old value replaced.
    """
    array[array == old_value] = new_value  # Replace the old value with the new value
    return array

def NetworkReconfiguration(grid, all=False, selected_configuration=None, value_all=False, value_configuration=False):
        """ Network reconfiguration : sets the line state to on/off depending on the arguments. The changes are applied to the DistributionNetwork structure, to the internal graph and to the original source of information (openDSS, GridCal, PandaPower)

        Args:
            all=False : if True all the lines will be set to the value in "value_all", if False the line status will remain as it is.
            selected_configuration=None : if different of None, the values in this list will be set to the value in "value_configuration"
            value_configuration=False
            value_all=False
            extraData=True : this atributte indicates if the functions :DistributionNetwork2Graph(), LinesOutService() and NetworkSummary() will be called
        Returns:
            None
        """
        if all == True:
            for line in grid.lines:
                line.active=value_all
            for trafo in grid.transformers2w:
                trafo.active=value_all
        if selected_configuration != None:
            for line in grid.lines:
                if line.idtag in selected_configuration:
                    line.active=value_configuration
            for trafo in grid.transformers2w:
                if trafo.idtag in selected_configuration:
                    trafo.active=value_configuration
                               # Function to check if a list belongs to any of the sublists within a list of lists

def list_belongs_to_list_of_lists(list_to_check, list_of_lists):
    """Checks if the given list is a subset of any list within a list of lists.

    Args:
        list_to_check: The list to check.
        list_of_lists: The list of lists to check against.

    Returns:
        True if the list is a subset of any list in the list of lists, otherwise False.
    """
    if any(all(item in sublist for item in list_to_check) for sublist in list_of_lists):
        return True
    else:
        return False

# Function to remove duplicates from a list of lists
def remove_duplicates(list_of_lists):
    """Removes duplicate items from each sublist in a list of lists.

    Args:
        list_of_lists: The list of lists from which to remove duplicates.

    Returns:
        The modified list of lists with duplicates removed.
    """
    unique_elements = []  # List to keep track of unique elements
    for sub_list in list_of_lists:
        unique_sub_list = []  # Initialize a unique sublist
        for item in sub_list:
            if item not in unique_elements:
                unique_sub_list.append(item)  # Add unique item to the sublist
                unique_elements.append(item)  # Track the unique item globally
        sub_list[:] = unique_sub_list  # Replace the original sublist with the unique one
    return list_of_lists  # Return the modified list of lists


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,  # Set the log level to DEBUG
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Set the log format
        datefmt='%Y-%m-%d %H:%M:%S'  # Set the date format
    )

    gridGC = FileOpen("D:\\15_Thesis-code\\DistributionNetwork_libraries\\NetworkExamples\\gridcal\\case33.gridcal").open()
    #Tie=['e97b6c0ecdd8456d8d5dc2c9cdbe4b01', '0a64b3f458e2483d9448104540a57da1', 'c685fb815f7f4331b230b2e164bcde2c', 'cf980a54675a4d52ba30436d24fbd529', 'c86eb3e7f72d4735b656c343c3d50061']
    #_, losses = GC_utils.GC_PowerFlow(gridGC, config=[])
    #print("losses=", losses)
    NetworkReconfiguration(gridGC, all=True, value_all=True)

    shuffle =False
    flagUsed = True

    solution = SearchLoopsLines(gridGC, flagUsed=flagUsed, shuffle=shuffle)


    print("solutions", solution[0])
    #for loop in solution:
    #    print(loop)