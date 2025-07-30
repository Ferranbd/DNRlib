import networkx as nx
import pandas as pd
import logging
import numpy as np
import random
import sys

from GC_Baran1989 import Baran1989
from GC_Jakus2020 import Jakus2020
from GC_Khalil_Gorpinich2012 import BPSO
from GC_Merlin1975 import Merlin1975
from GC_Salkuti2021 import Saltuki2021
from GC_MSTgreedy import MSTgreedy
from GC_Morton2000 import Morton2000
from GC_Taylor2012_pyomo import MICP_Pyomo as MICP_Pyomo

import GC_utils

import GridCalEngine.api as gce  # For interfacing with the GridCal API
from GridCalEngine.IO.file_handler import FileOpen, FileSave



class DistributionNetworkReconfiguration:
    '''
    Library of Distribution Network Reconfiguration algorithms
    Constraints : radial topology without unconnected nodes
    Available algorithms :
    - Heuristic : Baran (1989), Merlin (1975), Morton (2000), MSTgreedy, Salkuti (2021)
    - Metaheuristic : Jakus (2020), Khalil-Gorpinich (2012)
    - Mathematical : Jabr-Taylor (2012), Taylor-gurobi (2012), Taylor-pyomo 2012)

    The input is the grid
    The output is the list of disabled lines

    example:
    mst = sp.src(G)
    mst_Prim_edges = mst.Solve('Prim')
    '''
    def __init__(self, grid=None, verbose_logging=logging.INFO):
        '''
        to initialize the class requires the grid format (openDSS, GridCal or PandaPower) and the grid file name
        '''
        self.verbose_logging = verbose_logging
        self.grid = grid  #dn class
        self.NumPF = 0
        logging.getLogger('DistributionNetworkReconfiguration.py').debug(f"init DNR")

    def _Merlin1975(self, order=None):
        merlin = Merlin1975(grid=self.grid, verbose_logging=self.verbose_logging)
        DisabledLines = merlin.Solve()
        self.NumPF = merlin.NumPF
        return DisabledLines

    def _Baran1989(self, order=None, TieLines= None):
        baran = Baran1989(grid=self.grid, TieLines=TieLines, verbose_logging=self.verbose_logging)
        DisabledLines = baran.Solve(order=order)
        self.NumPF = baran.NumPF
        return DisabledLines

    def _Jakus2020(self, TieLines = None, PopulationSize=32, MutationProbability=0.02, PopulationSize_SBEA=None, ElitePopulation=None, fitness_ratio=1, loss_factor=0.08, Niter=8):
        self._jakus = Jakus2020(grid=self.grid, TieLines=TieLines, PopulationSize=PopulationSize, MutationProbability=MutationProbability,
                          ElitePopulation=ElitePopulation, Niter=Niter, fitness_ratio=fitness_ratio, loss_factor=loss_factor, verbose_logging=self.verbose_logging)
        #_jakus = Jakus2020(grid=self.grid,Niter=8, TieLines=TieLines, verbose_logging=self.verbose_logging)
        DisabledLines = self._jakus.Solve() 
        self.NumPF = self._jakus.NumPF
        return DisabledLines

    def _Khalil2012(self, TieLines=None, NumCandidates=10, maxiter=20):
        _bpso = BPSO(grid=self.grid, NumCandidates=NumCandidates, TieLines=TieLines, verbose_logging=self.verbose_logging)
        _bpso.config()
        DisabledLines = _bpso.Solve(wmax=1.2, wmin=0.9, beta=0.05, maxiter=maxiter, vmax=100, c1=2, c2=2)
        self.NumPF = _bpso.NumPF
        return DisabledLines

    def _Salkuti2021(self, TieLines=None):
        salkuti = Saltuki2021(grid=self.grid, TieLines=TieLines, verbose_logging=self.verbose_logging)
        DisabledLines = salkuti.Solve()
        self.NumPF = salkuti.NumPF
        return DisabledLines

    def _MSTgreedy(self, randomMST=False, algorithm=None, one=True, current_power=False):
        mstgreedy = MSTgreedy(grid=self.grid, verbose_logging=self.verbose_logging)
        DisabledLines = mstgreedy.Solve(randomMST=randomMST, algorithm=algorithm, one=one, current_power=current_power)
        self.NumPF = mstgreedy.NumPF
        return DisabledLines

    def _Taylor(self, algorithm="SOC", solver="IPOPT", bigM=1e8, Imax=0, vmin=0, vmax=1, TieLines=None):
        self.Vbase = self.grid.buses[0].Vnom
        self.Zbase = self.grid.Sbase / self.Vbase**2
        for idx, load in enumerate(self.grid.loads):
            self.grid.loads[idx].P = self.grid.loads[idx].P / self.grid.Sbase
            self.grid.loads[idx].Q = self.grid.loads[idx].Q / self.grid.Sbase
        for idx, line in enumerate(self.grid.lines):
            self.grid.lines[idx].R = self.grid.lines[idx].R / self.Zbase
            self.grid.lines[idx].X = self.grid.lines[idx].X / self.Zbase

        micp = MICP_Pyomo(grid=self.grid, NumTieLines=len(TieLines), algorithm=algorithm, bigM=bigM, Imax = Imax, vmin = vmin, vmax = vmax, verbose_logging=self.verbose_logging)
        convergence, MICPBestConfiguration = micp.Solve(solver=solver)

        for idx, load in enumerate(self.grid.loads):
            self.grid.loads[idx].P = self.grid.loads[idx].P * self.grid.Sbase
            self.grid.loads[idx].Q = self.grid.loads[idx].Q * self.grid.Sbase
        for idx, line in enumerate(self.grid.lines):
            self.grid.lines[idx].R = self.grid.lines[idx].R * self.Zbase
            self.grid.lines[idx].X = self.grid.lines[idx].X * self.Zbase
        return MICPBestConfiguration

    def _Morton2000(self, TieLines):
        GC_utils.NetworkReconfiguration(self.grid,all=True, selected_configuration=None, value_all=True, value_configuration=False)
        initial_loops = GC_utils.SearchLoopsLines(self.grid, flagUsed=False)
        print(initial_loops)
        self.morton = Morton2000(grid=self.grid, init_config=TieLines, verbose_logging=self.verbose_logging)
        disabled_lines = self.morton.Solve()
        self.NumPF = self.morton.NumPF
        return disabled_lines

    def Solve(self, method='Baran', *args, **kwargs):
        '''
        once the class is initialize, the only available function is the Solve function, which obtains the MST based
        on the algorithm given by the argument 'method', which is 'Kruskal' by default
        '''
        if method.lower() == 'merlin':
            return self._Merlin1975()
        elif method.lower() == 'baran':
            return self._Baran1989(TieLines=kwargs['TieLines'])
        elif method.lower() == 'jakus':
            if 'ElitePopulation' in kwargs:
                ElitePopulation = kwargs['ElitePopulation']
            else:
                ElitePopulation = int(round(kwargs['PopulationSize'] / 4, 0))
            return self._Jakus2020(PopulationSize=kwargs['PopulationSize'], MutationProbability=kwargs['MutationProbability'],
                    TieLines=kwargs['TieLines'], ElitePopulation=ElitePopulation, Niter=kwargs['Niter'], fitness_ratio=kwargs['fitness_ratio'], loss_factor=kwargs['loss_factor'])
        elif method.lower() == 'khalil':
            return self._Khalil2012(NumCandidates=kwargs['NumCandidates'], TieLines=[])
        elif method.lower() == 'salkuti':
            return self._Salkuti2021(TieLines=kwargs['TieLines'])
        elif method.lower() == 'mstgreedy':
            return self._MSTgreedy(randomMST=kwargs['randomMST'], algorithm=kwargs['algorithm'], one=kwargs['one'], current_power=kwargs['current_power'])
        elif method.lower() == 'taylor':
            return self._Taylor(algorithm=kwargs['algorithm'], solver=kwargs['solver'], bigM=kwargs['bigM'], Imax=kwargs['Imax'], vmin=kwargs['vmin'], vmax=kwargs['vmax'], TieLines=kwargs['TieLines'])
        elif method.lower() == 'morton':
            return self._Morton2000(TieLines=kwargs['TieLines'])



 

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
        import simbench as sb
        import pandapower as pp
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

    dnr = DistributionNetworkReconfiguration(gridGC, verbose_logging=logging.ERROR)

    try:
        dnr.NumPF=0
        disabled_lines = dnr.Solve(method="Merlin", TieLines=TieLinesID)
        _, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        print(f"Merlin: The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{dnr.NumPF} is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )
    except:
        print("falla Merlin")        

    try:
        dnr.NumPF=0
        disabled_lines = dnr.Solve(method="Baran", TieLines=TieLinesID)
        _, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        print(f"Baran: The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{dnr.NumPF} is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )
    except:
        print("falla Baran")        

    try:
        dnr.NumPF=0
        disabled_lines = dnr.Solve(method="Salkuti", TieLines=TieLinesID)
        _, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        print(f"Salkuti: The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{dnr.NumPF} is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )
    except:
        print("falla Salkuti")        

    try:
        dnr.NumPF=0
        disabled_lines = dnr.Solve(method="MSTgreedy", randomMST=False, one=False, current_power=True, algorithm="prim")
        _, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        print(f"MSTgreedy: The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{dnr.NumPF} is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )
    except:
        print("falla MSTGreedy")        
        
    try:
        dnr.NumPF=0
        disabled_lines = dnr.Solve(method="Khalil", NumCandidates=10)
        _, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        print(f"Khalil: The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{dnr.NumPF} is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )
    except:
        print("falla Khalil")

    try:
        dnr.NumPF=0
        disabled_lines = dnr.Solve(method="Jakus", MutationProbability=0.4, PopulationSize=16, Niter=20, ElitePopulation=2, fitness_ratio=1.0, loss_factor=0.08, TieLines=TieLinesID)
        _, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        print(f"Jakus: The new optimal configuration losses:{loss}, radiality:{radiality}, numPF:{dnr.NumPF} is {GC_utils.GC_Line_idtag2name_array(gridGC, disabled_lines)}" )
    except:
        print("falla Jakus")
    
    try:
        if case==1:
            disabled_lines_name = ['line 6', 'line 8', 'line 13', 'line 31', 'line 36']
            disabled_lines = GC_utils.GC_Line_Name2idtag_array(gridGC, disabled_lines_name)
        #disabled_lines = dnr.Solve(method="Morton", TieLines=TieLines)
        _, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        print("Morton: ", GC_utils.GC_Line_idtag2name_array(gridGC,disabled_lines), loss, radiality )   
    except:
        print("falla Morton")

