import numpy as np
import networkx as nx
import random
import sys
import os
import logging
import warnings
import pandas as pd
import math
import logging

import GC_utils 

import GridCalEngine.api as gce  # For interfacing with the GridCal API
from GridCalEngine.IO.file_handler import FileOpen, FileSave

from gurobipy import *
from Taylor2012_gurobi import MICP_Gurobi

warnings.filterwarnings('ignore')  # Ignore warnings during executionimport numpy as np

class MICP_QP(MICP_Gurobi):
    def __init__(self, grid, bigM=1e6, verbose_logging=logging.WARNING):
        logging.getLogger('micp.py').info("Init from MICP_Gurobi_QP")
        super().__init__(grid, bigM=bigM, verbose_logging=verbose_logging)
        self.Psubstation_i = self.model.addVars(self.SubstationNodes, lb=-GRB.INFINITY, name='Psubstation_i')
        self.Qsubstation_i = self.model.addVars(self.SubstationNodes, lb=-GRB.INFINITY, name='Qsubstation_i')
        self.P_ij = self.model.addVars(self.branch_ij, lb=-GRB.INFINITY, name='P_ij')  #Branch P
        self.Q_ij = self.model.addVars(self.branch_ij, lb=-GRB.INFINITY, name='Q_ij')  #Branch Q
        self.beta_iji = self.model.addVars(self.betas_iji, lb=0, ub=1, name='beta_iji')  #branch status
        self.alpha_ij = self.model.addVars(self.branch_ij, vtype=GRB.BINARY, name='alpha')  #branch flow

    def Solve(self,OutputFlag=1):
        logging.getLogger('micp.py').info("Solve from MICP_Gurobi_QP")
        obj = quicksum((self.P_ij[i, j] ** 2 + self.Q_ij[i, j] ** 2) * next(c for a, b, c in self.r_ij if a == i and b == j)
                for (i, j) in self.branch_ij)
        return super().Solve(obj,OutputFlag)

    def ContraintsDefinition(self):
        logging.getLogger('micp.py').info("ContraintsDefinition from MICP_Gurobi_QP")
        # #Distflow equations (5-6)
        self.model.addConstrs(
            (quicksum(self.P_ij[j,i] for j in self.BusIN[i]) - quicksum(self.P_ij[i,j] for j in self.BusOUT[i]) == self.Pload_i[i] for i in self.BusWithoutSubstation), name='DistFlow-5')
        self.model.addConstrs(
            (quicksum(self.Q_ij[j,i] for j in self.BusIN[i]) - quicksum(self.Q_ij[i,j] for j in self.BusOUT[i]) == self.Qload_i[i] for i in self.BusWithoutSubstation), name='DistFlow-6')
        # #Distflow equations (7-8)
        self.model.addConstrs(
            (quicksum(self.P_ij[i, j] for j in self.BusOUT[i]) == self.Psubstation_i[i] for i in self.SubstationNodes), name='DistFlow-7')
        self.model.addConstrs(
            (quicksum(self.Q_ij[i, j] for j in self.BusOUT[i]) == self.Qsubstation_i[i] for i in self.SubstationNodes), name='DistFlow-8')
        # #Distflow equations (9-10)
        self.model.addConstrs(((self.P_ij[i, j] <= self.bigM * self.beta_iji[i, j]) for (i, j) in self.branch_ij), name='DistFlow-9a')
        self.model.addConstrs(((self.P_ij[i, j] >= 0) for (i, j) in self.branch_ij), name='DistFlow-9b')
        self.model.addConstrs(((self.Q_ij[i, j] <= self.bigM * self.beta_iji[i, j]) for (i, j) in self.branch_ij), name='DistFlow-10a')
        self.model.addConstrs(((self.Q_ij[i, j] >= 0) for (i, j) in self.branch_ij), name='DistFlow-10b')
        ##condition 11 is implicit in the variable definition beta_ij>=0
        self.model.addConstrs(((self.beta_iji[i,j] == 0) for (i, j) in self.branch_ij if j in self.SubstationNodes), name='radiality-12')
        ##condition 13 do not apply, as it is considered that any branch can be disabled
        self.model.addConstrs((self.beta_iji[i, j] + self.beta_iji[j,i] == self.alpha_ij[i, j] for (i, j) in self.branch_ij), name='radiality-14')
        self.model.addConstrs((quicksum(self.beta_iji[j,i] for j in self.BusIN[i]) == 1 for i in self.BusWithoutSubstation),  name='radiality-15')
        ##condition 16 is implicit in the variable definition alpha_ij>=0


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,  # Set the log level to DEBUG
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Set the log format
        datefmt='%Y-%m-%d %H:%M:%S'  # Set the date format
    )

    # Create a network object
    gridGC = FileOpen("D:\\15_Thesis-code\\DistributionNetwork_libraries\\NetworkExamples\\gridcal\\case33.gridcal").open()
    TieLines = ['line 32','line 33','line 34','line 35','line 36']    
    # Create a network object
    #gridGC = FileOpen("D:\\15_Thesis-code\\DistributionNetwork_libraries\\NetworkExamples\\gridcal\\case69.gridcal").open()
    #TieLines = ['line 68','line 69','line 70','line 71','line 72']
    # Create an MSTgreedy object
    micp_qp = MICP_QP(gridGC, bigM=1e6,verbose_logging=logging.DEBUG)

    # Solve the Minimum Spanning Tree problem
    disabled_lines = micp_qp.Solve()

    # Print the list of disabled line indices
    _,loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
    radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
    print(disabled_lines, loss, radiality )





