import numpy as np
import random
import math
from gurobipy import *

import logging


class MICP_Gurobi:
    def __init__(self, grid, bigM=1e6, Imax=0, vmin=0, vmax=1, verbose_logging=logging.WARNING):
        logging.getLogger('micp.py').setLevel(verbose_logging)
        logging.getLogger('micp.py').info("init from MICP_Gurobi")
        #MICP parameters
        self.grid = grid
        # Configuration parameters
        self.bigM = bigM
        self.Imax = Imax
        # read and calculated bases
        self.basekV = grid.buses[0].Vnom
        self.baseMVA = grid.Sbase
        self.baseI = self.baseMVA / self.basekV
        self.baseZ = self.basekV ** 2 / self.baseMVA
        # Squared limits
        self.vmin = vmin ** 2
        self.vmax = vmax ** 2
        self.I_ijmax = Imax ** 2 / self.baseI ** 2

        self.SubstationNodes = [0]
        self.Buses = [self.grid.buses[i].name for i in range(len(self.grid.buses))]
        # line
        self.betas_iji = list(zip([self.grid.lines[i].bus_from.name for i in range(len(self.grid.lines))],[self.grid.lines[i].bus_to.name for i in range(len(self.grid.lines))]))+list(zip([self.grid.lines[i].bus_to.name for i in range(len(self.grid.lines))],[self.grid.lines[i].bus_from.name for i in range(len(self.grid.lines))]))
        self.branch_ij = list(zip([self.grid.lines[i].bus_from.name for i in range(len(self.grid.lines))],[self.grid.lines[i].bus_to.name for i in range(len(self.grid.lines))]))
        #impedances
        a=[self.grid.lines[i].bus_from.name for i in range(len(self.grid.lines))]
        b=[self.grid.lines[i].bus_to.name for i in range(len(self.grid.lines))]
        r=[self.grid.lines[i].R / self.baseZ for i in range(len(self.grid.lines))] 
        x=[self.grid.lines[i].X / self.baseZ for i in range(len(self.grid.lines))] 
        self.r_ij = list(zip(a,b,r))
        self.x_ij = list(zip(a,b,x))
        #powers
        a=[self.grid.loads[node].bus.name for node in range(len(self.grid.loads))]
        b=[self.grid.loads[node].P for node in range(len(self.grid.loads))]
        c=[self.grid.loads[node].Q for node in range(len(self.grid.loads))]
        self.Pload_i =dict(zip(a,b))
        self.Qload_i =dict(zip(a,c))
        # Nodes to where the node is sending the power
        self.BusIN = {node: [line.bus_to.name for line in self.grid.lines if line.bus_from.name == node] for node in self.Buses }
        # Nodes from where the node is receiving the power
        self.BusOUT = {node: [line.bus_from.name for line in self.grid.lines if line.bus_to.name == node] for node in self.Buses }
        # all nodes connected to a certain node
        self.BusALL = {node: [line.bus_to.name for line in self.grid.lines if line.bus_from.name == node]  +[line.bus_from.name for line in self.grid.lines if line.bus_to.name == node] for node in self.Buses}
        #list of branches with direction
        self.branchDirALL = [(key, val) for key, values in self.BusALL.items() for val in values]
        self.branchDirIN = [(key, val) for key, values in self.BusIN.items() for val in values]
        self.branchDirOUT = [(key, val) for key, values in self.BusOUT.items() for val in values]

        # root nodes
        self.BusWithoutSubstation = list(set(self.Buses) - set(self.SubstationNodes))

        self.model = Model('distflow')
        #self.model.Params.LogToConsole = 0


    def Solve(self,obj,OutputFlag=1):
        # Main loops
        logging.getLogger('micp.py').info("Solve from MICP_Gurobi")
        self.ContraintsDefinition()
        self.model.update()
        self.model.setObjective(obj, GRB.MINIMIZE)
        self.model.update()
        self.model.Params.OutputFlag = OutputFlag
        self.model.Params.NumericFocus = 3
        self.model.Params.DualReductions = 0
        self.model.optimize()
        return self.returnResult()

    def returnResult(self):
        logging.getLogger('micp.py').info("returnResult from MICP_Gurobi")
        DisabledLines = []
        if self.model.Status > GRB.OPTIMAL:
            return DisabledLines
        for line in self.branch_ij:
            lineij = tuple(line)
            lineji = lineij[::-1]

            status = round(self.beta_iji[lineij].x + self.beta_iji[lineji].x)
            if status == 0:
                DisabledLine = \
                    self.net.lines[((self.net.lines['Bus1'] == line[0]) & (self.net.lines['Bus2'] == line[1])) | (
                            (self.net.lines['Bus2'] == line[0]) & (self.net.lines['Bus1'] == line[1]))][
                        'id'].values[0]
                # print(line,DisabledLine)
                DisabledLines.append(DisabledLine)
        logging.getLogger('micp.py').debug(f"disabled lines:{DisabledLines}")
        # print("Disabled lines", DisabledLines)
        return DisabledLines











