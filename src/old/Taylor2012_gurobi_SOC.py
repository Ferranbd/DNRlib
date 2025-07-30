import numpy as np
import random
import math
from gurobipy import *
from Taylor2012_gurobi import MICP_Gurobi

import logging

class MICP_SOC(MICP_Gurobi):
    def __init__(self, net, bigM=1e6, vmin=0, vmax=1, verbose_logging=logging.WARNING):
        logging.getLogger('micp.py').info("Init from MICP_Gurobi_QP")
        super().__init__(net, bigM=bigM, vmin=vmin, vmax=vmax, verbose_logging=verbose_logging)
        self.Psubstation_i = self.model.addVars(self.SubstationNodes, lb=-GRB.INFINITY, name='Psubstation_i')
        self.Qsubstation_i = self.model.addVars(self.SubstationNodes, lb=-GRB.INFINITY, name='Qsubstation_i')
        self.P_ij = self.model.addVars(self.branch_ij, lb=-GRB.INFINITY, name='P_ij')  # Branch P
        self.Q_ij = self.model.addVars(self.branch_ij, lb=-GRB.INFINITY, name='Q_ij')  # Branch Q
        self.beta_iji = self.model.addVars(self.betas_iji, lb=0, ub=1, name='beta_iji')  # branch status
        #self.beta_iji = self.model.addVars(self.betas_iji, vtype=GRB.BINARY, name='beta_iji')  # branch status
        self.alpha_ij = self.model.addVars(self.branch_ij, vtype=GRB.BINARY, name='alpha')  # branch flow
        self.v_i = self.model.addVars(self.Buses, lb=self.vmin, ub=self.vmax, name='v_i')
        self.v2aux_i = self.model.addVars(self.Buses, lb=self.vmin, ub=self.vmax, name='v2aux_i')
        self.paux_i = self.model.addVars(self.Buses, lb=-GRB.INFINITY,  name='paux_i')
        self.qaux_i = self.model.addVars(self.Buses, lb=-GRB.INFINITY,  name='qaux_i')

    def Solve(self,OutputFlag=1):
        logging.getLogger('micp.py').info("Solve from MICP_Gurobi_QP")
        obj = quicksum((self.P_ij[i, j] ** 2 + self.Q_ij[i, j] ** 2) * next(c for a, b, c in self.r_ij if a == i and b == j)
                for (i, j) in self.branch_ij)
        return super().Solve(obj,OutputFlag)

    def ContraintsDefinition(self):
        logging.getLogger('micp.py').info("ContraintsDefinition from MICP_Gurobi_QP")
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


        #Distflow equations (21-22)
        self.model.addConstrs(
            (self.paux_i[i] == (quicksum(self.P_ij[j,i] for j in self.BusIN[i]) - quicksum(self.P_ij[i,j] for j in self.BusOUT[i]) - self.Pload_i[i]) for i in self.BusWithoutSubstation), name='DistFlow-21')
        self.model.addConstrs(
            (self.qaux_i[i] == (quicksum(self.Q_ij[j,i] for j in self.BusIN[i]) - quicksum(self.Q_ij[i,j] for j in self.BusOUT[i]) - self.Qload_i[i]) for i in self.BusWithoutSubstation), name='DistFlow-22')

        #Distflow equations (23-24)
        self.model.addConstrs( (  self.v2aux_i[i] <= self.v_i[j] + self.bigM * (1 - self.beta_iji[j,i])    for (i, j) in self.branchDirIN), name='DistFlow23')
        self.model.addConstrs( (  self.v2aux_i[i] >= self.v_i[j] - self.bigM * (1 - self.beta_iji[j,i])    for (i, j) in self.branchDirIN), name='DistFlow24')

        #Distflow equations (25-26)
        self.model.addConstrs(
            (  (next(c for a, b, c in self.r_iji if a == i and b == j) * (self.P_ij[j,i]**2 + self.Q_ij[j,i]**2) ) <= self.v2aux_i[i] * self.paux_i[i]
                                for (i, j) in self.branchDirIN), name='DistFlow25')
        self.model.addConstrs(
            (  (next(c for a, b, c in self.x_iji if a == i and b == j)  * (self.P_ij[j,i]**2 + self.Q_ij[j,i]**2) ) <= self.v2aux_i[i] * self.qaux_i[i]
                                for (i, j) in self.branchDirIN), name='DistFlow26')

        #Distflow equations (27-28)
        self.model.addConstrs(
            (self.v_i[i] <=     self.v_i[j]
                                - 2 * (next(c for a, b, c in self.r_iji if a == i and b == j)  * self.P_ij[j, i] + next(c for a, b, c in self.x_iji if a == i and b == j)  * self.Q_ij[j, i])
                                + self.bigM * (1 - self.beta_iji[j,i])
                                for (i, j) in self.branchDirIN), name='DistFlow27')
        self.model.addConstrs(
            (self.v_i[i] >=     self.v_i[j]
                                - 2 * (next(c for a, b, c in self.r_iji if a == i and b == j)  * self.P_ij[j, i] + next(c for a, b, c in self.x_iji if a == i and b == j)  * self.Q_ij[j, i])
                                - self.bigM * (1 - self.beta_iji[j,i])
                                for (i, j) in self.branchDirIN), name='DistFlow28')

        #Voltage in the substation connection is 1pu (29)
        self.model.addConstrs((self.v_i[i] == 1 for i in self.SubstationNodes),
                         name='DistFlow29')

