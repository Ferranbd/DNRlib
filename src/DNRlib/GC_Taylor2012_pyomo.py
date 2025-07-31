import numpy as np
import random
import math
import pyomo.environ as pyo

import logging
import GridCalEngine.api as gce  # For interfacing with the GridCal API
from GridCalEngine.IO.file_handler import FileOpen, FileSave
import GC_utils 

class MICP_Pyomo ():
    def __init__(self, grid, NumTieLines, algorithm="QP", bigM=1e6, Imax=0, vmin=0, vmax=1, verbose_logging=logging.WARNING):
        #MICP parameters
        self.verbose = verbose_logging
        logging.getLogger('micp.py').setLevel(self.verbose)
        logging.getLogger('micp.py').debug("init from MICP_Pyomo")
        # MICP parameters
        self.grid = grid
        self.method = algorithm
        self.NumTieLines = NumTieLines
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

        # buses
        self.BusesAll = [self.grid.buses[i].name for i in range(len(self.grid.buses))]
        self.BusesSubstations = [self.grid.buses[i].name for i in range(len(self.grid.buses)) if self.grid.buses[i].is_slack]
        self.BusesNoSubstations = [item for item in self.BusesAll if item not in self.BusesSubstations]
        # lines

        self.beta_iji = list((self.grid.lines[i].bus_from.name, self.grid.lines[i].bus_to.name) for i in range(len(self.grid.lines)))+list((self.grid.lines[i].bus_to.name, self.grid.lines[i].bus_from.name) for i in range(len(self.grid.lines)))
        self.Lines_ij = list((self.grid.lines[i].bus_from.name, self.grid.lines[i].bus_to.name) for i in range(len(self.grid.lines)))
        
        logging.getLogger('micp.py').debug(f"self.Lines_ij: {self.Lines_ij}")
 
        self.SubstationOutputLines = [tup for tup in self.Lines_ij if tup[0] in self.BusesSubstations]
        # Nodes to where the node is sending the power
        self.BusIN = {node: [line.bus_to.name for line in self.grid.lines if line.bus_from.name == node] for node in self.BusesAll }
        # Nodes from where the node is receiving the power
        self.BusOUT = {node: [line.bus_from.name for line in self.grid.lines if line.bus_to.name == node] for node in self.BusesAll }

        self.LinesDirIN = [(key, val) for key, values in self.BusIN.items() for val in values]

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

        #optimization model
        logging.getLogger('micp.py').setLevel(logging.WARNING)

        self.model = pyo.ConcreteModel(name="DNR")
        #model lines
        self.model.beta_iji = pyo.Var(self.beta_iji, initialize=0, domain=pyo.Binary)
        #self.model.beta_iji = pyo.Var(self.beta_iji, initialize=0, bounds=(0,1))
        self.model.alpha_ij = pyo.Var(self.Lines_ij, initialize=1, domain=pyo.Binary)
        #model powers
        self.model.Psubstation_i = pyo.Var(self.BusesSubstations, initialize=0)
        self.model.Qsubstation_i = pyo.Var(self.BusesSubstations, initialize=0)
        self.model.P_ij = pyo.Var(self.Lines_ij, initialize=0)
        self.model.Q_ij = pyo.Var(self.Lines_ij, initialize=0)
        self.model.Paux_i = pyo.Var(self.BusesAll, initialize=0)
        self.model.Qaux_i = pyo.Var(self.BusesAll, initialize=0)
        #model squared voltages
        self.model.V2_i = pyo.Var(self.BusesAll, initialize=0)
        self.model.V2aux_i = pyo.Var(self.BusesAll, initialize=0)

    #constraint definition
    # Distflow equations (5-6)
    def Constraint_5_Distflow(self, model, node):
        Pin=0
        Pout=0
        for i in self.BusIN[node]:
            Pin += self.model.P_ij[node,i]
        for i in self.BusOUT[node]:
            Pout += self.model.P_ij[i,node]

        return (Pin - Pout) == self.Pload_i[node]

    def Constraint_6_Distflow(self, model, node):
        Qin=0
        Qout=0
        for i in self.BusIN[node]:
            Qin += self.model.Q_ij[node,i]
        for i in self.BusOUT[node]:
            Qout += self.model.Q_ij[i,node]

        return (Qin - Qout) == self.Qload_i[node]

    # Distflow equations (7-8)
    def Constraint_7_Distflow(self, model, node):
        Pout = 0
        for i in self.BusOUT[node]:
            Pout += self.model.P_ij[node, i]

        return Pout == self.model.Psubstation_i[node]

    def Constraint_8_Distflow(self, model, node):
        Qout = 0
        for i in self.BusOUT[node]:
            Qout += self.model.Q_ij[node, i]

        return Qout == self.model.Qsubstation_i[node]

    # Distflow equations (9-10)

    def Constraint_9a_Distflow(self, model, line_i, line_j):
        (i,j)=(line_i, line_j)
        return (self.model.P_ij[i, j] <= self.bigM * self.model.beta_iji[i, j])

    def Constraint_9b_Distflow(self, model, line_i, line_j):
        (i, j) = (line_i, line_j)
        return (self.model.P_ij[i, j] >= 0)

    def Constraint_10a_Distflow(self, model, line_i, line_j):
        (i, j) = (line_i, line_j)
        return (self.model.Q_ij[i, j] <= self.bigM * self.model.beta_iji[i, j])

    def Constraint_10b_Distflow(self, model, line_i, line_j):
        (i, j) = (line_i, line_j)
        return (self.model.Q_ij[i, j] >= 0)


    # Radiality (12)
    def Constraint_12_Radiality(self, model, line_i, line_j):
        return self.model.beta_iji[line_j, line_i] == 0
    # Radiality (14)
    def Constraint_14_Radiality(self, model, line_i, line_j):
        return self.model.beta_iji[line_i, line_j] + self.model.beta_iji[line_j, line_i] == self.model.alpha_ij[line_i, line_j]
    # Radiality (15)
    def Constraint_15_Radiality(self, model, node):
        sumBeta=0
        for j in self.BusOUT[node]:
            sumBeta += self.model.beta_iji[j, node]
        return sumBeta == 1

    #Distflow equations (21-22)
    def Constraint_21_SOC(self, model, node):
        Pin = 0
        Pout = 0
        for i in self.BusIN[node]:
            Pin += self.model.P_ij[i, node]
        for i in self.BusOUT[node]:
            Pout += self.model.P_ij[node, i]

        return self.model.Paux_i[node] == - self.Pload_i[node] + (Pin - Pout)

    def Constraint_22_SOC(self, model, node):
        Qin = 0
        Qout = 0
        for i in self.BusIN[node]:
            Qin += self.model.Q_ij[i, node]
        for i in self.BusOUT[node]:
            Qout += self.model.Q_ij[node, i]

        return self.model.Qaux_i[node] == - self.Qload_i[node] + (Qin - Qout)
    def Constraint_23_SOC(self, model, line_i, line_j):
        return self.model.V2aux_i[line_i] <= self.model.V2_i[line_j] + self.bigM * (1 - self.model.beta_iji[line_j, line_i])

    def Constraint_24_SOC(self, model, line_i, line_j):
        return self.model.V2aux_i[line_i] >= self.model.V2_i[line_j] - self.bigM * (1 - self.model.beta_iji[line_j, line_i])

    def Constraint_25_SOC(self, model, line_i, line_j):
        r = next(c for a, b, c in self.r_ij if a == line_j and b == line_i)
        return r * (self.model.P_ij[line_j,line_i]**2 + self.model.Q_ij[line_j,line_i]**2) <= self.model.V2aux_i[line_i] * self.model.Paux_i[line_i]

    def Constraint_26_SOC(self, model, line_i, line_j):
        x = next(c for a, b, c in self.x_ij if a == line_j and b == line_i)
        return x * (self.model.P_ij[line_j,line_i]**2 + self.model.Q_ij[line_j,line_i]**2) <= self.model.V2aux_i[line_i] * self.model.Qaux_i[line_i]

    def Constraint_27_SOC(self, model, line_i, line_j):
        r = next(c for a, b, c in self.r_ij if a == line_j and b == line_i)
        x = next(c for a, b, c in self.x_ij if a == line_j and b == line_i)
        return self.model.V2_i[line_i] <= self.model.V2_i[line_j] - 2 * (r * self.model.P_ij[line_j, line_i] + x * self.model.Q_ij[line_j, line_i]) + self.bigM * (1-self.model.beta_iji[line_j, line_i])
    def Constraint_28_SOC(self, model, line_i, line_j):
        r = next(c for a, b, c in self.r_ij if a == line_j and b == line_i)
        x = next(c for a, b, c in self.x_ij if a == line_j and b == line_i)
        return self.model.V2_i[line_i] >= self.model.V2_i[line_j] - 2 * (r * self.model.P_ij[line_j, line_i] + x * self.model.Q_ij[line_j, line_i]) - self.bigM * (1 - self.model.beta_iji[line_j, line_i])

    def Constraint_29_SOC(self, model, node):
        return self.model.V2aux_i[node] >= self.model.V2_i[node]

    def TotalPower(self, model):
        power = 0
        for (i, j) in self.Lines_ij:
               power += (model.P_ij[i, j] ** 2 + model.Q_ij[i, j] ** 2) * next(c for a, b, c in self.r_ij if a == i and b == j)
        return power

    def ConstraintDefinition(self):
        if self.method=="QP":
            logging.getLogger('micp.py').debug("Constraints for QP solver")
            self.model.Constraint_5_Distflow = pyo.Constraint(self.BusesNoSubstations, rule=self.Constraint_5_Distflow)
            self.model.Constraint_6_Distflow = pyo.Constraint(self.BusesNoSubstations, rule=self.Constraint_6_Distflow)

        self.model.Constraint_7_Distflow = pyo.Constraint(self.BusesSubstations, rule=self.Constraint_7_Distflow)
        self.model.Constraint_8_Distflow = pyo.Constraint(self.BusesSubstations, rule=self.Constraint_8_Distflow)
        self.model.Constraint_9a_Distflow = pyo.Constraint(self.Lines_ij, rule=self.Constraint_9a_Distflow)
        self.model.Constraint_9b_Distflow = pyo.Constraint(self.Lines_ij, rule=self.Constraint_9b_Distflow)
        self.model.Constraint_10a_Distflow = pyo.Constraint(self.Lines_ij, rule=self.Constraint_10a_Distflow)
        self.model.Constraint_10b_Distflow = pyo.Constraint(self.Lines_ij, rule=self.Constraint_10b_Distflow)
        # condition 11 is implicit in the variable definition beta_ij>=0
        self.model.Constraint_12_Radiality = pyo.Constraint(self.SubstationOutputLines, rule=self.Constraint_12_Radiality)
        # condition 13 do not apply, as it is considered that any branch can be disabled
        self.model.Constraint_14_Radiality = pyo.Constraint(self.Lines_ij, rule=self.Constraint_14_Radiality)
        self.model.Constraint_15_Radiality = pyo.Constraint(self.BusesNoSubstations, rule=self.Constraint_15_Radiality)
        ##condition 16 is implicit in the variable definition alpha_ij>=0

        if self.method=="SOC":
            logging.getLogger('micp.py').debug("Constraints for SOC solver")
            self.model.Constraint_21_SOC = pyo.Constraint(self.BusesAll, rule=self.Constraint_21_SOC)
            self.model.Constraint_22_SOC = pyo.Constraint(self.BusesAll, rule=self.Constraint_22_SOC)
            self.model.Constraint_23_SOC = pyo.Constraint(self.Lines_ij, rule=self.Constraint_23_SOC)
            self.model.Constraint_24_SOC = pyo.Constraint(self.Lines_ij, rule=self.Constraint_24_SOC)
            self.model.Constraint_25_SOC = pyo.Constraint(self.LinesDirIN, rule=self.Constraint_25_SOC)
            self.model.Constraint_26_SOC = pyo.Constraint(self.LinesDirIN, rule=self.Constraint_26_SOC)
            self.model.Constraint_27_SOC = pyo.Constraint(self.LinesDirIN, rule=self.Constraint_27_SOC)
            self.model.Constraint_28_SOC = pyo.Constraint(self.LinesDirIN, rule=self.Constraint_28_SOC)
            self.model.Constraint_29_SOC = pyo.Constraint(self.BusesSubstations, rule=self.Constraint_29_SOC)
            pass


    def Solve(self, solver="ipopt"):
        # Main loops
        logging.getLogger('micp.py').setLevel(self.verbose)
        logging.getLogger('micp.py').info(f"Start solving by MICP by {self.method}")
        self.ConstraintDefinition()
        self.model.obj = pyo.Objective(rule=self.TotalPower, sense=pyo.minimize)
        if (solver=="ipopt"):
            opt = pyo.SolverFactory('ipopt', executable='C:/Users/ferra/AMPL/ipopt.exe')
        elif (solver == "gurobi"):
            opt = pyo.SolverFactory('cbc', executable='C:/Users/ferra/AMPL/gurobi.exe') # --solver-options="acceptable_tol=0.01"')
            #opt = pyo.SolverFactory('gurobi')
            opt.options['NonConvex'] = 2
        elif (solver=="couenne"):
            opt = pyo.SolverFactory('couenne', executable='C:/Users/ferra/AMPL/couenne.exe') # --solver-options="acceptable_tol=0.01"')
        elif (solver=="cbc"):
            opt = pyo.SolverFactory('cbc', executable='C:/Users/ferra/AMPL/cbc.exe') # --solver-options="acceptable_tol=0.01"')
        elif (solver=="ipopt2"):
            opt = pyo.SolverFactory('ipopt', executable='C:/Users/ferra/AMPL/ipopt.exe') # --solver-options="acceptable_tol=0.01"')
        elif (solver=="glpk"):
            opt = pyo.SolverFactory('glpk', executable='C:/Users/ferran.bohigas/glpk/w64/glpsol.exe') # --solver-options="acceptable_tol=0.01"')
        elif (solver=="appsi_highs"):
            opt = pyo.SolverFactory('appsi_highs')
        else:
            opt = pyo.SolverFactory('C:/Users/ferra/AMPL/gurobi')

        with open('model.txt', 'w') as f:
            self.model.pprint(f)
        self.result = opt.solve(self.model, tee=False)
        if (self.result.solver.termination_condition != 'infeasible'):
            return self.result,self.returnResult()
        else:
            logging.getLogger('micp.py').debug("Infeasible solution")
            return self.result,self.returnResult()

    def returnResult(self):
        logging.getLogger('micp.py').debug("returnResult from MICP_Gurobi")
#        TieLines = []
#        for branch in self.Lines_ij:
            #a = self.model.beta_iji[branch[0], branch[1]].value
            #b = self.model.beta_iji[branch[1], branch[0]].value
#            alpha = self.model.alpha_ij[branch[0], branch[1]].value
#            if alpha < 0.75:
#                TieLines.append(
#                    self.net.lines[(self.net.lines['Bus1'] == branch[0]) & (self.net.lines['Bus2'] == branch[1])].index[
#                        0])

        branches = []
        for branch in self.Lines_ij:
            a = self.model.beta_iji[branch[0], branch[1]].value
            b = self.model.beta_iji[branch[1], branch[0]].value
            alpha = self.model.alpha_ij[branch[0], branch[1]].value
            branches.append((branch[0], branch[1], alpha))
        branches= sorted(branches, key=lambda x: x[2])[:self.NumTieLines]

        #disabled_lines = [(el[0],el[1]) for el in b]
        disconnected_lines = []
        for disabled in branches:
            for line in self.grid.lines:
                if (line.bus_from.name == disabled[0]) and (line.bus_to.name == disabled[1]):
                    disconnected_lines.append(line.name)
                if (line.bus_from.name == disabled[1]) and (line.bus_to.name == disabled[0]):
                    disconnected_lines.append(line.name)
        disconnected_lines


        return disconnected_lines


if __name__ == '__main__':
    logging.basicConfig(
        level=logging.INFO,  # Set the log level to DEBUG
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Set the log format
        datefmt='%Y-%m-%d %H:%M:%S'  # Set the date format
    )

    # Create a network object
    caseName = "33 buses case based on Baran and Wu"
    gridGC33 = FileOpen("D:\\15_Thesis-code\\02_DistributionNetworkOperationFramework\\03_NetworkExamples\\gridcal\\case33bw.m").open()
    #TieLines33Name=['line 32','line 33','line 34','line 35','line 36']
    TieLines33Name=['21_8_1', '9_15_1', '12_22_1', '18_33_1','25_29_1']

    Tielines33ID = GC_utils.GC_Line_Name2idtag_array(gridGC33, TieLines33Name)
    GC_utils.NetworkReconfiguration(gridGC33, all=True, value_all=True, selected_configuration=[], value_configuration=False)

    micp = MICP_Pyomo(gridGC33, Tielines33ID, algorithm="QP", bigM=1e8, Imax=0, vmin=0, vmax=1, verbose_logging=logging.DEBUG)

    # Solve the Minimum Spanning Tree problem
    disabled_lines = micp.Solve(solver="ipopt")
    print("disabled:",disabled_lines[1])

