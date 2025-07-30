import networkx as nx
import pandas as pd
import logging
import numpy as np
import random
import os
import sys
import GridCalEngine.api as gce  # For interfacing with the GridCal API
from GridCalEngine.IO.file_handler import FileOpen, FileSave
import simbench as sb
import pandapower as pp
import time
import argparse

# Add the directory containing utility_script.py to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))

import GC_utils
import GC_PandaPowerImporter
import GC_DistributionNetworkReconfiguration

if __name__ == '__main__':
    print('Algorithms to find the optimal distribution network configuration')

    parser = argparse.ArgumentParser(
        description="arguments for the DNR summary example"
    )

    parser.add_argument("--case", "-c",   type=int,   default= 0, help="network to be solved, located in the root\\networks folder, coded as 0='16 buses', 1='33 buses', 2='69 buses', 3='118 buses', 4='1-HVMV-urban-2.203-0-no_sw by Simbench'"    )
    parser.add_argument("--verbose", "-v", type=int,  default= 0, help="verbose level : 0 (Default)=logging.ERROR, 1=logging.WARNING, 2=logging.DEBUG, 3=logging.INFO"    )

    # 3. Parse the arguments
    args = parser.parse_args()

    #parse logging level
    if args.verbose==0:
        verbose_logging = logging.ERROR
    elif args.verbose==1:
        verbose_logging = logging.WARNING
    elif args.verbose==2:
        verbose_logging = logging.DEBUG
    elif args.verbose==3:
        verbose_logging = logging.INFO
    else:
        verbose_logging=logging.ERROR

    logging.basicConfig(
        level=verbose_logging,  # Set the log level to DEBUG
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',  # Set the log format
        datefmt='%Y-%m-%d %H:%M:%S'  # Set the date format
    )

    #parse case and create the gridcal network object
    if args.case==0:
        caseName = "16 buses"
        gridGC = FileOpen("..\\networks\\case16ci.m").open()
        TieLinesName=['1_4_1', '1_8_1', '4_5_1']
    elif args.case==1:
        caseName = "33 buses"
        gridGC = FileOpen("..\\networks\\case33bw.m").open()
        TieLinesName=['21_8_1', '9_15_1', '12_22_1', '18_33_1','25_29_1']
    elif args.case==2:
        caseName = "69 buses"
        gridGC = FileOpen("..\\networks\\case69.gridcal").open()
        TieLinesName = ['line 57','line 10','line 69','line 14','line 19']
    elif args.case==3:
        caseName = "118 buses"
        gridGC = FileOpen("..\\networks\\case118.m").open()
        TieLinesName = ['1_2_1', '4_5_1', '3_12_1', '7_12_1', '13_15_1', '14_15_1', '12_16_1', '18_19_1', '19_20_1', '15_19_1', '29_31_1', '27_32_1', '15_33_1', '35_36_1', '39_40_1', '40_41_1', '40_42_1', '34_43_1', '46_47_1', '46_48_1', '47_49_1', '53_54_1', '54_55_1', '54_56_1', '55_56_1', '56_57_1', '56_58_1', '60_61_1', '60_62_1', '61_62_1', '63_64_1', '62_67_1', '65_68_1', '24_70_1', '24_72_1', '70_74_1', '70_75_1', '75_77_1', '78_79_1', '68_81_1', '84_85_1', '85_88_1', '90_91_1', '93_94_1', '82_96_1', '94_96_1', '94_100_1', '95_96_1', '96_97_1', '98_100_1', '99_100_1', '100_101_1', '103_104_1', '104_105_1', '105_106_1', '106_107_1', '108_109_1', '17_113_1', '114_115_1', '76_118_1', '30_17_1', '64_61_1']
    elif args.case==4:
        caseName = "1-HVMV-urban-2.203-0-no_sw"
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

    #init the system
    results = dict()
    fitness_ratio = 1
    TieLinesID=GC_utils.GC_Line_Name2idtag_array(gridGC, TieLinesName)  #convert the name of the line to the GridCal tagID
    dnr = GC_DistributionNetworkReconfiguration.DistributionNetworkReconfiguration(gridGC, verbose_logging=verbose_logging)     #creates the DistributionNetworkReconfiguration object

    try:
        #calculates the losses, radiality and voltages for the initial case, using the case tielines
        print("starting DNR by original configuration", caseName)
        dnr.NumPF=0     #number of Power Flow executions
        start = time.time() #init the calculation time
        disabled_lines = TieLinesID
        PFresult, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        Vmax = PFresult.results.get_bus_df().Vm.max()
        Vmin = PFresult.results.get_bus_df().Vm.min()
        runtime = time.time() - start
        print(f"TieLines: The new optimal configuration losses:{loss:.3f}, radiality:{radiality}, fitness_ratio=1, numPF:{dnr.NumPF}, vmax:{Vmax:.3f}, vmin:{Vmin:.3f}, time={runtime:.3f} ") 
        #results are stored in dictionary
        results['original '+caseName] = [caseName, 'TieLines',['TieLines'],loss,Vmax,Vmin,runtime, fitness_ratio, radiality, dnr.NumPF, TieLinesName]
    except:
        print("Original configuration fails")   

    try:
        #calculates the result by using the Merlin 1975 paper implementation
        print("starting DNR by Merlin", caseName)
        dnr.NumPF=0
        start = time.time()
        disabled_lines = dnr.Solve(method="Merlin", TieLines=TieLinesID)
        PFresult, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        Vmax = PFresult.results.get_bus_df().Vm.max()
        Vmin = PFresult.results.get_bus_df().Vm.min()
        runtime = time.time() - start
        print(f"Merlin: The new optimal configuration losses:{loss:.3f}, radiality:{radiality}, numPF:{dnr.NumPF}, vmax:{Vmax:.3f}, vmin:{Vmin:.3f}, time={runtime:.3f}   ") 
        results['Merlin'+caseName] = [caseName, 'Merlin',['TieLines'],loss,Vmax,Vmin,runtime, fitness_ratio, radiality, dnr.NumPF, disabled_lines]
    except:
        print("Merlin method has failed")        

    try:
        #calculates the result by using the Baran 1989 paper implementation
        print("Baran", caseName)
        dnr.NumPF=0
        start = time.time()
        disabled_lines = dnr.Solve(method="Baran", TieLines=TieLinesID)
        PFresult, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        Vmax = PFresult.results.get_bus_df().Vm.max()
        Vmin = PFresult.results.get_bus_df().Vm.min()
        runtime = time.time() - start
        print(f"Baran: The new optimal configuration losses:{loss:.3f}, radiality:{radiality}, numPF:{dnr.NumPF}, vmax:{Vmax:.3f}, vmin:{Vmin:.3f}, time={runtime:.3f}   ") 
        results['Baran'+caseName] = [caseName, 'Baran', ['TieLines'],loss,Vmax,Vmin,runtime, fitness_ratio, radiality, dnr.NumPF, disabled_lines]
    except:
        print("Baran method has failed")        

    try:
        #calculates the result by using the Salkuti 2021 paper implementation
        print("starting DNR by Salkuti", caseName)
        dnr.NumPF=0
        start = time.time()
        disabled_lines = dnr.Solve(method="Salkuti", TieLines=TieLinesID)
        PFresult, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        Vmax = PFresult.results.get_bus_df().Vm.max()
        Vmin = PFresult.results.get_bus_df().Vm.min()
        runtime = time.time() - start
        print(f"Salkuti: The new optimal configuration losses:{loss:.3f}, radiality:{radiality}, numPF:{dnr.NumPF}, vmax:{Vmax:.3f}, vmin:{Vmin:.3f}, time={runtime:.3f}   ") 
        results['Salkuti'+caseName] = [caseName, 'Salkuti',['TieLines'],loss,Vmax,Vmin,runtime, fitness_ratio, radiality, dnr.NumPF, disabled_lines]
    except:
        print("Salkuti method has failed")        

    try:
        #calculates the result by using the Montoya 2012 (MST greedy) paper implementation
        print("starting DNR by Montoya/MSTGreedy", caseName)
        dnr.NumPF=0
        start = time.time()
        disabled_lines = dnr.Solve(method="MSTgreedy", randomMST=False, one=False, current_power=True, algorithm="prim")
        PFresult, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        Vmax = PFresult.results.get_bus_df().Vm.max()
        Vmin = PFresult.results.get_bus_df().Vm.min()
        runtime = time.time() - start
        print(f"Montoya/MSTgreedy: The new optimal configuration losses:{loss:.3f}, radiality:{radiality}, numPF:{dnr.NumPF}, vmax:{Vmax:.3f}, vmin:{Vmin:.3f}, time={runtime:.3f}   ") 
        results['Montoya/MSTgreedy'+caseName] = [caseName, 'MSTgreedy', [False,False,True,'prim'],loss,Vmax,Vmin,runtime, fitness_ratio, radiality, dnr.NumPF, disabled_lines]
    except:
        print("MSTGreedy method has failed")        
        
    try:
        #calculates the result by using the Khalil 2012 paper implementation
        print("starting DNR by Khalil", caseName)
        if caseName=="118 buses":
            print("Khalil method do not obtain a solution for the 118-bus case, as it is not able to find the initial set of candidates, due to the loop selection")
        else:
            dnr.NumPF=0
            start = time.time()
            disabled_lines = dnr.Solve(method="Khalil", NumCandidates=10, fitness_ratio=1, loss_factor=0.02)
            PFresult, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
            radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
            Vmax = PFresult.results.get_bus_df().Vm.max()
            Vmin = PFresult.results.get_bus_df().Vm.min()
            runtime = time.time() - start
            print(f"Khalil: The new optimal configuration losses:{loss:.3f}, radiality:{radiality}, numPF:{dnr.NumPF}, vmax:{Vmax:.3f}, vmin:{Vmin:.3f}, time={runtime:.3f}   ")
            results['Khalil'+caseName] = [caseName, 'Khalil', [10,1,0.2],loss,Vmax,Vmin,runtime, fitness_ratio, radiality, dnr.NumPF, disabled_lines]
    except:
        print("Khalil method has failed")

    try:
        #calculates the result by using the Jakus 2020 paper implementation
        print("starting DNR by Jakus", caseName, fitness_ratio)
        dnr.NumPF=0
        start = time.time()
        disabled_lines = dnr.Solve(method="Jakus", MutationProbability=0.4, PopulationSize=16, Niter=20, ElitePopulation=2, TieLines=TieLinesID, fitness_ratio=fitness_ratio, loss_factor=0.02)
        PFresult, loss = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines)
        radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
        Vmax = PFresult.results.get_bus_df().Vm.max()
        Vmin = PFresult.results.get_bus_df().Vm.min()
        runtime = time.time() - start
        print(f"Jakus: The new optimal configuration losses:{loss:.3f}, radiality:{radiality}, numPF:{dnr.NumPF}, vmax:{Vmax:.3f}, vmin:{Vmin:.3f}, time={runtime:.3f}   ") 
        results['Jakus'+caseName+str(fitness_ratio)] = [caseName,'Jakus' ,[0.4,16,20,2,'TieLines',1,0.02],loss,Vmax,Vmin,runtime, fitness_ratio, radiality, dnr.NumPF, disabled_lines]
    except:
        print("Jakus"+caseName+str(fitness_ratio)+" method has failed")
        

    #the results stored in the dictionary are saved into an Excel file
    pd.DataFrame.from_dict(results, orient='index', columns=['case', 'Method','parameters', 'Loss','vmax','vmin','time','fitness_ratio','Radiality', 'NumPF', 'DisabledLines']).to_excel('.\\data\\results.xlsx')