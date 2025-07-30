import networkx as nx
import pandas as pd
import logging
import numpy as np
import random
import time
import os
import sys
from datetime import date
import matplotlib.pyplot as plt
import GridCalEngine.api as gce  # For interfacing with the GridCal API
from GridCalEngine.IO.file_handler import FileOpen, FileSave
import simbench as sb
import pandapower as pp
import json
# Add the directory containing utility_script.py to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'src')))
import argparse
import GC_utils
import GC_PandaPowerImporter
from GC_DistributionNetworkReconfiguration import *

import warnings
warnings.filterwarnings("ignore")

def save_dict_progressively(data, filename):
    with open(filename, 'w') as f:
        json.dump(data, f)

def load_dict(filename):
    with open(filename, 'r') as f:
        return json.load(f)

# define a function to apply absolute values
def apply_absolute_values(gridPP, absolute_values_dict, case_or_time_step):
    for elm_param in absolute_values_dict.keys():
        if absolute_values_dict[elm_param].shape[1]:
            elm = elm_param[0]
            param = elm_param[1]
            gridPP[elm].loc[:, param] = absolute_values_dict[elm_param].loc[case_or_time_step]

if __name__ == '__main__':


    parser = argparse.ArgumentParser(
        description="arguments for timeseries example with Simbench 1-HVMV-urban-2.203-0-no_sw"
    )

    parser.add_argument("--firstday", "-f",     type=int,   default= 90, help="Staring day, usually 0, 90,180,270"    )
    parser.add_argument("--days", "-d",         type=int,   default= 7, help="Number of days in the timeseries"    )
    parser.add_argument("--sampling", "-s",     type=int,   default= 4, help="It is use one of every n samples (i.e if sample_coef=4 uses one of each 4 samples = 1 sample/hour)"    )
    parser.add_argument("--verbose", "-v",      type=int,   default= 0, help="verbose level : 0 (Default)=logging.ERROR, 1=logging.WARNING, 2=logging.DEBUG, 3=logging.INFO"    )

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

    first_day = args.firstday
    num_days = args.days
    sample_coef = args.sampling  


    print('Algorithms to find the optimal distribution network configuration')


    #sb_code1 = "1-MVLV-rural-1.108-0-no_sw"
    #sb_code1 = "1-MV-rural--1-sw"
    #sb_code1 = "1-MVLV-rural-all-1-no_sw"
    sb_code1 = "1-HVMV-urban-2.203-0-no_sw"

    print(f"loading the example network: {sb_code1}")

    gridPP = sb.get_simbench_net(sb_code1)
    #delete duplicate lines and trafos, and switches, to avoid issues with GC conversion
    gridPP.switch.drop([232,234,236,238,240, 242,244,246], inplace=True)
    gridPP.trafo.drop([1,3,4], inplace=True)
    gridPP.line.drop(set([123,226,139,140,151,161,166,170,173,178,180,186,187,188,208,223,225,123,226,227,232,228,229,230,231,227,232,233]), inplace=True)
    gridPP.ext_grid.at[0,'name']="grid_ext"
    gridPP.line['in_service'] = True
    # set trafo tap position so that no voltage limits are violated
    #gridPP.trafo.tap_pos = 1
    #run pandapower powerflow
    pp.runpp(gridPP)
    #import the network to gridcal
    gridGC = GC_PandaPowerImporter.PP2GC(gridPP)

    #for the "1-HVMV-urban-2.203-0-no_sw" the following lines are choosen as tie_lines to achieve a radial fully connected network
    TieLinesName=['1_2_1', '1_24_1', '1_36_1', '1_47_1', '51_52_1', '1_60_1', '1_74_1', '1_85_1', '117_181_1', '171_117_1', '117_125_1', '127_164_1', '121_188_1', '146_147_1', '171_181_1', '116_196_1', '116_154_1']
    TieLinesName=['196_165_1', '121_188_1', '127_164_1', '146_147_1', '147_179_1', '147_151_1', '171_181_1', '171_152_1', '116_154_1', '108_109_1', '1_24_1', '1_36_1', '1_85_1', '9_10_1', '51_52_1', '112_113_1', '71_72_1']# calculate absolute profiles
    profiles = sb.get_absolute_values(gridPP, profiles_instead_of_study_cases=True)
    # check that all needed profiles existent
    assert not sb.profiles_are_missing(gridPP)


    # Define the start date (January 1st of the year)
    start_date = '2024-01-01 00:00:00'
    # Create a date range with 15-minute intervals for a year
    timestamps = pd.date_range(start=start_date, periods=96 * 366, freq='15T')
    total_n_steps = len(profiles[('load', 'q_mvar')])
    hours = total_n_steps / 4
    days = hours / 24
    samples_per_day = 4*24
    print(f"A total of {total_n_steps} samples which means {hours} hours and {days} days, from {timestamps[0]} to {timestamps[-1]}")

    TimeRange = int(num_days * samples_per_day/sample_coef)
    initial_sample = samples_per_day*first_day
    print(f"Starting at sample : {initial_sample}, the analysis will range a total of {TimeRange} samples, which corresponds to {num_days} days, taking 1 of each {sample_coef} samples")

    # time range calculation
    time_steps = range(initial_sample,initial_sample+TimeRange+1)  #96*10)
    print(time_steps,time_steps[0], time_steps[-1], time_steps[-1] - time_steps[0])

    previousTieLinesBaranName = TieLinesName
    previousTieLinesMerlinName = TieLinesName
    previousTieLinesJakusName = TieLinesName
    previousTieLinesSalkutiName = TieLinesName

    # run the time series and store results into a DataFrame
    # two options
    # - the method uses the initial TieLines as TieLines
    # - the method uses the solution from the previuos iteration as TieLines
    results = {}
    for time_step in time_steps:
        #charges the new profile and converts into GC grid
        apply_absolute_values(gridPP, profiles, time_step)
        gridGC = GC_PandaPowerImporter.PP2GC(gridPP)
        TieLinesID=GC_utils.GC_Line_Name2idtag_array(gridGC, TieLinesName)
        _, loss_init = GC_utils.GC_PowerFlow(gridGC, config=TieLinesID)
        print(f"    initial powerflow: {loss_init}")
        dnr = DistributionNetworkReconfiguration(gridGC, verbose_logging=logging.ERROR)
        #as the GC grid changes, the idtags change as well, so we need to recalculate the TieLinesID
        previousTieLinesBaranID = GC_utils.GC_Line_Name2idtag_array(gridGC, previousTieLinesBaranName)
        previousTieLinesMerlinID = GC_utils.GC_Line_Name2idtag_array(gridGC, previousTieLinesMerlinName)
        previousTieLinesJakusID = GC_utils.GC_Line_Name2idtag_array(gridGC, previousTieLinesJakusName)
        previousTieLinesSalkutiID = GC_utils.GC_Line_Name2idtag_array(gridGC, previousTieLinesSalkutiName)
        #Without optimization
        _, loss_init = GC_utils.GC_PowerFlow(gridGC, config=TieLinesID)
        radiality_init = GC_utils.CheckRadialConnectedNetwork(gridGC)
        ##With MST method 
        disabled_lines_MST_ID = dnr.Solve(method="MSTgreedy", randomMST=False, one=False, current_power=True, algorithm="prim")
        _, loss_MST = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines_MST_ID)
        radiality_MST = GC_utils.CheckRadialConnectedNetwork(gridGC)
        ##With Khalil method 
        disabled_lines_PSO_ID = dnr.Solve(method="Khalil", NumCandidates=10, fitness_ratio=1, loss_factor=0.08)
        _, loss_PSO = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines_PSO_ID)
        radiality_PSO = GC_utils.CheckRadialConnectedNetwork(gridGC)
        ##With Jakus method and TieLinesID
        disabled_lines_Jakus_TL_ID = dnr.Solve(method="Jakus", MutationProbability=0.4, PopulationSize=16, Niter=20, ElitePopulation=4, TieLines=TieLinesID, fitness_ratio=0.5, loss_factor=0.08)
        _, loss_Jakus_TL = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines_Jakus_TL_ID)
        radiality_Jakus_TL = GC_utils.CheckRadialConnectedNetwork(gridGC)
        ##with Jakus method and previous solution as TieLines
        disabled_lines_Jakus_OW_ID = dnr.Solve(method="Jakus", MutationProbability=0.4, PopulationSize=16, Niter=20, ElitePopulation=4, TieLines=previousTieLinesJakusID, fitness_ratio=0.2, loss_factor=0.08)
        _, loss_Jakus_OW = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines_Jakus_OW_ID)
        radiality_Jakus_OW = GC_utils.CheckRadialConnectedNetwork(gridGC)
        previousTieLinesJakusName = GC_utils.GC_Line_idtag2name_array(dnr.grid,disabled_lines_Jakus_OW_ID)
        ##With Salkuti method and TieLinesID
        disabled_lines_Salkuti_TL_ID = dnr.Solve(method="Salkuti", TieLines=TieLinesID)
        _, loss_Salkuti_TL = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines_Salkuti_TL_ID)
        radiality_Salkuti_TL = GC_utils.CheckRadialConnectedNetwork(gridGC)
        ##With Salkuti method and previous solution as TieLines
        disabled_lines_Salkuti_OW_ID = dnr.Solve(method="Salkuti", TieLines=previousTieLinesSalkutiID)
        _, loss_Salkuti_OW = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines_Salkuti_OW_ID)
        radiality_Salkuti_OW = GC_utils.CheckRadialConnectedNetwork(gridGC)
        previousTieLinesSalkutiName = GC_utils.GC_Line_idtag2name_array(dnr.grid,disabled_lines_Salkuti_OW_ID)
        #With Baran method and TieLinesID
        disabled_lines_Baran_TL_ID = dnr.Solve(method="Baran", TieLines=TieLinesID)
        _, loss_Baran_TL = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines_Baran_TL_ID)
        radiality_Baran_TL = GC_utils.CheckRadialConnectedNetwork(gridGC)
        #With Baran method and previous solution as TieLines
        disabled_lines_Baran_OW_ID = dnr.Solve(method="Baran", TieLines=previousTieLinesBaranID)
        _, loss_Baran_OW = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines_Baran_OW_ID)
        radiality_Baran_OW = GC_utils.CheckRadialConnectedNetwork(gridGC)
        previousTieLinesBaranName = GC_utils.GC_Line_idtag2name_array(dnr.grid,disabled_lines_Baran_OW_ID)
        ##With Merlin method and TieLinesID
        disabled_lines_Merlin_TL_ID = dnr.Solve(method="Merlin", TieLines=TieLinesID)
        _, loss_Merlin_TL = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines_Merlin_TL_ID)
        radiality_Merlin_TL = GC_utils.CheckRadialConnectedNetwork(gridGC)
        ##With Merlin method and previous solution as TieLines
        disabled_lines_Merlin_OW_ID = dnr.Solve(method="Merlin", TieLines=TieLinesID)
        _, loss_Merlin_OW = GC_utils.GC_PowerFlow(gridGC, config=disabled_lines_Merlin_OW_ID)
        radiality_Merlin_OW = GC_utils.CheckRadialConnectedNetwork(gridGC)
        previousTieLinesMerlinName = GC_utils.GC_Line_idtag2name_array(dnr.grid,disabled_lines_Merlin_OW_ID)
        
        results[time_step] = {
                        'load':        np.array([gridGC.loads[i].P for i in range(len(gridGC.loads))]).sum(),
                        'init':        (loss_init,radiality_init),
                        'PSO':         (loss_PSO,radiality_PSO),
                        'MST':         (loss_MST,radiality_MST),
                        'Jakus_TL':    (loss_Jakus_TL,radiality_Jakus_TL),
                        'Jakus_OW':    (loss_Jakus_OW,radiality_Jakus_OW),
                        'Salkuti_TL':  (loss_Salkuti_TL,radiality_Salkuti_TL),
                        'Salkuti_OW':  (loss_Salkuti_OW,radiality_Salkuti_OW),
                        'Baran_TL':    (loss_Baran_TL,radiality_Baran_TL),   
                        'Baran_OW':    (loss_Baran_OW,radiality_Baran_OW),   
                        'Merlin_TL':   (loss_Merlin_TL,radiality_Merlin_TL),
                        'Merlin_OW':   (loss_Merlin_OW,radiality_Merlin_OW),

                        'Jakus_TL_Solultion':    GC_utils.GC_Line_idtag2name_array(dnr.grid,disabled_lines_Jakus_TL_ID),
                        'Jakus_OW_Solultion':    previousTieLinesJakusName,
                        'Salkuti_TL_Solultion':  GC_utils.GC_Line_idtag2name_array(dnr.grid,disabled_lines_Salkuti_TL_ID),
                        'Salkuti_OW_Solultion':  previousTieLinesSalkutiName,
                        'Baran_TL_Solultion':    GC_utils.GC_Line_idtag2name_array(dnr.grid,disabled_lines_Baran_TL_ID),
                        'Baran_OW_Solultion':    previousTieLinesBaranName,
                        'Merlin_TL_Solultion':   GC_utils.GC_Line_idtag2name_array(dnr.grid,disabled_lines_Merlin_TL_ID),
                        'Merlin_OW_Solultion':   previousTieLinesMerlinName,
                        }       
        
        save_dict_progressively(results, '.\\data\\results_Simbech_timestep.json')

        print(f"    timestep: {time_step}/{time_steps[-1]}  ")    








