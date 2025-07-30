## Distribution Network Library

##paper_results Folders

This folder contains the files used to obtain the results of the paper associated to this library:
- Results_DNR_00_Summary : executes all the methods for the selected case
    python Results_DNR_00_Summary.py -c 2 -v 3 #executes all the reconfiguration methods for the 69-buses case (c=2) with verbose level = logging.INFO (v=3)
- Results_DNR_Pareto : executes several iterations of a set of fitness ratios ([0,0.01,0.1,0.3,0.5,0.7,0.9,0.99,1]) on the selected case using the Jakus's method
    python Results_DNR_Pareto.py -v 0 -c 1 -i 10     #executes the pareto analysis, for case 33-buses case (c=1), averaging 10 iterations (i=10) and uses verbose level = logging.ERROR (v=0)
- Results_DNR_VoltageProfile : stores the voltage profile of the selected case with all the available methods
    python Results_DNR_VoltageProfile.py -c 1 -v 0      #executes the voltage profile analysis for the 118-buses case (c=3) and uses verbose leve = logging.DEGUG (v=2)
- Results_DNR_Simbench.py : stores the timeseries result for the 1-HVMV-urban-2.203-0-no_sw Simbench case, using several methods and configurations for those methods. It is possible to select the number of days which will be analyzed and how many samples are used, one of each 'n' samples
    python Results_DNR_Simbench.py -f 0 -d 1 -s 4 -v 0  #executes the time serie analysis, starting at day 0, analyzing 1 day and taking one of every 4 samples (as the file has measurements every 15min, 1 of every 4 it means a measure per hour)

Results are stored in the .\\paper_results\\data (xlsx,json,csv) and .\\paper_results\\results (pdf,png)