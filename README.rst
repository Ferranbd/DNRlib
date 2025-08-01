Distribution Network Reconfiguration Library
============================================

The main class/file in this library is DistributionNetworkReconfiguration, located in GC_DistributionNetworkReconfiguration.py, which manages the Distribution Network Reconfiguration, by importing
the initial network using GridCal simulator (accepts Matpower .m files) and calling the different methods with a single function DistributionNetworkReconfiguration.Solve()

The Distribution Network Reconfiguration library implements the following papers:

- "Network reconfiguration in distribution systems for loss reduction and load" by Baran and Wu, 1989 (GC_Baran1989.py) : https://ieeexplore-ieee-org.recursos.biblioteca.upc.edu/document/25627
- "Optimal Reconfiguration of Distribution Networks Using Hybrid Heuristic-Genetic Algorithm" by Jakus et al., 2020 (GC_Jakus2020.py) : https://www.mdpi.com/1996-1073/13/7/1544
- "Reconfiguration for loss reduction of distribution systems using selective particle swarm optimization" by Khalil and Gorpinich, 2012 (GC_Khalil_Gorpinich_2012.py) : http://www.ijmse.org/Volume3/Issue6/paper4.pdf
- "Search for a minimal loss operating spanning tree in an urban power distribution system" by Merlin and Back, 1975 (GC_Merlin1975.py) : https://pscc-central.epfl.ch/repo/papers/1975/18.pdf
- "An efficient brute-force solution to the network reconfiguration problem" by Morton and Mareels, 2000 (GC_Morton2000.py) : https://ieeexplore-ieee-org.recursos.biblioteca.upc.edu/document/871365/
- "A minimal spanning tree algorithm for distribution networks configuration" by Montoya and Ramirez, 2012 (MSTgreedy.py) : https://ieeexplore-ieee-org.recursos.biblioteca.upc.edu/document/6344718
- "An effective network reconfiguration approach of radial distribution system for loss minimization and voltage profile improvement" et al. by Salkuti et al., 2021 (GC_Salkuti2012.py)
- "Convex Models of Distribution System Reconfiguration" by Taylor and Hover, 2012 (GC_Taylor2012_pyomo.py)

There are also two auxiliary files :

- GC_PandaPowerImporter.py : with functions to import PandaPower cases to GridCal, which is necessary for the use of the Simbench timeseries
- GC_utils.py : basic functions used in more than one file

Installation
------------

The module can be installed via pip with the command::

- clone or copy the project in your computer and execute "pip install ." from the project root folder
- pip install DNRlib   (to be finished)


Tests [to be implemented]
-----
There are unit tests for each of the algorithms::

python test/test_DistributionNetworkReconfiguration.py

Dependencies [to be updated]
------------

Python packages dependencies :

- networkx-3.5
- numpy-2.3.2
- pandas-2.2.6
- matplotlib-3.10.3
- gridcalengine-5.3.53
- pandapower-3.1.2
- simbench-1.6.1
- pyomo-6.9.2
- ipykernel-6.30.0
- seaborn-0.13.2 (only to visualize the paper results)
- sphinx-8.2.3 (to regenerate the documentation)
- sphinx_rtd_theme-3.0.2 (to regenerate the documentation)
- rst-to-myst[sphinx]-0.4.0 (to regenerate the documentation)
- myst_parser-4.0.1 (to regenerate the documentation)
- sphinx_autodoc_typehints==3.2.0 (to regenerate the documentation)
- numpydoc==1.9.0 (to regenerate the documentation)


Written using Python 3.13.0

The exact requirements can be found in the `requirements.txt` file.

Folders
-------
- docs
- examples
- networks : network examples, in .m format, copied from matpower github repository
- paper_results : python scripts used to obtain the results shown in the paper
- src : sourced code in python for all the methods and classes
- test : test units (TBD)

Example
-------

.. code-block:: python
   :linenos:
   :caption: A simple code

   from GridCalEngine.IO.file_handler import FileOpen

   gridGC = FileOpen("case69.m").open()
   dnr = DistributionNetworkReconfiguration(gridGC)
   radiality = GC_utils.CheckRadialConnectedNetwork(gridGC)
   disabled_lines = dnr.Solve(method="Khalil", NumCandidates=10)

Documentation
-------------

Detailed functions documentation can be found in https://dnrlib.readthedocs.io/en/latest/
