# Installation

We assume the code is stored in some root folder:

```
/path/to/code/
```

## Dependencies

**Required**
* Python 3, Version 3.11.2
* SageMath version 9.5

**Optional**
* tqmd, Version: 4.65
* adjustText, Version: 1.1.1 (https://github.com/Phlya/adjustText)


```
pip3 install adjustText
pip3 install tqdm 
```

## Build alfk 
```
cd /path/to/code/LatticeHeuristics/alfk
python3 build_alfk.py
```


# Cost Estimation 
```commandline
usage: main_interface.py [-h] [-parameter-set PARAMETER_SET] [-y Y] [-z Z] [-stepsize-z STEPSIZE_Z]
                         [-circuits CIRCUITS] [-max-depth MAX_DEPTH] [-bound BOUND] [-cache]
                         [-m-lower M_LOWER] [-m M] [-df DF] [-qd QD] [-wq WQ]
                         [-nodesbranching NODESBRANCHING] [-forcedf FORCEDF]

options:
  -h, --help            show this help message and exit
  -parameter-set PARAMETER_SET
                        One of the list [kyber512, kyber768, kyber1024] (default: None)
  -y Y                  Upper bound for QRACM/ log-Bound on number of combined nodes for virtual tree.
                        (default: 64)
  -z Z                  Upper bound for log-Jensen's z. (default: 64)
  -stepsize-z STEPSIZE_Z
                        Step size for Jensen's z. (default: 1)
  -circuits CIRCUITS    One of the list [Query, Minimal] (default: None)
  -max-depth MAX_DEPTH  One of the list [40, 64, 96, 9999, -1], where 9999 is unbounded, -1 is trivial
                        attack (default: None)
  -bound BOUND          Uses upper/lower bound for subtree sizes, [LBUB, UBUB, LBLB] (default: LBUB)
  -cache                Cache intermediate computation. (default: False)
  -m-lower M_LOWER      Dev: Log of the minimal number of bases considered for extreme pruning. Pruning
                        Radi only support [64]. (default: 64)
  -m M                  Dev: Upper bound for log-Bound on number of combined enumeration trees. Pruning
                        Radi only support [64]. (default: 64)
  -df DF                Bound constant in number of calls to DetectMV in FindMV. (default: 1)
  -qd QD                Bound constant in number of calls to QPE in DetectMV. (default: 1)
  -wq WQ                Bound constant in number of calls to operator W in QPE. (default: 1)
  -nodesbranching NODESBRANCHING
                        Bound on number of nodes checked on each level of the DFS. (default: 1)
  -forcedf FORCEDF      Force a specific value for DF (overwrites -df). (default: 1)

```

# Generating cost estimation results. 
```
cd /path/to/code/Cost_Estimation

python3 main_interface.py -bound UBUB -cache
python3 main_interface.py -bound LBUB -cache
python3 main_interface.py -bound LBLB -cache

cd /path/to/code/Cost_Estimation/costResults

python3 criticalToLaTeX_Bounds.py -all

cd /path/to/code/Cost_Estimation/costResults/TeX/

pdflatex Tables.tex

```
The default values correspond to the lower bounds for the constant factors. 
in ```main_interface.py``` you can uncomment line 241 to get the respective "beyond lower bounds" from the paper, 
or provide them via the command line parameters. 

This computes all tables and plots from the paper (https://eprint.iacr.org/2023/1423) of Section 5, Appendix G, H, saves intermediate results to `costData` and final results are written to `costResults`.
Additional information can be found via
```
cd /path/to/code/Cost_Estimation

python3 main_interface.py --help 
```


## Results

From the paper:
* Pattern for filenames is: costResults/CIRCUIT/BLOCKSIZE/CONSTANTS_CIRCUIT_BOUND_LESSSIMPLE
    * CIRCUIT = {Query, Minimal}
    * BLOCKSIZE = {406, 623, 873}
    * BOUND = {LBUB, LBLB, UBUB}
* `/path/to/code/Cost_Estimation/costResults/Query` contains estimations with instantiation of quantum operator W as in Section 4.1
* `/path/to/code/Cost_Estimation/costResults/Minimal` contains estimations with instantiation of quantum operator W as in Section 4.2
* The graphs the manuscript can be found in the folder with the corresponding blocksize:
    * for Kyber512 in folder `406`
    * for Kyber768 in folder `623`
    * for Kyber1024 in folder `873`
* Tables with crossover-points vs cost of Grover on AES are printed to console, and can also be found in `*.critical` files 




### Section 5 
* Table 5, printed to console with data in
    * `/path/to/code/Cost_Estimation/costResults/Query/DF=QD=WQ=1_Query_LBUB_LESSIMPLE.critical`
    * `/path/to/code/Cost_Estimation/costResults/Query/DF=QD=WQ=1_Query_UBUB_LESSIMPLE.critical`
    * `/path/to/code/Cost_Estimation/costResults/Query/DF=QD=WQ=1_Query_LBLB_LESSIMPLE.critical`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/DF=QD=WQ=1_Minimal_LBUB_LESSIMPLE.critical`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/DF=QD=WQ=1_Minimal_UBUB_LESSIMPLE.critical`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/DF=QD=WQ=1_Minimal_LBLB_LESSIMPLE.critical`
* Figure 5 in
    * `/path/to/code/costResults/Query/873/DF=QD=WQ=1_Query_40_873_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/costResults/Query/873/DF=QD=WQ=1_Query_96_873_LBUB_LESSIMPLE.pdf`


### Appendix G.2 
* Table 9 and 10, printed to console with data in 
    * `/path/to/code/Cost_Estimation/costResults/Query/DF=QD=WQ=1_Query_LBUB_LESSIMPLE.critical`
    * `/path/to/code/Cost_Estimation/costResults/Query/DF=QD=WQ=1_Query_UBUB_LESSIMPLE.critical`
    * `/path/to/code/Cost_Estimation/costResults/Query/DF=QD=WQ=1_Query_LBLB_LESSIMPLE.critical`
* Figure 21 in 
    * `/path/to/code/Cost_Estimation/costResults/Query/873/DF=QD=WQ=1_Query_64_873_LBUB_LESSIMPLE.pdf`
* Figure 22 in 
    * `/path/to/code/Cost_Estimation/costResults/Query/873/DF=QD=WQ=1_Query_9999_873_LBUB_LESSIMPLE.pdf`
* Figure 23 in 
    * `/path/to/code/Cost_Estimation/costResults/Query/406/DF=QD=WQ=1_Query_40_406_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Query/406/DF=QD=WQ=1_Query_64_406_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Query/406/DF=QD=WQ=1_Query_96_406_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Query/623/DF=QD=WQ=1_Query_40_623_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Query/623/DF=QD=WQ=1_Query_64_623_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Query/623/DF=QD=WQ=1_Query_96_623_LBUB_LESSIMPLE.pdf`

* Figure 24 in 
    * `/path/to/code/Cost_Estimation/costResults/Minimal/406/DF=QD=WQ=1_Query_9999_406_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/623/DF=QD=WQ=1_Query_9999_623_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/873/DF=QD=WQ=1_Query_9999_873_LBUB_LESSIMPLE.pdf`

* Figure 25 in 
    * `/path/to/code/Cost_Estimation/costResults/Minimal/406/DF=QD=WQ=1_Query_9999_406_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/623/DF=QD=WQ=1_Query_9999_623_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/873/DF=QD=WQ=1_Query_9999_873_LBUB_LESSIMPLE.pdf`

* Figure 26 in 
    * `/path/to/code/Cost_Estimation/costResults/Minimal/406/DF=QD=WQ=1_Query_40_406_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/406/DF=QD=WQ=1_Query_64_406_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/406/DF=QD=WQ=1_Query_96_406_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/623/DF=QD=WQ=1_Query_40_623_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/623/DF=QD=WQ=1_Query_64_623_LBUB_LESSIMPLE.pdf`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/623/DF=QD=WQ=1_Query_96_623_LBUB_LESSIMPLE.pdf`

### Appendix H.2 
* Table 11, 12, 13, printed to console with data in
    * `/path/to/code/Cost_Estimation/costResults/Query/DF=nLogC_e=20_b=64_Query_LBUB_LESSIMPLE.critical`
    * `/path/to/code/Cost_Estimation/costResults/Query/DF=nLogC_e=20_b=64_Query_UBUB_LESSIMPLE.critical`
    * `/path/to/code/Cost_Estimation/costResults/Query/DF=nLogC_e=20_b=64_Query_LBLB_LESSIMPLE.critical`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/DF=nLogC_e=20_b=64_Minimal_LBUB_LESSIMPLE.critical`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/DF=nLogC_e=20_b=64_Minimal_UBUB_LESSIMPLE.critical`
    * `/path/to/code/Cost_Estimation/costResults/Minimal/DF=nLogC_e=20_b=64_Minimal_LBLB_LESSIMPLE.critical`

Specific Files:
* `*.critical` contain the intersection of lines for varying security relevant parameters
* `*.pdf` are the corresponding plots 
* `*.costs` contain the detailed costs of each point plotted in the corresponding `*.pdf` 


# Code Overview 

* main_interface.py: User interface, calls methods in `CostEstimation.py`
* Circuit.py: Code for instantiation of quantum operator W as in Section 4
* CostEstimation.py: Code for costing loop over the parameters m, y, z and levels k 
* CostFunctions.py: Code for quantum and classical GCost and T-Depth.
* plot_interface.py: Matplotlib code to produce figures from estimation results
* Tools.py: Data caching and formatted printing 
* TreeHeuristics.py: Acts as an interface to `LatticeHeuristics/enumTrees.py` and holds functions corresponding to H_k^M, N_kh^M and S_kh^M for the combined enumertion tree. 
* istarmap.py: https://stackoverflow.com/questions/57354700/starmap-combined-with-tqdm
* configSecurity.py: Holds the predefined parameters for Kyber ans AES 