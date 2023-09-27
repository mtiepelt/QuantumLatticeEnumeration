# Installation

We assume the code is stores in some root folder:

```
/path/to/code/
```

## Dependencies

* Python 3, tqmd 


## Build alfk 

```
cd /path/to/code/LatticeHeuristics/alfk
python3 build_alfk.py
```

# Generating cost estimation


## Section 3

Naive depth of QPE(W) using the LWE Estimator (page 12).

```
cd /path/to/code/LWEEstimator 

sage our_kyber_estimates.py
```



## Section 5

```
cd /path/to/code/Cost_Estimation

python3 main_interface.py [-num-cores X]
```
This computes all results from the paper, saves intermediate results to `costData` and final results are written to `costResults`.
Additional informations can be found via

```
cd /path/to/code/Cost_Estimation

python3 main_interface.py --help
```


## Results

From the paper:

* ...  `/path/to/code/costResults/Query` contains estimations with instantiation of quantum opererator W as in Section 4.1
* ... `/path/to/code/costResults/Minimal` contains estimations with instantiation of quantum opererator W as in Section 4.2
* The graphs of Figure 4 in Section 5 and Figure 12,13,14,15 of th appendix can be found in the folder with the corresponding blocksize:
    * for Kyber512 in folder `406`
    * for Kyber768 in folder `623`
    * for Kyber1024 in folder `873`
* Table 6 and 9 (Summary of values for Jensen's gap at corssover points) is printed to console, and can also be found in `*.critical` files 


Specific Files:
* `*.critical` contain the intersection of lines for varying security relevant parameters
* `*.pdf` are the corresponding plots 
* `*.costs` contain the detailed costs of each point plotted in the corresponding `*.pdf` 


# Code Overview 

* Circuit.py: Code for instantiation of quantum operator W as in Section 5
* CostEstimation.py: Code for costing loop over the parameters m, y, z and levels k 
* CostFunctions.py: Code for quantum and classical GCost and T-Depth
* plot_interface.py: Matplotlib code to produce figures from estimation results
* Tools.py: Data caching and formatted printing 
* TreeHeuristics.py: Acts as an interface to `LatticeHeuristics/enumTrees.py` and holds functions corresponding to H_k^M, N_kh^M and S_kh^M for the combined enumertion tree 