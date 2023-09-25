# Experiments on the number of children

In this directory we release the code used to run and report the esperiments on:
- The number of children of a node in the enumeration tree, C(g) (App. B.2)
- The number of descendants of a node in the enumeration tree, W_{k,h}(g) (App. B.2)
- The average number of nodes in a subtree, to measure the multiplicative Jensen's Gap (App. C)
- The unbiased sample variance of nodes in a subtree (App G.3)

## Preparing to run the experiments

The code used to run the experiments is based on fplll 5.4.2.
As such, the [code dependencies](https://github.com/fplll/fplll/tree/5.4.2#installation-from-source) of fplll should be installed in the system.

However, for the sake of ease of use, our Makefile will take care of fetching fplll, patching it and compiling it.

Downloading, patching, and locally installing fplll:
```
make fetch_fplll
make patch_fplll
make first_compilation
```

In order to analyse and plot the results building the alfk module located in `LatticeHeuristics/alfk` is required (see the readme file in that directory).
Furhtermore, `sagemath` and GNU `parallel` are required.

## Running and plotting experiments for the number of children and descendants of a node (App. B.2)

Prepare directory structure for experiments:
```
make prepare_experiments
```

Set the instances to run (meaning LWE and BKZ properties) by editing `subtree_experiments`. 

Run experiments:
```
make run_subtree_experiments
```

The resulting statistics are contained in `stats`.
The plots for the number of children can be found in `plots/<instance>/<tour>-<block index>-0.pdf`.
The plots for the number of descendants can be found in `subtree-plots/<instance>/<tour>-<block index>-0.{png,pdf}`.

If you already run the experiments and only want to regenerate the plots, run
```
make plot_subtree_experiments
```

## Running and plotting experiments for the Jensen's gap (App C and G.3)

Prepare directory structure for experiments:
```
make prepare_experiments
```

Set the instances to run (meaning LWE and BKZ properties) by editing `jensen_experiments`. 

Run experiments:
```
make run_jensen_experiments
```

The resulting statistics are contained in `pruned/stats`.
The plots for App C can be found in `pruned/jensen-plots/<instance>/<tour>.{png,pdf}`.
The plots for App G.3 can be found in `pruned/jensen-plots/<instance>/<tour>-incl-sigma.{png,pdf}`.

If you already run the experiments and only want to regenerate the plots, run
```
make plot_jensen_experiments
```

### Where to find the figures

After running all experiments, the figures should have been reproduced. We list where to find each in particular. We note that tour and vectors are counted from 0, so "tour 10" will be saved as 9, and so on.
- Fig 5(a): `/fplll_experiments/plots/n72-q97-sd1.0-m87-beta40-tours15-ghbound1.1-seed1337/9-79-0.pdf`
- Fig 5(b): `/fplll_experiments/plots/n72-q97-sd1.0-m87-beta40-tours15-ghbound1.1-seed1337/14-69-0.pdf`
- Fig 6(a): `/fplll_experiments/subtree-plots/n72-q97-sd1.0-m87-beta40-tours15-ghbound1.1-seed1337/10/9-79-0.pdf`
- Fig 6(b): `/fplll_experiments/subtree-plots/n72-q97-sd1.0-m87-beta40-tours15-ghbound1.1-seed1337/20/9-79-0.pdf`
- Fig 7(a): `/fplll_experiments/subtree-plots/n72-q97-sd1.0-m87-beta40-tours15-ghbound1.1-seed1337/10/9-89-0.pdf`
- Fig 7(b): `/fplll_experiments/subtree-plots/n72-q97-sd1.0-m87-beta40-tours15-ghbound1.1-seed1337/20/9-89-0.pdf`
- Fig 8(a): `/fplll_experiments/subtree-plots/n72-q97-sd1.0-m87-beta40-tours15-ghbound1.1-seed1337/10/14-79-0.pdf`
- Fig 8(b): `/fplll_experiments/subtree-plots/n72-q97-sd1.0-m87-beta40-tours15-ghbound1.1-seed1337/20/14-79-0.pdf`
- Fig 9(a): `/fplll_experiments/pruned/jensen-plots/n72-q97-sd1.0-m87-beta50-tours15-ghbound1.1-seed1337/20/9.pdf`
- Fig 9(b): `/fplll_experiments/pruned/jensen-plots/n72-q97-sd1.0-m87-beta60-tours15-ghbound1.1-seed1337/20/9.pdf`
- Fig 9(c): `/fplll_experiments/pruned/jensen-plots/n72-q97-sd1.0-m87-beta70-tours10-ghbound1.1-seed1337/20/9.pdf`
- Fig 10(a): `/fplll_experiments/pruned/jensen-plots/n72-q97-sd1.0-m87-beta50-tours15-ghbound1.1-seed1337/30/9.pdf`
- Fig 10(b): `/fplll_experiments/pruned/jensen-plots/n72-q97-sd1.0-m87-beta60-tours15-ghbound1.1-seed1337/30/9.pdf`
- Fig 10(c): `/fplll_experiments/pruned/jensen-plots/n72-q97-sd1.0-m87-beta70-tours10-ghbound1.1-seed1337/30/9.pdf`
- Fig 19(a): `/fplll_experiments/pruned/jensen-plots/n72-q97-sd1.0-m87-beta50-tours15-ghbound1.1-seed1337/20/9-incl-sigma.pdf`
- Fig 19(b): `/fplll_experiments/pruned/jensen-plots/n72-q97-sd1.0-m87-beta60-tours15-ghbound1.1-seed1337/20/9-incl-sigma.pdf`
- Fig 19(c): `/fplll_experiments/pruned/jensen-plots/n72-q97-sd1.0-m87-beta70-tours10-ghbound1.1-seed1337/20/9-incl-sigma.pdf`
- Fig 20(a): `/fplll_experiments/pruned/jensen-plots/n72-q97-sd1.0-m87-beta50-tours15-ghbound1.1-seed1337/30/9-incl-sigma.pdf`
- Fig 20(b): `/fplll_experiments/pruned/jensen-plots/n72-q97-sd1.0-m87-beta60-tours15-ghbound1.1-seed1337/30/9-incl-sigma.pdf`
- Fig 20(c): `/fplll_experiments/pruned/jensen-plots/n72-q97-sd1.0-m87-beta70-tours10-ghbound1.1-seed1337/30/9-incl-sigma.pdf`


## What every file does

- Makefile: scrips used to prepare the environment and run the experiments
- bkz.py: python wrapper to BKZ with non-pruned enumeration
- diff/*.patch: changes to the fplll library required to run the experiments (see details section below)
- enum_ub_exp.sh: runs a single experiment consisting of generating an LWE challenge, passing it to fplll, then plotting statistics about it
- enum_subtree_jensen.sh: script generating LWE instances, running experiments and generating plots for App. C and G
- enum_subtree_size.sh: script generating LWE instances, running experiments and generating plots for App. B
- gen_patch.sh: script that automatically generates the `diff/` directory from edits to the fplll source
- jensen_experiments: list of experiment parameters for the experiments in App. C and G
- lll_sim.py: simulates LLL rediced basis profile
- lwe.py: generates an LWE problem instance and outputs it in a format that fplll can read
- plot_children_predictions.py: generates plots for the distribution of children (App. B)
- plot_jensen_gap.py: generates plots for the average and unbiased sample variance of tree sizes (App. C and G)
- plot_subtree_predictions.py: generates plots for the distribution of descendants (App. B)
- se93_strat.py: simple preprocessing and pruning strategy for BKZ with no prning. Crucially, it sets the pruning radii!
- settings.py: contains basic parameters for the scripts
- subtree_experiments: list of experiment parameters for the experiments in App. B
- tree.py: given a dump of an enumeration tree from fplll, it extracts statistics about it and draws a plot
- utilities.py: utilities used by files above

## Changes made to fplll

The experiments run on a modified version of fplll.
Patches can be found under `diff/`.

Base: fplll-5.4.2

Modified files:
- fplll/main.h, fplll/main.cpp: added flags for a full enumeration tree or just the measured averages and sample variances for the number of children
- defs.h, bkz.h, bkz_param.h: generic code for storing trees or stats
- bkz.cpp: code that creates requested dump files, then collects the data/statistics at every SVP call of every tour and appends it to the dump files in a separate entry
- tests/test_bkz.cpp, tests/test_svp.cpp: reduced the number of test for respective functions, to call them only on a simpler input case (and generate then less output)
- enum/enumerate.h: removes support to external enumerators, to force basic (modified) enumeration routine; adds handles to get statistics
- enum/enumerate.cpp: prints out debug information about what routine fplll is executing, and some of the statistics collected, creates/destroys tree structure before/after enumeration is run
- enum/enumerate_base.h: defines tree structure class for data collection, prints some debug information about which routine is being executed
- enum/enumerate_base.cpp: where enumeration actually happens: defines functions for transforming tree data into JSON, populates/creates tree representation of the enumeration tree, computes statistics about number of children per node, implements methods that return JSON representations of tree and stats, prints some debug information about which routine is being executed
- fplll/svpcvp.cpp: contains the api for SVP solving; our changes here are only logging that could be useful during debug
ging

## Contributors

Fernando Virdia