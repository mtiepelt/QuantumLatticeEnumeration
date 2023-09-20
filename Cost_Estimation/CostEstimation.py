# !/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Anonymous"
__version__ = "0.0"
__status__ = "Research"

import os

import TreeHeuristics as TH
import CostFunctions as CF

import Tools
import sys
import plot_interface
import math
sys.path.append("../LatticeHeuristics")
import enumTrees
import Circuits

from dataclasses import asdict
from collections import defaultdict

from multiprocessing import Pool
import multiprocessing


class CostEstimation:
    """
    Cost estimation following section 3, 4 and 5.
    """

    def __init__(self, n, q, boundCoefficient, const_reps_QPE, const_reps_W, num_nodes_branching_dfs, force_DF, lattice_class=enumTrees.CilinderPruningLowerBound, bound='upper', use_cached=True, num_cores=2):
        """

        :param n:  BKZ blocksize
        :param q: LWE modulus
        :param boundCoefficient: Bound on node degree for FMV (cf. Section 3.1)
        :param boundCoefficient:    Bound C on the node degree used for enumeration of the tree (cf. Section 3.1).
        :param const_reps_QPE:  Constant 1/b for number of repetitions in QPE (cf. Corollary 1)
        :param const_reps_W: Constant epsilon for number of repetition of QPE (cf. Section 3.1)
        :param num_nodes_branching_dfs: Number of nodes checked for branchi
        :param force_DF: Number of DMV calls in FMV (cf. Section 3.1)
        :param lattice_class: Pruning strategy of lattice heuristics.
        :param bound: Upper or Lower bounds for prunung startegy.
        :param use_cached: Flag if results are cached to files.
        :param num_cores: Number of cores for multiprocessing.
        """
        self.n = n
        self.modulus_q = q
        self.boundCoefficient = boundCoefficient

        self.const_reps_QPE = const_reps_QPE
        self.const_reps_W = const_reps_W
        self.num_nodes_branching_dfs = num_nodes_branching_dfs
        self.force_DF = force_DF

        self.bound = bound
        self.lattice_class = lattice_class

        self.use_cached = use_cached
        self.num_cores = num_cores

    def get_all_k_ascending(self, max_depth, prevZ_k_to_cost, k_start, k_lowest, costF, log_m, log_Y, log_z, circuit, verb=False):
        """
        Find levels k in enumeration tree such that quantum circuit rooted at level k is below MAXDEPTH.

        :param max_depth: Max depth constraint for quantum circuit.
        :param prevZ_k_to_cost: Optimal level k of previous value for Jensen's Gap z.
        :param k_start: Level k to start estimating cost.
        :param k_lowest: Lowest possible k to estimate cost to.
        :param costF: Cost class object.
        :param log_m:  Log Number of combined, randomized bases for extreme pruning.
        :param log_Y:  Log Bound on number of combined nodes for a virtual tree on level k.
        :param log_z: Log Jensen's Gap z.
        :param circuit: Circuit for instantiation of quantum operator W (cf. Section 3.1)
        :param verb:
        :return:
    """
        assert k_start < self.n, f" k_start {k_start} must be smaller than n:={self.n}, No quantum circuit for single layer!"
        assert k_start >= k_lowest, f" Expected k_start > k_lowest := {k_lowest}, was k_start={k_start}"

        # Keep only cheapest
        log_enumeration_classical = costF.log_enumeration_cost_classical()
        cheapest = Tools.Cost()
        cheapest.set_parameters(log_m, 0, log_z, self.n + 1)
        cheapest.set_gcost(log_enumeration_classical, 0)
        cheapest.set_misc(0, 0)

        """
        Optimization to exaustive search:
            From z --> z+1
              Depth decreases by factor sqrt(2)/2,
            Increasing y --> y+1 Increases depth again by a factor of 2,
              This decreases GCost by factor of at most sqrt(2)/2
    
            Thus after increasing z, also increasing y can lower cost by at most 1/2 * sqrt(2)/2
              But if classical cost of z is 2* larger than gcost of z, then this does not result in a lower total cost
              Therefore, in that case, we can skip all previous levels k and just try to further decrease k
        """
        for k in reversed(sorted(prevZ_k_to_cost.keys())):
            if prevZ_k_to_cost[k].log_gcost_classical > (prevZ_k_to_cost[k].log_gcost_quantum + (log_Y - prevZ_k_to_cost[k].log_y)):
                k_start = k
            else:
                break
        k_to_cost = {}

        for k in range(k_start, max(k_lowest-1, -1), -1):
            h = self.n + 1 - k

            log_y_lower = 0
            if k in prevZ_k_to_cost.keys():
                log_y_lower = prevZ_k_to_cost[k].log_y

            # Depth of smallest circuit
            log_qw_depth = costF.log_depth_total(k=k, h=h, log_y=0, log_z=log_z, circuit=circuit)

            # Number of combined nodes limited by number of total nodes
            log_Y_upper = min(log_Y, math.ceil(costF.treeH.log_Hk_M(k, self.n, nbases=2**log_m, pprime=1)))

            if log_qw_depth < max_depth:
                # Find largest y that fits into MD, which is lowest total cost
                log_y = self.binary_search_y(log_z, k, circuit, costF, max_depth, log_y_upper=log_Y_upper, log_y_lower=log_y_lower)

                log_gcost_classical = costF.log_gcost_classical(k)
                log_gcost_quantum = costF.log_gcost_quantum_enumeration(k, h, log_y=log_y, log_z=log_z, circuit=circuit, verb=verb)
                log_Hk_M = costF.treeH.log_Hk_M(k=k, n=self.n, nbases=2 ** log_m, pprime=1)
                qracm = min(log_Hk_M, log_y)

                # Cost
                cost = Tools.Cost()
                cost.set_parameters(log_m, log_y, log_z, k)
                cost.set_gcost(log_gcost_classical, log_gcost_quantum)
                cost.set_misc(qracm, log_qw_depth)

                k_to_cost[k] = cost
            else:
                return k_to_cost

        return k_to_cost

    def pool_func_get_costs(self, proc_id, data, log_Y, max_depth, return_dict):
        """
        Entry point for multiprocessing.
        :param proc_id: Process id.
        :param data: Process data input.
        :param log_Y: Log Bound on number of combined nodes for a virtual tree on level k.
        :param max_depth: Max depth constraint for quantum circuit.
        :param return_dict:
        :return:
        """
        log_z = data['log_z']
        log_m = data['log_m']
        circuit = data['circuit']
        treeH = data['treeH']
        k_start = self.n - 1  #data['k_start']
        prevZ_k_to_cost = data['prevZ_k_to_cost']
        costF = CF.CostFunctions(treeH, n=self.n, bound_coefficient=self.boundCoefficient, nbases=2 ** log_m, pprime=1,
                                 const_reps_QPE=self.const_reps_QPE, const_reps_W=self.const_reps_W, num_nodes_branching_dfs=self.num_nodes_branching_dfs, force_DF=self.force_DF)
        k_to_cost = self.get_all_k_ascending(max_depth, prevZ_k_to_cost, k_start=k_start, k_lowest=0, costF=costF, log_m=log_m, log_Y=log_Y, log_z=log_z, circuit=circuit)
        return_dict[proc_id] = k_to_cost
        return


    def binary_search_y(self, log_z, k, circuit, costF, max_depth, log_y_upper, log_y_lower):
        """
        Binary search for largest possible number of combined nodes 2^y, such that quantum circuit is within MaxDepth:
            Increasing y always lowers overall cost, as it moves cost under the square-root.

        :param log_z: Log Jensen's Gap z
        :param k: Level k of enumeration tree.
        :param circuit: Circuit instantiatinfg quantum operator W (cf. Section 5)
        :param costF: Cost function object.
        :param max_depth: Max depth constraint for quantum circuit.
        :param log_y_upper: Upper bound on number of combined nodes.
        :param log_y_lower: Lower bound on number of combined nodes.
        :return:
        """
        if log_y_lower == log_y_upper:
            return log_y_lower

        h = self.n - k + 1
        log_y = log_y_lower
        while log_y_lower <= log_y_upper:
            log_y = math.ceil((log_y_lower + log_y_upper) / 2)

            log_qw_depth = costF.log_depth_total(k=k, h=h, log_y=log_y, log_z=log_z, circuit=circuit)

            if log_qw_depth <= max_depth:
                if log_y_lower == log_y_upper:
                    return log_y
                log_y_lower = log_y
            else:
                log_y_upper = log_y - 1
        return log_y

    def optimal_cost_is_x(self, log_costs_x, log_costs_y):
        """
        Comparison of lists of costs. All log costs in list are summed

        :param log_costs_x:
        :param log_costs_y:
        :return: Returns true, if sum of costs in first listr (x) is smaller.
        """
        log_total_x = log_costs_x[0]
        for x in log_costs_x[1:]:
            log_total_x = Tools.add_logs(log_total_x, x)

        log_total_y = log_costs_y[0]
        for y in log_costs_y[1:]:
            log_total_y = Tools.add_logs(log_total_y, y)

        if log_total_x < log_total_y:
            return True
        else:
            return False

    def tqdm_body(self, z_m_k_to_cost, log_z, step_size_z, log_M, log_M_lower, log_Y, step_size, circuit, max_depth):
        """
        Body for TQDM progress bar.
        :param z_m_k_to_cost: Dictionary mapping [Jenzen's Z] --> [#Combined Bases] --> [Level k] --> Cost
        :param log_M:   Log Bound on number of combined, randomized bases for extreme pruning.
        :param log_Y:   Log Bound on number of combined nodes for a virtual tree on level k.
        :param step_size:   Step-size for looping over m,y.
        :param log_Z:   Log bound for Jensen's Gap z.
        :param step_size_Z: Step-size for loopzing over z.
        :param log_M_lower: Lower log bound for Bound on number of combined, randomized bases.
        :param circuit: Circuit instantiatinfg quantum operator W (cf. Section 5)
        :param max_depth: Max depth constraint for quantum circuit.
        :return:
        """
        manager = multiprocessing.Manager()
        return_dict = manager.dict()

        treeH = TH.TreeHeuristics(lattice_class=self.lattice_class, lattice_bound=self.bound)
        proc_id = 0
        pool_inputs = []

        DBG_cnt_saving = 0
        for log_m in range(log_M_lower, log_M + 1, step_size):
            proc_id += 1

            data = {
                'log_z': log_z,
                'log_m': log_m,
                'treeH': treeH,
                'circuit': circuit,
                'prevZ_k_to_cost': {}
            }

            if len(z_m_k_to_cost) > 0:
                prev_log_z = log_z - step_size_z  # Could be improved by checking for all previous values of z
                if prev_log_z in z_m_k_to_cost.keys():
                    if log_m in z_m_k_to_cost[prev_log_z].keys():
                        data['prevZ_k_to_cost'] = z_m_k_to_cost[prev_log_z][log_m]

            pool_inputs.append([proc_id, data, log_Y, max_depth, return_dict])

        # Run Pool
        pool_size = min(self.num_cores, multiprocessing.cpu_count())

        if pool_size > 0:
            pool = Pool(pool_size)
            # pool_inputs.reverse() # depending on the order in which pool_input is populated, reversing it can help the overall runtime
            results = pool.starmap(self.pool_func_get_costs, pool_inputs)
            pool.close()
            pool.join()

        m_k_to_cost = {}
        # Extract data relative to z
        for proc_id, k_to_cost in return_dict.items():
            m_k_to_cost[next(iter(k_to_cost.values())).log_m] = k_to_cost

        if self.use_cached:
            costFiles = Tools.DataFiles(max_depth, self.n, self.boundCoefficient, self.modulus_q, self.const_reps_QPE, self.const_reps_W, self.num_nodes_branching_dfs, self.force_DF)
            costFiles.write_data_cost_file(m_k_to_cost, circuit, log_M, log_Y, step_size, log_z)

        return m_k_to_cost

    def get_cost_all(self, max_depth, circuit, log_M, log_Y, step_size, log_Z, step_size_z, log_M_lower):
        """
        Get cost for all combinations of parameters.
        :param circuit: Circuit instantiatinfg quantum operator W (cf. Section 5)
        :param max_depth: Max depth constraint for quantum circuit.
        :param log_M:   Log bound on number of combined, randomized bases for extreme pruning.
        :param log_Y:   Log bound on number of combined nodes for a virtual tree on level k.
        :param step_size: Stepsize for number of combined bases M.
        :param log_Z: Log Bound on Jensen's Gap z.
        :param step_size_z: Step-size for looping over z.
        :param log_M_lower: Lower log bound on number of combined, randomized bases for extreme pruning.
        :return: z_m_k_to_cost: Dictionary mapping [Jenzen's Z] --> [#Combined Bases] --> [Level k] --> Cost
        """
        # Results
        z_m_k_to_cost = {}
        log_jensen_z_range = range(0, log_Z+1, step_size_z)

        # Process cached data
        if self.use_cached == True:
            costFiles = Tools.DataFiles(max_depth, self.n, self.boundCoefficient, self.modulus_q, self.const_reps_QPE, self.const_reps_W, self.num_nodes_branching_dfs, self.force_DF)
            z_m_k_to_cost = costFiles.read_cached_data_from_file(circuit, log_M, log_Y, step_size, log_Z, step_size_z)

            log_jensen_z_range = [z for z in log_jensen_z_range if z not in z_m_k_to_cost.keys()]

        # Estimate rest
        if len(log_jensen_z_range) > 0:
            try:
                import tqdm
                tqdm_present = True
            except:
                tqdm_present = False

            if tqdm_present:
                for count, log_z in enumerate(tqdm.tqdm(log_jensen_z_range)):
                    m_k_to_cost = self.tqdm_body(z_m_k_to_cost, log_z, step_size_z, log_M, log_M_lower, log_Y, step_size, circuit, max_depth)
                    z_m_k_to_cost[log_z] = m_k_to_cost

            else:
                for count, log_z in enumerate(log_jensen_z_range):
                    m_k_to_cost = self.tqdm_body(z_m_k_to_cost, log_z, step_size_z, log_M, log_M_lower, log_Y, step_size, circuit, max_depth)
                    z_m_k_to_cost[log_z] = m_k_to_cost

        return z_m_k_to_cost

    def find_cheapest_among_m(self, z_m_k_to_cost, log_M):
        """
        Find cheapest cost among dictionary.
        :param z_m_k_to_cost:
        :param log_M: Log bound on number of combined, randomized bases for extreme pruning.
        :return:
        """
        # Find cheaptest m
        treeH = TH.TreeHeuristics(lattice_class=self.lattice_class, lattice_bound=self.bound)

        z_to_cheapest = {}
        for log_z in z_m_k_to_cost.keys():
            costF = CF.CostFunctions(treeH, n=self.n, bound_coefficient=0, nbases=2 ** log_M, pprime=1, const_reps_QPE=0, const_reps_W=0, num_nodes_branching_dfs=0, force_DF=0)
            cheapest = Tools.Cost()
            cheapest.set_parameters(log_m=log_M, log_y=0, log_z=log_z, k=self.n)
            cheapest.set_gcost(costF.log_enumeration_cost_classical(), 0)

            z_to_cheapest[log_z] = cheapest

            for log_m in z_m_k_to_cost[log_z].keys():
                for k in z_m_k_to_cost[log_z][log_m].keys():
                    cost = z_m_k_to_cost[log_z][log_m][k]
                    log_costs_x = [cost.log_gcost_quantum, cost.log_gcost_classical, cost.qracm]
                    log_costs_y = [cheapest.log_gcost_quantum, cheapest.log_gcost_classical, cheapest.qracm]

                    x_is_cheaper = self.optimal_cost_is_x(log_costs_x, log_costs_y)
                    if x_is_cheaper:
                        z_to_cheapest[log_z].overwrite(cost)

        return z_to_cheapest

    def estimate_quantum_enumeration(self, max_depth, circuit, log_M, log_Y, step_size, log_Z, step_size_Z, log_M_lower):
        """
        Run cost estimation and find lowest overall cost.
        :param circuit: Circuit instantiatinfg quantum operator W (cf. Section 5)
        :param max_depth: Max depth constraint for quantum circuit.
        :param log_M:   Log Bound on number of combined, randomized bases for extreme pruning.
        :param log_Y:   Log Bound on number of combined nodes for a virtual tree on level k.
        :param step_size:   Step-size for looping over m,y.
        :param log_Z:   Log bound for Jensen's Gap z.
        :param step_size_Z: Step-size for loopzing over z.
        :param log_M_lower: Lower log bound on number of combined, randomized bases for extreme pruning.
        :return: Dictionary [Jensen's Gap z] --> cost
        """
        if max_depth == 0:
            return self.cost_no_max_depth(circuit, log_M, step_size, log_Z, step_size_Z)
        else:
            z_m_k_to_cost = self.get_cost_all(max_depth=max_depth, circuit=circuit, log_M=log_M, log_Y=log_Y, step_size=step_size, log_Z=log_Z, step_size_z=step_size_Z,
                                              log_M_lower=log_M_lower)
            z_to_cheapest = self.find_cheapest_among_m(z_m_k_to_cost, log_M)

            return z_to_cheapest

    def cost_no_max_depth(self, circuit, log_M, step_size, log_Z, step_size_Z):
        """
        Run cost estimation and find lowest overall cost for no limit on maxDepth.
        :param circuit: Circuit instantiatinfg quantum operator W (cf. Section 5)
        :param max_depth: Max depth constraint for quantum circuit.
        :param log_M:   Log Bound on number of combined, randomized bases for extreme pruning.
        :param step_size:   Step-size for looping over m,y.
        :param log_Z:   Log bound for Jensen's Gap z.
        :param step_size_Z: Step-size for loopzing over z.
        :return:
        """
        treeH = TH.TreeHeuristics(lattice_class=self.lattice_class, lattice_bound=self.bound)
        # costF = CF.CostFunctions(treeH, n=self.n, bound_coefficient=self.boundCoefficient, nbases=2 ** log_M, pprime=1,
        #                          const_reps_QPE=self.const_reps_QPE, const_reps_W=self.const_reps_W, num_nodes_branching_dfs=self.num_nodes_branching_dfs, force_DF=self.force_DF)
        # log_cost_classical = costF.log_enumeration_cost_classical()

        z_to_cost = {}
        for z in range(0, log_Z+1, 1):
            cost = Tools.Cost()
            cost.set_parameters(log_m=log_M, log_y=0, log_z=z, k=0)
            costF = CF.CostFunctions(treeH, n=self.n, bound_coefficient=self.boundCoefficient, nbases=2 ** log_M, pprime=1,
                                     const_reps_QPE=self.const_reps_QPE, const_reps_W=self.const_reps_W, num_nodes_branching_dfs=self.num_nodes_branching_dfs, force_DF=self.force_DF)

            log_gcost_quantum = costF.log_gcost_quantum_enumeration(k=0, h=self.n + 1, log_y=0, log_z=z, circuit=circuit)


            log_t_depth = costF.log_depth_total(k=0, h=self.n + 1, log_y=0, log_z=z, circuit=circuit)

            cost.set_gcost(0, log_gcost_quantum)
            cost.set_misc(0,log_t_depth)
            z_to_cost[z] = cost
        return z_to_cost