# !/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Marcel Tiepelt"
__version__ = "1.0"
__status__ = "Research Code"

import TreeHeuristics as TH
import CostFunctions as CF

import Tools
import math

import multiprocessing as mpp

class CostEstimation:
    """
    Abstract class for cost estimation
    """
    def __init__(self, n, q, constants, log_M, log_Y, log_Z, bound, use_cached, step_size=1, infty_const=9999):
        """
        :param n:  BKZ blocksize
        :param q: LWE modulus
        :param constants: Dictionary with constants (cf. Section 3.1)
        :param log_M: Log of randomized bases
        :param log_Y: Log of combined nodes
        :param log_Z: Log of Jensen's gap
        :param bound: Upper or Lower bounds for prunung startegy.
        :param use_cached: Flag if results are cached to files.
        :param step_size:   Step-size of variable z representing the multiplicate Jensen's gap.
        :param infty_const: Constant denoting infinity for max depth; represents an exponent.
        """

        self.n = n
        self.modulus_q = q
        self.step_size = step_size

        self.log_M = log_M
        self.log_Y = log_Y
        self.log_Z = log_Z

        self.constants = constants
        self.bound = bound
        self.use_cached = use_cached

        self.infty_const = infty_const
        self.infty_trivial = -1

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
        :return: log_y
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

    def cheapest_cost(self, z_m_k_to_cost, maxdepths):
        """
        Finds the lowest Cost() object with the lowest combined quantum-classical cost over all z,m,k for individual max-depths.
        :param z_m_k_to_cost:   Dictionary mapping z,m,k to cost object
        :param maxdepths: MaxDepth restriction.
        :return: dict [md] --> [log_z] --> Cost()
        """

        md_z_to_cheapest = dict.fromkeys(maxdepths, {})

        for md in maxdepths:
            for log_z in z_m_k_to_cost.keys():
                costF = CF.CostFunctions(n=self.n, nbases=2**self.log_M, pprime=1, constants=self.constants, bound=self.bound)
                cheapest = Tools.Cost()
                cheapest.set_parameters(log_m=self.log_M, log_y=0, log_z=log_z, k=self.n)
                cheapest.set_gcost(costF.log_enumeration_cost_classical(), 0)
                md_z_to_cheapest[md][log_z] = cheapest

                for log_m in z_m_k_to_cost[log_z].keys():
                    for k in z_m_k_to_cost[log_z][log_m].keys():
                        cost = z_m_k_to_cost[log_z][log_m][k]

                        log_total_cheapest = Tools.add_logs(cheapest.log_gcost_classical, cheapest.log_gcost_quantum)
                        log_total = Tools.add_logs(cost.log_gcost_classical, cost.log_gcost_quantum)

                        if log_total < log_total_cheapest:
                            md_z_to_cheapest[md][log_z].overwrite(cost)
        return md_z_to_cheapest

class CostEstimationSimple(CostEstimation):
    """
    Cost estimation following Section 3, 4 and 5.
    Simple, but slow costing loop.
    """

    def __init__(self, n, q, constants, log_M, log_Y, log_Z, bound, use_cached, step_size=1):
        super().__init__(n, q, constants, log_M, log_Y, log_Z, bound, use_cached, step_size)

    def simpleInterface(self, circuit, maxdepths):
        """
        Interface for cost estimation. Manages multiprocessing over loop of Jensen's z.
        :param circuit: Circuit instance.
        :param maxdepths: List of possible maxdepths.
        :return: Dict [md] --> [z] --> Cost()
        """
        # Costing Loop
        max_depth = max(maxdepths)

        manager = mpp.Manager()
        z_m_k_costs = manager.dict()

        pool_inputs = []
        for log_z in range(0, self.log_Z + 1, self.step_size):
            costF = CF.CostFunctions(n=self.n, nbases=2 ** self.log_M, pprime=1, constants=self.constants, bound=self.bound)
            pool_inputs.append([z_m_k_costs, costF, log_z, self.log_M, self.log_Y, circuit, max_depth])

        pool_size = min(len(pool_inputs), int((0.9 * mpp.cpu_count()) - 2))
        print(f"Using: {pool_size} cores on {len(pool_inputs)} inputs.")

        if pool_size > 0:
            try:
                import tqdm
                tqdm_present = True
            except:
                tqdm_present = False

            pool = mpp.Pool(pool_size)
            if tqdm_present:
                for _ in tqdm.tqdm(pool.istarmap(self.loop_k_y, pool_inputs), total=len(pool_inputs)):
                    pass
            else:
                r = pool.starmap(self.loop_k_y, pool_inputs)
            pool.close()
            pool.join()

        # Find cheapest for each maxdepth and z
        md_z_to_cheapest = self.cheapest_cost(z_m_k_costs, maxdepths)

        return md_z_to_cheapest

    def loop_k_y(self, z_m_k_costs, costF, log_z, log_M, log_Y, circuit, max_depth):
        """
            Computes Cost() for all levels k, as long as TDepth(quantum) <= maxdepth
        """
        m_k_costs = {log_M: {}}
        for k in range(self.n - 1, -1, -1):
            h = self.n + 1 - k

            # Number of combined nodes limited by number of total nodes
            if self.bound == 'LBUB' or self.bound == 'LBLB':
                log_Hk_M = TH.log_Hk_M(k, self.n, nbases=2 ** log_M, pprime=1, bound='LB')
            else:  # 'UBUB'
                log_Hk_M = TH.log_Hk_M(k, self.n, nbases=2 ** log_M, pprime=1, bound='UB')

            log_Y_upper = min(log_Y, math.floor(log_Hk_M))
            log_y = self.binary_search_y(log_z, k, circuit, costF, max_depth, log_y_upper=log_Y_upper, log_y_lower=0)
            log_qw_depth = costF.log_depth_total(k=k, h=h, log_y=log_y, log_z=log_z, circuit=circuit)

            # If depth of circuit too large, any circuit for y' > y will also be to large
            if max_depth > 0 and log_qw_depth > max_depth:
                # If depth of circuit for single node too large, any circuit for k' < k will also be too large
                if log_y == 0:
                    z_m_k_costs[log_z] = m_k_costs
                    break
                else:
                    continue

            log_gcost_classical = costF.log_gcost_classical(k)
            log_gcost_quantum = costF.log_gcost_quantum_enumeration(k=k, h=h, log_y=log_y, log_z=log_z, circuit=circuit)

            qracm = min(log_Hk_M, log_y)

            cost = Tools.Cost()
            cost.set_parameters(log_M, log_y, log_z, k)
            cost.set_gcost(log_gcost_classical, log_gcost_quantum)
            cost.set_misc(qracm, log_qw_depth)
            m_k_costs[log_M][k] = cost

        z_m_k_costs[log_z] = m_k_costs

class CostEstimationFast(CostEstimation):
    """
    Cost estimation following section 3, 4 and 5.

    Complicated, but fast costing loop.
    """
    def __init__(self, n, q, constants, log_M, log_Y, log_Z, bound, use_cached, step_size=1):
        super().__init__(n, q, constants, log_M, log_Y, log_Z, bound, use_cached, step_size)

    def get_all_k_ascending(self, max_depth, prevZ_k_to_cost, k_start, k_lowest, costF, log_m, log_Y, log_z, circuit, verb=False):
        """
        Find levels k in enumeration tree (starting at lowest level) such that TDepth(quantum)<= maxdpeth.

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
        :return: [k] --> Cost()
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
            if self.bound == 'LBUB' or self.bound == 'LBLB':
                log_Hk_M = TH.log_Hk_M(k, self.n, nbases=2 ** log_m, pprime=1, bound='LB')
            else: # 'UBUB'
                log_Hk_M = TH.log_Hk_M(k, self.n, nbases=2 ** log_m, pprime=1, bound='UB')

            # ARTEFACT
            # CilinderPruningUpperBound
            # Example: enumTrees.CilinderPruningUpperBound.log_Hk(k=1, n=n, nbases=2**log_m, pprime=1) --> -inf
            if log_Hk_M == float('-inf'):
                log_Hk_M = -9999

            log_Y_upper = min(log_Y, math.floor(log_Hk_M))

            if log_qw_depth <= max_depth:
                # Find largest y that fits into MD, which is lowest total cost
                log_y = self.binary_search_y(log_z, k, circuit, costF, max_depth, log_y_upper=log_Y_upper, log_y_lower=log_y_lower)

                log_qw_depth = costF.log_depth_total(k=k, h=h, log_y=log_y, log_z=log_z, circuit=circuit)

                log_gcost_classical = costF.log_gcost_classical(k)
                log_gcost_quantum = costF.log_gcost_quantum_enumeration(k=k, h=h, log_y=log_y, log_z=log_z, circuit=circuit, verb=verb)

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

    def estimate_quantum_enumeration(self, max_depth, circuit, log_M_lower):
        """
        Find lowest Cost() for all values z,m,k
        :param max_depth: Max depth constraint for quantum circuit.
        :param circuit: Circuit for instantiation of quantum operator W (cf. Section 3.1)
        :param log_M_lower:  Lower bound for log number of combined, randomized bases for extreme pruning.
        :return: [Jensen's Gap z] --> cost()
        """
        # Results
        z_m_k_to_cost = {}
        log_jensen_z_range = range(0, self.log_Z + 1, self.step_size)

        # Process cached data
        if self.use_cached == True:
            costFiles = Tools.DataFiles(self.n, self.modulus_q, max_depth, self.constants, self.bound)
            z_m_k_to_cost = costFiles.read_cached_data_from_file(circuit, self.log_M, self.log_Y, self.step_size, self.log_Z, self.step_size)
            log_jensen_z_range = [z for z in log_jensen_z_range if z not in z_m_k_to_cost.keys()]

        # Estimate rest
        for count, log_z in enumerate(log_jensen_z_range):
            m_k_to_cost = {}
            for log_m in range(log_M_lower, self.log_M + 1, self.step_size):

                prevZ_k_to_cost = {}

                if len(z_m_k_to_cost) > 0:
                    prev_log_z = log_z - self.step_size  # Could be improved by checking for all previous values of z
                    if prev_log_z in z_m_k_to_cost.keys():
                        if log_m in z_m_k_to_cost[prev_log_z].keys():
                            prevZ_k_to_cost = z_m_k_to_cost[prev_log_z][log_m]

                k_start = self.n - 1
                costF = CF.CostFunctions(n=self.n, nbases=2 ** log_m, pprime=1, constants=self.constants, bound=self.bound)

                # TO GET CONSISTENT RESULTS WITH CRYPTO24 Submission
                # START
                if max_depth == self.infty_trivial:
                    log_qw_depth = costF.log_depth_total(k=0, h=self.n + 1, log_y=0, log_z=log_z, circuit=circuit)
                    log_gcost_quantum = costF.log_gcost_quantum_enumeration(k=0, h=self.n + 1, log_y=0, log_z=log_z, circuit=circuit)

                    # Cost
                    cost = Tools.Cost()
                    cost.set_parameters(log_m, 0, log_z, k=0)
                    cost.set_gcost(0, log_gcost_quantum)
                    cost.set_misc(0, log_qw_depth)
                    m_k_to_cost[log_m] = {0: cost}
                    continue
                # END

                k_to_cost = self.get_all_k_ascending(max_depth, prevZ_k_to_cost, k_start=k_start, k_lowest=0, costF=costF, log_m=log_m, log_Y=self.log_Y, log_z=log_z, circuit=circuit)

                m_k_to_cost[log_m] = k_to_cost

            if self.use_cached:
                costFiles = Tools.DataFiles(self.n, self.modulus_q, max_depth, self.constants, self.bound)
                costFiles.write_data_cost_file(m_k_to_cost, circuit, self.log_M, self.log_Y, self.step_size, log_z)

            z_m_k_to_cost[log_z] = m_k_to_cost
        # Estimation done

        md_z_to_cheapest = self.cheapest_cost(z_m_k_to_cost, [max_depth])

        return md_z_to_cheapest[max_depth]
