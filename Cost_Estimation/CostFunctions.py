#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Anonymous"
__version__ = "0.0"
__status__ = "Research"

import Tools
import math
import TreeHeuristics
from functools import lru_cache

class CostFunctions:
    """
        Cost function for classical and quantum enumeration cost according to Section 3,4 and 5
    """
    def __init__(self, treeH, n, bound_coefficient, nbases, pprime, num_nodes_branching_dfs, const_reps_QPE, const_reps_W, force_DF = -1):
        """

        :param treeH:
        :param n: BKZ blocksize.
        :param bound_coefficient:   Bound on node degree for DMV and FMV (cf. Section 3.1)
        :param nbases: Number of combined random bases for extreme pruning
        :param pprime: Success probability of shortest vector to be in enumeration tree
        :param num_nodes_branching_dfs: Bound on number of nodes in binary search for DMV and FMV (cf. Section 3.1)
        :param const_reps_QPE: Constant epsilon for number of repetitions of QPE (cf. Appendix F)
        :param const_reps_W: Constant 1/b for number of applications of quantum operator W in QPE (cf. Appendix F)
        :param force_DF: Number of DMV per FMV (cf. Section 3.1)
        """
        self.n = n
        self.bound_coefficient = bound_coefficient

        self.num_nodes_branching_dfs = num_nodes_branching_dfs
        self.force_DF = force_DF

        self.const_reps_QPE = const_reps_QPE
        self.const_reps_W = const_reps_W

        self.nbases = nbases
        self.pprime = pprime

        self.treeH = treeH

    def log_num_combined_nodes(self, k, log_y):
        """
        Return number of combined nodes in virtual root at level k.
        :param k: Level k
        :param log_y: Log bound on number of combined nodes.
        :return:
        """
        if k == 0:
            return 0
        log_Hk_M = self.treeH.log_Hk_M(k=k, n=self.n, nbases=self.nbases, pprime=self.pprime)
        return min(log_Hk_M, log_y)


    def log_virtual_treesize(self, k, h, log_num_comb_nodes_on_k):
        """
        Treesize of virtual tree rooted at level k
        :param k: Level k
        :param h: Height of virtual tree
        :param log_num_comb_nodes_on_k: Number of combined nodes on level k
        :return:
        """
        # Subtree + subtree-root
        log_single_subtreesize = Tools.add_logs(self.treeH.log_avg_N_kh_M(k=k, h=h, n=self.n, nbases=self.nbases, pprime=self.pprime), 0)

        # Virtual tree without virtual root
        log_virtual_treesize_no_root = log_num_comb_nodes_on_k + log_single_subtreesize

        # Virtual tree + root
        log_virtual_treesize = Tools.add_logs(0, log_virtual_treesize_no_root)

        return log_virtual_treesize


    def log_reps_W_kh(self, k, h, log_num_comb_nodes, log_z, verb=False):
        """
        Number of applications of quantum operator W in QPE (cf. Section 3.1)
        :param k: Level k
        :param h: Height of subtree
        :param log_num_comb_nodes: Number of combined nodes on level k
        :param log_z: Log Jensen's z
        :param verb:
        :return:
        """
        b = self.const_reps_W
        # Treesize
        log_virt_treesize = self.log_virtual_treesize(k, h, log_num_comb_nodes)

        # W/Q: b * 1/2^z * sqrt(treesize * h)
        log_sqrt = (log_virt_treesize + Tools.doLog2(h)) / 2
        log_num_repetitions_W = Tools.doLog2(b) - log_z + log_sqrt

        if verb:
            print(f" (log_repetitions_operator_ra) log_num_repetitions_W {round(float(log_num_repetitions_W), 2)} "
                  f"= log(b) {round(math.log(b, 2))} - log_jensen_z {log_z} + log_treesize/2 {round(float(log_virt_treesize) / 2, 6)} + log_h/2 {math.log(h) / 2}")
        return log_num_repetitions_W

    """
        --------------------
                G-COST
        --------------------
    """
    def log_gcost_quantum_enumeration(self, k, h, log_y, log_z, circuit, verb=False):
        """
        GCost of enumerating subtree from level k of level h
        :param k: Level k
        :param h: Height of tree
        :param log_y: Log bound on number of combined nodes.
        :param log_z: Log Jensen's gap
        :param circuit: Instantiation for quantum operator W (cf. Section 5)
        :param verb:
        :return:
        """
        if h == 0:
            return 0

        # Number of combined nodes on level k: log_num_comb_nodes <= H_k_M
        log_num_comb_nodes = self.log_num_combined_nodes(k, log_y)

        # Number of "classical" subtrees
        log_Hk_M = self.treeH.log_Hk_M(k=k, n=self.n, nbases=self.nbases, pprime=self.pprime)
        log_num_subtree = log_Hk_M - log_num_comb_nodes

        # Classical binary tree: h_prime [Eq. 53, Our Work]
        h_prime = h * math.floor(Tools.doLog2(self.bound_coefficient))

        # D/F
        log_expected_dmw_in_fmv = math.log(max(math.ceil(self.num_nodes_branching_dfs * h_prime), 1), 2)

        # Q/D
        if self.force_DF != -1:
            log_expected_qw_in_dmv = math.log(max(math.ceil(self.const_reps_QPE * self.force_DF), 1), 2)
        else:
            log_expected_qw_in_dmv = math.log(max(math.ceil(self.const_reps_QPE * log_expected_dmw_in_fmv), 1), 2)

        # W/Q
        log_reps_W = self.log_reps_W_kh(k=k, h=h, log_num_comb_nodes=log_num_comb_nodes, log_z=log_z, verb=False)

        log_gcost_W = Tools.doLog2(circuit.get_operator_G_cost(k, h, log_num_its_qpe=log_reps_W))

        # Total
        log_gost_enumeration = log_num_subtree + log_expected_dmw_in_fmv + log_expected_qw_in_dmv + log_reps_W + log_gcost_W

        if verb:
            print(f" (log_gcost_quantum_enumeration) k={k}, log_num_combined_nodes={round(float(log_y), 2)}, log_jensen_z {log_z} -->"
                  f" log_gost_enumeration {round(float(log_gost_enumeration), 3)} "
                  f"= log_num_subtree {round(float(log_num_subtree), 3)} + log_expected_dmw_in_fmv {round(float(log_expected_dmw_in_fmv), 3)} "
                  f"+ log_reps_W {round(float(log_reps_W), 3)} + log_gcost_W {round(float(log_gcost_W), 3)}")

        return log_gost_enumeration

    def log_gcost_classical(self, k):
        """
        Expected cost of classical enumeration up to level k.
        :param k:
        :return:
        """
        if k == 0:
            return 0
        else:
            # Half tree because of symmetry + root
            return Tools.add_logs(self.treeH.log_avg_N_kh_M(k=0, h=k, n=self.n, nbases=self.nbases, pprime=self.pprime) - 1, 0)

    def log_enumeration_cost_classical(self):
        """
        Expected cost of classical enumeration of complete tree.
        :return:
        """
        # Half tree because of symmetry
        return Tools.add_logs(self.treeH.log_avg_N_kh_M(k=0, h=self.n+1, n=self.n, nbases=self.nbases, pprime=self.pprime) - 1, 0)


    def log_depth_total(self, k, h, log_y, log_z, circuit, verb=False):
        """
        Log depth of circuit for quantumly detecting markex vertex in subtree.
        :param k: Level k
        :param h:height of subtree
        :param log_y: Log bound on number of combined nodes on level k
        :param log_z: Jensen's z
        :param circuit: Instantiation for quantum operator W (cf. Section 5)
        :param verb:
        :return:
        """
        assert k >= 0, f"Level expected >= 0, was {k}! -- log_nbases {math.log(self.nbases, 2)}, log_pprime {math.log(self.pprime, 2)}"

        if h == 0:
            return 0

        log_num_comb_nodes = self.log_num_combined_nodes(k, log_y)

        # Log Depth QW
        log_operator_W_depth = Tools.doLog2(circuit.get_resource_count(k, h, gateid=circuit.idTdepth))
        log_repetitions_W = self.log_reps_W_kh(k, h, log_num_comb_nodes, log_z, verb=verb)

        log_depth_qc = log_repetitions_W + log_operator_W_depth

        if verb:
            print(f" (log_depth_total) log_depth_qc {round(float(log_depth_qc), 2)} "
                  f"= 2**(log_operator_W_depth {round(float(log_operator_W_depth), 2)} + log_repetitions_W {round(float(log_repetitions_W), 2)})")

        return log_depth_qc
