#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Marcel Tiepelt"
__version__ = "1.0"
__status__ = "Research"

import Tools
import math
import TreeHeuristics as TH

class CostFunctions:
    """
        Cost function for classical and quantum enumeration cost according to Section 3,4 and 5
    """
    def __init__(self, n, nbases, pprime, constants, bound):
        """

        :param n: BKZ blocksize.
        :param nbases: Number of combined random bases for extreme pruning
        :param pprime: Success probability of shortest vector to be in enumeration tree
        :param constants: Dict of constants 'boundCoefficient', 'num_nodes_branching_dfs', 'force_DF', 'const_reps_QPE', 'const_reps_W' determining bound on DetectMV and FindMV
        :param bound: Bounds on estimation of subtree size, 'LBLB', 'UBLB', 'UBUB'
        """

        self.n = n
        self.bound_coefficient = constants['boundCoefficient']

        self.num_nodes_branching_dfs = constants['num_nodes_branching_dfs']
        self.force_DF = constants['force_DF']

        self.const_reps_QPE = constants['const_reps_QPE']
        self.const_reps_W = constants['const_reps_W']

        self.nbases = nbases
        self.pprime = pprime

        self.bound = bound

    def log_num_combined_nodes(self, k, log_y):
        """
        Return number of combined nodes in virtual root at level k.
        :param k: Level k
        :param log_y: Log bound on number of combined nodes.
        :return: log_num_combined_nodes
        """
        if k == 0:
            return 0

        if self.bound == 'LBUB' or self.bound == 'LBLB':
            log_Hk_M = TH.log_Hk_M(k=k, n=self.n, nbases=self.nbases, pprime=self.pprime, bound='LB')
        else:  # 'UBUB'
            log_Hk_M = TH.log_Hk_M(k=k, n=self.n, nbases=self.nbases, pprime=self.pprime, bound='UB')

        return min(log_Hk_M, log_y)

    def log_reps_W_kh(self, k, h, log_num_comb_nodes, log_z, verb=False):
        """
        Number of applications of quantum operator W in QPE (cf. Section 3.1)
        :param k: Level k
        :param h: Height of subtree
        :param log_num_comb_nodes: Number of combined nodes on level k
        :param log_z: Log Jensen's z
        :param verb:
        :return: log_reps_W_kh
        """
        b = self.const_reps_W

        # Treesize: Subtree + subtree-root

        # Nodes of subtree +subtree root; -1 from expected value of Conjecture 3
        log_single_subtreesize = Tools.add_logs(TH.log_avg_N_kh_M(k=k, h=h, n=self.n, nbases=self.nbases, pprime=self.pprime, bound=self.bound) - 1, 0)

        # Virtual tree without virtual root
        log_virtual_treesize_no_root = log_num_comb_nodes + log_single_subtreesize

        # Virtual tree + root
        log_virtual_treesize = Tools.add_logs(0, log_virtual_treesize_no_root)

        # W/Q: b * 1/2^z * sqrt(treesize * h)
        log_sqrt = (log_virtual_treesize + Tools.doLog2(h)) / 2
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
        :return: log_gost_enumeration
        """
        if h == 0:
            return 0

        # Number of combined nodes on level k: log_num_comb_nodes <= H_k_M
        log_num_comb_nodes = self.log_num_combined_nodes(k, log_y)

        # Number of "classical" subtrees
        if self.bound == 'LBUB' or self.bound == 'LBLB':
            log_Hk_M = TH.log_Hk_M(k=k, n=self.n, nbases=self.nbases, pprime=self.pprime, bound='LB')
        else:  # 'UBUB'
            log_Hk_M = TH.log_Hk_M(k=k, n=self.n, nbases=self.nbases, pprime=self.pprime, bound='UB')

        # - 1 for symmetry of the tree. This are classical repetitions of virtual subtrees.
        log_num_subtree = (log_Hk_M - log_num_comb_nodes) - 1

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

        log_gcost_W = Tools.doLog2(circuit.g_cost(k, h, log_num_its_qpe=log_reps_W))

        # Total
        log_gost_enumeration = log_num_subtree + log_expected_dmw_in_fmv + log_expected_qw_in_dmv + log_reps_W + log_gcost_W

        if verb:
            print(f" (log_gcost_quantum_enumeration) k={k}, log_num_combined_nodes={round(float(log_y), 2)}, log_jensen_z {log_z} -->"
                  f" log_gost_enumeration {round(float(log_gost_enumeration), 3)} = \n"
                  f"   " + f"log_num_subtree {round(float(log_num_subtree), 3)} \n"
                  f"  +" + f"log_expected_dmw_in_fmv {round(float(log_expected_dmw_in_fmv), 3)} \n"
                  f"  +" + f"log_reps_W {round(float(log_reps_W), 3)} \n"
                  f"  +" + f"log_gcost_W {round(float(log_gcost_W), 3)}")

        return log_gost_enumeration

    def log_gcost_classical(self, k):
        """
        Expected cost of classical enumeration up to level k.
        :param k: Level to enumerate to.
        :return: log of classical + quantum GCost
        """
        if k == 0:
            return 0
        else:
            # Half tree because of symmetry + root
            return Tools.add_logs(TH.log_avg_N_kh_M(k=0, h=k, n=self.n, nbases=self.nbases, pprime=self.pprime, bound=self.bound) - 1, 0)

    def log_enumeration_cost_classical(self):
        """
        Expected cost of classical enumeration of complete tree.
        :return: log of classical GCost
        """
        return self.log_gcost_classical(k=self.n+1)


    def log_depth_total(self, k, h, log_y, log_z, circuit, verb=False):
        """
        Log depth of circuit for quantumly detecting markex vertex in subtree.
        :param k: Level k
        :param h:height of subtree
        :param log_y: Log bound on number of combined nodes on level k
        :param log_z: Jensen's z
        :param circuit: Instantiation for quantum operator W (cf. Section 5)
        :param verb:
        :return: log of TDepth(quantum)
        """
        assert k >= 0, f"Level expected >= 0, was {k}! -- log_nbases {math.log(self.nbases, 2)}, log_pprime {math.log(self.pprime, 2)}"

        if h == 0:
            return 0

        log_num_comb_nodes = self.log_num_combined_nodes(k, log_y)

        # Log Depth QW
        log_operator_W_depth = Tools.doLog2(circuit.t_depth(k, h, gateid=circuit.idTdepth))
        log_repetitions_W = self.log_reps_W_kh(k, h, log_num_comb_nodes, log_z, verb=verb)
        log_depth_qc = log_repetitions_W + log_operator_W_depth

        return log_depth_qc
