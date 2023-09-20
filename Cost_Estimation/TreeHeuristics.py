#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Anonymous"
__version__ = "0.0"
__status__ = "Research"

import Tools
import math
from functools import lru_cache

import sys
sys.path.append("../LatticeHeuristics")
import enumTrees


class TreeHeuristics:
    """
        Interface for enumTrees.
    """
    def __init__(self, lattice_class=enumTrees.CilinderPruningLowerBound, lattice_bound='upper'):
        """
        :param lattice_class: Pruning strategy for lattice heuristics.
        :param lattice_bound: Bounds parameter on lattice heuristics.
        """
        self.bound = lattice_bound
        self.lattice_class = lattice_class

    #@lru_cache(maxsize=None)
    def log_Hk_M(self, k, n, nbases, pprime, verb=False):
        """
        Number of nodes on a level of the combined enumeration tree.
        :param k: Level k
        :param n: BKZ blocksize
        :param nbases: Number of combined, randomized bases
        :param pprime: Success probability for shortest vector to be in tree
        :param verb:
        :return:
        """

        assert k >= 0, f"Level: {k} cannot be smaller than 0!"
        assert k <= n + 2, f"Level: {k} cannot be larger than n+2: {n+2}"

        if k == 0:
            # Only root node
            if verb:
                print(f" . . (log_Hk)  root only, log_Hk_M: {0}")
            return Tools.doLog2(1)

        if k == 1:
            # Only top tree
            if verb:
                print(f" . . (log_Hk)  top tree only, log_Hk_M: {math.log(nbases, 2)}")
            return Tools.doLog2(nbases)
        else:
            # Top tree * Enum tree
            x = Tools.doNumber(
                self.lattice_class.log_Hk(k=k-1, n=n, nbases=nbases, pprime=pprime)
            )
            if verb:
                print(f" . . (log_Hk_M) enum + top tree, k {k}, n {n}, log_nbases {round(math.log(nbases, 2))}, bound {self.bound}, "
                      f"---> (on k={k}) log_Hk_M {round(float(x + Tools.doLog2(nbases)), 2)} = log_H (on k={k-1}) {round(float(x), 3)} + log(nbases) {math.log(nbases, 2)}")
            return x + Tools.doLog2(nbases)

    #@lru_cache(maxsize=None)
    def log_avg_Skh_M(self, k, h, n, nbases, pprime, verb=False):
        """
        Number of nodes on level h of a subtree rooted at leevel k of the combined enumeration tree.
        :param k: Level k
        :param h: Level of subtree
        :param n: BKZ blocksize
        :param nbases: Number of combined, randomized bases
        :param pprime: Success probability for shortest vector to be in tree
        :return:
        """
        assert k >= 0, f"Level {k} cannot be smaller than 0!  -- log_nbases {math.log(nbases, 2)}, log_pprime {math.log(pprime, 2)}"
        assert k <= n + 1, f"Level {k} cannot be larger than n + 2 {n + 2}!  -- log_nbases {math.log(nbases, 2)}, log_pprime {math.log(pprime, 2)}"
        assert h > 0, f"Depth {h} has to be at least 1!  -- log_nbases {math.log(nbases, 2)}, log_pprime {math.log(pprime, 2)}"

        x = self.log_Hk_M(k + h, n, nbases, pprime, verb=False) - self.log_Hk_M(k, n, nbases, pprime, verb=False)
        if verb:
            print(
                f" . . (log_avg_Skh_M) k {k}, h {h}, n {n}, log_n {round(math.log(nbases, 2))}, --> log_avg_Skh_M: {round(float(x), 3)} = {self.log_Hk_M(k + h, n, nbases, pprime, verb)} - {self.log_Hk_M(k, n, nbases, pprime, verb)}")
        return x

    #@lru_cache(maxsize=None)
    def log_avg_N_kh_M(self, k, h, n, nbases, pprime, verb=False):
        """
        Number of nodes in the subtree rooted at level k of height h of the combined enumeration tree without the root node on level k.
        :param k: Level k
        :param h: Level of subtree
        :param n: BKZ blocksize
        :param nbases: Number of combined, randomized bases
        :param pprime: Success probability for shortest vector to be in tree
        :param verb:
        :return:
        """
        assert k >= 0, f"Level {k} cannot be smaller than 0!"
        assert k <= n + 1, f"Level {k} cannot be larger than n+2: {n + 2}"
        assert h > 0, f"Depth {h} has to be at least 1!"
        assert k + h <= n + 1, f"Level {k} + D {h} cannot be larger than n+2: {n + 2}"

        log_sum_x = self.log_avg_Skh_M(k, 1, n, nbases, pprime, verb)
        for i in range(2, h+1):
            x = self.log_avg_Skh_M(k, i, n, nbases, pprime, verb)
            log_sum_x = Tools.add_logs(x, log_sum_x)

        if verb:
            print(f" . . (log_avg_N_kh_M) k: {k}, h: {h}, n: {n}, log_nbases: {round(math.log(nbases, 2))}, bound: {self.bound} "
                       f"--> log_avg_N_kh_M {round(float(log_sum_x), 15)}")

        return log_sum_x
