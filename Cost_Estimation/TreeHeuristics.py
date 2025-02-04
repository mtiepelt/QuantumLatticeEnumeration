#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Marcel Tiepelt"
__version__ = "1.0"
__status__ = "Research"

import Tools
import math

import sys
sys.path.append("../LatticeHeuristics")
import enumTrees

#@lru_cache(maxsize=None)
def log_Hk_M(k, n, nbases, pprime, bound, verb=False):
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
        return Tools.doLog2(1)

    if k == 1:
        # Only top tree
        return Tools.doLog2(nbases)
    else:
        # Top tree * Enum tree
        if bound == "LB":
            x = Tools.doNumber(
                enumTrees.CilinderPruningLowerBound.log_Hk(k=k-1, n=n, nbases=nbases, pprime=pprime)
            )
        elif bound == "UB":
            x = Tools.doNumber(
                enumTrees.CilinderPruningUpperBound.log_Hk(k=k-1, n=n, nbases=nbases, pprime=pprime)
            )
        else:
            raise ValueError(f"log_Hk_M bound can only be LB or UB, was {bound}")

        return x + Tools.doLog2(nbases)


#@lru_cache(maxsize=None)
def log_avg_Skh_M(k, h, n, nbases, pprime, bound, verb=False):
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

    if bound == "UBUB":
        log_Hk_M_kph = log_Hk_M(k + h, n, nbases, pprime, bound="UB", verb=False)
        log_Hk_M_k = log_Hk_M(k, n, nbases, pprime, bound="UB", verb=False)
    elif bound == "LBUB":
        log_Hk_M_kph = log_Hk_M(k + h, n, nbases, pprime, bound="LB", verb=False)
        log_Hk_M_k = log_Hk_M(k, n, nbases, pprime, bound="UB", verb=False)
    elif bound == "LBLB":
        log_Hk_M_kph = log_Hk_M(k + h, n, nbases, pprime, bound="LB", verb=False)
        log_Hk_M_k = log_Hk_M(k, n, nbases, pprime, bound="LB", verb=False)
    else:
        raise ValueError("lattice_bound can only be UBUB, LBLB, or LBUB")

    x = log_Hk_M_kph - log_Hk_M_k
    return x

#@lru_cache(maxsize=None)
def log_avg_N_kh_M(k, h, n, nbases, pprime, bound, verb=False):
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

    log_sum_x = log_avg_Skh_M(k, 1, n, nbases, pprime, bound, verb)
    for i in range(2, h+1):
        x = log_avg_Skh_M(k, i, n, nbases, pprime, bound, verb)
        log_sum_x = Tools.add_logs(x, log_sum_x)
    return log_sum_x
