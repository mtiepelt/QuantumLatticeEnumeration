#!/usr/bi /python3
# -*- coding: utf-8 -*-

__author__ = "Marcel Tiepelt"
__version__ = "1.0"
__status__ = "Research"

import math

class Circuit:
    """
    Abstract class for circuit instantiations of quantum operator W (cf. Section 4)
    """
    def __init__(self, desc):
        self.desc = desc
        self.gate = 'gate'
        pass

    def t_depth(self, k, h, **kwargs):
        pass
    def individual_gate_count(self, k, h, **kwargs):
        pass
    def g_cost(self, k, h, **kwargs):
        pass

class QueryCircuit(Circuit):
    """
        'Black-box'/ Query class for circuit instantiations of quantum operator W (cf. Section 4.1).
    """
    idT = 'gate'
    idTdepth = 'gate'
    def __init__(self, desc):
        super().__init__(desc)

    def individual_gate_count(self, k, h, **kwargs):
        """
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """
        assert h > 0, f"No circuits for tree of height 0."
        return {self.gate : 1}

    def t_depth(self, k, h, **kwargs):
        """
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """
        assert h > 0, f"No circuits for tree of height 0."
        return 1

    def g_cost(self, k, h, **kwargs):
        """
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """
        assert h > 0, f"No circuits for tree of height 0."
        return 1

class MinimalCircuit(Circuit):
    """
        Minimal circuit instantiations of quantum operator W (cf. Section 4.2).
    """
    idT = 'gate'
    idTdepth = 'gate'

    def __init__(self, desc, q=3329):
        super().__init__(desc)
        self.q = q

    def t_depth(self, k, h, **kwargs):
        """
        Depth of Minimal Circuit according to Section 4.2, Table 1, Table 2
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """
        assert h > 0, f"No circuits for tree of height 0."

        b = math.ceil(math.log(self.q, 2))

        # Double precision for floating point operations according to
        # Hermans, J., Schneider, M., Buchmann, J., Vercauteren, F., Preneel, B.: Parallel shortest lattice vector enumeration on graphics cards.
        bound_prec = 53

        def depth_mult(bitsize):
            return math.log(bitsize, 2)**2

        def depth_add(bitsize):
            return math.log(bitsize, 2)

        def depth_cmp(bitsize):
            return math.log(bitsize, 2)

        x_0 = min(bound_prec, b)
        mult_project = depth_mult(x_0)

        x_1 = bound_prec + b
        add_norm = depth_add(x_1) * math.log(h, 2)

        x_2 = x_1 + math.log(h, 2)
        square_norm = depth_mult(x_2)

        x_3 = 2 * x_2
        add_norm_sum = depth_add(x_3) * math.log(h, 2)

        x_4 = x_3 + math.log(h, 2)
        cmp_norm = depth_cmp(x_4)

        return mult_project + add_norm + square_norm + add_norm_sum + cmp_norm

    def individual_gate_count(self, k, h, **kwargs):
        """
        Count of Minimal Circuit for each possible gate is lower bounded by depth.
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """
        assert h > 0, f"No circuits for tree of height 0."
        return {self.gate : self.t_depth(k, h)}

    def g_cost(self, k, h, **kwargs):
        """
        GCost of Minimal Circuit according to Section 4.2, Table 1, Table 2
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """
        assert h > 0, f"No circuits for tree of height 0."

        b = math.ceil(math.log(self.q, 2))

        # Double precision for floating point operations according to
        # Hermans, J., Schneider, M., Buchmann, J., Vercauteren, F., Preneel, B.: Parallel shortest lattice vector enumeration on graphics cards.
        bound_prec = 53

        def gcost_mult(bitsize):
            return bitsize * math.log(bitsize, 2) * math.log(math.log(bitsize,2), 2)

        def gcost_add(bitsize):
            return bitsize

        def gcost_cmp(bitsize):
            return bitsize

        x_0 = min(bound_prec, b)
        mult_project = gcost_mult(x_0) * h**2

        x_1 = bound_prec + b
        add_norm = gcost_add(x_1) * h**2

        x_2 = x_1 + math.log(h, 2)
        square_norm = gcost_mult(x_2) * h

        x_3 = 2 * x_2
        add_norm_sum = gcost_add(x_3) * h

        x_4 = x_3 + math.log(h, 2)
        cmp_norm = gcost_cmp(x_4)

        return mult_project + add_norm + square_norm + add_norm_sum + cmp_norm