#!/usr/bi /python3
# -*- coding: utf-8 -*-

__author__ = "Anonymous"
__version__ = "0.0"
__status__ = "Research"

import math

class Circuit:
    """
    Abstract class for circuit instantiations of quantum operator W (cf. Section 4)
    """
    def __init__(self, desc):
        self.desc = desc
        pass

    def get_resource_count(self, k, h, **kwargs):
        pass
    def get_operator_gate_count(self, k, h, **kwargs):
        pass
    def get_operator_G_cost(self, k, h, **kwargs):
        pass

class QueryCircuit(Circuit):
    """
        'Black-box'/ Query class for circuit instantiations of quantum operator W (cf. Section 4.1).
    """
    idT = 'gate'
    idTdepth = 'gate'
    def __init__(self, desc):
        super().__init__(desc)

    def get_resource_count(self, k, h, **kwargs):
        """
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """
        if h == 0:
            return 0
        return 1

    def get_operator_gate_count(self, k, h, **kwargs):
        """
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """
        gate_costs = {}
        if h == 0:
            gate_costs['gate'] = 0
        else:
            assert h > 1, f"No quantum circuits for single layer of tree. Expected: h > 1, was h = {h}."
            gate_costs['gate'] = 1
        return gate_costs

    def get_operator_G_cost(self, k, h, **kwargs):
        """
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """
        if h == 0:
            return 0
        assert h > 1, f"No quantum circuits for single layer of tree. Expected: h > 1, was h = {h}."
        return self.get_operator_gate_count(k, h)[self.idT]


class MinimalCircuit(Circuit):
    """
        Minimal circuit instantiations of quantum operator W (cf. Section 4.2).
    """
    idT = 'gate'
    idTdepth = 'gate'

    def __init__(self, desc, q=3329):
        super().__init__(desc)
        self.q = q

    def get_resource_count(self, k, h, **kwargs):
        """
        Depth of Minimal Circuit according to Section 4.2, Table 1, Table 2
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """

        if h == 0:
            return 0
        assert h > 1, f"No quantum circuits for single layer of tree. Expected: h > 1, was h = {h}."

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

    def get_operator_gate_count(self, k, h, **kwargs):
        """
        Depth of Minimal Circuit for each possible gate.
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """
        gate_costs = {}
        if h == 0:
            gate_costs[self.idT] = 0
        else:
            assert h > 1, f"No quantum circuits for single layer of tree. Expected: h > 1, was h = {h}."
            gate_costs[self.idT] = self.get_resource_count(k, h)

        return gate_costs

    def get_operator_G_cost(self, k, h, **kwargs):
        """
        GCost of Minimal Circuit according to Section 4.2, Table 1, Table 2
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """
        """
            G-cost: Sum of Clifford, CNOT, T, R, M
        """
        if h == 0:
            return 0
        assert h > 1, f"No quantum circuits for single layer of tree. Expected: h > 1, was h = {h}."

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



class QSharpCircuit(Circuit):
    """
        QSharp circuit instantiations of quantum operator W (cf. -- NOT IN PAPER --).
    """
    idWidth = 'initial width'
    idxWidth = 'extra width'

    idTdepth = 'T depth'
    idCnot = 'CNOT count'
    idT = 'T count'
    idCliffs = '1-qubit Clifford count'

    idM = 'M count'
    idR = 'R count'

    idfactor = 'comment'
    idOperation = 'operation'

    def __init__(self, beta, max_bound_below, fn_dir, desc, q=3329):
        # ORDER IS IMPORTANT!
        super().__init__(desc)
        self.header = [QSharpCircuit.idCnot, QSharpCircuit.idCliffs, QSharpCircuit.idT,
                       QSharpCircuit.idR, QSharpCircuit.idM, QSharpCircuit.idTdepth,
                       QSharpCircuit.idWidth, QSharpCircuit.idxWidth]

        # [blocksize] --> [op] --> [gate] --> number
        self.blocksize_op_gate_cost = {}
        self.parse_all_qsharp_files_factor(q, beta, fn_dir=fn_dir)

        # [blocksize] --> [maximal bound on node degree below]
        self.max_bound_below = max_bound_below
        self.q = q

    def parse_all_qsharp_files_factor(self, q, beta, fn_dir, fn_prefix='qsharp'):
        """
            Reads QSharp Resources estiamation files to dictionary.
            Write to self.blocksize_op_gate_cost: blocksize_op_gate_cost[block_size] --> [gate] --> number
        :param q:   LWE modulus.
        :param beta: BKZ blocksize.
        :param fn_dir:  Target directory of qsharp estimate files.
        :param fn_prefix: Prefix for qsharp estimate files.
        :return:
        """
        rank_not_found = []
        beta_top = beta + 1

        for b in range(2, beta_top + 1):
            filename = f"{fn_dir}/{fn_prefix}D_{q}_{b}.Y"
            try:
                self.blocksize_op_gate_cost[b] = self.parse_qsharp_file_factor(filename)
            except FileNotFoundError:
                print(f" (!) File {filename} expected, but not found!")
                rank_not_found.append(b)

        if len(self.blocksize_op_gate_cost.keys()) + 1 < beta_top:
            print(" (!-estc) Largest rank for Q# estimation is " + str(
                len(self.blocksize_op_gate_cost.keys())) + f", expected {beta_top}.")
            if rank_not_found:
                print(f" (!-estc) Ranks {rank_not_found} for Q# estimation not found.")

        # def print_op_cost(b_g_c):
        #     print("-" * 20)
        #     for b,g in b_g_c.items():
        #         print(f"beta {b}")
        #         print(f"   {', '.join([str(k) + ' ' + str(v) for k,v in g.items()])}")
        #         print("-" * 5)

        # print_op_cost(self.blocksize_op_gate_cost)

    def parse_qsharp_file_factor(self, filename):
        """
        Parse QSharp resource estimator files.
                op_gate_counts[operation] --> [gateid] --> number
        :param filename:
        :return: Dictionary with operations, gateid and number of gates.
        """
        file = open(filename, "r")
        op_gate_counts = {}

        for line in file.readlines():
            if not line.strip():
                continue

            values = [x.strip() for x in line.strip().split(',') if x.strip()]
            # assert len(values) == (len(self.header) + 2), "q-sharp outputs have wrong format, length " + str(len(values)) + ", but expected " + str(len(header))

            # Header line
            if values[0] == QSharpCircuit.idOperation:
                continue

            op_gate_counts[values[0]] = dict.fromkeys(self.header, 0)

            # Counts != Width
            for i in range(1, 9):
                op_gate_counts[values[0]][self.header[i - 1]] = float(values[i])  # first entry operation name is ignored

        return op_gate_counts

    def scaling_QW_gates_dict(self, block_size, number_qaa_combinations, number_qaa_iterations, root_level):
        """
            Assume Gate Depth Optimization/ Parallel execution of gates on independent states
                - do_setup, * 2 for uncompute
                - do_proj_mult, * 4 for uncompute
                                * number_qaa_iterations for QAA
                                * number_qaa_combinations for QAA
                                * block_size**2 to account for multiplication of complete vector
                                + 1 to account for projection on first level of operator
                - do_proj_add,  * 2 for uncompute
                                * number_qaa_iterations for QAA
                                * number_qaa_combinations for QAA
                                * block_size**2 to account for multiplication of complete vector
                                + 1 to account for projection on first level of operator
                - do_norm,      * 2 for uncompute
                                * number_qaa_iterations for QAA
                                * number_qaa_combinations for QAA
                                * block_size to account for multiplication of complete vector
                                + 1 to account for projection on first level of operator
                - do_norm_misc, * 2 for uncompute
                                * number_qaa_iterations for QAA
                                * number_qaa_combinations for QAA
                                + 1 to account for projection on first level of operator
                - qaa_misc, * 2 for uncompute
                                * number_qaa_iterations for QAA
                - do_state_init,    * 2 for uncompute
                                    * number_qaa_combinations
                                    * number_qaa_iterations
                                    +
                                    * 2 for uncompute
                                    * 2**y  for preparation of first level
                                    * h for preparation of state spanning previous coefficients
                - do_qaa_misc,  * 2 for uncompute
                                * number_qaa_combinations
                                * block_size for every coefficient

        :param block_size: BKZ blocksize
        :param number_qaa_combinations: Number of QAA combined to reduce failure probability.
        :param number_qaa_iterations:  Number of iterations within QAA.
        :param root_level: Level k.
        :return:
        """
        number_qaa_combinations = math.ceil(number_qaa_combinations)
        number_qaa_iterations = math.ceil(number_qaa_iterations)
        qsharp_op_ra_factor = {
            'QOpEstimates.Enumeration.do_setup': 2,
            'QOpEstimates.Enumeration.do_proj_mult': 2 * (number_qaa_iterations + 1) * number_qaa_combinations * block_size**2,
            'QOpEstimates.Enumeration.do_proj_add': 2 * (number_qaa_iterations + 1) * number_qaa_combinations * block_size,
            'QOpEstimates.Enumeration.do_proj_misc': 2 * (number_qaa_iterations + 1),
            'QOpEstimates.Enumeration.do_norm': 2 * (number_qaa_iterations + 1) * number_qaa_combinations * block_size,
            'QOpEstimates.Enumeration.do_norm_misc': 2 * (number_qaa_iterations + 1) * number_qaa_combinations,
            'QOpEstimates.Enumeration.do_qaa_misc': 2 * number_qaa_iterations * number_qaa_combinations,
            'QOpEstimates.Enumeration.do_qaa_state_init': 2 * number_qaa_iterations * number_qaa_combinations * self.max_bound_below[
                root_level],
            'QOpEstimates.Enumeration.doL6': 1,
        }
        return qsharp_op_ra_factor

    def get_resource_count(self, k, h, **kwargs):
        """
        Depth of QSharp Circuit according to Q# Resource Estimation
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param kwargs:
        :return:
        """
        if h == 0:
            return 0
        assert h > 1, f"No quantum circuits for single layer of tree. Expected: h > 1, was h = {h}."

        verb = False
        if 'verb' in kwargs.keys():
            verb = kwargs['verb']

        gateid = self.idTdepth
        if 'gateid' in kwargs.keys():
            gateid = kwargs['gateid']

        block_size = h
        root_level = k
        assert block_size > 0, f"block_size {block_size} has to be > 0!"
        # Dummy partition. QC exists only for depth at least 2, but partition may have single layers. Upper bound with half cost for 2.
        cost_factor = 1
        if block_size == 1:
            block_size = 2
            cost_factor = 0.5
        # If root node bound is maximum
        if root_level == 0:
            root_level = 1

        def num_qaa_repetitions(self, k):
            """
            Number of iterations within QAA (cf. -- NOT IN PAPER --)
            :param self:
            :param k: Level of enumeration tree.
            :return:
            """
            C = self.q
            D = self.max_bound_below[k]
            if D == 0:
                return 0
            r_nominator = C - (C / D) + math.sqrt(C) - math.sqrt(C / D)
            r_denominator = 2 * (math.sqrt(C) - math.sqrt(C / D))
            r = r_nominator / r_denominator
            return r

        def num_qaa_combinations(self, k, log_num_repetitions_W, verb=False):
            """
            Number of QAA combination (cf. -- NOT IN PAPER --)
            :param self:
            :param k: Level of enumeration tree.
            :param log_num_repetitions_W: Upper bound on number of application of operator W within QPE
            :param verb:
            :return:
            """
            if self.max_bound_below[k] == 0:
                return 0

            if log_num_repetitions_W == 0:
                return 0

            C = self.q
            D = self.max_bound_below[k]

            first_log = math.log(0.01 / 2 ** log_num_repetitions_W, 2)
            second_log = math.log(1 / 2 + 1 / (4.84 * (1 - math.sqrt(1 / C) * (math.sqrt(1 - D / C)))), 2)
            num_qaa_combinations = first_log / second_log
            if verb:
                print(f" . . (get_num_qaa_combinations) num_qaa_combinations {round(num_qaa_combinations, 2)}")

            return num_qaa_combinations

        number_qaa_iterations = self.num_qaa_repetitions(root_level) + math.log(cost_factor, 2)

        number_qaa_combinations = 1
        if 'number_qaa_combinations' in kwargs.keys():
            number_qaa_combinations = kwargs['number_qaa_combinations']

        qsharp_op_ra_factor = self.scaling_QW_gates_dict(block_size, number_qaa_combinations=number_qaa_combinations, number_qaa_iterations=number_qaa_iterations, root_level=root_level)
        gate_count = 0

        for opid, v in qsharp_op_ra_factor.items():
            gate_count += (self.blocksize_op_gate_cost[block_size][opid][gateid] * v)

        return gate_count * cost_factor

    def get_circuit_count(self, k, h, gateid, log_num_its_qpe, **kwargs): #block_size, gate_id, root_level, log_N_kD_M, verb=False):
        """
        Depth of QSharp Circuit according to Q# Resource Estimation for gateid
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :param gateid: Gate id
        :param log_num_its_qpe:  Upper bound on number of application of operator W within QPE
        :param kwargs:
        :return:
        """
        verb = False
        if 'verb' in kwargs.keys():
            verb = kwargs['verb']
        if h == 0:
            return 0
        assert h > 1, f"No quantum circuits for single layer of tree. Expected: h > 1, was h = {h}."

        number_qaa_combinations = self.num_qaa_combinations(k, log_num_its_qpe, verb=verb)

        # Only affects width and gates
        if gateid == self.idTdepth:
            number_qaa_combinations = 1

        gate_count = self.get_resource_count(k, h, gateid=gateid, number_qaa_combinations=number_qaa_combinations)

        if verb:
            print(f" . . . (get_circuit_depth) {self.desc} number_qaa_combinations {number_qaa_combinations}, gates ~ {gate_count[self.idT]}")

        return gate_count


    def get_operator_gate_count(self, k, h, **kwargs):
        """
        Sum of gate counts of QSharp Circuit according to Q# Resource Estimation for all gateids
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :return:
        """
        log_num_its_qpe = 0

        if 'log_num_its_qpe' in kwargs.keys():
            log_num_its_qpe = kwargs['log_num_its_qpe']
        if 'verb' in kwargs.keys():
            verb = kwargs['verb']

        # Combines QAA's
        number_qaa_combinations = self.num_qaa_combinations(k, log_num_its_qpe, verb=verb)

        # Only affects width and gates
        if gateid == self.idTdepth:
            number_qaa_combinations = 1

        gate_costs = {}
        for gate in self.header[:6]:
            if h == 0:
                gate_costs[gate] = 0
            else:
                gate_costs[gate] = self.get_resource_count(k, h, gateid=gateid, number_qaa_combinations=number_qaa_combinations)

        return gate_costs

    def get_operator_G_cost(self, k, h, **kwargs):
        """
        GCost of QSharp Circuit according to Q# Resource Estimation
        :param k: Level of enumeration tree.
        :param h: Height of subtree.
        :return:
        """
        if h == 0:
            return 0
        assert h > 1, f"No quantum circuits for single layer of tree. Expected: h > 1, was h = {h}."
        log_num_its_qpe = kwargs['log_num_its_qpe']
        gate_costs = self.get_operator_gate_count(k, h, log_num_its_qpe=log_num_its_qpe)
        g_cost = 0
        for gate in self.header[:5]:
            g_cost += gate_costs[gate]
        return g_cost
