#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Marcel Tiepelt"
__version__ = "1.0"
__status__ = "Research"

from sage.all import RealField
RRR = RealField(200)

from functools import lru_cache

import math
import numpy as np
import Circuits
import CostFunctions as CF

import os
import json

from dataclasses import dataclass, asdict
from collections import defaultdict

# Load expected security values from configuration file
from configSecurity import aes_expected_security, kyber_expected_security


def get_intersection_x(x_coords, y_total, y_critical):
    """
    Find smallest x_coord where y_totalcost <= y_crititcal
    :param x_coords: X coordinates
    :param y_total: Y coordinates
    :param y_critical: Coordinates to check for intersection with with Y coordinates
    :return: x_coords of intersection, {'>>', '<=', '<'}
    """

    # result: 1 <==> total <= critical
    # argwhere: array of Indices of non-zero elements
    # flatten: gives indices in list => last is position of first element s.t. total <= larger
    result = []
    for i, j in zip(y_total, y_critical):
        if i <= j:
            result.append(1)
        else:
            result.append(0)
    all_idz = np.argwhere(result).flatten()


    if len(all_idz) > 0:
        # Intersection exists: First indice where total <= critical
        idz = all_idz[0]

        if y_total[idz] == y_critical[idz]:
            # Exact intersection, for z >= idz, all y_total[z] <= y_critial
            return x_coords[idz], '<='
        elif y_total[idz] < y_critical[idz]:
            # Fractional intersection, for z >= idz, y_total[z] < y_critial
            return x_coords[idz], '<'
        # This cannot happen
        else:
            assert False, f"y_total[idz]={y_total[idz]} < y_critical[idz]={y_critical[idz]}, numpy error!"
    else:
        # All y_total[x] are above y_critial, for any z, y_total[z] > y_critical
        return x_coords[-1], '>>'


def get_intersections(n, max_depth, log_M, z_to_cost, constants, bound):
    """
    Find intersection for
        1. Expected cost of classical enumeration
		2. Quasi-Sqrt(classical cost)
		3. Target security of Kyber
		4. Expected cost of Grover on AES
    :param n: BKZ blocksize
    :param max_depth: MaxDepth constraint for quantum circuit as used in cost estimation
    :param log_M: Log bound on number of combindes randomized enumeration bases
    :param z_to_cost: Dictionary [Jensen's z] --> cost
    :return: Dictionary of all intersections
    """
    x_coords = sorted(list(z_to_cost.keys()))
    y_total_cost = [z_to_cost[z].log_gcost_total() for z in x_coords]

    intersections = {
        'aes_grover' : [-1,'none'],
        'target': [-1, 'none'],
        'quasi-quadratic': [-1, 'none']
    }

    # Intersection with Expected Cost for AES
    z_value, total_INEQ_crit = get_intersection_x(x_coords, y_total_cost, [aes_expected_security[n][max_depth]] * len(x_coords))
    intersections['aes_grover'] = [z_value, total_INEQ_crit, y_total_cost[x_coords.index(z_value)],]

    # Intersection with Expected Cost for Security Target
    z_value, total_INEQ_crit = get_intersection_x(x_coords, y_total_cost, [kyber_expected_security[n]] * len(x_coords))
    intersections['target'] = [z_value, total_INEQ_crit, y_total_cost[x_coords.index(z_value)], ]

    # Intersection with Quasi-Quadratic Quantum Speedup
    costF =CF.CostFunctions(n=n, nbases=2**log_M, pprime=1, constants=constants, bound=bound)
    log_cost_classical = costF.log_enumeration_cost_classical()

    z_value, total_INEQ_crit = get_intersection_x(x_coords, y_total_cost, [(log_cost_classical + math.log(n, 2)) / 2] * len(x_coords))
    intersections['quasi-quadratic'] = [z_value, total_INEQ_crit, y_total_cost[x_coords.index(z_value)]]
    return intersections

@dataclass
class Cost:
    """
    Class for keeping track of cost of quantum enumeration for a given set of parameters.
    """
    #
    name: str = 'NONE'
    #
    k: int = -1
    #
    log_m: int = -1
    log_y: int = -1
    log_z: int = -1
    #
    log_gcost_classical: int = 99999999
    log_gcost_quantum: int = 99999999
    #
    #
    qracm: int = -1
    log_t_depth: int = -1

    def overwrite(self, otherCost):
        self.name = otherCost.name
        self.set_parameters(otherCost.log_m, otherCost.log_y, otherCost.log_z, otherCost.k)
        self.set_gcost(otherCost.log_gcost_classical, otherCost.log_gcost_quantum)
        self.set_misc(otherCost.qracm, otherCost.log_t_depth)

    def log_gcost_total(self):
        return add_logs(self.log_gcost_classical, self.log_gcost_quantum)

    def set_parameters(self, log_m, log_y, log_z, k):
        self.k = k
        self.log_m = log_m
        self.log_y = log_y
        self.log_z = log_z

    def set_gcost(self, log_gcost_classical, log_gcost_quantum):
        self.log_gcost_classical = log_gcost_classical
        self.log_gcost_quantum = log_gcost_quantum

    def set_misc(self, qracm, log_t_depth):
        self.log_t_depth = log_t_depth
        self.qracm = qracm


class CostFiles:
    """
        Output exact costs for every Jensen's z to file.
    """
    def __init__(self, q, constants, bound):
        """

        :param q: LWE Modulus.
        :param constants: Dict of constants 'boundCoefficient', 'num_nodes_branching_dfs', 'force_DF', 'const_reps_QPE', 'const_reps_W' determining bound on DetectMV and FindMV
        :param bound: Bounds on estimation of subtree size, 'LBLB', 'UBLB', 'UBUB'
        """

        self.bound_coef = constants['boundCoefficient']
        self.modulus_q = q
        self.const_reps_W = constants['const_reps_W']
        self.const_reps_QPE = constants['const_reps_QPE']
        self.num_nodes_branching_dfs = constants['num_nodes_branching_dfs']
        self.force_DF = constants['force_DF']
        self.bound = bound

    def get_parameter_str(self):
        """
        Unique identifier for instance based on constant and modulus.
        :return: unique string
        """
        return f"C{self.bound_coef}_Q{self.modulus_q}_b{self.const_reps_W}_e{self.const_reps_QPE}_branch{self.num_nodes_branching_dfs}_DF{self.force_DF}_b{self.bound}"

class DataFiles(CostFiles):
    """
        Data cache for cost estimation
    """
    def __init__(self, n, q, max_depth, constants, bound):
        self.n = n
        self.max_depth = max_depth
        super().__init__(q, constants, bound)

    def __get_dir(self, circuit, cache_dir='costData'):
        """
        Get shared directory for all cache files.
        :param circuit: Circuit instantiation object as used in cost esteiamtion.
        :param cache_dir: Target directory.
        :return: string
        """
        return cache_dir + '/' + str(self.n) + '/' + str(circuit.desc)

    def __get_filename(self, log_M, log_Y, step_size, log_z, circuit):
        """
        Get cache data filename
        :param log_M: Log bound on number of combined bases
        :param log_Y: Log bound on number of combined nodes on level k
        :param step_size: Stepsize for M
        :param log_z: Log bound on Jensen's z
        :param circuit: Circuit instantiation.
        :return: string
        """
        return f"{self.n}_{self.max_depth}_{circuit.desc}_M{log_M}_Y{log_Y}_ss{step_size}_Z{log_z}_{super(DataFiles, self).get_parameter_str()}.data"

    def write_data_cost_to_file(self, z_to_cost_new, circuit, log_M, log_Y, step_size):
        """
        Write cost for all Jensen's gap z values to file.
        :param z_to_cost_new: Dictionary [Jensen z] --> cost
        :param circuit: Circuit instantiation.
        :param log_M: Log bound on number of combined bases
        :param log_Y: Log bound on number of combined nodes on level k
        :param step_size: Stepsize for M
        """
        for key, value in z_to_cost_new.items():
            self.write_data_cost_file(value, circuit, log_M, log_Y, step_size, key)

    def write_data_cost_file(self, m_k_to_cost, circuit, log_M, log_Y, step_size, log_z):
        """
        Write cost data for one Jensen's gap z of a unique instance to a file.
        :param m_k_to_cost: Dictionary [#combined bases] --> [Level k] --> cost
        :param circuit: Circuit instantiation.
        :param log_M: Log bound on number of combined bases
        :param log_Y: Log bound on number of combined nodes on level k
        :param step_size: Stepsize for M
        :param log_z: Log Jensen's z
        """
        full_dir = self.__get_dir(circuit)
        filename = self.__get_filename(log_M, log_Y, step_size, log_z, circuit)

        os.system('mkdir -p ' + full_dir)
        os.system('touch ' + full_dir + '/PLACEHOLDER')
        file = open(f"{full_dir}/{filename}", 'w')

        file.write('# Comments\n')
        file.write(f'\n# data\n')

        for m, k_to_cost in m_k_to_cost.items():
            for k, cost in k_to_cost.items():
                full_data = asdict(cost)
                data = json.dumps(full_data)
                file.write(data + '\n\n')

    def read_cached_data_from_file(self, circuit, log_M, log_Y, step_size, log_Z, step_size_z):
        """
        Read all cached data from files
        :param circuit: Circuit instantiation.
        :param log_M: Log bound on number of combined bases
        :param log_Y: Log bound on number of combined nodes on level k
        :param step_size: Stepsize for M
        :param log_z: Log Jensen's z
        :param step_size_z: Stepsize for Jensen's z
        :return: [z] --> [m] --> [k] --> Cost()
        """
        full_dir = self.__get_dir(circuit)
        z_m_k_cost = {}

        for log_z in range(0, log_Z+1, step_size_z):
            filename = self.__get_filename(log_M, log_Y, step_size, log_z, circuit)
            exists = os.path.isfile(f"{full_dir}/{filename}")

            if exists:
                z_m_k_cost[log_z] = self.__read_data_cost_from_file(full_dir, filename, log_M, log_Y, step_size, read_step_size=0, log_M_lower=0)
            else:
                for ss in range(1, step_size):
                    if step_size % ss == 0:
                        filename = self.__get_filename(log_M, log_Y, ss, log_z, circuit)
                        exists = os.path.isfile(f"{full_dir}/{filename}")
                        if exists:
                            z_m_k_cost[log_z] = self.__read_data_cost_from_file(full_dir, filename, log_M, log_Y, step_size, read_step_size=ss, log_M_lower=0)

        return z_m_k_cost

    def __read_data_cost_from_file(self, full_dir, filename, log_M, log_Y, step_size, read_step_size=0, log_M_lower=0):
        """
        Reads cost estimation data from a specific file.

        :param full_dir: Target directory
        :param filename: target filename
        :param log_M: Log bound on number of combined bases
        :param log_Y: Log bound on number of combined nodes on level k
        :param step_size: Stepsize for M
        :param log_z: Log Jensen's z
        :param read_step_size: Stepsize for m to read from file
        :param log_M_lower: Lower bound on number of combined bases
        :return: [m] --> [k] --> Cost()
        """
        os.system('mkdir -p ' + full_dir)
        os.system('touch ' + full_dir + '/PLACEHOLDER')
        file = open(f"{full_dir}/{filename}", 'r')

        m_k_cost = defaultdict(lambda: defaultdict(None))

        if read_step_size == 0:
            read_step_size = step_size

        relevant_m = list(range(log_M_lower, log_M + 1, read_step_size))

        for line in file:
            if line[0] == '#' or not line.strip():
                continue
            cost_dict = json.loads(line)
            if cost_dict is not None:
                cost_dict['k'] = int(cost_dict['k'])

                if cost_dict['log_m'] in relevant_m:
                    cost = Cost()
                    cost.set_parameters(log_m = cost_dict['log_m'], log_y = cost_dict['log_y'], log_z = cost_dict['log_z'], k = cost_dict['k'])
                    cost.set_gcost(log_gcost_classical=cost_dict['log_gcost_classical'], log_gcost_quantum=cost_dict['log_gcost_quantum'])
                    cost.set_misc(qracm=cost_dict['qracm'], log_t_depth=cost_dict['log_t_depth'])
                    m_k_cost[cost_dict['log_m']][cost_dict['k']] = cost
        return m_k_cost

class ResultFiles(CostFiles):
    """
        Output result for cost estimation
    """
    def __init__(self, q, constants, bound, results_dir='costResults'):
        """
        :param q: LWE modulus
        :param constants: Dict of constants 'boundCoefficient', 'num_nodes_branching_dfs', 'force_DF', 'const_reps_QPE', 'const_reps_W' determining bound on DetectMV and FindMV
        :param bound: Bounds on estimation of subtree size, 'LBLB', 'UBLB', 'UBUB'
        :param results_dir:
        """

        self.results_dir=results_dir
        self.bound = bound
        self.constants = constants
        super().__init__(q, constants, bound)

    def prepare_plot_dir(self, n, circuit, results_dir='costResults'):
        """
        Prepare dorectory for plots
        :param n: BKZ blocksize
        :param circuit: Circuit instantiation
        :param results_dir: Target directory for output files
        :return: string
        """
        full_dir = self.results_dir + '/' + str(circuit.desc) + '/' + str(n) + '/'
        os.system('mkdir -p ' + full_dir)
        os.system('touch ' + full_dir + '/PLACEHOLDER')
        return full_dir

    def prepare_critical_dir(self, n, circuit, results_dir='costResults'):
        """
        Prepare dorectory for plots
        :param n: BKZ blocksize
        :param circuit: Circuit instantiation
        :param results_dir: Target directory for output files
        :return: string
        """
        full_dir = self.results_dir + '/' + str(circuit.desc) + '/'
        os.system('mkdir -p ' + full_dir)
        os.system('touch ' + full_dir + '/PLACEHOLDER')
        return full_dir

    def get_plot_filename(self, circuit, n, max_depth, prefix, postfix):
        """
        (Instance-Unique) name for plot files.
        :param n: BKZ blocksize
        :param circuit: Circuit instantiation
        :param max_depth: MaxDepth constraint
        :param prefix: Filename prefix.
        :param postfix: Filename postfix.
        :return: string
        """
        strings = [prefix, circuit.desc, str(max_depth), str(n), self.bound, postfix]
        return f"{'_'.join(filter(None, strings))}"

    def write_table_critical(self, z_to_cost, intersections, n, max_depth, circuit, prefix, postfix, print_header=False):
        """
        Output of crititical values (ie. intersection with security bounds), eqaivalent to
            Table 6 in Section 6
            Table 9 in Appendix G
        :param n: BKZ blocksize
        :param max_depth: MaxDepth constraint on circuit depth
        :param circuit: Circuit instantiation
        :param z_to_cost: Results Dictionary [Jensen's z] --> cost
        :param log_M: Log bound on number of combined bases
        :param prefix: Prefix for filename
        :param postfix: Postfix for filenames
        :param print_header: Flag to print header of tables
        """
        full_dir = self.prepare_critical_dir(n, circuit)
        filename = f"{'_'.join(filter(None, [prefix, circuit.desc, self.bound, postfix]))}.critical"

        header = True
        if os.path.isfile(f"{full_dir}/{filename}"):
            header = False
        file = open(f"{full_dir}/{filename}", 'a')

        def line(circ_desc, n, md, z_aes_sec, z_quasi_quadratic, z_exp_sec):
            return '{z: >10};{n: >4};{m: >8};{a: >37};{q: >37};{e: >37}'.format(z=circ_desc,n=n, m=md, a=f"{z_aes_sec}", e=f"{z_exp_sec}", q=f"{z_quasi_quadratic}")

        l = line('Circuit', 'n', 'MaxDepth', 'Grover on AES', 'Quasi-Sqrt', '{128,192,256}', )
        l2 = line('----------', '----', '--------', '-------------', '-------------', '----------')

        if header:
            file.write(l + '\n')
            file.write(l2 + '\n')

        if print_header == True:
            print()
            print(l, flush=True)
            print(l2, flush=True)

        def format_intersection(z_value, inequality, k_value, y_value, total_g_cost):
            if inequality == '<=':
                return 'z >= ' + str(z_value) + ', k <= ' + str(k_value) + ', y = ' + str(y_value) + ', tgcost >= ' + str(math.floor(total_g_cost))
            elif inequality == '<':
                return 'z >= ' + str(z_value) + ', k <= ' + str(k_value) + ', y = ' + str(y_value) + ', tgcost = ' + str(math.floor(total_g_cost))
            else:
                return 'z <= ' + str(z_value) + ', k = - ' + ', y = ' + str(y_value) + ', tgcost >= ' + str(math.floor(total_g_cost))

        target_str = format_intersection(z_value=intersections['target'][0], inequality=intersections['target'][1], k_value=z_to_cost[intersections['target'][0]].k, y_value=z_to_cost[intersections['target'][0]].log_y, total_g_cost = intersections['target'][2])
        aes_grover_str = format_intersection(z_value = intersections['aes_grover'][0], inequality=intersections['aes_grover'][1], k_value=z_to_cost[intersections['aes_grover'][0]].k, y_value=z_to_cost[intersections['aes_grover'][0]].log_y, total_g_cost = intersections['aes_grover'][2])
        quasi_quadratic_str = format_intersection(z_value = intersections['quasi-quadratic'][0], inequality=intersections['quasi-quadratic'][1], k_value=z_to_cost[intersections['quasi-quadratic'][0]].k, y_value=z_to_cost[intersections['quasi-quadratic'][0]].log_y, total_g_cost = intersections['quasi-quadratic'][2])

        l = line(circuit.desc, n, max_depth,  aes_grover_str, quasi_quadratic_str, target_str)
        file.write(l + '\n')

        print(l)

        if n == 873:
            l = line('----------', '----', '--------', '-------------', '-------------', '----------')
            file.write(l + '\n')


    def write_gatecost_to_file(self, file, dim, circuit, op_Ra_gate_costs, verb=False):
        """
        Output of gatecost for quantum operator W.
        :param file: Target file
        :param dim: Dimension for quantum operator W
        :param circuit: Circuit instantiation
        :param op_Ra_gate_costs: Dictionary of gate costs
        :param verb:
        """
        log_op_Ra_gate_costs = {}
        for gateid, count in op_Ra_gate_costs.items():
            if count > 0:
                log_op_Ra_gate_costs[gateid] = round(float(doLog2(count)), 2)
            else:
                log_op_Ra_gate_costs[gateid] = 'Zero'

        def line(dim, Gcost, T_depth, T, CNOT, Cliffs, M, R):
            x = '{y: >11},{z:>7},{b: >7},{a: >7},{c: >7},{d: >7},{e: >7},{f: >7}'.format(y=dim,z=Gcost, b=T_depth, a=T, c=CNOT, d=Cliffs, e=M, f=R)
            return x

        if isinstance(circuit, Circuits.QSharpCircuit):
            head = line( 'Lattice-Dim','T-depth', 'GCost', 'T', 'CNOT', '1-Cliff', 'M', 'R')
            gcost = op_Ra_gate_costs[circuit.idT] + op_Ra_gate_costs[circuit.idCnot] + op_Ra_gate_costs[circuit.idCliffs] + op_Ra_gate_costs[circuit.idM] + op_Ra_gate_costs[circuit.idR]
            if gcost > 0:
                log_gcost = round(doLog2(gcost), 2)
            else:
                log_gcost = "Zero"

            l = line(dim, log_op_Ra_gate_costs[circuit.idTdepth],
                     log_gcost,
                     (log_op_Ra_gate_costs[circuit.idT]),
                     (log_op_Ra_gate_costs[circuit.idCnot]),
                     (log_op_Ra_gate_costs[circuit.idCliffs]),
                     (log_op_Ra_gate_costs[circuit.idM]),
                     (log_op_Ra_gate_costs[circuit.idR]))

        else:
            head = line('Lattice-Dim', 'gate', '-,' '-', '-', '-', '-', '-')

            gcost = op_Ra_gate_costs[circuit.idT]
            if gcost > 0:
                log_gcost = round(doLog2(gcost), 2)
            else:
                log_gcost = "Zero"

            l = line(log_gcost, round(doLog2(op_Ra_gate_costs[circuit.idT]), 2), '-','-', '-', '-', '-', '-')

        file.write(head + '\n')
        file.write(l)
        if verb==True:
            print(head)
            print(l)

    def write_table_costs(self, z_to_cost, n, max_depth, circuit, prefix, postfix):
        """
        Write table with costs for every Jensen'z of result to file.
        :param n:  BKZ blocksize
        :param max_depth: MaxDepth constraint.
        :param circuit: Circuit instantiation
        :param z_to_cost: Result Dictionary [Jensen's z] --> cost
        :param prefix: Filename prefix
        :param postfix: Filename postfix
        """
        full_dir = self.prepare_plot_dir(n, circuit)
        filename = self.get_plot_filename(circuit, n, max_depth, prefix, postfix)

        file = open(f"{full_dir}/{filename}.costs", 'w')
        def line(log_z, log_m, log_y, num_k, log_classical, log_quantum, log_qracm, log_t_depth):
            x = '{z: >3},{m: >3},{y: >3},{k: >3},{c: >11},{q: >11},{r: >7},{t: >6}'.format(z=log_z, m=log_m, y=log_y, k=num_k, c=log_classical, q=log_quantum, r=log_qracm,
                                                                                                   t=log_t_depth)
            return x

        file.write(line('z', 'm', 'y', 'k', 'classical', 'quantum', 'qracm', 't-dep'))
        file.write('\n')

        z_range = z_to_cost.keys()
        for z in sorted(z_range):
            k = z_to_cost[z].k

            if k < n + 1:
                log_m = z_to_cost[z].log_m
            cost = z_to_cost[z]


            file.write(line(z, cost.log_m, cost.log_y, cost.k, round(float(cost.log_gcost_classical), 2),
                            round(float(cost.log_gcost_quantum), 2), round(float(cost.qracm), 1), round(float(cost.log_t_depth), 2)))
            file.write('\n')


def add_logs(log_x, log_y):
    """
    Global function for addition of numbers given as logarithms.
    :param log_x: log of x
    :param log_y: log of y
    :return: log_2(2^x + 2^y)
    """
    return doNumber((RRR(2**log_x) + RRR(2**log_y)).log2())
    #return doLog2(2**log_x + 2**log_y)

def doLog2(x):
    """
    Global function for logarithm with base 2.
    :param x:
    :return:
    """
    #return RRR(x).log2()
    return math.log(x, 2)

@lru_cache(maxsize=None)
def doNumber(x):
    """
    Global function for converting number to float.
    :param x:
    :return:
    """
    #return RRR(x)
    return float(x)
