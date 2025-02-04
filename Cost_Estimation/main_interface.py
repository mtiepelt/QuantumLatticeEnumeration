#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Marcel Tiepelt"
__version__ = "1.0"
__status__ = "Research"
import istarmap  # import to apply patch

import CostEstimation
import Tools
import Circuits
import plots
import multiprocessing as mpp

# Import parameters of Kyber
from configSecurity import Kyber

def get_circuits(circ, **kwargs):
    """
    Instantiation of quantum operator W according to section 4
    :param circ: 'Query' (Section 4.1) or 'Minimal' (Section 4.2)
    :param kwargs:
    :return:
    """
    if circ == 'Query':
        return Circuits.QueryCircuit('Query')
    elif circ == 'Minimal':
        if 'q' in kwargs.keys():
            q = kwargs['q']
        else:
            print(f"Minimal circuit requires the modulus q")
            exit(0)

        return Circuits.MinimalCircuit('Minimal', q=q)
    else:
        assert True, f"Circuit {circ} not implemented!"
        exit(0)


def optimized_pool_func(ins_para, cost_classical_quantum, key, md, circuit, log_M, log_M_lower, return_dict):
    """
    Pool function for multiprocessing of optimized_costing_loop
    :param ins_para: Instance parameters.
    :param cost_classical_quantum: Cost Estimation instance.
    :param key: Identifier for parameter set.
    :param md: Max Depth.
    :param circuit: Circuit instance.
    :param log_M: Bound on number of randomized, combined bases.
    :param log_M_lower: Lower bound on number of randomized, combined bases.
    :param return_dict: MPP dictionary for results.
    """
    z_to_cheapest = cost_classical_quantum.estimate_quantum_enumeration(md, circuit, log_M_lower)
    intersections = Tools.get_intersections(cost_classical_quantum.n, md, log_M, z_to_cheapest, constants=cost_classical_quantum.constants, bound=ins_para['bound'])
    return_dict[(circuit.desc, md, key)] = [z_to_cheapest,intersections]

def optimized_costing_loop(ins_para, cl_para, constants, use_cached=False):
    """
    Loop over all combinations of circuits, max-depths and kyber parameters.

    :param ins_para: Instance parameters.
    :param cl_para: Costing loop parameters
    :param constants: Quantum walk constants
    :param use_cached: Flag determining if intermediate results are cached
    :return: Dictionary mapping key of instnace parameters to results.
    """
    # Computing
    manager = mpp.Manager()
    return_dict = manager.dict()

    pool_inputs = []
    for circ in ins_para['circuits']:
        for md in ins_para['lst_md']:
            for key in ins_para['kyber_sets']:
                # Choose which circuit models to include here
                circuit = get_circuits(circ, n=Kyber[key]['n'], beta=Kyber[key]['beta'], q=Kyber[key]['q'])
                fast_cost_classical_quantum = CostEstimation.CostEstimationFast(n=Kyber[key]['beta'], q=Kyber[key]['q'], constants=constants, log_M=cl_para['log_M'], log_Y=cl_para['log_Y'], log_Z=cl_para['log_Z'], bound=ins_para['bound'], step_size=cl_para['step_size_z'], use_cached=use_cached)
                pool_inputs.append([ins_para, fast_cost_classical_quantum, key, md, circuit, cl_para['log_M'], cl_para['log_M_lower'], return_dict])
            # End Kyber
        # End MaxDepths

    # POOL IT
    pool_size = min(len(pool_inputs), int((0.9 * mpp.cpu_count()) - 2))
    print(f"Using: {pool_size} cores on {len(pool_inputs)} inputs.")

    pool = mpp.Pool(pool_size)
    try:
        import tqdm
        tqdm_present = True
    except:
        tqdm_present = False

    if tqdm_present:
        for _ in tqdm.tqdm(pool.istarmap(optimized_pool_func, pool_inputs), total=len(pool_inputs)):
            pass
    else:
        results = pool.starmap(optimized_pool_func, pool_inputs)

    pool.close()
    pool.join()
    return return_dict

def plot_and_print(ins_para, cl_para, constants, return_dict, prefix, postfix):
    """
    :param ins_para: Instance parameters.
    :param cl_para: Costing loop parameters
    :param constants: Quantum walk constants
    :return: Dictionary mapping key of instnace parameters to results.
    :param prefix: File prefix.
    :param postfix: File postfix.
    :return:
    """
    print_header = True
    # Ploting and Writing Results
    for circ in ins_para['circuits']:
        for md in ins_para['lst_md']:
            for key in ins_para['kyber_sets']:
                n = Kyber[key]['n']
                beta = Kyber[key]['beta']
                q = Kyber[key]['q']
                circuit = get_circuits(circ, n=n, beta=beta, q=q)

                if (circuit.desc, md, key) in return_dict.keys():
                    z_to_cheapest = return_dict[(circuit.desc, md, key)][0]
                    intersections = return_dict[(circuit.desc, md, key)][1]
                else:
                    continue

                resultsFiles = Tools.ResultFiles(q, bound=ins_para['bound'], constants=constants)
                resultsFiles.write_table_critical(z_to_cheapest, intersections, beta, md, circuit, prefix=prefix, postfix=postfix, print_header=print_header)
                resultsFiles.write_table_costs(z_to_cheapest, beta, md, circuit, prefix=prefix, postfix=postfix)
                print_header = False

                plots.plot_cost_estimation_z_to_g_cost(resultsFiles.prepare_plot_dir(beta, circuit),
                                                                resultsFiles.get_plot_filename(circuit, beta, md, prefix, postfix),
                                                                z_to_cheapest, intersections, md, beta, cl_para['log_M'], constants, bound=ins_para['bound'])
            # End Kyber
            print('-' * 10)
        # End MaxDepths
        print('=' * 20)

def print_cost_header(ins_para, constants, prefix, use_cached):
    """
    Output header of current cost estimation.
    """
    print("=" * 50)
    print(' ' * 3 + f"Quantum Enumeration Cost Estimation for")
    print(' ' * 8 + f" {ins_para['kyber_sets']} over operator W as {ins_para['circuits']} and for all MaxDepth {ins_para['lst_md']} using ")
    print(' ' * 8 + f"... quantum walk constants: C={constants['boundCoefficient']}, epsilon={constants['const_reps_QPE']}, b={constants['const_reps_W']}, "
                     f"#nodes/branching on DFS = {constants['num_nodes_branching_dfs']}, D/F={constants['force_DF']}")
    print(' ' * 8 + f"... lattice estimation bound: {ins_para['bound']}")
    print(' ' * 8 + f"... caching results {use_cached} and writing files to ./costResults/ with prefix {prefix}")
    print("=" * 50)

def main():
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-parameter-set', required=False, type=str, help="One of the list [kyber512, kyber768, kyber1024]")
    parser.add_argument('-y', required=False, type=int, default=64, help="Upper bound for QRACM/ log-Bound on number of combined nodes for virtual tree.")
    #
    parser.add_argument('-z', required=False, type=int, default=64, help="Upper bound for log-Jensen's z.")
    parser.add_argument('-stepsize-z', required=False, type=int, default=1, help="Step size for Jensen's z.")
    #
    parser.add_argument('-circuits', required=False, type=str, help="One of the list [Query, Minimal]")
    #
    parser.add_argument('-max-depth', required=False, type=int, help="One of the list [40, 64, 96, 9999, -1], where 9999 is unbounded, -1 is trivial attack")
    #
    parser.add_argument('-bound', required=False, type=str, default="LBUB", help="Uses upper/lower bound for subtree sizes, [LBUB, UBUB, LBLB]")
    parser.add_argument('-cache', required=False, action='store_true', help="Cache intermediate computation.")

    # Development parameters
    parser.add_argument('-m-lower', required=False, type=int, default=64, help="Dev: Log of the minimal number of bases considered for extreme pruning. Pruning Radi only support [64].")
    parser.add_argument('-m', required=False, type=int, default=64, help="Dev: Upper bound for log-Bound on number of combined enumeration trees. Pruning Radi only support [64].")

    parser.add_argument('-df', required=False, type=int, default=1, help="Bound constant in number of calls to DetectMV in FindMV.")
    parser.add_argument('-qd', required=False, type=int, default=1, help="Bound constant in number of calls to QPE in DetectMV.")
    parser.add_argument('-wq', required=False, type=int, default=1, help="Bound constant in number of calls to operator W in QPE.")
    parser.add_argument('-nodesbranching', required=False, type=int, default=1, help="Bound on number of nodes checked on each level of the DFS.")
    parser.add_argument('-forcedf', required=False, type=int, default=1, help="Force a specific value for DF (overwrites -df).")

    args = parser.parse_args()

    # Kyber Parameters
    kyber_sets = ['kyber512', 'kyber768', 'kyber1024']
    if args.parameter_set is not None:
        if args.parameter_set not in kyber_sets:
            parser.print_help()
            exit()
        else:
            kyber_sets = [args.parameter_set]

    # Max Depth
    lst_md = [40, 64, 96, 9999, -1]
    if args.max_depth is not None:
        if args.max_depth not in lst_md:
            parser.print_help()
            exit()
        else:
            lst_md = [args.max_depth]

    # Suppported Circuit Models
    circuits = ['Query', 'Minimal']
    if args.circuits is not None:
        if args.circuits not in circuits:
            parser.print_help()
            exit()
        else:
            circuits = [args.circuits]

    ins_para = {
        'kyber_sets': kyber_sets,
        'circuits': circuits,
        'lst_md': lst_md,
        'bound': args.bound,
    }

    cl_para = {
        'log_M_lower' : args.m_lower,
        'log_M': args.m,
        'log_Y': args.y,
        'log_Z': args.z,
        'step_size_z': args.stepsize_z
    }

    # ==========================================================================================
    #   BOUND ON QW
    # ==========================================================================================
    # D/F = (x * n) Log C, Q/D, W/Q = b * sqrt()
    prefix = f"DF={args.df}_QD={args.qd}_WQ={args.wq}"
    constants = {
        'boundCoefficient': args.df,  # =: C, Lower bounds D/F
        'const_reps_QPE': args.qd,  # Lower bounds Q/D
        'const_reps_W': args.wq,  # =: b, Lower bound W/Q
        'num_nodes_branching_dfs': args.nodesbranching,  # =: x, Lower bound on number of nodes per level of binary search, for D/F
        'force_DF': args.forcedf  # Force DF to a specific value
    }

    # ==========================================================================================
    #   EXAMPLE BOUNDS ON QW for https://eprint.iacr.org/2023/1423.pdf
    # ==========================================================================================
    prefix, constants = paperBounds()

    # ==========================================================================================
    #   EXAMPLE BOUNDS ON QW for Dissertation
    # ==========================================================================================
    #prefix, constants = dissBounds()

    print_cost_header(ins_para, constants, prefix, args.cache)

    results = optimized_costing_loop(ins_para, cl_para, constants, use_cached=args.cache)
    plot_and_print(ins_para, cl_para, constants, results, prefix, 'LESSSIMPLE')

def paperBounds():
    # ==========================================================================================
    #   PRACTICAL BOUND ON QW
    # ==========================================================================================
    # D/F = (x * n) Log C, Q/D = 20, W/Q = b * sqrt()
    prefix = 'DF=nLogC_e=20_b=64'
    constants = {
        'boundCoefficient': 2,  # =: C, Lower bounds D/F
        'const_reps_QPE': 20,  # Lower bounds Q/D
        'const_reps_W': 64,  # := b Lower bound W/Q
        'num_nodes_branching_dfs': 1,  # Lower bound on number of nodes per level of binary search, for D/F
        'force_DF': -1  # =: x, Force DF to a specific value, -1 means its ignored
    }
    return prefix, constants

def dissBounds():
    # ==========================================================================================
    #   PRACTICAL BOUND ON QW
    # ==========================================================================================
    # D/F = (x * n) Log C, Q/D =20, W/Q = b * sqrt()
    prefix = 'DF=nLogC_e=20_b=64'
    constants = {
        'boundCoefficient': 1,  # =: Cm Lower bounds D/F
        'const_reps_QPE': 1,  # Lower bounds Q/D
        'const_reps_W': 64,  # =: b, Lower bound W/Q
        'num_nodes_branching_dfs': 1,  # Lower bound on number of nodes per level of binary search, for D/F
        'force_DF': 1  # =: x,  Force DF to a specific value
    }

    return prefix, constants


# Main body
if __name__ == '__main__':
    main()

