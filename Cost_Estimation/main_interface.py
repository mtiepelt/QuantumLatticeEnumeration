#!/usr/bin/python3
# -*- coding: utf-8 -*-

__author__ = "Anonymous"
__version__ = "0.0"
__status__ = "Research"

import CostEstimation
import Tools
import Circuits
import TreeHeuristics as TH
import plot_interface
import os
import CostFunctions as CF

# Hard-coded directories for QSharp Circuit
# Resource Estimates
qsharp_dir = '../QuantumLatticeEnum/qsharpEstimates/'
#  Bounds on number of children in subtree
all_bounddegree = {
    512: '../LatticeHeuristics/children/children_Kyber512_512_381.X',
    768: '../LatticeHeuristics/children/children_Kyber768_768_623.X',
    1024: '../LatticeHeuristics/children/children_Kyber1024_1024_873.X'
}

# Kyber Parameters and BKZ blocksizes from
# Aono, Y., Nguyen, P.Q., Seito, T., Shikata, J.: Lower bounds on lattice enu meration with extreme pruning. pp. 608–637 (2018). https://doi.org/10.1007/978-3-319-96881-0_21
# Eq. (16)
Kyber = {
        'kyber512' :
            {
                'n' : 512,
                'q' : 3329,
                'beta' : 406
            },
        'kyber768' :
            {
                'n' : 768,
                'q' : 3329,
                'beta' : 623,
            },
        'kyber1024' :
            {
                'n' : 1024,
                'q' : 3329,
                'beta' : 873
            }
    }

def get_circuits(circ, **kwargs):
    """
    Instantiation of quantum operator W according to section 4
    :param circ: 'Query' (Section 4.1) or 'Minimal' (Section 4.2)
    :param kwargs:
    :return:
    """
    circuits = []
    if circ == 'Query':
        return Circuits.QueryCircuit('Query')
    elif circ == 'Minimal':
        return Circuits.MinimalCircuit('Minimal')

    elif circ == 'QSharp':
        n = 512
        beta = 406
        if 'n' in kwargs.keys():
            n = kwargs['n']

        if 'beta' in kwargs.keys():
            beta = kwargs['beta']
        n = 1024

        fn_bounddegree = all_bounddegree[n]
        qsharp_dir = '../QuantumLatticeEnum/qsharpEstimates/'
        max_bound_below = Tools.parse_bound_node_degree_file(beta, fn_bounddegree)
        return Circuits.QSharpCircuit(beta=beta, max_bound_below=max_bound_below, fn_dir=qsharp_dir, desc='QSharp')
    else:
        assert True, f"Circuit {circ} not implemented!"
        exit(0)


def print_log_cost_operator_W(h, n, beta, log_z=1, log_m=1, log_y=1, boundCoefficient=1, const_reps_QPE=1, const_reps_W=1, num_nodes_branching_dfs=1, force_DF=0, verb=True):
    """
    Print a txt-table with example size for cost of operator W.
    :return:
    """
    circuit = get_circuits(n=n, beta=beta)

    full_dir = 'costResults/Gates'
    os.system('mkdir -p ' + full_dir)
    os.system('touch ' + full_dir + '/PLACEHOLDER')
    filename = 'Gates_dim-' + str(h) + '.gates'
    file = open(f"{full_dir}/{filename}", 'w')


    treeH = TH.TreeHeuristics()
    costF = CF.CostFunctions(treeH, n=h, bound_coefficient=boundCoefficient, nbases=2 ** log_m, pprime=1, const_reps_QPE=const_reps_QPE, const_reps_W=const_reps_W,
                             num_nodes_branching_dfs=num_nodes_branching_dfs)

    # Operator Gates
    log_num_its_qpe = costF.log_reps_W_kh(0, h, log_y, log_z)
    op_Ra_gate_costs = circuit.get_operator_gate_count(0, h, log_num_its_qpe=log_num_its_qpe)

    if verb == True:
        print("Estimated Gate Cost of Operator W (cf. Table 13)")

    file.write("-" * 20)
    file.write(f"\nLog-Cost of operator RA for log_m {log_m}, log_y {log_y}, log_z {log_z}  of lattice dim h {h}  \n")
    resultsFiles = Tools.ResultFiles(boundCoefficient, const_reps_QPE, const_reps_W, False)
    resultsFiles.write_gatecost_to_file(file, h, circuit, op_Ra_gate_costs, verb=verb)

def loop_cost_estimation(kyber_sets, lst_md, circs,
                         log_M, log_Y, step_size, log_Z, step_size_Z, log_M_lower,
                         num_cores, use_cached,
                         prefix, postfix,
                         boundCoefficient, const_reps_QPE, const_reps_W, num_nodes_branching_dfs, force_DF):
    """
    Loop over all combinations of circuits, max-depths and kyber parameters.
    :param kyber_sets:  List of Kyber paramater keys.
    :param lst_md:  List of MaxDepths.
    :param circs: List of circuit identifiers.
    :param log_M:   Log Bound on number of combined, randomized bases for extreme pruning.
    :param log_Y:   Log Bound on number of combined nodes for a virtual tree on level k.
    :param step_size:   Step-size for looping over m,y.
    :param log_Z:   Log bound for Jensen's Gap z.
    :param step_size_Z: Step-size for loopzing over z.
    :param log_M_lower: Lower log bound for Bound on number of combined, randomized bases.
    :param num_cores:   Number of cores for multi-processing.
    :param use_cached: Flag indicating if results are cached to files.
    :param prefix: Prefix for written files.
    :param postfix: Postfix for writting files.
    :param boundCoefficient:    Bound C on the node degree used for enumeration of the tree (cf. Section 3.1).
    :param const_reps_QPE:  Constant 1/b for number of repetitions in QPE (cf. Corollary 1)
    :param const_reps_W: Constant epsilon for number of repetition of QPE (cf. Section 3.1)
    :param num_nodes_branching_dfs: Number of nodes checked for branching in classical binary search within FMV.
    :param force_DF: Number of DMV calls in FMV (cf. Section 3.1)
    :return:
    """
    print_header = True
    for circ in circs:
        for md in lst_md:
            for key in kyber_sets:
                # Choose which circuit models to include here
                circuit = get_circuits(circ, n=Kyber[key]['n'], beta=Kyber[key]['beta'])

                cost_classical_quantum = CostEstimation.CostEstimation(n=Kyber[key]['beta'], q=Kyber[key]['q'], boundCoefficient=boundCoefficient,
                                                                             const_reps_QPE=const_reps_QPE, const_reps_W=const_reps_W, num_nodes_branching_dfs=num_nodes_branching_dfs, force_DF=force_DF,
                                                                             use_cached=use_cached, num_cores=num_cores)

                z_to_cheapest = cost_classical_quantum.estimate_quantum_enumeration(md, circuit, log_M, log_Y, step_size, log_Z, step_size_Z, log_M_lower)


                resultsFiles = Tools.ResultFiles(boundCoefficient, Kyber[key]['q'], const_reps_QPE, const_reps_W, num_nodes_branching_dfs, force_DF)
                resultsFiles.write_table_critical(Kyber[key]['beta'], md, circuit, z_to_cheapest, log_M, prefix=prefix, postfix=postfix, print_header=print_header)
                resultsFiles.write_table_costs(Kyber[key]['beta'], md, circuit, z_to_cheapest, prefix=prefix, postfix=postfix)
                print_header = False

                plot_interface.plot_cost_estimation_z_to_g_cost(resultsFiles.prepare_plot_dir(Kyber[key]['beta'], circuit), resultsFiles.get_plot_filename(circuit, Kyber[key]['beta'], md, prefix, postfix), z_to_cheapest, md, Kyber[key]['beta'], log_M, circuit, prefix, postfix=postfix)

            # End Kyber
            print('-' * 10)
        # End MaxDepths
        print('=' * 20)

def print_cost_header(prefix, boundCoefficient, const_reps_QPE, const_reps_W, DF, num_nodes_branching_dfs):
    """
    Output header of current cost estimation.
    :param prefix: Prefix for written files.
    :param postfix: Postfix for writting files.
    :param boundCoefficient:    Bound C on the node degree used for enumeration of the tree (cf. Section 3.1).
    :param const_reps_QPE:  Constant 1/b for number of repetitions in QPE (cf. Corollary 1)
    :param const_reps_W: Constant epsilon for number of repetition of QPE (cf. Section 3.1)
    :param num_nodes_branching_dfs: Number of nodes checked for branching in classical binary search within FMV.
    :param DF: Number of DMV calls in FMV (cf. Section 3.1)
    :return:
    """
    print("=" * 30)
    if DF == -1:
        print(f"  Cost Estimation with C={boundCoefficient}, epsilon={const_reps_QPE}, b={const_reps_W}, #nodes per branching on classical binary search = {num_nodes_branching_dfs}")
    else:
        print(f"  Cost Estimation with C={boundCoefficient}, epsilon={const_reps_QPE}, b={const_reps_W}, #nodes per branching on classical binary search = {num_nodes_branching_dfs}, D/F={DF}")
    print(f"  /w prefix {prefix}")
    print("=" * 30)


def lower_bound(kyber_sets, lst_md, circuits, log_M=64, log_Y=64, step_size=1, log_Z=64, step_size_z=1, num_cores=4, use_cached=True):
    """
    Cost estimation for lower bound (cf. Section 3.1).

    :param kyber_sets:  List of Kyber paramater keys.
    :param lst_md:  List of MaxDepths.
    :param circs: List of circuit identifiers.
    :param log_M:   Log Bound on number of combined, randomized bases for extreme pruning.
    :param log_Y:   Log Bound on number of combined nodes for a virtual tree on level k.
    :param step_size:   Step-size for looping over m,y.
    :param log_Z:   Log bound for Jensen's Gap z.
    :param step_size_Z: Step-size for loopzing over z.
    :param log_M_lower: Lower log bound for Bound on number of combined, randomized bases.
    :param num_cores:   Number of cores for multi-processing.
    :param use_cached: Flag indicating if results are cached to files.
    :return:
    """
    prefix = 'DF=QD=WQ=1'
    boundCoefficient = 1    # Lower bounds D/F
    const_reps_QPE = 1      # Lower bounds Q/D
    const_reps_W = 1    # Lower bound W/Q
    num_nodes_branching_dfs = 1 # Lower bound on number of nodes per level of binary search, for D/F
    force_DF = 1        # Force DF to a specific value

    print_cost_header(prefix, boundCoefficient, const_reps_QPE, const_reps_W, force_DF, num_nodes_branching_dfs)

    # Result for Lower Bounds on D/F, Q/D, W/Q
    # Plots Section X and Table 2
    loop_cost_estimation(kyber_sets, lst_md, circuits,
                         log_M, log_Y, step_size, log_Z, step_size_z, log_M_lower=0,
                         num_cores=num_cores, use_cached=use_cached,
                         prefix=prefix, postfix='',
                         boundCoefficient=boundCoefficient, const_reps_QPE=const_reps_QPE, const_reps_W=const_reps_W, num_nodes_branching_dfs=num_nodes_branching_dfs, force_DF=force_DF)


def non_lower_bound(kyber_sets, lst_md, circuits, log_M, log_Y, step_size, log_Z, step_size_z, num_cores, use_cached):
    """
    Cost estimation for non-lower bound (cf. Appendix F).
    :param kyber_sets:  List of Kyber paramater keys.
    :param lst_md:  List of MaxDepths.
    :param circs: List of circuit identifiers.
    :param log_M:   Log Bound on number of combined, randomized bases for extreme pruning.
    :param log_Y:   Log Bound on number of combined nodes for a virtual tree on level k.
    :param step_size:   Step-size for looping over m,y.
    :param log_Z:   Log bound for Jensen's Gap z.
    :param step_size_Z: Step-size for loopzing over z.
    :param log_M_lower: Lower log bound for Bound on number of combined, randomized bases.
    :param num_cores:   Number of cores for multi-processing.
    :param use_cached: Flag indicating if results are cached to files.
    :return:
    """
    # ------------------------------------------------
    # CASE 2
    # D/F = Q/D = 1, W/Q = 256 * sqrt()
    prefix = 'DF=1_e=20_b=64'
    boundCoefficient = 1   # ==> Q/D will be = const_reps_QPE
    const_reps_QPE = 20  # Lower bounds Q/D
    const_reps_W = 64  # Lower bound W/Q
    num_nodes_branching_dfs = 1
    force_DF = -1
    print_cost_header(prefix, boundCoefficient, const_reps_QPE, const_reps_W, force_DF, num_nodes_branching_dfs)

    # Plots Section X and Table 2
    loop_cost_estimation(kyber_sets, lst_md, circuits,
                         log_M, log_Y, step_size, log_Z, step_size_z, log_M_lower=0,
                         num_cores=num_cores, use_cached=use_cached,
                         prefix=prefix, postfix='',
                         boundCoefficient=boundCoefficient, const_reps_QPE=const_reps_QPE, const_reps_W=const_reps_W, num_nodes_branching_dfs=num_nodes_branching_dfs, force_DF=force_DF)

def main():
    import argparse
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-parameter-set', required=False, type=str, help="One of the list [kyber512, kyber768, kyber1024]")
    parser.add_argument('-m', required=False, type=int, default=64, help="Upper bound for log-Bound on number of combined enumeration trees.")
    parser.add_argument('-y', required=False, type=int, default=64, help="Upper bound for log-Bound on number of combined nodes for virtual tree.")
    parser.add_argument('-stepsize', required=False, type=int, default=1, help="Step size for combined nodes and trees.")
    parser.add_argument('-z', required=False, type=int, default=64, help="Upper bound for log-Jensen's z.")
    parser.add_argument('-stepsize-z', required=False, type=int, default=1, help="Step size for Jensen's z.")
    parser.add_argument('-circuits', required=False, type=str, help="One of the list [Query, Minimal]")
    #
    parser.add_argument('-max-depth', required=False, type=int, help="One of the list [40, 64, 96]")
    #
    parser.add_argument('-num-cores', required=False, type=int, default=2, help="Number of individual cores used.")

    args = parser.parse_args()

    #print(f"args.parameter_set {args.parameter_set}")

    # Kyber Parameters
    kyber_sets = ['kyber512', 'kyber768', 'kyber1024']
    if args.parameter_set is not None:
        if args.parameter_set not in kyber_sets:
            parser.print_help()
            exit()
        else:
            kyber_sets = [args.parameter_set]

    # Max Depth
    lst_md = [40, 64, 96, 0]
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

    print(f"Quantum Enumeration Cost Estimation")
    print(f" Output Dir: ./costResults")
    lower_bound(kyber_sets=kyber_sets, lst_md=lst_md, circuits=circuits, log_M=args.m, log_Y=args.y, step_size=args.stepsize, log_Z=args.z, step_size_z=args.stepsize_z,
                        num_cores=args.num_cores, use_cached=True)
    non_lower_bound(kyber_sets=kyber_sets, lst_md=lst_md, circuits=circuits, log_M=args.m, log_Y=args.y, step_size=args.stepsize, log_Z=args.z, step_size_z=args.stepsize_z,
                    num_cores=args.num_cores, use_cached=True)

# Main body
if __name__ == '__main__':
    main()
