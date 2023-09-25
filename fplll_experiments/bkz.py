from sage.all import matrix
import fpylll
from fpylll import BKZ, Enumeration, EnumerationError
from fpylll.fplll.bkz_param import BKZParam
from fpylll.algorithms.bkz2 import BKZReduction
from fpylll.tools.bkz_stats import BKZTreeTracer, Accumulator, pretty_dict, dummy_tracer
from fpylll.tools.quality import basis_quality
import time
from time import process_time
from collections import OrderedDict
from utilities import sage_to_fplll_matrix, fplll_to_sage_matrix
from settings import FPLLL_PATH
from enumeration import enum_cost_for_svp


class MyBKZTreeTracer(BKZTreeTracer):
    """
    """
    def exit(self, **kwds):  # noqa, shut up linter about this function being too complex
        """
        By default CPU and wall time are recorded.  More information is recorded for "enumeration"
        and "tour" labels.  When the label is a tour then the status is printed if verbosity > 0.
        """
        node = self.current
        label = node.label

        node.data["cputime"] += process_time()
        node.data["walltime"] += time.time()

        if kwds.get("dump_gso", False):
            node.data["r"] = node.data.get("r", []) + [self.instance.M.r()]

        if label == "enumeration":
            full = kwds.get("full", True)
            if full:
                try:
                    # for ind in range(-1, kwds['block_size']):
                    #     print(kwds["enum_obj"].get_nodes(ind))
                    # print("--------")
                    # print(enum_nodes)
                    # print(node)
                    node.data["#enum"] = Accumulator(kwds["enum_obj"].get_nodes(), repr="sum") + node.data.get("#enum", None)  # noqa
                except KeyError:
                    pass
                try:
                    node.data["%"] = Accumulator(kwds["probability"], repr="avg") + node.data.get("%", None)
                except KeyError:
                    pass
            try:
                # for ind in range(-1, kwds['block_size']):
                #     print(kwds["enum_obj"].get_nodes(ind))
                # print("--------")
                nodes_list = [kwds["enum_obj"].get_nodes(ind) + 1 for ind in range(kwds['block_size'])]
                # node.data["nodes"] = [kwds["enum_obj"].get_nodes(ind) + 1 for ind in range(kwds['block_size'])]
                node.data["nodes"] = [(kwds["kappa"], nodes_list)] + (node.data.get("nodes", None) if node.data.get("nodes", None) else []) # noqa
            except KeyError:
                pass

        if label[0] == "tour":
            data = basis_quality(self.instance.M)
            for k, v in data.items():
                if k == "/":
                    node.data[k] = Accumulator(v, repr="max")
                else:
                    node.data[k] = Accumulator(v, repr="min")

        if self.verbosity and label[0] == "tour":
            report = OrderedDict()
            report["i"] = label[1]
            report["cputime"] = node["cputime"]
            report["walltime"] = node["walltime"]
            try:
                report["preproc"] = node.find("preprocessing", True)["cputime"]
            except KeyError:
                pass
            try:
                report["svp"] = node.find("enumeration", True)["cputime"]
            except KeyError:
                pass
            report["#enum"] = node.sum("#enum")
            report["lll"] = node.sum("cputime", label="lll")
            try:
                report["pruner"] = node.find("pruner", True)["cputime"]
            except KeyError:
                pass
            report["r_0"] = node["r_0"]
            report["/"] = node["/"]

            print(pretty_dict(report))

        self.current = self.current.parent


class MyBKZReduction(BKZReduction):

    def __call__(self, params, min_row=0, max_row=-1):
        """Run the BKZ algorithm with parameters `param`.

        :param params: BKZ parameters
        :param min_row: start processing in this row
        :param max_row: stop processing in this row (exclusive)

        """
        try:
            label = params["name"]
        except KeyError:
            label = "bkz"
        tracer = MyBKZTreeTracer(self, root_label=label, verbosity=params.flags & BKZ.VERBOSE, start_clocks=True)

        if params.flags & BKZ.AUTO_ABORT:
            auto_abort = BKZ.AutoAbort(self.M, self.A.nrows)

        cputime_start = process_time()

        with tracer.context("lll"):
            self.lll_obj()

        i = 0
        while True:
            with tracer.context("tour", i, dump_gso=params.flags & BKZ.DUMP_GSO):
                clean = self.tour(params, min_row, max_row, tracer)
            i += 1
            if clean or params.block_size >= self.A.nrows:
                break
            if (params.flags & BKZ.AUTO_ABORT) and auto_abort.test_abort():
                break
            if (params.flags & BKZ.MAX_LOOPS) and i >= params.max_loops:
                break
            if (params.flags & BKZ.MAX_TIME) and process_time() - cputime_start >= params.max_time:
                break

        tracer.exit()
        self.trace = tracer.trace
        return clean

    def svp_reduction(self, kappa, block_size, params, tracer=dummy_tracer):
        """

        :param kappa:
        :param block_size:
        :param params:
        :param tracer:

        """

        self.lll_obj.size_reduction(0, kappa+1)
        old_first, old_first_expo = self.M.get_r_exp(kappa, kappa)

        remaining_probability, rerandomize = 1.0, False

        while remaining_probability > 1. - params.min_success_probability:
            with tracer.context("preprocessing"):
                if rerandomize:
                    with tracer.context("randomization"):
                        self.randomize_block(kappa+1, kappa+block_size,
                                             density=params.rerandomization_density, tracer=tracer)
                with tracer.context("reduction"):
                    self.svp_preprocessing(kappa, block_size, params, tracer=tracer)

            with tracer.context("pruner"):
                radius, re, pruning = self.get_pruning(kappa, block_size, params, tracer)

            try:
                enum_obj = Enumeration(self.M)
                with tracer.context("enumeration",
                                    enum_obj=enum_obj,
                                    probability=pruning.expectation,
                                    full=block_size==params.block_size,
                                    kappa=kappa,
                                    block_size=block_size):
                    max_dist, solution = enum_obj.enumerate(kappa, kappa + block_size, radius, re,
                                                            pruning=pruning.coefficients)[0]
                with tracer.context("postprocessing"):
                    self.svp_postprocessing(kappa, block_size, solution, tracer=tracer)
                rerandomize = False

            except EnumerationError:
                rerandomize = True

            remaining_probability *= (1 - pruning.expectation)

        self.lll_obj.size_reduction(0, kappa+1)
        new_first, new_first_expo = self.M.get_r_exp(kappa, kappa)

        clean = old_first <= new_first * 2**(new_first_expo - old_first_expo)
        return clean


def run_fpylll_BKZ(L, b, max_tours,strategies=fpylll.BKZ.DEFAULT_STRATEGY, lwe=None, verbose=False):
    """Set up and run Algorithm 2 from the paper, recording detailed statistics.
    :param L:           lattice basis
    :param b:           BKZ block size
    :param max_tours:   max number of BKZ tours
    :returns:           the BKZ object and the tracer containing statistics

    TEST:
        >>> basis = matrix([
        >>>     [4, 3, 7],
        >>>     [3, 0, 1],
        >>>     [3, 5, 3],
        >>> ])
        >>> reduced_basis = run_fpylll_BKZ(basis, 2, 5)
        >>> assert reduced_basis == matrix([
        >>>     [ 3, 0,  1],
        >>>     [ 2, 2, -3],
        >>>     [ 0, 5,  2]
        >>> ])
    """

    # preprocess basis for low precision, else after GSO has float_type
    BC = fpylll.IntegerMatrix.from_matrix(L)
    fpylll.LLL.reduction(BC)
    BC_GSO = fpylll.GSO.Mat(BC, float_type="d")
    BC_GSO.update_gso()

    profile = [BC_GSO.get_r(i,i) for i in range(BC.nrows)]

    # # get lattice volume
    # vol = sqrt(prod([RR(BC_GSO.get_r(i, i)) for i in range(len(L))]))

    # set up BKZ
    params_fplll = BKZParam(block_size=b,
                            strategies=strategies,
                            flags=0
                            | (int(fpylll.BKZ.VERBOSE) * int(verbose))
                            | fpylll.BKZ.AUTO_ABORT
                            # | fpylll.BKZ.GH_BND
                            | fpylll.BKZ.MAX_LOOPS,
                            # ,
                            max_loops=max_tours)
    if verbose:
        print("strategies:")
        for strat in params_fplll.strategies:
            print(strat)
        print()
    bkz = MyBKZReduction(BC_GSO)
    bkz(params_fplll)

    if verbose:
        tour = 0
        for node in bkz.trace.children[1:]:
            tour += 1
            # if tour != 10:
            #     continue
            assert(node.label[0] == 'tour')
            enum = list(filter(lambda o: 'enumeration' in o.label, node.children))[0]
            nodes = enum.data['nodes']
            # print(f"tour {node.label[1]}, nodes {nodes}")
            print()
            for index in range(len(nodes)):
                sim_nodes = enum_cost_for_svp(tour, index+1, None, {'beta': b, 'ghfactor': 1.}, initial_profile=profile, verbose=False)
                measured_nodes = list(filter(lambda x: x[0] == index, nodes))[0][1]
                # print(sim_nodes)
                # print(measured_nodes)
                print(f"tour {node.label[1]}, index %02d" % index, end="  ")
                print([float("%.1f" % (measured_nodes[x]/sim_nodes[x])) for x in range(len(measured_nodes))])
                # print()
            print()
            print()
            print()
            # exit(0)

    M = matrix(nrows=BC_GSO.B.nrows, ncols=BC_GSO.B.ncols)
    BC_GSO.B.to_matrix(M)

    return M


def run_fplll_BKZ(basis, b, max_tours, fplll_path=FPLLL_PATH+"fplll", tree_filename=None, gso_filename=None, tree_stats_filename=None, gh_bound=1.1):
    """Set up and run Algorithm 2 from the paper, recording detailed statistics.
    :param L:           lattice basis
    :param b:           BKZ block size
    :param max_tours:   max number of BKZ tours
    :returns:           the BKZ object and the tracer containing statistics

    TEST:
        >>> basis = matrix([
        >>>     [4, 3, 7],
        >>>     [3, 0, 1],
        >>>     [3, 5, 3],
        >>> ])
        >>> reduced_basis = run_fplll_BKZ(basis, 2, 5, "./fplll/out/bin/fplll", "test_tree.json", "test_gso.json")
        >>> assert reduced_basis == matrix([
        >>>     [ 3, 0,  1],
        >>>     [ 2, 2, -3],
        >>>     [ 0, 5,  2]
        >>> ])
    """

    basis_s = sage_to_fplll_matrix(basis)
    import subprocess
    command = [
        fplll_path,
        "-a", "bkz",
        "-b", str(b),
        "-bkzautoabort",
        "-bkzmaxloops", str(max_tours),
        "-bkzghbound", str(gh_bound)
    ]
    if tree_filename:
        command += ["-bkzdumpenumtrees", tree_filename]

    if gso_filename:
        command += ["-bkzdumpgso", gso_filename]

    if tree_stats_filename:
        command += ["-bkzdumpenumtreestats", tree_stats_filename]

    process = subprocess.run(command, input=bytes(basis_s, encoding="ascii"), capture_output=True)
    reduced_basis_s = process.stdout.decode(encoding='ascii')
    reduced_basis = fplll_to_sage_matrix(reduced_basis_s)
    return reduced_basis


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', type=int, default=20, help="BKZ block size")
    parser.add_argument('-maxtours', type=int, default=20, help="BKZ maximum number of tours")
    parser.add_argument('-fpylll', action='store_true', help='Run the standard fpylll BKZ')
    parser.add_argument('-tree-filename', type=str, default=None)
    parser.add_argument('-tree-stats-filename', type=str, default=None)
    parser.add_argument('-SE93', action='store_true', help='disable pruning and preprocessing; -fpylll only')
    parser.add_argument('-verbose', action='store_true')
    parser.add_argument('files', metavar='FILE', nargs='*', help='files to read, if empty, stdin is used')
    args = parser.parse_args()

    import fileinput
    basis_s = "".join(fileinput.input(files=args.files if len(args.files) > 0 else ('-', )))
    basis = fplll_to_sage_matrix(basis_s)

    if args.fpylll:
        if args.SE93:
            reduced_basis = run_fpylll_BKZ(basis, args.b, args.maxtours, strategies=None, verbose=args.verbose)
        else:
            reduced_basis = run_fpylll_BKZ(basis, args.b, args.maxtours, verbose=args.verbose)
    else:
        reduced_basis = run_fplll_BKZ(basis, args.b, args.maxtours, FPLLL_PATH + "fplll", tree_filename=args.tree_filename, tree_stats_filename=args.tree_stats_filename)

    print(reduced_basis)
