from sage.all import save, sqrt, log, prod, RR, pi, gamma, e, RealField, cached_function, copy
from bkz_simulators.sim import Sim, LLLProfile


RRR = RealField(1000)
RRRR = RealField(10000)
@cached_function
def stirling_gamma(x):
    x = RRR(x)
    _gamma = ((2*pi*(x-1))**RRR(.5)) * (((x-1)/e)**(x-1))
    return _gamma

@cached_function
def gaussian_heuristic(v, n):
    """
        :param v:   lattice volume
        :param n:   lattice rank
    """
    from math import pi, sqrt, e
    v = RRR(v)
    n = RRR(n)
    pi_part = ((pi * n)**(1/(2*n)))
    sqrt_part = sqrt(n/(2 * pi * e))
    v_part = (v ** (1/n))
    print(f"v {v}, pi_part {pi_part}, sqrt_part {sqrt_part}, v_part {v_part}")
    return pi_part * sqrt_part * v_part


def lwe_gh(lwe):
    n, q, m = lwe['n'], lwe['q'], lwe['m']
    n, q, m = RRR(n), RRR(q), RRR(m)
    vol = q**m
    dim = n+m+1
    return gaussian_heuristic(vol, dim)


def profile_gh(profile):
    # profile contains square norms
    dim = len(profile)
    vol = sqrt(prod(map(RRR, profile)))
    print(f"vol {vol}")
    return gaussian_heuristic(vol, dim)

@cached_function
def vol_ball(dim, rad):
    dim, rad = RRR(dim), RRR(rad)
    rad_part = rad**dim
    pi_part = pi**(dim/2)
    gamma_part = stirling_gamma(dim/2 + 1)
    # print(f"vol_ball rad {rad_part}, pi {pi_part}, gamma {gamma_part}")
    vol = rad_part * pi_part / gamma_part
    return vol


def deltaf(beta):
    """
    # Taken from https://github.com/malb/lattice-estimator/blob/f38155c863c8a1fb8da83a10664685d3496fe219/estimator/reduction.py#L12
    Compute δ from block size β without enforcing β ∈ ZZ.
    δ for β ≤ 40 were computed as follows:
    ```
    # -*- coding: utf-8 -*-
    from fpylll import BKZ, IntegerMatrix
    from multiprocessing import Pool
    from sage.all import mean, sqrt, exp, log, cputime
    d, trials = 320, 32
    def f((A, beta)):
        par = BKZ.Param(block_size=beta, strategies=BKZ.DEFAULT_STRATEGY, flags=BKZ.AUTO_ABORT)
        q = A[-1, -1]
        d = A.nrows
        t = cputime()
        A = BKZ.reduction(A, par, float_type="dd")
        t = cputime(t)
        return t, exp(log(A[0].norm()/sqrt(q).n())/d)
    if __name__ == '__main__':
        for beta in (5, 10, 15, 20, 25, 28, 30, 35, 40):
            delta = []
            t = []
            i = 0
                while i < trials:
                threads = int(open("delta.nthreads").read()) # make sure this file exists
                pool = Pool(threads)
                A = [(IntegerMatrix.random(d, "qary", beta=d//2, bits=50), beta) for j in range(threads)]
                for (t_, delta_) in pool.imap_unordered(f, A):
                    t.append(t_)
                    delta.append(delta_)
                i += threads
                print u"β: %2d, δ_0: %.5f, time: %5.1fs, (%2d,%2d)"%(beta, mean(delta), mean(t), i, threads)
            print
    ```
    """
    small = (
        (2, 1.02190),  # noqa
        (5, 1.01862),  # noqa
        (10, 1.01616),
        (15, 1.01485),
        (20, 1.01420),
        (25, 1.01342),
        (28, 1.01331),
        (40, 1.01295),
    )

    if beta <= 2:
        return RR(1.0219)
    elif beta < 40:
        for i in range(1, len(small)):
            if small[i][0] > beta:
                return RR(small[i - 1][1])
    elif beta == 40:
        return RR(small[-1][1])
    else:
        return RR(beta / (2 * pi * e) * (pi * beta) ** (1 / beta)) ** (1 / (2 * (beta - 1)))


def gsa_alpha(beta):
    delta = RRR(deltaf(beta))
    return delta ** -2   


def enum_cost_for_svp(tour, index, lwe, bkz):
    """
        :param gamma:   R = gamma * GH(lambda_1) is the enumeration radius
    """
    if tour <= 0:
        raise ValueError("Need -tour >= 1")
    if index <= 0:
        raise ValueError("Need -index >= 1")

    n, q, m = lwe['n'], lwe['q'], lwe['m']
    # run LLL
    sim = Sim(LLLProfile(n, q, m))

    # first simulate completion of previous tours 
    if tour > 1:
        sim(bkz['beta'], tour-1)
        profile = prev_tour_profile = sim.profile

    # # atempt simulating partial completion of current loop
    # sim(bkz['beta'], 1)
    # final_cur_tour_profile = sim.profile
    # cur_partial_prof = final_cur_tour_profile[:index-1] + prev_tour_profile[index-1:]
    # # NOTE: this does not work, the total volume is off.
    # # This difference may matter in practice because after previous SVP call, LLL is called on the basis,
    # # meaning some change may occur to further-down blocks. Will have to ignore

    n = min(bkz['beta'], len(profile) - index + 1)
    gamma = bkz['ghfactor']
    GH_lambda_1_block = profile_gh(profile[index-1:index-1+n])
    R = RRR(gamma * GH_lambda_1_block)
    print(f"R {R}")
    print(f"GH {GH_lambda_1_block}")
    print(f"gamma {gamma}")
    # exit(0)
    print(f"eumerating a block of rank {n}")
    b1_norm = profile[index-1]**RRR(.5)
    print(f"b1 {b1_norm}")
    alpha = gsa_alpha(bkz['beta'])
    # H_k = lambda k: .5 * (gamma**n) * (b1_norm ** (n-k)) * (alpha ** (.5*(n-k-1)*(n-k))) # wrong
    # H_k = lambda k: .5 * vol_ball(k, R) * (b1_norm ** (-k)) * (alpha ** (.5 * (k-1)*(k-2*n)))
    def H_k(k):
        k = RRR(k)
        volb = vol_ball(k, R)
        b1_part = (RRR(b1_norm) ** RRR(-k))
        # print("b1_norm", b1_norm)
        # print("-k", -k)
        # print("b1^k", RRR(b1_norm) ** RRR(k))
        # print("b1^-k", RRR(b1_norm) ** -RRR(k))
        # print("1/b1^k", RRR(1)/(RRR(b1_norm) ** RRR(k)))
        alpha_part = (alpha ** (.5 * (k-1)*(k-2*n)))
        # if (b1_part == 0.0):
        #     print(f"k {k}, volb {volb}, b1 {b1_part}, alpha {alpha_part}")
        #     print()
        #     exreturnit(0)
        _cost = RRRR(volb) * RRRR(b1_part) * RRRR(alpha_part) / RRRR(2)
        # print(f"k {k}, _cost {_cost}, volb {volb}, b1 {b1_part}, alpha {alpha_part}")
        return _cost
    level_nodes = {}
    for k in range(1, n+1):
        level_nodes[k] = H_k(k)
    # exit(0)

    cost = sum((level_nodes[k] for k in level_nodes.keys()))
    
    print(f"Cost of enumerating a block of rank {n} starting at index {index} on tour {tour} of BKZ-{bkz['beta']}")
    print("k\tlog(Hk,2)")
    for k in level_nodes.keys():
        print(f"{k}\t%.3f" % RRRR(log(level_nodes[k],2)))
    print(f"tot:\t%.3f" % RRRR(log(cost,2)))

    return level_nodes

    #print(f"csvp sieve\t%.3f" % (.292 * bkz['beta']))
    #chen13 = lambda beta: RR(0.270188776350190 * beta * log(beta) - 1.0192050451318417 * beta + 16.10253135200765 + log(100, 2))
    #print(f"chen13\t%.3f" % chen13(n))
    #print(f"q^(n^2/8)\t%.3f" % log(alpha**(-n**2/2),2))

    #from sage.all import line, save
    #g = line([(k, RRR(log(level_nodes[k],2))) for k in sorted(level_nodes.keys())])
    #save(g, "diamond.png")


def gen_cost_estimates(tour, index, lwe, bkz, estimates_dir = './estimations/nodes/'):
    """
        Gets estimates and writes them to file
    """
    import os
    level_nodes = enum_cost_for_svp(tour, index, lwe, bkz)

    filename = f"{estimates_dir}/nodes_{bkz['beta']}-{lwe['q']}.X"
    os.system('mkdir -p '+ estimates_dir)

    f = open(f"{filename}", "w")
    f.write(', '.join([f"%.3f" % RRRR(log(v,2)) for v in level_nodes.values()]))


def interface(n = 72, q = 97, sd = .1, m =87, beta=30, tours = 20, ghfactor=1.1, tour=None, index=None):

    if tour == None:
        tour = tours - 1 
    if index == None:
        index = 1

    lwe = { 'n': n, 'q': q, 'sd': sd, 'm': m }
    bkz = { 'beta': beta, 'tours': tours, 'ghfactor': ghfactor}

    for i in range(8, beta+1):
        lwe = { 'n': n, 'q': q, 'sd': sd, 'm': m }
        bkz = { 'beta': i, 'tours': tours, 'ghfactor': ghfactor}
        gen_cost_estimates(tour, index, lwe, bkz)

    return enum_cost_for_svp(tour, index, lwe, bkz)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', type=int, default=72, help="LWE secret dimension")
    parser.add_argument('-q', type=int, default=97, help="LWE modulo")
    parser.add_argument('-sd', type=float, default=1., help="LWE discrete gaussian standard deviation")
    parser.add_argument('-m', type=int, default=87, help="LWE samples")
    parser.add_argument('-beta', type=int, default=30, help="BKZ block build/cythonized/sage/structure/element.c:35999size")
    parser.add_argument('-tours', type=int, default=20, help="Max BKZ tours")
    parser.add_argument('-ghfactor', type=float, default=1.1, help="BKZ GH factor")
    parser.add_argument('-tour', type=int, help="count from 1")
    parser.add_argument('-index', type=int, help="count from 1")
    args = parser.parse_args()

    interface(args.n, args.q, args.sd, args.m, args.beta, args.tours, args.ghfactor, args.tour, args.index)


def sims(lwe, bkz):
    n, q, sd, m = lwe['n'], lwe['q'], lwe['sd'], lwe['m']
    beta, tours, ghfactor = bkz['beta'], bkz['tours'], bkz['ghfactor']
    print(n, q, sd, m, beta, tours, ghfactor)

    lll_prof = LLLProfile(n, q, m)
    sim = Sim(lll_prof)
    save(sim.plot(), "test_lll.png")
    sim(beta, tours)
    save(sim.plot(), "test_bkz.png")
