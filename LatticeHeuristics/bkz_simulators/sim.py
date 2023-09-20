from matplotlib import legend_handler
from sympy import reduced
import bkz_simulators.CN11 as CN11# this is identical to fpylll.tools.bkz_simulator
import bkz_simulators.BSW18 as BSW18
from fpylll import BKZ
from sage.all import line, log, e, RR
from lll_sim import lll_simulator


def deltaf(beta):
    """
    Taken from https://github.com/malb/lattice-estimator/blob/f38155c863c8a1fb8da83a10664685d3496fe219/estimator/reduction.py#L12

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
        from sage.all import pi, e
        return RR(beta / (2 * pi * e) * (pi * beta) ** (1 / beta)) ** (1 / (2 * (beta - 1)))


def gsa_alpha(beta):
    delta = RR(deltaf(beta))
    return delta ** -2   


def LLLProfile(n, q, m, nu=1, embedding="baigal", use_gsa=False):
    """
    Returns simulated LLL profile with Z-shape due to q-vectors being at the
    beginning of the input basis.

    :returns:   list of squared norms
    """
    if embedding == "kannan":
        _dim = m+1
        _k = m-n
        scale = 1
        if nu != 1:
            raise ValueError("ν != 1 makes sense only using baigal.")
    elif embedding == "baigal":
        _dim = n+m+1
        _k = m
        scale = nu
        assert(nu > 0)

    if use_gsa:
        log_delta = float(log(1.021900))
        log_vol = float(_k * log(q) + n * log(scale))
        log_alpha = 2 * _dim * log_delta/(1-_dim)
        # log_alpha = -2 * float(log(delta))
        log_bi = lambda i: (i-1) * log_alpha + log_vol/_dim + _dim * log_delta # i from 1 to _dim
        # vvol = sum([log_bi(i+1) for i in range(_dim)])
        # print("original logvol", log_vol)
        # print("recomputed lvol", vvol)
        return [e**(2*log_bi(i+1)) for i in range(_dim)]

    return lll_simulator(_dim, _k, q, scale=scale)


class Sim:
    """
    Class to simulate BKZ reduction and variants on random q-ary lattices.

    TESTS:
        >>> import bkz_simulators.CN11 as CN11
        >>> from fpylll import BKZ
        >>> from bkz_simulators.sim import Sim, LLLProfile
        >>> n, q, m = 50, 2**10, 50
        >>> beta, tours = 40, 16
        >>> lll_prof = LLLProfile(n, q, m)
        >>> r1, l1 = CN11.simulate(lll_prof, BKZ.Param(block_size=beta, max_loops=tours))
        >>> sim = Sim(lll_prof)
        >>> sim(beta, 8)
        >>> len(sim.tours)
        1
        >>> sim.tours[0] == (beta, 8)
        True
        >>> sim(beta, tours-8)
        >>> len(sim.tours)
        2
        >>> sim.tours[1] == (beta, tours-8)
        True
        >>> sim.profile == r1
        True
        >>> sim(40, 2)
        >>> len(sim.tours)
        3
        >>> sim.tours[2] == (40, 2)
        True
    """
    def __init__(self, initial_profile):
        self.profile = initial_profile[::]
        self.tours = []

    def __call__(self, beta, tours, sim="CN11", prng_seed=0xdeadbeef):
        if sim == "CN11":
            simulate = CN11.simulate
        elif sim == "BSW18":
            simulate = lambda prof, pars: BSW18.simulate(prof, pars, prng_seed=prng_seed)
        elif sim == "GSA":
            def gsa_sim(prof, pars):
                r = list(map(lambda x: log(x, 2) / 2.0, prof))
                n = len(r)
                log_vol = sum(r)
                beta = pars.block_size
                alpha = gsa_alpha(beta)
                log_b1 = -(n-1) * log(alpha, 2) / 2 + log_vol/n
                log_prof = [log_b1]
                for j in range(1, n):
                    i = j+1
                    log_prof.append((i-1) * log(alpha, 2) + log_b1)
                reduced_prof = list(map(lambda x: 2**(2*x), log_prof))
                return reduced_prof, pars.max_loops
            simulate = gsa_sim
        pars = BKZ.Param(block_size=beta, max_loops=tours)
        l, r = simulate(self.profile, pars)
        self.tours.append((beta, tours))
        self.profile = l

    def plot(self, base=2, label="log profile"):
        g = line(
                [(i, log(self.profile[i], base)/2) for i in range(len(self.profile))],
                axes_labels = ['$i$','$\\log_{%s}\\|b_{i}^*\\|$' % ("" if base == e else base)],
                legend_label=label
            )
        return g

