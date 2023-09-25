from sage.all import vector, matrix, identity_matrix, ZZ, QQ, RR, GF
from sage.all import sqrt, floor, prod, log
from sage.all import randint, seed
from sage.crypto.lwe import LWE, DiscreteGaussianDistributionIntegerSampler
from utilities import balance, sage_to_fplll_matrix
import fpylll


def genLWEInstance(n,
                   q,
                   sd,
                   m,
                   float_type="d",
                   mpfr_precision=None,
                   embedding="baigal",
                   secret_dist="noise",
                   nu=1,
                   moot=False,
                   use_suboptimal_embedding=False):
    """Generate lattices from LWE instances using Kannan's embedding.
    :param n:                   secret dimension
    :param q:                   lwe modulo
    :param sd:                  standard deviation
    :param m:                   number of lwe samples ('baigal' only uses m-n)
    :param float_type:          floating point type
    :param mpfr_precision:      floating point precision (if using mpfr)
    :param embedding:           "kannan" or "baigal"
    :param nu:                  scaling factor for "baigal" embedding
    :param moot:                if true, c is uniform rather than As+e
    :param use_suboptimal_embedding: if False, put the q vectors on top of basis
    :returns:                   the lwe generator, the samples, the lattice and
                                its volume
    """

    # generate LWE instance
    lwe = LWE(n=n, q=q, secret_dist=secret_dist,
                D=DiscreteGaussianDistributionIntegerSampler(sd))

    # get m different LWE samples
    samples = [lwe() for i in range(m)]

    A = matrix(a for a, _ in samples)
    if moot:
        C = matrix(randint(0, q-1) for _ in range(len(samples)))
    else:
        C = matrix(c for _, c in samples)

    # print(lwe._LWE__s)
    # print(lwe._LWE__e)
    # print(vector(GF(q), C[0]) == vector(GF(q), A * lwe._LWE__s + vector(lwe._LWE__e)))
    evec = vector(GF(q), C[0]) - vector(GF(q), A * lwe._LWE__s)

    if embedding == "kannan":
        # generate kannan's embedding lattice
        AT = A.T.echelon_form()
        qs = matrix(ZZ, m-n, n).augment(q*identity_matrix(m-n))
        if use_suboptimal_embedding:
            B = AT.change_ring(ZZ).stack(qs)
        else:
            B = qs.stack(AT.change_ring(ZZ))
        # embed the ciphertext to the lattice, so that error vector
        # becomes the (most likely unique) SVP in the lattice
        BC = B.stack(matrix(C).change_ring(ZZ))
        BC = BC.augment(matrix(m+1, 1))
        BC[-1, -1] = max(floor(sd), 1)
    elif embedding == "baigal":
        # generate scaled Bai-and-Galbraith's embedding lattice
        assert(nu > 0)
        nu_rat = QQ(round(nu*100)/100)
        nu_num, nu_denom = nu_rat.numerator(), nu_rat.denominator()
        if nu_denom != 1:
            print(
                f"WARNING: due to fractional ν, output lengths are scaled by {nu_denom}")
        AT = (nu_num*identity_matrix(n)).augment(nu_denom*A.change_ring(ZZ).T)
        qs = matrix(ZZ, m, n).augment(nu_denom*q*identity_matrix(m))
        if use_suboptimal_embedding:
            B = AT.change_ring(ZZ).stack(qs)
        else:
            B = qs.stack(AT.change_ring(ZZ))
        # embed the ciphertext to the lattice, so that error vector
        # becomes the (most likely unique) SVP in the lattice
        BC = B.stack(matrix(1, n).augment(nu_denom*matrix(C).change_ring(ZZ)))
        BC = BC.augment(matrix(m+n+1, 1))
        BC[-1, -1] = nu_denom*max(floor(sd), 1)
    else:
        raise ValueError("embedding can only be 'kannan' or 'baigal'")

    # preprocess basis for low precision, else after GSO has float_type
    BC = fpylll.IntegerMatrix.from_matrix(BC)
    if float_type == "d":
        fpylll.LLL.reduction(BC)

    # set floating point precision
    if float_type == "mpfr":
        _ = fpylll.FPLLL.set_precision(mpfr_precision)

    BC_GSO = fpylll.GSO.Mat(BC, float_type=float_type)

    BC_GSO.update_gso()

    if float_type != "d":
        lll = fpylll.LLL.Reduction(BC_GSO)
        lll()

    # get lattice volume
    vol = sqrt(prod([RR(BC_GSO.get_r(i, i)) for i in range(n+m+1)]))

    return (lwe, samples, A, C, BC_GSO, vol, evec)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-randseed', type=int, default=0xdeadbeef, help="PRNG seed")
    parser.add_argument('-n', type=int, default=72, help="LWE secret dimension")
    parser.add_argument('-q', type=int, default=97, help="LWE modulo")
    parser.add_argument('-sd', type=float, default=1., help="LWE discrete gaussian standard deviation")
    parser.add_argument('-m', type=int, default=87, help="LWE samples")
    parser.add_argument('-printsol', action="store_true", help="Print secret, error and target vectors")
    args = parser.parse_args()

    n, q, sd, m = args.n, args.q, args.sd, args.m
    with seed(args.randseed):
        lwe, samples, A, C, L, vol, evec = genLWEInstance(n, q, sd, m)
        if args.printsol:
            import json
            sol = {
                's': list(map(int, balance(lwe._LWE__s))),
                'e': list(map(int, balance(evec))),
                'sv': list(map(int, balance(vector(GF(q), list(lwe._LWE__s) + list(evec) + [max(floor(sd), 1)]))))
            }
            print(json.dumps(sol, indent=2))
        else:
            M = matrix(nrows=L.B.nrows, ncols=L.B.ncols)
            L.B.to_matrix(M)
            basis_s = sage_to_fplll_matrix(M)
            print(basis_s)