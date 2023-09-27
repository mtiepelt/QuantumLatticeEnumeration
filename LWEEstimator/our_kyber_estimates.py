from sage.all import load, sqrt, RR, log, ZZ, pi
try:
    load('estimator.py')
except:
    load('../OurEstimator/estimator.py')
    pass
import sys
sys.path.append("../LatticeHeuristics")
import enumTrees


# try adding a joblib.Memory decorator
try:
    from joblib import Memory
    memory = Memory("cachedir", verbose=0)
except:
    joblib_available = False
else:
    joblib_available = True


def conditional_decorator(dec, condition):
    # from https://stackoverflow.com/a/10724898
    def decorator(func):
        if not condition:
            # Return the function unchanged, not decorated.
            return func
        return dec(func)
    return decorator


# @conditional_decorator(memory.cache, joblib_available)
def classical_cost(n, sd, q, pruning="cilinder", verbose=False, **kwargs):

    if pruning == "no":
        PruningClass = enumTrees.NoPruning
    elif pruning == "linear":
        PruningClass = enumTrees.LinearPruning
    elif pruning == "cilinder":
        PruningClass = enumTrees.CilinderPruningLowerBound

    asympt = enumTrees.fit_curve(lambda n: PruningClass.log_cost(n, **kwargs))[0]

    if verbose:
        print("asymptotic total cost", asympt)

    cost_model = lambda beta, d, B: ZZ(2)**RR(asympt['nl2']*beta*log(beta,2) + asympt['n']*beta + asympt['l2']*log(beta,2) + asympt['c'] + (0 if 'nbases' not in kwargs else log(kwargs['nbases'], 2)))

    alpha = sqrt(2*pi)*sd/RR(q)
    secret_distribution = "normal"
    success_probability = 0.99
    m = 2*n

    return primal_usvp(n, alpha, q, secret_distribution=secret_distribution, m=m, success_probability=success_probability, reduction_cost_model=cost_model)


if __name__ == "__main__":
    q = 3329
    pprime = 1
    nbases = 2**64
    for n in [512, 768, 1024]:
        if n == 512:
            sd = sqrt(3/2)
        else:
            sd = 1
        cost = classical_cost(n, sd, q, pruning="cilinder", pprime=pprime, nbases=nbases)
        qcost = (log(cost['rop'], 2) + log(cost['beta'],2)).n()/2.
        print()
        print(f"Kyber {n}")
        print()
        print("cost assuming cilinder pruning")
        print(cost)
        print()
        print('naive cost approximation as in "Depth of QPE(W)" (log2)')
        print(qcost)
        print()
        print("----")
