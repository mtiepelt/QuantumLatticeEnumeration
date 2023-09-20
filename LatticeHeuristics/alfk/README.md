## Dependencies

```
sudo apt install libgsl-dev python3-cffi
```

## Build

```
python3 build_alfk.py
```

## Run

From LatticeHeuristics:
 
```python
from alfk import _alfk

alpha = 1/16777216
k = 540
n = 623

def _higher_prec_alpha_k(alpha, k, n, max_iterations=64, convergence_err=1e-10):
    """ Note, this is not GSA's alpha, but rather alpha as in Corollary 2
        and Eq. 16 of [ANSS18]. 
    """
    a = k/2
    b = 1 + (n-k)/2
    return _alfk.lib.my_gsl_cdf_beta_Pinv(alpha, a, b, max_iterations, convergence_err)

print(_higher_prec_alpha_k(alpha, k, n)) # returns nan
print(_higher_prec_alpha_k(alpha, k, n, convergence_err=1e-11)) # returns 0.7422...
```