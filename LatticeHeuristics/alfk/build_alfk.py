from cffi import FFI
ffibuilder = FFI()

ffibuilder.cdef("double my_gsl_cdf_beta_Pinv (const double P, const double a, const double b, const unsigned int max_iterations, const double convergence_err);")

ffibuilder.set_source("_alfk",  # name of the output C extension
"""
    #include "alfk.h"
""",
    sources=['alfk.c'],   # includes pi.c as additional sources
    libraries=['m', 'gsl'])    # on Unix, link with the math library

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)