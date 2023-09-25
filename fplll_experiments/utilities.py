from sage.all import sage_eval, matrix
from sage.all import parent, ZZ


def balance(e, q=None):
    """ Return a representation of `e` with elements balanced between `-q/2` and `q/2`
    :param e: a vector, polynomial or scalar
    :param q: optional modulus, if not present this function tries to recover it from `e`
    :returns: a vector, polynomial or scalar over/in the integers
    """
    try:
        p = parent(e).change_ring(ZZ)
        return p([balance(e_, q=q) for e_ in e])
    except (TypeError, AttributeError):
        if q is None:
            try:
                q = parent(e).order()
            except AttributeError:
                q = parent(e).base_ring().order()
        e = ZZ(e)
        e = e % q
        return ZZ(e-q) if e > q//2 else ZZ(e)


def fplll_to_sage_matrix(mat_s):
    """
    TEST:
        >>> mat_s = "[[4 3 7]\n[3 0 1]\n[3 5 3]\n]\n"
        >>> mat = fplll_to_sage_matrix(mat_s)
        >>> assert mat == matrix([
        >>>     [4, 3, 7],
        >>>     [3, 0, 1],
        >>>     [3, 5, 3],
        >>> ])
    """
    L = []
    for line in mat_s.split('\n'):
        row_s = line.replace("[[","[").replace("]]","]").replace(" ", ", ")
        if row_s not in ["", "]"]:
            L.append(sage_eval(row_s))
    return matrix(L)


def sage_to_fplll_matrix(mat):
    """
    TEST:
        >>> mat = matrix([
        >>>     [4, 3, 7],
        >>>     [3, 0, 1],
        >>>     [3, 5, 3],
        >>> ])
        >>> mat_s = sage_to_fplll_matrix(mat)
        >>> assert mat_s == "[[4 3 7]\n[3 0 1]\n[3 5 3]\n]\n"
    """
    mat_s = "["
    nrows = mat.nrows()
    for r in range(nrows):
        row = mat[r]
        row_s = str(row).replace('(', '[').replace(')', ']').replace(",", "") + "\n"
        mat_s += row_s
    mat_s += "]\n"
    return mat_s
