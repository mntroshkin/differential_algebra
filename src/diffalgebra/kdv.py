from .diff_ring import FuncGenerator, DifferentialPolynomial

def KdV_hierarchy(variable: FuncGenerator, order: int = 1) -> DifferentialPolynomial:
    if order < 0:
        raise ValueError
    K = variable[1]
    for i in range(order):
        L = K.integral()
        K = L.diff(3) + 4 * variable * L.diff() + 2 * variable[1] * L
    return K
