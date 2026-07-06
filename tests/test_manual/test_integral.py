from fractions import Fraction
import pytest
import diffalgebra as da


def test_integral_partial():
    R = da.ConstantRing(constants=["t"])
    t = R.gen("t")
    f = t ** 2
    assert f.integrate(t) == Fraction(1, 3) * t ** 3

def test_integral_total():
    A = da.DifferentialRing(functions=["u"])
    u = A.gen("u")
    assert u[1].integrate() == u
    assert (2 * u * u[1]).integrate() == u ** 2


A = da.DifferentialRing(functions=["u", "v"])
u, v = A.gens()
@pytest.mark.parametrize("f", [u, u * u[1] ** 2, u * u[1] ** 2 + u[4], u[1] + v])
def test_integral_does_not_exist(f: da.DifferentialPolynomial):
    with pytest.raises(da.exceptions.NonIntegrableError):
        F = f.integrate()