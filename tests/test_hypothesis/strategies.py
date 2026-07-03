from hypothesis import strategies as st
import diffalgebra as da
from diffalgebra.constant_ring import Monomial
from diffalgebra.diff_ring import DiffMonomial, DiffFactor


@st.composite
def polynomial(draw, ring: da.ConstantRing, max_terms: int = 5) -> da.ConstantPolynomial:
    gen_count = ring._gen_count

    terms = draw(st.lists(st.builds(Monomial, 
                                    exponents=st.lists(st.integers(min_value=0, max_value=5),
                                                       min_size=gen_count, max_size=gen_count).map(tuple), 
                                    coefficient=st.integers(min_value=-5, max_value=5)), 
                            max_size=max_terms))
    return da.ConstantPolynomial(ring, terms)


@st.composite
def diff_polynomial(draw, ring: da.DifferentialRing, 
                    max_terms: int = 5, max_nonlinearity: int = 2) -> da.DifferentialPolynomial:
    gens = ring.gens()
    gen_count = len(gens)

    terms = draw(st.lists(st.builds(DiffMonomial, factors=st.lists(st.builds(DiffFactor, 
                                                                             gen_id=st.integers(min_value=0, max_value=gen_count - 1),
                                                                             derivative=st.integers(min_value=0, max_value=3),
                                                                             power=st.integers(min_value=0, max_value=3)),
                                                                    max_size=max_nonlinearity).map(tuple),
                                                    coefficient=st.integers(min_value=-5, max_value=5)),
                        max_size=max_terms))
    return da.DifferentialPolynomial(ring, terms)
