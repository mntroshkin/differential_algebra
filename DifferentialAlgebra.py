from PolynomialAlgebra import GeneralAlgebra, GeneralMonomial, GeneralPolynomial, Algebra, Polynomial


class DiffAlgebra(GeneralAlgebra):
    _empty_exponent = tuple()
    _default_exponent = (0,)

    def __init__(self, variable_identifiers, coefficient_algebra=None):
        super().__init__(variable_identifiers)
        if not coefficient_algebra:
            coefficient_algebra = self.CoefficientAlgebraType([])
        if not isinstance(coefficient_algebra, self.CoefficientAlgebraType):
            raise TypeError
        self.coefficient_algebra = coefficient_algebra

    def is_coefficient(self, expression):
        return self.coefficient_algebra.is_element(expression)

    def normalize_coefficient(self, coefficient):
        return self.CoefficientType(self.coefficient_algebra, coefficient)

    def normalize_exponent(self, exponent):
        if not isinstance(exponent, tuple):
            raise TypeError
        for diff_degree in exponent:
            if not isinstance(diff_degree, int):
                raise TypeError
            if diff_degree < 0:
                raise ValueError
        return tuple(sorted(exponent))


class DiffMonomial(GeneralMonomial):
    AlgebraType = DiffAlgebra

    def __init__(self, algebra, exponents, coefficient=1):
        super().__init__(algebra, exponents, coefficient)

    def __str__(self):
        factors = []
        if self.coefficient != 1 or self.max_diff_degree() == -1:
            factors.append(self.coefficient.str_parenthesis())
        for var_id in self.algebra.get_all_ids():
            var_name = self.algebra.variables[var_id].representation
            for p in range(0, self.max_diff_degree() + 1):
                p_count = self.exponents[var_id].count(p)
                if p_count > 0:
                    if p == 0:
                        if p_count == 1:
                            factors.append(f"{var_name}")
                        else:
                            factors.append(f"{var_name}^{p_count}")
                    else:
                        if p_count == 1:
                            factors.append(f"{var_name}_{p}")
                        else:
                            factors.append(f"({var_name}_{p})^{p_count}")
        return '*'.join(factors)

    def max_diff_degree(self):
        return max([max(self.exponents[var_id] + (-1,)) for var_id in self.algebra.get_all_ids()])


class DiffPolynomial(GeneralPolynomial):
    def __init__(self, algebra, argument):
        super().__init__(algebra, argument)

    def max_diff_degree(self):
        return max([m.max_diff_degree() for m in self.monomials])

    def diff_x(self, n=1):
        if n == 0:
            return self
        if n == 1:
            summands = []
            for var_id in self.algebra.get_all_ids():
                for i in range(0, self.max_diff_degree() + 1):
                    dvar_i = self.algebra.get_monomial({var_id : (i+1, )})
                    summands.append(dvar_i * self.diff_wrt_variable(var_id, i))
            return sum(summands)
        else:
            return self.diff_x(n - 1).diff_x()

    def diff_wrt_variable(self, var_id, degree):
        summands = []
        for monomial in self.monomials:
            var_count = monomial.exponents[var_id].count(degree)
            if var_count != 0:
                var_index = monomial.exponents[var_id].index(degree)
                new_coefficient = var_count * monomial.coefficient
                var_degrees = monomial.exponents[var_id][:var_index] + monomial.exponents[var_id][var_index + 1:]
                new_factors = {var_id2 : (monomial.exponents[var_id2] if var_id2 != var_id else var_degrees) for var_id2 in self.algebra.get_all_ids()}
                summands.append(DiffMonomial(self.algebra, new_factors, new_coefficient))
        return DiffPolynomial(self.algebra, summands)

    def restrict_to_zero(self):
        summands = []
        for m in self.monomials:
            if sum([m.exponents[var_id].count(0) for var_id in self.algebra.get_all_ids()]) == 0:
                summands.append(m)
        return DiffPolynomial(self.algebra, summands)


DiffAlgebra.CoefficientAlgebraType = Algebra
DiffAlgebra.CoefficientType = Polynomial
DiffAlgebra.MonomialType = DiffMonomial
DiffAlgebra.PolynomialType = DiffPolynomial


class Homomorphism:
    source: DiffAlgebra
    target: DiffAlgebra
    components: dict

    def __init__(self, source, target, components):
        if not isinstance(source, DiffAlgebra) or not isinstance(target, DiffAlgebra):
            raise TypeError
        self.source = source
        self.target = target
        self.components = dict()
        if not isinstance(components, dict):
            raise TypeError
        if sorted(components.keys()) != sorted(source.get_all_ids()):
            raise KeyError
        for var_id in source.get_all_ids():
            if not target.is_element(components[var_id]):
                raise TypeError
            self.components[var_id] = DiffPolynomial(target, components[var_id])

    def apply(self, polynomial):
        if not self.source.is_element(polynomial):
            raise TypeError
        polynomial = DiffPolynomial(self.source, polynomial)
        summands = []
        for monomial in polynomial.monomials:
            product = DiffPolynomial(self.target, monomial.coefficient)
            for var_id in self.source.get_all_ids():
                exponents = monomial.exponents[var_id]
                for diff_degree in exponents:
                    product *= self.components[var_id].diff_x(diff_degree)
            summands.append(product)
        return sum(summands)

    def __mul__(self, other):
        if not isinstance(other, Homomorphism):
            raise TypeError
        if not self.source == other.target:
            raise TypeError
        prod_components = {var_id : self.apply(other.components[var_id]) for var_id in other.source.get_all_ids()}
        return Homomorphism(other.source, self.target, prod_components)


class EvolutionaryOperator:
    def __init__(self, algebra, components):
        if not isinstance(algebra, DiffAlgebra):
            raise TypeError
        self.algebra = algebra
        self.components = dict()
        if not isinstance(components, dict):
            raise TypeError
        if sorted(components.keys()) != sorted(algebra.get_all_ids()):
            raise KeyError
        for var_id in algebra.get_all_ids():
            if not algebra.is_element(components[var_id]):
                raise TypeError
            self.components[var_id] = DiffPolynomial(algebra, components[var_id])

    def apply(self, polynomial) -> DiffPolynomial:
        if not self.algebra.is_element(polynomial):
            raise TypeError
        polynomial = DiffPolynomial(self.algebra, polynomial)
        summands = []
        for var_id in self.algebra.get_all_ids():
            for j in range(polynomial.max_diff_degree() + 1):
                summands.append(self.components[var_id].diff_x(j) * polynomial.diff_wrt_variable(var_id, j))
        return sum(summands)

    def __matmul__(self, other):
        if not isinstance(other, EvolutionaryOperator):
            raise TypeError
        if not self.algebra == other.algebra:
            raise TypeError
        commutator_components = dict()
        for var_id in self.algebra.get_all_ids():
            var = self.algebra.get_variable(var_id)
            commutator_components[var_id] = self.apply(other.apply(var)) - other.apply(self.apply(var))
        return EvolutionaryOperator(self.algebra, commutator_components)