from fractions import Fraction
from typing import Self


class DiffAlgebra:
    def __init__(self, variable_identifiers):
        self.variables = dict()
        for var_id in variable_identifiers:
            if var_id in self.variables.keys():
                raise KeyError("Variable identifiers must be unique within each algebra")
            self.variables[var_id] = Variable(self, var_id)

    def is_element(self, expression):
        if isinstance(expression, int) or isinstance(expression, Fraction):
            return True
        if isinstance(expression, DiffPolynomial) and expression.algebra == self:
            return True
        return False

    def partitions_normal_form(self, partitions):
        factors = {var_id : tuple() for var_id in self.variables.keys()}
        if isinstance(partitions, dict):
            for var_id in partitions.keys():
                if not var_id in self.variables.keys():
                    raise KeyError("Variable identifier not found in the algebra")
                if not (isinstance(partitions[var_id], list) or isinstance(partitions[var_id], tuple)):
                    raise TypeError
                for degree in partitions[var_id]:
                    if not isinstance(degree, int):
                        raise TypeError
                    if degree < 0:
                        raise ValueError
                factors[var_id] = tuple(sorted(partitions[var_id]))
        return factors

    def get_variable(self, var_id):
        if not var_id in self.variables.keys():
            raise KeyError("Variable identifier not found in the algebra")
        return self.variables[var_id].to_polynomial()

    def get_monomial(self, partitions):
        return DiffMonomial(self, partitions).to_polynomial()

    def get_all_ids(self):
        return tuple([var_id for var_id in self.variables.keys()])


class Variable:
    def __init__(self, algebra, identifier, representation=None):
        representation = representation or identifier
        if not isinstance(algebra, DiffAlgebra):
            raise TypeError
        if not isinstance(identifier, str) or not isinstance(representation, str):
            raise TypeError
        if identifier == "" or representation == "":
            raise ValueError
        self.algebra = algebra
        self.identifier = identifier
        self.representation = representation or identifier

    def set_representation(self, new_rep):
        if not isinstance(new_rep, str):
            raise TypeError
        if not new_rep:
            raise ValueError("Representation string must be non-empty.")
        self.representation = new_rep

    def diff_monomial(self, partition):
        return DiffMonomial(self.algebra, {self.identifier : partition}).to_polynomial()

    def to_polynomial(self):
        return self.diff_monomial([0])


class DiffMonomial:
    algebra: DiffAlgebra
    factors: dict
    coefficient: int | Fraction

    def __init__(self, algebra, partitions, coefficient=1):
        self.algebra = algebra
        self.factors = {var_id : tuple() for var_id in self.algebra.variables.keys()}
        if isinstance(partitions, dict):
            self.factors = self.algebra.partitions_normal_form(partitions)
            self.coefficient = Fraction(coefficient)
        elif isinstance(partitions, int) or isinstance(partitions, Fraction):
            self.coefficient = Fraction(partitions * coefficient)
        else:
            raise TypeError

    def factors_tuple(self):
        return tuple(zip(self.factors.keys(), self.factors.values()))

    def __str__(self):
        if self.coefficient != 1 or self.max_degree() == -1:
            coef_str = str(self.coefficient) + '*'
        else:
            coef_str = ''
        degree_str = []
        for var_id in self.algebra.get_all_ids():
            var_name = self.algebra.variables[var_id].representation
            for p in range(0, self.max_degree() + 1):
                p_count = self.factors[var_id].count(p)
                if p_count > 0:
                    if p == 0:
                        if p_count == 1:
                            degree_str.append(f"{var_name}")
                        else:
                            degree_str.append(f"{var_name}^{p_count}")
                    else:
                        if p_count == 1:
                            degree_str.append(f"{var_name}_{p}")
                        else:
                            degree_str.append(f"({var_name}_{p})^{p_count}")
        return coef_str + '*'.join(degree_str)

    def to_polynomial(self):
        return DiffPolynomial(self.algebra, self)

    def max_degree(self):
        return max([max(self.factors[var_id] + (-1,)) for var_id in self.algebra.get_all_ids()])


class DiffPolynomial:
    algebra: DiffAlgebra
    monomials: list[DiffMonomial]

    def __init__(self, algebra, argument):
        self.algebra = algebra
        if isinstance(argument, list):
            self.monomials = []
            for monomial in argument:
                if not isinstance(monomial, DiffMonomial):
                    raise TypeError
            argument.sort(key=lambda monomial: monomial.factors_tuple())
            i = 0
            while i < len(argument):
                factors_tuple = argument[i].factors_tuple()
                factors_dict = argument[i].factors
                coefficient_sum = 0
                while i < len(argument) and argument[i].factors_tuple() == factors_tuple:
                    coefficient_sum += argument[i].coefficient
                    i += 1
                if coefficient_sum != 0:
                    self.monomials.append(DiffMonomial(algebra, factors_dict, coefficient_sum))
        elif isinstance(argument, DiffPolynomial):
            self.monomials = argument.monomials[:]
        elif isinstance(argument, DiffMonomial):
            self.monomials = [argument]
        elif isinstance(argument, Variable):
            self.monomials = [argument.diff_monomial([0])]
        elif isinstance(argument, Fraction) or isinstance(argument, int):
            self.monomials = [DiffMonomial(algebra, argument)]

    def __str__(self):
        if len(self.monomials) == 0:
            return '0'
        else:
            return ' + '.join(map(str, self.monomials))

    def __add__(self, other) -> Self:
        if not self.algebra.is_element(other):
            return NotImplemented
        other = DiffPolynomial(self.algebra, other)
        return DiffPolynomial(self.algebra, self.monomials + other.monomials)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + other * (-1)

    def __rsub__(self, other):
        return self * (-1) + other

    def __mul__(self, other) -> Self:
        if not self.algebra.is_element(other):
            return NotImplemented
        other = DiffPolynomial(self.algebra, other)
        product_monomials = []
        for m1 in self.monomials:
            for m2 in other.monomials:
                product_factors = {var_id: m1.factors[var_id] + m2.factors[var_id] for var_id in self.algebra.get_all_ids()}
                product_coefficient = m1.coefficient * m2.coefficient
                product_monomial = DiffMonomial(self.algebra, product_factors, product_coefficient)
                product_monomials.append(product_monomial)
        return DiffPolynomial(self.algebra, product_monomials)

    def __rmul__(self, other):
        return self * other

    def max_degree(self):
        return max([m.max_degree() for m in self.monomials])

    def diff_x(self, n=1) -> Self:
        if n == 0:
            return self
        if n == 1:
            summands = []
            for var_id in self.algebra.get_all_ids():
                for i in range(0, self.max_degree() + 1):
                    dvar_i = self.algebra.get_monomial({var_id : (i+1, )})
                    summands.append(dvar_i * self.diff_wrt_variable(var_id, i))
            return sum(summands)
        else:
            return self.diff_x(n - 1).diff_x()

    def diff_wrt_variable(self, var_id, degree) -> Self:
        summands = []
        for monomial in self.monomials:
            var_count = monomial.factors[var_id].count(degree)
            if var_count != 0:
                var_index = monomial.factors[var_id].index(degree)
                new_coefficient = var_count * monomial.coefficient
                var_degrees = monomial.factors[var_id][:var_index] + monomial.factors[var_id][var_index + 1:]
                new_factors = {var_id2 : (monomial.factors[var_id2] if var_id2 != var_id else var_degrees) for var_id2 in self.algebra.get_all_ids()}
                summands.append(DiffMonomial(self.algebra, new_factors, new_coefficient))
        return DiffPolynomial(self.algebra, summands)

    def evol_field(self, field: list[Self]) -> Self:
        summands = []
        for (var_id, component) in zip(self.algebra.get_all_ids(), field):
            summand = sum([component.diff_x(j) * self.diff_wrt_variable(var_id, j) for j in range(self.max_degree() + 1)])
            summands.append(summand)
        return sum(summands)

    def restrict_to_zero(self) -> Self:
        summands = []
        for m in self.monomials:
            if sum([m.factors[var_id].count(0) for var_id in self.algebra.get_all_ids()]) == 0:
                summands.append(m)
        return DiffPolynomial(self.algebra, summands)

    def coefficient(self, factors) -> Fraction:
        factors = self.algebra.partitions_normal_form(factors)
        for m in self.monomials:
            if m.factors == factors:
                return m.coefficient
        return Fraction(0)

class Mapping:
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
        polynomial = DiffPolynomial(self.target, polynomial)
        summands = []
        for monomial in polynomial.monomials:
            product = DiffPolynomial(self.target, monomial.coefficient)
            for var_id in self.source.get_all_ids():
                degrees = monomial.factors[var_id]
                for part in degrees:
                    product *= self.components[var_id].diff_x(part)
            summands.append(product)
        return sum(summands)

    def __mul__(self, other):
        if not isinstance(other, Mapping):
            raise TypeError
        if not self.source == other.target:
            raise TypeError
        prod_components = {var_id : self.apply(other.components[var_id]) for var_id in other.source.get_all_ids()}
        return Mapping(other.source, self.target, prod_components)


class Derivation:
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
            for j in range(polynomial.max_degree() + 1):
                summands.append(self.components[var_id].diff_x(j) * polynomial.diff_wrt_variable(var_id, j))
        return sum(summands)

    def __matmul__(self, other: Self) -> Self:
        if not isinstance(other, Derivation):
            raise TypeError
        if not self.algebra == other.algebra:
            raise TypeError
        commutator_components = dict()
        for var_id in self.algebra.get_all_ids():
            var = self.algebra.get_variable(var_id)
            commutator_components[var_id] = self.apply(other.apply(var)) - other.apply(self.apply(var))
        return Derivation(self.algebra, commutator_components)
