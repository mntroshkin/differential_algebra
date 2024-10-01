from DifferentialAlgebra import DiffAlgebra, DiffPolynomial
from typing import Self


class GrassmannDiffAlgebra:
    coefficient_algebra: DiffAlgebra
    variables: dict

    def __init__(self, coefficient_algebra, variable_identifiers):
        if not isinstance(coefficient_algebra, DiffAlgebra):
            raise TypeError
        self.coefficient_algebra = coefficient_algebra
        self.variables = dict()
        for var_id in variable_identifiers:
            if var_id in self.variables.keys():
                raise KeyError("Variable identifiers must be unique within each algebra")
            self.variables[var_id] = GrassmannVariable(self, var_id)

    def partitions_normal_form_sgn(self, partitions):
        factors = {var_id : tuple() for var_id in self.variables.keys()}
        total_sign = 1
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
                factor, sign = signed_sort(partitions[var_id])
                total_sign *= sign
                factors[var_id] = tuple(factor)
        if total_sign == 0:
            factors = {var_id: tuple() for var_id in self.variables.keys()}
        return factors, total_sign

    def is_element(self, expression):
        if self.coefficient_algebra.is_element(expression):
            return True
        if isinstance(expression, GrassmannDiffPolynomial) and expression.algebra == self:
            return True
        return False

    def get_variable(self, var_id):
        if not var_id in self.variables.keys():
            raise KeyError("Variable identifier not found in the algebra")
        return self.variables[var_id].to_polynomial()

    def get_monomial(self, partitions):
        return GrassmannDiffMonomial(self, partitions).to_polynomial()

    def get_all_ids(self):
        return tuple([var_id for var_id in self.variables.keys()])


class GrassmannVariable:
    def __init__(self, algebra, identifier, representation=None):
        representation = representation or identifier
        if not isinstance(algebra, GrassmannDiffAlgebra):
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
        return GrassmannDiffMonomial(self.algebra, {self.identifier: partition}).to_polynomial()

    def to_polynomial(self):
        return self.diff_monomial([0])


class GrassmannDiffMonomial:
    algebra: GrassmannDiffAlgebra
    factors: dict
    coefficient: DiffPolynomial

    def __init__(self, algebra, partitions, coefficient=1):
        if not isinstance(algebra, GrassmannDiffAlgebra):
            raise TypeError
        if not algebra.coefficient_algebra.is_element(coefficient):
            raise TypeError
        self.algebra = algebra
        self.coefficient = 1
        self.factors = {var_id : tuple() for var_id in self.algebra.variables.keys()}
        if isinstance(partitions, dict):
            factors, sign = self.algebra.partitions_normal_form_sgn(partitions)
            coefficient *= sign
            if coefficient != 0:
                self.factors = factors
            self.coefficient = DiffPolynomial(self.algebra.coefficient_algebra, coefficient)
        elif self.algebra.coefficient_algebra.is_element(partitions):
            self.coefficient = DiffPolynomial(self.algebra.coefficient_algebra, partitions * coefficient)
        else:
            raise TypeError

    def factors_tuple(self):
        return tuple(zip(self.factors.keys(), self.factors.values()))

    def __str__(self):
        if self.coefficient != 1 or self.factors == -1:
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
                        degree_str.append(f"{var_name}")
                    else:
                        degree_str.append(f"{var_name}_{p}")
        return coef_str + '*'.join(degree_str)

    def to_polynomial(self):
        return GrassmannDiffPolynomial(self.algebra, [self])

    def max_degree(self):
        return max([max(self.factors[var_id] + (-1,)) for var_id in self.algebra.get_all_ids()])


class GrassmannDiffPolynomial:
    algebra: GrassmannDiffAlgebra
    monomials: list[GrassmannDiffMonomial]

    def __init__(self, algebra, argument):
        if not isinstance(algebra, GrassmannDiffAlgebra):
            raise TypeError
        self.algebra = algebra
        if isinstance(argument, list):
            self.monomials = []
            for monomial in argument:
                if not isinstance(monomial, GrassmannDiffMonomial):
                    raise TypeError
                if monomial.algebra != algebra:
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
                    self.monomials.append(GrassmannDiffMonomial(algebra, factors_dict, coefficient_sum))
        elif isinstance(argument, GrassmannDiffPolynomial):
            if self.algebra != algebra:
                raise TypeError
            self.monomials = argument.monomials[:]
        elif algebra.coefficient_algebra.is_element(argument):
            if argument != 0:
                self.monomials = [GrassmannDiffMonomial(algebra, argument)]
            else:
                self.monomials = []
        else:
            raise TypeError

    def __str__(self):
        if len(self.monomials) == 0:
            return '0'
        else:
            return ' + '.join(map(str, self.monomials))

    def __add__(self, other) -> Self:
        if not self.algebra.is_element(other):
            return NotImplemented
        other = GrassmannDiffPolynomial(self.algebra, other)
        return GrassmannDiffPolynomial(self.algebra, self.monomials + other.monomials)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        return self + other * (-1)

    def __rsub__(self, other):
        return self * (-1) + other

    def __mul__(self, other) -> Self:
        if not self.algebra.is_element(other):
            return NotImplemented
        other = GrassmannDiffPolynomial(self.algebra, other)
        product_monomials = []
        for m1 in self.monomials:
            for m2 in other.monomials:
                product_factors = {var_id: m1.factors[var_id] + m2.factors[var_id] for var_id in
                                   self.algebra.get_all_ids()}
                product_coefficient = m1.coefficient * m2.coefficient
                product_monomial = GrassmannDiffMonomial(self.algebra, product_factors, product_coefficient)
                product_monomials.append(product_monomial)
        return GrassmannDiffPolynomial(self.algebra, product_monomials)

    def __rmul__(self, other):
        return self * other

    def max_degree(self):
        return max([m.max_degree() for m in self.monomials])

    def diff_x(self, n=1) -> Self:
        if n == 0:
            return self
        if n == 1:
            summands = []
            for var_id in self.algebra.coefficient_algebra.get_all_ids():
                for monomial in self.monomials:
                    coefficient = monomial.coefficient.diff_x()
                    factors = monomial.factors
                    new_monomial = GrassmannDiffMonomial(self.algebra, factors, coefficient)
                    summands.append(GrassmannDiffPolynomial(self.algebra, [new_monomial]))
            for var_id in self.algebra.get_all_ids():
                for i in range(0, self.max_degree() + 1):
                    dvar_i = self.algebra.get_monomial({var_id: (i + 1,)})
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
                new_coefficient = var_count * monomial.coefficient * (-1) ** var_index
                var_degrees = monomial.factors[var_id][:var_index] + monomial.factors[var_id][var_index + 1:]
                new_factors = {var_id2: (monomial.factors[var_id2] if var_id2 != var_id else var_degrees) for var_id2 in
                               self.algebra.get_all_ids()}
                summands.append(GrassmannDiffMonomial(self.algebra, new_factors, new_coefficient))
        return GrassmannDiffPolynomial(self.algebra, summands)


def signed_sort(lst):
    lst = list(lst)
    sign = 1
    for i in range(len(lst)):
        for j in range(i + 1, len(lst)):
            if lst[i] > lst[j]:
                sign *= -1
                lst[i], lst[j] = lst[j], lst[i]
    for i in range(len(lst) - 1):
        if lst[i] == lst[i + 1]:
            sign = 0
    return tuple(lst), sign
