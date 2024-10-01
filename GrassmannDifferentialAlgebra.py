from DifferentialAlgebra import DiffAlgebra, DiffPolynomial


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
                factors[var_id] = factor
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
            self.factors = factors
            self.coefficient = DiffPolynomial(self.algebra.coefficient_algebra, coefficient * sign)
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
        return DiffPolynomial(self.algebra, self)

    def max_degree(self):
        return max([max(self.factors[var_id] + (-1,)) for var_id in self.algebra.get_all_ids()])


class GrassmannDiffPolynomial:
    pass


def signed_sort(lst):
    sign = 1
    for i in range(len(lst)):
        for j in range(i + 1, len(lst)):
            if lst[i] > lst[j]:
                sign *= -1
                lst[i], lst[j] = lst[j], lst[i]
    for i in range(len(lst) - 1):
        if lst[i] == lst[i + 1]:
            sign = 0
    return lst, sign
