from fractions import Fraction
from typing import Sequence, Optional, NamedTuple

from .exceptions import SymbolNameError
from .exceptions import IncompatibleRingsError, InvalidGeneratorError


type Constant = int | Fraction | ConstantPolynomial
type Rational = int | Fraction

class Monomial(NamedTuple):
    exponents: tuple[int, ...]
    coefficient: Rational


def _multiply_terms(left: Monomial, right: Monomial) -> Monomial:
    exponents = tuple(exp_left + exp_right for exp_left, exp_right in zip(left.exponents, right.exponents))
    coefficient = left.coefficient * right.coefficient
    return Monomial(exponents, coefficient)

def _partial_termwise(term: Monomial, gen_id: int) -> Monomial:
    exponent = term.exponents[gen_id]
    if exponent == 0:
        return Monomial(tuple(0 for var in term.exponents), 0)
    new_exponents = tuple(exp - 1 if i == gen_id else exp for i, exp in enumerate(term.exponents))
    new_coefficient = term.coefficient * exponent
    return Monomial(new_exponents, new_coefficient)

def _integrate_termwise(term: Monomial, gen_id: int) -> Monomial:
    factor = Fraction(1, term.exponents[gen_id] + 1)
    new_exponents = tuple(exp + 1 if i == gen_id else exp for i, exp in enumerate(term.exponents))
    new_coefficient = term.coefficient * factor
    return Monomial(new_exponents, new_coefficient)

def _str_termwise(ring: ConstantRing, term: Monomial) -> str:
    exponents, coefficient = term
    if exponents == tuple([0] * ring._gen_count):
        return str(coefficient)
    
    factors = []
    for i, gen_name in enumerate(ring._gen_names):
        if exponents[i] == 1:
            factors.append(gen_name)
        if exponents[i] > 1:
            factors.append(f"{gen_name}^{exponents[i]}")
    
    if coefficient == 1:
        return "*".join(factors)
    elif coefficient == -1:
        return "-" + "*".join(factors)
    else:
        return str(coefficient) + "*" + "*".join(factors)

def _normalize_terms(terms: Sequence[Monomial]) -> tuple[Monomial, ...]:
    terms = sorted(terms, reverse=True)
    normalized_terms: list[Monomial] = []
    i = 0
    while i < len(terms):
        exponents = terms[i].exponents
        coefficient_sum = 0
        while i < len(terms) and terms[i].exponents == exponents:
            coefficient_sum += terms[i].coefficient
            i += 1
        if coefficient_sum != 0:
            normalized_terms.append(Monomial(exponents, coefficient_sum))
    return tuple(normalized_terms)


class ConstantPolynomial:
    _ring: ConstantRing
    _terms: tuple[Monomial, ...]

    def __init__(self, ring: ConstantRing, terms: Sequence[Monomial]):
        self._ring = ring
        self._terms = _normalize_terms(terms)

    def __add__(self, other):
        if not isinstance(other, int | Fraction | ConstantPolynomial):
            return NotImplemented
        if isinstance(other, ConstantPolynomial) and self._ring != other._ring:
            raise IncompatibleRingsError(f"Cannot add {self} and {other} from incompatible rings {self._ring} and {other._ring}")
        if isinstance(other, int | Fraction):
            other = self._ring.promote(other)
        return ConstantPolynomial(ring=self._ring, terms = self._terms + other._terms)

    def __radd__(self, other):
        return self + other
    
    def __mul__(self, other):
        if not isinstance(other, int | Fraction | ConstantPolynomial):
            return NotImplemented
        if isinstance(other, ConstantPolynomial) and self._ring != other._ring:
            raise IncompatibleRingsError(f"Cannot multiply {self} and {other} from incompatible rings {self._ring} and {other._ring}")
        if isinstance(other, int | Fraction):
            other = self._ring.promote(other)
        product_terms = tuple(_multiply_terms(left, right) for left in self._terms for right in other._terms)
        return ConstantPolynomial(ring=self._ring, terms=product_terms)
        
    def __rmul__(self, other):
        return self * other
    
    def __neg__(self):
        return self * (-1)
    
    def __sub__(self, other):
        return self + other * (-1)
    
    def __rsub__(self, other):
        return self * (-1) + other
    
    def __pow__(self, other) -> ConstantPolynomial:
        if not isinstance(other, int):
            raise TypeError(f"Exponent {other} must be an non-negative integer")
        if other < 0:
            raise ValueError(f"Exponent {other} must be a non-negative integer")
        if other == 0:
            return self._ring.promote(1)
        if other % 2 == 0:
            sqrt = self ** (other // 2)
            return sqrt * sqrt
        else:
            return self * self ** (other - 1)

    def d(self, var: ConstantGenerator, order: int = 1) -> ConstantPolynomial:
        if not isinstance(order, int):
            raise TypeError(f"Derivative order {order} must be an non-negative integer")
        if order < 0:
            raise ValueError(f"Derivative order {order} must be an non-negative integer")
        if not isinstance(var, ConstantGenerator):
            raise InvalidGeneratorError(f"{var} is not a valid ring generator")
        if self._ring != var._ring:
            raise IncompatibleRingsError(f"Generator {var} belongs to ring {var._ring}, "
                                 f"but operation is performed on element {self} of {self._ring}")
        if order == 0:
            return self
        if order == 1:
            return ConstantPolynomial(ring=self._ring,
                                    terms = tuple(_partial_termwise(term, var._gen_id) for term in self._terms))
        else:
            return self.d(var, order - 1).d(var)
    
    def integrate(self, var: ConstantGenerator) -> ConstantPolynomial:
        if not isinstance(var, ConstantGenerator):
            raise InvalidGeneratorError(f"{var} is not a valid ring generator")
        if self._ring != var._ring:
            raise IncompatibleRingsError(f"Generator {var} belongs to ring {var._ring}, "
                                 f"but operation is performed on element {self} of {self._ring}")
        return ConstantPolynomial(ring=self._ring,
                                  terms = tuple(_integrate_termwise(term, var._gen_id) for term in self._terms))

    def __eq__(self, other):
        if self._ring.is_element(other):
            if not isinstance(other, ConstantPolynomial):
                other = self._ring.promote(other)
            return self._terms == other._terms
        return False

    def __str__(self) -> str:
        if len(self._terms) == 0:
            return "0"
        terms_str = [_str_termwise(self._ring, term) for term in self._terms]
        for i in range(1, len(self._terms)):
            if terms_str[i][0] != "-":
                terms_str[i] = "+" + terms_str[i]
        return ''.join(terms_str)
    
    def parenthesis_str(self) -> str:
        return f"({str(self)})" if len(self._terms) > 1 else str(self)
    
    def __lt__(self, other: Constant) -> bool:
        if not isinstance(other, int | Fraction | ConstantPolynomial):
            return NotImplemented
        if isinstance(other, ConstantPolynomial) and self._ring != other._ring:
            raise IncompatibleRingsError(f"Cannot compare {self} and {other} from incompatible rings {self._ring} and {other._ring}")
        other = self._ring.promote(other)
        return self._terms < other._terms
    
    def __hash__(self) -> int:
        return hash((id(self._ring), self._terms))


def _partial_for_const(expression: Constant, var: ConstantGenerator, order: int = 1) -> Constant:
    if not isinstance(order, int):
        raise TypeError(f"Derivative order {order} must be an non-negative integer")
    if order < 0:
        raise ValueError(f"Derivative order {order} must be an non-negative integer")
    if not isinstance(var, ConstantGenerator):
        raise InvalidGeneratorError(f"{var} is not a valid ring generator")
    
    if order == 0:
        return expression
    else:
        if isinstance(expression, (int, Fraction)):
            return 0
        return expression.d(var, order)
    

class ConstantGenerator(ConstantPolynomial):
    _gen_id: int

    def __init__(self, ring: ConstantRing, gen_id: int):
        self._gen_id = gen_id
        exponents = tuple(1 if i == gen_id else 0 for i in range(ring._gen_count))
        super().__init__(ring, terms=[Monomial(exponents, 1)])


class ConstantRing:
    _name: str
    _description: str
    _gen_count: int
    _gen_names: tuple[str]

    def __init__(self, constants: Sequence[str], ring_name: Optional[str] = None):
        gen_names = []
        for gen_name in constants:
            if not gen_name:
                raise SymbolNameError(f"Empty generator names are not allowed")
            if gen_name in gen_names:
                raise SymbolNameError(f"Repeating generator names are not allowed")
            gen_names.append(gen_name)

        self._gen_names = tuple(gen_names)
        self._gen_count = len(self._gen_names)

        ring_description = f"QQ[{', '.join(self._gen_names)}]"

        self._name = str(ring_name)
        self._description_brief = ring_description
        self._description = f"{ring_name} = {ring_description}" if ring_name else ring_description

    def __str__(self) -> str:
        return self._description

    def is_element(self, expression) -> bool:
        if isinstance(expression, (int, Fraction)):
            return True
        elif isinstance(expression, ConstantPolynomial) and expression._ring == self:
            return True
        else:
            return False
        
    def is_generator(self, variable) -> bool:
        if isinstance(variable, ConstantGenerator) and variable._ring == self:
            return True
        else:
            return False

    def promote(self, expression) -> ConstantPolynomial:
        if not self.is_element(expression):
            raise TypeError(f"Expression {expression} is not an element of {self._name} and can't be promoted")
        if isinstance(expression, ConstantPolynomial):
            return expression
        elif isinstance(expression, (int, Fraction)):
            exponents = tuple(0 for gen in self._gen_names)
            return ConstantPolynomial(ring=self, terms=[Monomial(exponents, expression)])
        raise TypeError
    
    def gen(self, gen_name: str) -> ConstantGenerator:
        if not gen_name in self._gen_names:
            raise SymbolNameError(f"{gen_name} is a not a generator of {self}")
        else:
            gen_id = self._gen_names.index(gen_name)
            return ConstantGenerator(self, gen_id)
    
    def gens(self, *gen_names: str) -> tuple[ConstantGenerator, ...]:
        if len(gen_names) == 0:
            gen_names = self._gen_names
        return tuple(self.gen(gen_name) for gen_name in gen_names)

class Rationals(ConstantRing):
    def __init__(self):
        super().__init__(constants=[], ring_name="QQ")
        self._name = "QQ"
        self._description = "QQ"
        self._description_brief = "QQ"

QQ = Rationals()