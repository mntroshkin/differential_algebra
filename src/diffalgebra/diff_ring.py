from typing import Optional, Sequence, NamedTuple
from fractions import Fraction

from .exceptions import RingMismatchError, SymbolNameError
from .exceptions import InvalidGeneratorError, IncompatibleRingsError
from .constant_ring import ConstantRing, ConstantPolynomial, ConstantGenerator, QQ, Constant, _partial_for_const


type Expression = int | Fraction | ConstantPolynomial | DifferentialPolynomial
type AnyGenerator = ConstantGenerator | FuncGeneratorDerivative

class DiffFactor(NamedTuple):
    gen_id: int
    derivative: int
    power: int

class DiffMonomial(NamedTuple):
    factors: tuple[DiffFactor, ...]
    coefficient: Constant


def _normalize_factors(factors: Sequence[DiffFactor]) -> tuple[DiffFactor, ...]:
    factors = sorted(factors)
    normalized_factors: list[DiffFactor] = []
    i = 0
    while i < len(factors):
        gen_id, derivative, power = factors[i]
        total_power = 0
        while i < len(factors) and factors[i].gen_id == gen_id and factors[i].derivative == derivative:
            total_power += factors[i].power
            i += 1
        if total_power != 0:
            normalized_factors.append(DiffFactor(gen_id, derivative, total_power))
    return tuple(normalized_factors)


def _normalize_monomial(monomial: DiffMonomial) -> DiffMonomial:
    return DiffMonomial(_normalize_factors(monomial.factors), monomial.coefficient)


def _differential_degree(monomial: DiffMonomial) -> int:
    degree = 0
    for gen_id, derivative, power in monomial.factors:
        degree += derivative * power
    return degree


def _normalize_terms(terms: Sequence[DiffMonomial]) -> tuple[DiffMonomial, ...]:
    terms = sorted(map(_normalize_monomial, terms))
    terms.sort(key=_differential_degree)
    normalized_terms: list[DiffMonomial] = []
    i = 0
    while i < len(terms):
        factors = terms[i].factors
        coefficient_sum = 0
        while i < len(terms) and terms[i].factors == factors:
            coefficient_sum += terms[i].coefficient
            i += 1
        if coefficient_sum != 0:
            normalized_terms.append(DiffMonomial(factors, coefficient_sum))
    return tuple(normalized_terms)


def _multiply_terms(left: DiffMonomial, right: DiffMonomial) -> DiffMonomial:
    factors = left.factors + right.factors
    coefficient = left.coefficient * right.coefficient
    return DiffMonomial(factors, coefficient)


def _diff_termwise(term: DiffMonomial) -> list[DiffMonomial]:
    diff_terms: list[DiffMonomial] = []
    for gen_id, derivative, power in term.factors:
        extra_factors = (DiffFactor(gen_id, derivative, - 1), DiffFactor(gen_id, derivative + 1, 1))
        factors = term.factors + extra_factors
        diff_terms.append(DiffMonomial(factors, term.coefficient * power))
    return diff_terms


def _exponent_in_term(term: DiffMonomial, gen_id: int, derivative: int) -> int:
    for factor in term.factors:
        if factor.gen_id == gen_id and factor.derivative == derivative:
            return factor.power
    return 0


def _partial_termwise(term: DiffMonomial, gen_id: int, derivative: int) -> DiffMonomial:
    for factor in term.factors:
        if factor.gen_id == gen_id and factor.derivative == derivative:
            extra_factor = (DiffFactor(gen_id, derivative, -1), )
            factors = term.factors + extra_factor
            return DiffMonomial(factors, term.coefficient * factor.power)
    return DiffMonomial(tuple(), 0)
    

def _str_termwise(ring: DifferentialRing, term: DiffMonomial) -> str:
    if term.factors == tuple():
        return str(term.coefficient)
    
    str_factors = []
    for gen_id, derivative, power in term.factors:
        func_name = ring._gen_names[gen_id]
        if power == 1:
            if derivative == 0:
                str_factors.append(func_name)
            elif derivative <= 3:
                str_factors.append(f"{func_name}_{"x" * derivative}")
            else:
                str_factors.append(f"{func_name}_{derivative}")
        else:
            if derivative == 0:
                str_factors.append(f"{func_name}^{power}")
            elif derivative <= 3:
                str_factors.append(f"({func_name}_{"x" * derivative})^{power}")
            else:
                str_factors.append(f"({func_name}_{derivative})^{power}")
    
    if term.coefficient == 1:
        return "*".join(str_factors)
    elif term.coefficient == -1:
        return "-" + "*".join(str_factors)
    else:
        if isinstance(term.coefficient, ConstantPolynomial):
            coefficient_str = term.coefficient.parenthesis_str()
        else:
            coefficient_str = str(term.coefficient)
        return coefficient_str + "*" + "*".join(str_factors)
    

class DifferentialPolynomial:
    _ring: DifferentialRing
    _terms: tuple[DiffMonomial, ...]

    def __init__(self, ring: DifferentialRing, terms: Sequence[DiffMonomial]):
        self._ring = ring
        self._terms = _normalize_terms(terms)

    def __add__(self, other):
        if isinstance(other, DifferentialPolynomial) and other._ring != self._ring:
            raise RingMismatchError
        if isinstance(other, ConstantPolynomial) and other._ring != self._ring._base_ring:
            raise RingMismatchError
        if self._ring.is_element(other):
            if not isinstance(other, DifferentialPolynomial):
                other = self._ring.promote(other)
            return DifferentialPolynomial(self._ring, terms=self._terms + other._terms)
        else:
            return NotImplemented
        
    def __radd__(self, other):
        return self + other
    
    def __mul__(self, other):
        if isinstance(other, DifferentialPolynomial) and other._ring != self._ring:
            raise RingMismatchError
        if isinstance(other, ConstantPolynomial) and other._ring != self._ring._base_ring:
            raise RingMismatchError
        if self._ring.is_element(other):
            if not isinstance(other, DifferentialPolynomial):
                other = self._ring.promote(other)
            product_terms = [_multiply_terms(left, right) for left in self._terms for right in other._terms]
            return DifferentialPolynomial(ring=self._ring, terms=product_terms)     
        else:
            return NotImplemented
        
    def __rmul__(self, other):
        return self * other
    
    def __neg__(self):
        return self * (-1)
    
    def __sub__(self, other):
        return self + other * (-1)
    
    def __rsub__(self, other):
        return self * (-1) + other
    
    def __pow__(self, other: int) -> DifferentialPolynomial:
        if not isinstance(other, int):
            raise TypeError
        if other < 0:
            raise ValueError
        if other == 0:
            return self._ring.promote(1)
        if other % 2 == 0:
            sqrt = self ** (other // 2)
            return sqrt * sqrt
        else:
            return self * self ** (other - 1)
            
    def __eq__(self, other):
        if self._ring.is_element(other):
            if not isinstance(other, DifferentialPolynomial):
                other = self._ring.promote(other)
            return self._terms == other._terms
        return False
    
    def __hash__(self) -> int:
        return hash((id(self._ring), self._terms))

    def __str__(self) -> str:
        if len(self._terms) == 0:
            return "0"
        terms_str = [_str_termwise(self._ring, term) for term in self._terms]
        for i in range(1, len(self._terms)):
            if terms_str[i][0] != "-":
                terms_str[i] = "+" + terms_str[i]
        return ''.join(terms_str)

    def __repr__(self) -> str:
        return str(self)
    
    def diff(self, order: int = 1) -> DifferentialPolynomial:
        if order < 0:
            raise ValueError
        if order == 0:
            return self
        if order == 1:
            return DifferentialPolynomial(ring=self._ring, 
                                          terms=[diff_term 
                                                 for term in self._terms 
                                                 for diff_term in _diff_termwise(term)])
        else:
            return self.diff(order - 1).diff()
        
    def d(self, gen: AnyGenerator, order: int = 1) -> DifferentialPolynomial:
        if not isinstance(order, int):
            raise TypeError(f"Derivative order {order} must be an non-negative integer")
        if order < 0:
            raise ValueError(f"Derivative order {order} must be an non-negative integer")
        if not isinstance(gen, ConstantGenerator | FuncGeneratorDerivative):
            raise InvalidGeneratorError(f"{gen} is not a valid ring generator")
        if order == 0:
            return self
        if order == 1:
            if isinstance(gen, ConstantGenerator):
                if self._ring._base_ring != gen._ring:
                    raise IncompatibleRingsError(f"Generator {gen} belongs to ring {gen._ring}, "
                                                 f"but operation is performed on element {self} of {self._ring}")
                new_terms = [DiffMonomial(factors, _partial_for_const(coefficient, gen)) 
                             for (factors, coefficient) in self._terms]
                return DifferentialPolynomial(ring=self._ring, terms=new_terms)
            if isinstance(gen, FuncGeneratorDerivative):
                if self._ring != gen._ring:
                    raise IncompatibleRingsError(f"Generator {gen} belongs to ring {gen._ring}, "
                                                 f"but operation is performed on element {self} of {self._ring}")
                gen_id = gen._gen_id
                derivative_order = gen._derivative
                return DifferentialPolynomial(ring=self._ring,
                                    terms = [_partial_termwise(term, gen_id, derivative_order) 
                                             for term in self._terms])
        return self.d(gen, order - 1).d(gen)

    def _highest_derivative(self, gen: FuncGenerator) -> int:
        highest_derivative = -1
        for term in self._terms:
            factors = filter(lambda factor: factor.gen_id == gen._gen_id, term.factors)
            highest_for_term = max((factor.derivative for factor in factors), default=-1)
            highest_derivative = max(highest_derivative, highest_for_term)
        return highest_derivative
    
    def delta(self, gen: FuncGenerator) -> DifferentialPolynomial:
        if not isinstance(gen, FuncGenerator):
            raise InvalidGeneratorError(f"{gen} is not a valid ring generator")
        if not self._ring == gen._ring:
            raise IncompatibleRingsError(f"Generator {gen} belongs to ring {gen._ring}, "
                                        f"but operation is performed on element {self} of {self._ring}")
        d = self._highest_derivative(gen)
        result_terms: list[DiffMonomial] = []
        for i in range(d + 1):
            summand = (-1) ** i * self.d(gen.diff(i)).diff(order=i)
            result_terms.extend(self._ring.promote(summand)._terms)
        return DifferentialPolynomial(ring=self._ring, terms=result_terms)

    def coefficient(self, monomial: Expression) -> Constant:
        monomial = self._ring.promote(monomial)
        if len(monomial._terms) != 1:
            raise ValueError
        given_term = monomial._terms[0]
        for term in self._terms:
            if term.factors == given_term.factors:
                return term.coefficient
        return 0
    
    def integral(self) -> DifferentialPolynomial | None:
        if self == 0:
            return self._ring.promote(0)
        if self.coefficient(1) != 0:
            return
        for var_id, var in enumerate(self._ring.gens()):
            if self._highest_derivative(var) != -1:
                new_integrand = self._ring.promote(0)
                integral_found = self._ring.promote(0)
                k = self._highest_derivative(var)
                if k == 0:
                    return
                for term in self._terms:
                    monomial, coefficient = term
                    j = _exponent_in_term(term, var_id, k)
                    if j == 0:
                        new_integrand += DifferentialPolynomial(self._ring, terms=[term])
                    elif j > 1:
                        return
                    else:
                        j = _exponent_in_term(term, var_id, k - 1)
                        new_factor = DiffFactor(var_id, k - 1, j + 1)
                        const = Fraction(1, j + 1)
                        filtered_factors = tuple(filter(lambda factor: factor.gen_id != var_id or factor.derivative < k - 1, monomial))
                        h = DifferentialPolynomial(self._ring, terms=[DiffMonomial(filtered_factors, coefficient)])
                        new_factor = DifferentialPolynomial(self._ring, terms=[DiffMonomial(factors=(new_factor, ), coefficient=const)])
                        integral_found += h * new_factor
                        new_integrand -= h.diff() * new_factor
                integral_remaining = new_integrand.integral()
                if integral_remaining is None:
                    return
                return integral_found + integral_remaining
                        

def total_derivative(expression: Expression, order: int = 1) -> Expression:
    if order < 0:
        raise ValueError
    if order == 0:
        return expression

    if isinstance(expression, (int, Fraction, ConstantPolynomial)):
        return 0
    else:
        return expression.diff(order)

def partial_derivative(expression: Expression,
                       var: FuncGenerator,
                       order: int = 1) -> Expression:
    if isinstance(expression, (int, Fraction)):
        return 0
    if isinstance(expression, ConstantPolynomial) and isinstance(var, ConstantGenerator):
        return expression.d(var, order)
    if isinstance(expression, DifferentialPolynomial) and isinstance(var, (ConstantGenerator, FuncGeneratorDerivative)):
        return expression.d(var, order)
    raise TypeError


class FuncGeneratorDerivative(DifferentialPolynomial):
    _gen_id: int
    _derivative: int

    def __init__(self, ring: DifferentialRing, gen_id: int, derivative: int = 0):
        factors = (DiffFactor(gen_id, derivative, 1), )
        super().__init__(ring, terms=[DiffMonomial(factors, 1)])

        self._gen_id = gen_id
        self._derivative = derivative
    
    def diff(self, order: int = 1) -> FuncGeneratorDerivative:
        if not isinstance(order, int):
            raise TypeError("Derivative order must be a non-negative integers")
        if order < 0:
            raise ValueError("Derivative order must be a non-negative integers")
        else:
            return FuncGeneratorDerivative(self._ring, self._gen_id, self._derivative + order)

    def __getitem__(self, key) -> FuncGeneratorDerivative:
        return self.diff(order=key)
    

class FuncGenerator(FuncGeneratorDerivative):
    _gen_id: int
    _derivative = 0
    
    def __init__(self, ring: DifferentialRing, gen_id: int):
        super().__init__(ring, gen_id, derivative=0)


class DifferentialRing:
    _name: str
    _base_ring: ConstantRing
    _gen_names: tuple[str]
    _gen_count: int

    def __init__(self, functions: Sequence[str], base_ring: ConstantRing = QQ, ring_name: Optional[str] = None):
        self._base_ring = base_ring
        
        gen_names = []
        for gen_name in functions:
            if not gen_name:
                raise SymbolNameError(f"Empty function names are not allowed")
            if gen_name in gen_names:
                raise SymbolNameError(f"Repeating function names are not allowed")
            if gen_name in self._base_ring._gen_names:
                raise SymbolNameError(f"Function name coincides with a constant name")
            gen_names.append(gen_name)

        self._gen_names = tuple(gen_names)
        self._gen_count = len(self._gen_names)
        
        ring_definition = f"{self._base_ring._description_brief}{{{ ', '.join(self._gen_names) }}}"

        self._name = str(ring_name)
        self._description_brief = ring_definition
        self._description = f"{ring_name} = {ring_definition}" if ring_name else ring_definition

    def __str__(self) -> str:
        return self._description

    def __repr__(self) -> str:
        return f"{str(self)}: differential ring with id={id(self)}" 

    def is_element(self, expression) -> bool:
        if self._base_ring.is_element(expression):
            return True
        if isinstance(expression, DifferentialPolynomial) and expression._ring == self:
            return True
        return False
    
    def promote(self, expression) -> DifferentialPolynomial:
        if not self.is_element(expression):
            raise TypeError(f"Expression {expression} is not an element of {self._name} and can't be promoted")
        if isinstance(expression, DifferentialPolynomial):
            return expression
        if isinstance(expression, (int, Fraction, ConstantPolynomial)):
            return DifferentialPolynomial(ring=self, terms=[DiffMonomial(tuple(), expression)])
        raise TypeError
    
    def constants(self, *const_names: str) -> tuple[ConstantGenerator, ...]:
        return self._base_ring.gens(*const_names)
    
    def gen(self, gen_name: str) -> FuncGenerator:
        if not gen_name in self._gen_names:
            raise SymbolNameError(f"{gen_name} is a not a generator of {self}")
        else:
            gen_id = self._gen_names.index(gen_name)
            return FuncGenerator(self, gen_id)
    
    def gens(self, *gen_names: str) -> tuple[FuncGenerator, ...]:
        if len(gen_names) == 0:
            gen_names = self._gen_names
        return tuple(self.gen(gen_name) for gen_name in gen_names)
