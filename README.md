# `differential_algebra` package  

`differential_algebra` is a symbolic computer algebra package created to manipulate expressions in both commutative and anticommutative free differential algebras. It is specifically designed to explore problems in the theory of integrable hierarchies of partial differential equations (PDEs) and supports computations for both Hamiltonian and non-Hamiltonian evolutionary PDEs in one spatial variable.

## Mathematical background  

A *differential* $R$-*algebra* $(A, \partial)$ is a (graded) commutative algebra $A$ over some commutative ring $R$ of constants, equipped with an $R$-linear derivation $\partial: A \to A$. In our package, $R$ is chosen to be the polynomial ring $R = \mathbb{Q}[x_1, \dots, x_n]$, representing constants in the algebra.

A *homomorphism* of differential $R$-algebras $f: (A, \partial_A) \to (B, \partial_B)$ is a homomorphism of the underlying $R$-algebras such that it commutes with the derivations: $f \circ \partial_A = \partial_B \circ f$.

### Free Differential Algebras

A *free differential* $R$-*algebra* $\mathrm{fda}(R, U)$ over $R$ with a set of generators $U = (u^1, \dots, u^N)$ is defined by the following universal property: for every differential $R$-algebra $(A, \partial)$ and for every set map $\phi: U \to A$, there exists a unique differential $R$-algebra homomorphism $\mathrm{fda}(U) \to A$ extending $\phi$. 

We treat the generators $u^i$ as functions of a spatial variable $x$, and the derivation $\partial$ is understood as the spatial derivative $\partial_x := \frac{\partial}{\partial x}$. 

Concretely, the free differential algebra $\mathrm{fda}(R, U)$ is isomorphic to a polynomial algebra $R[u^i_j]$, where $i = 1, \dots, N$ and $j \geq 0$, and the symbol $u^i_j$ represents the $j$-th derivative of $u^i$ with respect to $x$, i.e., $u^i_j := \frac{\partial^j u^i}{\partial x^j}$. The derivation acts on these symbols as: $\partial u^i_j = u^i_{j+1}$

### Evolutionary Operators

An *evolutionary operator* on a differential $R$-algebra $(A, \partial)$ is an $R$-linear derivation $D: A \to A$ which commutes with $\partial$, i.e., $D(\partial a) = \partial(Da)$ for all $a \in A$. The evolutionary operator represents time-evolution in a system of PDEs and is uniquely determined by its action on the generators $u^i$:

$$
D(u^i) = F^i(u, \partial u, \dots)
$$

where $F^i$ is a polynomial in generators $u^j$ and their derivatives. The operator $D$ algebraically corresponds to the differential equation:

$$
\frac{\partial u^i}{\partial t} = D(u^i)
$$

The commutator of two evolutionary operators $[D, D'] := D \circ D' - D' \circ D$ is itself an evolutionary operator. Commutativity of two evolutionary operators $[D, D'] = 0$ implies the corresponding systems of PDEs are compatible, meaning they admit a common solution.

## Usage example

### 1. Importing classes from the submodules

```python
from PolynomialAlgebra import Algebra
from DifferentialAlgebra import DiffAlgebra, EvolutionaryOperator
```

### 2. Creating algebras and obtaining variables

We first create a commutative polynomial algebra $R$ over $\mathbb{Q}$ in one variable `t`, and a differential algebra $A$ over $R$ with two generators `u` and `v`:

```python
R = Algebra(["t"])  
A = DiffAlgebra(["u", "v"], R)
```

Next, we obtain the variables `t`, `u`, and `v` as algebraic expressions from their respective algebra objects:

```python
t = R.get_variable("t")  
u = A.get_variable("u")  
v = A.get_variable("v")
```

The resulting expressions can be used in algebraic operations.

### 3. Applying the derivation $\partial_x$

We can apply the differential operator $\partial_x$ (spatial derivative) to the variables. The method `.diff_x()` computes the first derivative, while passing an integer argument applies higher-order derivatives.

```python
ux = u.diff_x()        # First derivative: ∂u/∂x
vxx = v.diff_x(2)      # Second derivative: ∂²v/∂x²
```

These derivatives are automatically represented as `u_1` (for $\partial_x u$) and `v_2` (for $\partial_x^2 v$).

### 4. Displaying expressions

Let's print some expressions involving derivatives and combinations of terms:

```python
print(ux)                          # Output: u_1 (first derivative of u)
print(vxx)                         # Output: v_2 (second derivative of v)
print((ux - vxx) * (ux + vxx))     # Output: -1*(v_2)^2 + (u_1)^2
print((1 - t) * ux * (t + 1) * v)  # Output: (1 + -1*t^2)*u_1*v
```

The expressions are automatically simplified and displayed in a user-friendly format.


### 5. Extracting coefficients from an expression

The method `.coefficient()` is used for extracting coefficients of individual monomials with given exponents from an ordinary or a differential polynomial.

Coefficients of ordinary polynomials are rational numbers:

```python
P = (1 + t) * (1 - t) * (2 + 3 * t)    # P = 2 + 3*t + -2*t^2 + -3*t^3
print(P.coefficient(t * t))            # Output: -2
```

Coefficients of differential polynomials are polynomials in the base algebra:

```python
F = (1 + t + t * ux) * (2 - t + (3 + t) * vxx)
print(F.coefficient(ux * vxx))         # Output: 3*t + t^2
print(F.coefficient(u))                # Output: 0
```

The method .coefficient() accepts a Polynomial expression that contains exactly one monomial as its argument and returns the coefficient of a term with matching exponents. If the argument is a polynomial containing zero or more than one monomial, a ValueError will be raised.

```python
print(F.coefficient(0))            # raises a ValueError, because `0` has zero monomials
print(F.coefficient(u + v))        # raises a ValueError, because the argument has more than one monomial
```

### 6. Defining an evolutionary operator

An *evolutionary operator* is defined by its action on the generators of the algebra. We specify how the operator `D` acts on `u` and `v`:

```python
D = EvolutionaryOperator(A, {"u": ux, "v": u * vxx})
```

This defines $D(u) = \partial_x u$ and $D(v) = u \cdot \partial_x^2 v$. Now, we can apply this operator to more complex expressions:

```python
print(D.apply(vxx + ux * ux))  # Output: u*v_4 + 2*u_1*v_3 + 2*u_1*u_2 + u_2*v_2
```

### 7. Commutators of evolutionary operators

We now create two differential operators $D_1$ and $D_2$ that act on a variable `w` in another differential algebra. The commutator of $D_1$ and $D_2$ is computed to check whether they commute:

```python
B = DiffAlgebra(["w"], R)  
w = B.get_variable("w")  
D1 = EvolutionaryOperator(B, {"w": (w * w).diff_x()})  # D1(w) = ∂(w^2)
D2 = EvolutionaryOperator(B, {"w": (w * w * w).diff_x()})  # D2(w) = ∂(w^3)
```

To compute the commutator $C = [D_1, D_2]$, we use the `@` operator, which is overloaded to represent the commutator:

```python
C = D1 @ D2          # C = [D1, D2]
print(C.apply(w))    # Output: 0 (the operators commute)
```

Since the output is 0, this confirms that the two operators commute, and the corresponding differential equations are compatible.