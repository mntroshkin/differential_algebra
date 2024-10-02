# `differential-algebra` package

The `differential-algebra` package is created to manipulate expressions in commutative and anticommutative free  differential algebras.

## Mathematical background

A *differential algebra* $(A, \partial)$ is a (graded) commutative algebra $A$ together with its derivation $\partial$. For the purposes of our computations, we always take $\mathbb{Q}$ as the base field for $A$.

A homomorphism of differential algebras $f: (A, \partial_A) \to (B, \partial_B)$ is a homomorphism of underlying algebras such that $f \circ \partial_A = \partial_B \circ f$.

A *free differential algebra* $\mathrm{fda}(U)$ on a set of generators $U = \{u^1, \dots, u^N\}$ is defined by its universal property: for every differential algebra $(A, \partial)$ and every map of sets $\phi: U \to A$, there exists a unique differential algebra homomorphism $\mathrm{fda}(U) \to A$ extending $\phi$.

Concretely, we have $\mathrm{fda}(u^1, \dots, u^N) = \mathbb{Q}[u^i_{(j)}]$ where $i=1\dots N$, $j\geq 0$. The differential $\partial$ acts on these generators by $\partial u^i_{(j)} = u^i_{(j+1)}$.

We understand free differential variables $u^i$ as symbols for functions $u^i(x)$, and identify $\partial$ with $\partial_x := \frac{\partial}{\partial x}$.

## Usage example

Importing the module:
```python
import DifferentialAlgebra as da
```
Creating an object representing a commutative free differential algebra in two variables:
```python
A = da.DiffAlgebra(["u", "v"])
```
Obtaining variables as an expressions, which can be used in operations, from the algebra object:
```python
u = A.get_variable("u")
v = A.get_variable("v")
```
Applying the derivation $\partial_x$ to the variables:
```python
ux = u.diff_x()
vxx = v.diff_x(2)
```
Display format of expressions:
```python
print(ux)                       // u_1
print(vxx)                      // v_2
print((ux - vxx) * (ux + vxx))  // -1(v_2)^2 + (u_1)^2
```
Creating an evolutionary operator by defining it on the generators, and applying it to an expression:
```python
phi = da.Derivation(A, {"u" : ux, "v" : u * vxx})
print(phi.apply(vxx + ux * ux)) // u*v_4 + 2u_1*v_3 + 2u_1*u_2 + u_2*v_2
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbMTQwOTY5MzY5MCwxNzMxOTkwNzE1XX0=
-->