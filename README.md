# `differential_algebra` package  
  
`differential_algebra` is a small symbolic computer algebra package, created to manipulate expressions in commutative and anticommutative free differential algebras.

It was created to explore problems in the theory of integrable hierarchies of partial differential equations; the intended scope is to cover both non-Hamiltonian and Hamiltonian evolutionary PDEs in one spatial variable. 
  
## Mathematical background  
  
A *differential* $R$-*algebra* $(A, \partial)$ is a (graded) commutative algebra $A$ over some commutative ring $R$ of constants, together with an $R$-linear derivation $\partial: A \to A$. For the purposes of our computations, we take some polynomial algebra $R = \mathbb{Q}[x_1, \dots, x_n]$ as the ring of constants.  
  
A homomorphism of differential $R$-algebras $f: (A, \partial_A) \to (B, \partial_B)$ is a homomorphism of underlying $R$-algebras such that $f \circ \partial_A = \partial_B \circ f$.  
  
A *free differential* $R$-*algebra* $\mathrm{fda}(R, U)$ over $R$ on a set of generators $U = (u^1, \dots, u^N)$ is defined by its universal property: for every differential $R$-algebra $(A, \partial)$ and for every map of sets $\phi: U \to A$, there exists a unique differential $R$-algebra homomorphism $\mathrm{fda}(U) \to A$ extending $\phi$.  

We understand free differential variables $u^i$ as symbols for functions $u^i(x)$, and identify $\partial$ with $\partial_x := \frac{\partial}{\partial x}$.  
  
Concretely, the algebra $\mathrm{fda}(R, U)$ is isomorphic to the polynomial algebra $R[u^i_j]$, $i=1 \dots N$, $j \geq 0$ where the symbol $u^i_j$ represents $\partial^j u^i = \frac{\partial^j u^i}{\partial x}$. Thus, the derivation $\partial$ acts on these symbols as $\partial u^i_{j} = u^i_{j+1}$.

An *evolutionary operator*, on a differential $R$-algebra $(A, \partial)$ is an $R$-linear derivation $D: A \to A$ which commutes with $\partial$. An evolutionary operator $D$ is uniquely defined by its values on the generators $D(u^i)$. It algebraically represents a system of differential equation $\left( \frac{\partial u^i}{\partial t} = D(u_i) \right)$.

The commutator $[D, D']$ of two evolutionary operators is again an evolutionary operator. Two differential operators $D$ and $D'$ commute if and only if the corresponding systems of differential equations $\left( \frac{\partial u^i}{\partial t} = D(u_i) \right)$ and $\left( \frac{\partial u^i}{\partial s} = D'(u_i) \right)$ are compatible.

## Usage example  
  
Importing the classes from the submodules:

```python  
from PolynomialAlgebra import Algebra
from DifferentialAlgebra import DiffAlgebra, EvolutionaryOperator  
```  
Creating objects representing a polynomial algebra $R$ in one variable and a differential algebra $A$ over $R$ in two variables:  
```python  
R = Algebra(["t"])  
A = DiffAlgebra(["u", "v"], R)  
```  
Obtaining variables as expressions, which can be used in operations, from the algebra object:  
```python  
t = R.get_variable("t")  
u = A.get_variable("u")  
v = A.get_variable("v")  
```  
Applying the derivation $\partial_x$, or its powers, to the variables:  
```python  
ux = u.diff_x()  
vxx = v.diff_x(2)  
```  
Display format of expressions:  
```python  
print(ux)                          # u_1  
print(vxx)                         # v_2  
print((ux - vxx) * (ux + vxx))     # -1*(v_2)^2 + (u_1)^2  
print((1 - t) * ux * (t + 1) * v)  # (1 + -1*t^2)*u_1*v  
```  
Creating an evolutionary operator by defining it on the generators, and applying it to an expression:  
```python  
D = EvolutionaryOperator(A, {"u" : ux, "v" : u * vxx})  
print(D.apply(vxx + ux * ux))      # u*v_4 + 2*u_1*v_3 + 2*u_1*u_2 + u_2*v_2  
```  
Checking that two evolutionary operators from the Riemann-Hopf hierarchy commute:  
```python  
B = DiffAlgebra(["w"], R)  
w = B.get_variable("w")  
D1 = EvolutionaryOperator(B, {"w" : (w * w).diff_x()})  
D2 = EvolutionaryOperator(B, {"w" : (w * w * w).diff_x()})  
C = D1 @ D2          # notation for the commutator C = [D1, D2]  
print(C.apply(w))    # 0  
```