
# `differential_algebra` package  
  
The `differential_algebra` package is created to manipulate expressions in commutative and anticommutative free differential algebras.  
  
## Mathematical background  
  
A *differential* $R$-*algebra* $(A, \partial)$ is a (graded) commutative algebra $A$ over some commutative ring $R$ of constants, together with a derivation $\partial$ such that $\partial(R) \equiv 0$. For the purposes of our computations, we take some polynomial algebra $R = \mathbb{Q}[x_1, \dots, x_n]$ as the ring of constants.  
  
A homomorphism of differential $R$-algebras $f: (A, \partial_A) \to (B, \partial_B)$ is a homomorphism of underlying $R$-algebras such that $f \circ \partial_A = \partial_B \circ f$.  
  
A *free differential* $R$-*algebra* $\mathrm{fda}(R, U)$ over $R$ on a set of generators $U = (u^1, \dots, u^N)$ is defined by its universal property: for every differential $R$-algebra $(A, \partial)$ and for every map of sets $\phi: U \to A$, there exists a unique differential $R$-algebra homomorphism $\mathrm{fda}(U) \to A$ extending $\phi$.  
  
Concretely, we have $\mathrm{fda}(u^1, \dots, u^N) = R[u^i_{(j)}]$ where $i=1\dots N$, $j\geq 0$. The differential $\partial$ acts on these generators by $\partial u^i_{(j)} = u^i_{(j+1)}$.  
  
We understand free differential variables $u^i$ as symbols for functions $u^i(x)$, and identify $\partial$ with $\partial_x := \frac{\partial}{\partial x}$.  

A *derivation*, or an *evolutionary operator*, on a differential $R$-algebra $(A, \partial)$ is an $R$-linear map $D: A \to A$ such that $D$ commutes with $\partial$. The commutator $[D, D']$ of two evolutionary operators is again an evolutionary operator.

## Usage example  
  
Importing the module:  
```python  
from PolynomialAlgebra import Algebra  
from DifferentialAlgebra import DiffAlgebra, Derivation  
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
D = Derivation(A, {"u" : ux, "v" : u * vxx})  
print(D.apply(vxx + ux * ux))      # u*v_4 + 2*u_1*v_3 + 2*u_1*u_2 + u_2*v_2  
```  
Checking that two evolutionary operators from the Riemann-Hopf hierarchy commute:  
```python  
B = DiffAlgebra(["w"], R)  
w = B.get_variable("w")  
D1 = Derivation(B, {"w" : (w * w).diff_x()})  
D2 = Derivation(B, {"w" : (w * w * w).diff_x()})  
C = D1 @ D2           # notation for the commutator tau = [phi, psi]  
print(tau.apply(w))  # 0  
```
<!--stackedit_data:
eyJoaXN0b3J5IjpbLTc3NzY0MDE5MV19
-->