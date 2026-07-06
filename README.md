# diffalgebra

A symbolic differential algebra library, created for thesis calculations and as an exercise in computer algebra implementation.

**Disclaimer:** This is a learning project, not a production-grade library.
For serious research-level tools, I recommend checking out the links at <https://gdeq.org/Category:Software>.

## Background and motivation

During my thesis research, I was working with PDEs of the form:

$$u_t = F(u, u_x, u_{xx}, \dots),$$

for example, the Korteweg-de Vries (KdV) equation $u_t = 6 u u_x + u_{xxx}$.

To manipulate expressions such as the right-hand side of those equations (*differential polynomials*), I wrote some ad-hoc Python scripts. 
Later, I expanded and consolidated them into the first version of this library, as a personal exploration in the design of computer symbolic algebra tools.

## Features

- Perform arithmetic with unknown symbolic functions of one spatial variable and their derivatives;
- Detect total derivatives and integrate them;
- Compute substitutions via homomorphisms between differential rings;
- Check whether two differential equations of the form $u_t = F(u, u_x, \dots)$ and $u_s = G(u, u_x, \dots)$ are compatible; same for systems of PDEs in several functions $(u(x), v(x), \dots)$.

## Example: the KdV equation and its first higher symmetry

We initialize a differential ring $R$ in one symbolic function $u$ and obtain its generator:

```python
>>> from diffalgebra import DifferentialRing
>>> R = DifferentialRing(functions=["u"])
>>> u = R.gen("u")
```

We enter the right-hand side of the KdV equation $F = 6uu_x + u_{xxx}$ by hand:

```python
>>> F = 6*u*u.diff() + u.diff(order=3)
>>> F
6*u*u_x+u_xxx
```

We obtain the right-hand side of the second equation from KdV hierarchy from the library function:
```python
>>> from diffalgebra import KdV_hierachy
>>> G = KdV_hierarchy(variable=u, order=2)
>>> G
30*u^2*u_x+10*u*u_xxx+20*u_x*u_xx+u_5
```

We create two derivations $\frac{\partial}{\partial t}$ and $\frac{\partial}{\partial s}$ encoding the first two KdV equations $\frac{\partial u}{\partial t} = F(u, u_x, \dots)$ and $\frac{\partial u}{\partial s} = G(u, u_x, \dots)$:

```python
>>> from diffalgebra import EvolutionOperator
>>> d_dt = EvolutionOperator(ring=R, mapping={u: F})
>>> d_ds = EvolutionOperator(ring=R, mapping={u: G})
```

We can evaluate those derivations on any elements of $R$, for example, we compute $\frac{\partial}{\partial t} (u^2)$ - check that it satisfies the identity $\frac{\partial}{\partial t} (u^2) = 2uu_t$!

```python
>>> d_dt(u ** 2)
12*u^2*u_x+2*u*u_xxx
```

Finally, we check that $\left[\frac{\partial}{\partial t}, \frac{\partial}{\partial s}\right] (u) := \frac{\partial}{\partial t}\left( \frac{\partial u}{\partial s} \right) - \frac{\partial}{\partial s}\left( \frac{\partial u}{\partial t} \right) = 0$, that is, mixed partial derivatives commute, and the two equations are compatible:

```python
>>> d_dts = d_dt(d_ds(u))
>>> d_dst = d_ds(d_dt(u))
>>> d_dts - d_dst
0
```

## Limitations

- Only supports symbolic functions of one independent variable and polynomials in their derivatives;
- Rational and transcendental functions are not implemented;
- No integration with other computer algebra systems;
- No performance optimizations for large expressions.

## Installation and running tests

Clone the repository with
```bash
git clone https://github.com/mntroshkin/diffalgebra.git
```

Navigate to the root directory of the project and install with
```bash
pip install -e .
```

Run all tests with
```bash
pytest
```

Note: `hypothesis`-powered property tests might take >10s to run in total, since they check many random examples in each test.


## Potential refinements

**Memoization of derivatives.** Currently, differential ring homomorphisms recompute derivatives naively on each occurrence. Memoizing higher derivatives on first retrieval could be the natural next step; however, given the actual scale of my thesis calculations, I chose not to implement it.

## References

* Dickey L. A. *Soliton Equations and Hamiltonian Systems.*