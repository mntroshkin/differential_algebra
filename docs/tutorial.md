# diffalgebra - Tutorial

## Basics

The workflow of `diffalgebra` begins with defining a differential ring. We create two differential rings, in one and two functional generators respectively. 
We will use these two rings throughout this tutorial.

```python
>>> from diffalgebra import DifferentialRing
>>> R = DifferentialRing(functions=["u"])
>>> S = DifferentialRing(functions=["v", "w"])
```

Then we obtain functional generators via methods `.gen()` and `.gens()` of `DifferentialRing`.

```python
>>> u = R.gen("u")  # returns one specified generator
>>> v, w = S.gens() # returns the tuple of all ring generators
```

Ring generators can be used in arithmetic expressions:
```python
>>> (1 + u) * (2 + u)
2+3*u+u^2
```

Method `.diff()` returns derivatives of a differential polynomial:

```python
>>> (v * w).diff()
v*w_x+v_x*w
>>> (v * w).diff(order=2)   # higher derivative order as optional argument
v*w_xx+2*v_x*w_x+v_xx*w
```

For derivatives of ring generators, square brackets syntax is also available:

```python
>>> u * u[1] + u[3]
u*u_x+u_xxx
```

## Evolution operators

### Mathematical background

For a differential ring $R = \mathbb{Q}\{u\}$, an evolution operator $\phi = \frac{\partial}{\partial t}$ encodes a differential equation $\frac{\partial u}{\partial t} = F(u, u_x, u_{xx}, \dots)$, where $F$ is any differential polynomial.
Operator $\phi$ maps a differential polynomial $f$ to $\frac{\partial f}{\partial t}$ under this equation. 

Similarly, for ring with a tuple of generators $\mathbf{u} = (u_1, \dots, u_1)$, we can define a system of equation and the corresponding operator $\Phi := \frac{\partial}{\partial t}$ by $\frac{\partial u_i}{\partial t} = F_i(\mathbf{u}, \mathbf{u}_x, \mathbf{u}_{xx}, \dots)$. Here $(F_1, \dots, F_n)$ is an arbitrary tuple of differential polynomials.

### Implementation

To create an evolution operator, we pass a dictionary `{u_i: F_i}` mapping generators to their images:

```python
>>> from diffalgebra import EvolutionOperator
>>> phi = EvolutionOperator(ring=R, mapping={u: u * u[1]})
```

By definition, $\phi$ is a ring derivation which commutes with $\frac{\partial}{\partial x}$, which defines it on any differential polynomial:

```python
>>> phi(u ** 2)
2*u^2*u_x
>>> phi(u[1])
u*u_xx+(u_x)^2
```

Similarly for several functional variables:
```python
>>> Phi = EvolutionOperator(ring=S, mapping={v: v * v[1], w: v * w[1]})
>>> Phi(v * w)
v*v_x*w+v^2*w_x
```

## Example: the KdV equation and its first higher symmetry

The right-hand side of the KdV hierarchy equations is available via library functions:

```python
>>> from diffalgebra import KdV_hierachy
>>> F = KdV_hierarchy(variable=u)   # with no extra arguments, returns the original KdV equation
>>> F
6*u*u_x + u_xxx
>>> G = KdV_hierarchy(variable=u, order=2)  # higher-order symmetry of the KdV equation
>>> G
30*u^2*u_x+10*u*u_xxx+20*u_x*u_xx+u_5
```

We create two derivations $\frac{\partial}{\partial t}$ and $\frac{\partial}{\partial s}$ encoding the first two KdV equations $\frac{\partial u}{\partial t} = F(u, u_x, \dots)$ and $\frac{\partial u}{\partial s} = G(u, u_x, \dots)$:

```python
>>> from diffalgebra import EvolutionOperator
>>> d_dt = EvolutionOperator(ring=R, mapping={u: F})
>>> d_ds = EvolutionOperator(ring=R, mapping={u: G})
```

Finally, we check that $\left[\frac{\partial}{\partial t}, \frac{\partial}{\partial s}\right] (u) := \frac{\partial}{\partial t}\left( \frac{\partial u}{\partial s} \right) - \frac{\partial}{\partial s}\left( \frac{\partial u}{\partial t} \right) = 0$, that is, mixed partial derivatives commute, and the two equations are compatible:

```python
>>> d_dts = d_dt(d_ds(u))
>>> d_dst = d_ds(d_dt(u))
>>> d_dts - d_dst
0
```
