from diffalgebra import DifferentialRing, KdV_hierarchy, EvolutionOperator

R = DifferentialRing(functions=["u"])
u = R.gen("u")

F = 6*u*u.diff() + u.diff(order=3)
G = KdV_hierarchy(variable=u, order=2)

print(f"F = {F}")
print(f"G = {G}")

d_dt = EvolutionOperator(ring=R, mapping={u: F})
d_ds = EvolutionOperator(ring=R, mapping={u: G})

print(f"∂/∂t(u ** 2) = {d_dt(u ** 2)}")

d_dts = d_dt(d_ds(u))
d_dst = d_ds(d_dt(u))
print(f"[∂/∂t, ∂/∂s](u) = {d_dts - d_dst}")