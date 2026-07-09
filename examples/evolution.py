from diffalgebra import DifferentialRing, EvolutionOperator
R = DifferentialRing(functions=["u"])
S = DifferentialRing(functions=["v", "w"])
u = R.gen("u")
v, w = S.gens()

phi = EvolutionOperator(ring=R, mapping={u: u * u[1]})
print(phi(u ** 2))
print(phi(u[1]))

Phi = EvolutionOperator(ring=S, mapping={v: v * v[1], w: v * w[1]})
print(Phi(v * w))