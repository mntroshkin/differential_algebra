from diffalgebra import DifferentialRing
R = DifferentialRing(functions=["u"])
S = DifferentialRing(functions=["v", "w"])

u = R.gen("u")  # returns one specified generator
v, w = S.gens() # returns the tuple of all ring generators

print((1 + u) * (2 + u))        # compute arithmetic expressions
print((v * w).diff())           # and derivatives
print((v * w).diff(order=2))    # higher order derivative

print(u * u[1] + u[3])          # shorthand syntax for derivatives of generators