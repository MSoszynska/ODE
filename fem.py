from __future__ import print_function
from fenics import*
from parameters import*

# Define mesh
mesh = UnitIntervalMesh(N)

# Define function space
W = FiniteElement('Lagrange', mesh.ufl_cell(), 1)
V = FunctionSpace(mesh, MixedElement([W, W, W]))

# Define boundary
def boundary(x, on_boundary):
    return near(x[0], 0.0) and on_boundary
bc = DirichletBC(V, Constant((1.0, 1.0, 0.0)), boundary)

# Define trial and test functions 
U = TrialFunction(V)
(u, v, w) = split(U)
Phi = TestFunction(V)
(phi, psi, xi) = split(Phi)

# Define variational problem
F = grad(u)[0]*phi*dx + gamma*u*phi*dx - beta*v*v*phi*dx + \
    grad(v)[0]*psi*dx + mu*w*psi*dx - alpha*u*psi*dx + \
    grad(w)[0]*xi*dx - v*xi*dx

# Solve the problem
U_ = Function(V)
F = action(F, U_)
J = derivative(F, U_, U)
problem = NonlinearVariationalProblem(F, U_, bc, J)
solver = NonlinearVariationalSolver(problem)
solver.solve()
(u, v, w) = U_.split(U_)

# Save solution
u_pvd = File('solutions/u_f.pvd')
u.rename('u_f', 'Fluid')
u_pvd << u
v_pvd = File('solutions/v_s.pvd')
v.rename('v_s', 'Solid')
v_pvd << v
w_pvd = File('solutions/u_s.pvd')
w.rename('u_s', 'Solid')
w_pvd << w
