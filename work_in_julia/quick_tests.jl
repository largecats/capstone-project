#############################################################################
# Date created: 2018/10/06
# Description: Quick tests along the way.
#############################################################################
using SymPy


t = symbols("t")
symL = SymLinearDifferentialOperator([t+1 t+1 t+1], (0,1), t)
L = LinearDifferentialOperator([t->t+1, t->t+1, t->t+1], (0,1), symL)
symDerivMatrix = get_symPDerivMatrix(L; substitute = true)
symDerivMatrix = get_symPDerivMatrix(L; substitute = false)
pStringMatrix = get_pStringMatrix(L)
u, v = SymFunction("u")(t), SymFunction("v")(t)
symUVForm = get_symUvForm(L, u, v; substitute = true)
symUVForm = get_symUvForm(L, u, v; substitute = false)
pDerivMatrix = [t->t+1 t->t; t->t+1 t->t]
B = get_B(L; pDerivMatrix = pDerivMatrix)
symB = get_B(L; symbolic = true)
symB = get_B(L; symbolic = true, substitute = false)
BHat = get_BHat(L, B)
symXi = get_symXi(L; substitute = true, xDef = t^2+2)
evaluate_xi(L, 1, symXi)
xi = [t->t^2+2; t->2t]
evaluate_xi(L, 1, xi)
U = VectorBoundaryForm([1 2; 3 4], [4 3; 2 1])
get_boundaryCondition(L, U, symXi)
get_boundaryCondition(L, U, xi)

t = symbols("t")
symL = SymLinearDifferentialOperator([t+1 t+1 t+1], (0,1), t)
L = LinearDifferentialOperator([t->t+1, t->t+1, t->t+1], (0,1), symL)
symDerivMatrix = get_symPDerivMatrix(L; substitute = true)
symDerivMatrix = get_symPDerivMatrix(L; substitute = false)
pStringMatrix = get_pStringMatrix(L)
u, v = SymFunction("u")(t), SymFunction("v")(t)
symUVForm = get_symUvForm(L, u, v; substitute = true)
symUVForm = get_symUvForm(L, u, v; substitute = false)
symB = get_symB(L; substitute = true)
symB = get_symB(L; substitute = false)
pDerivMatrix = [t->t+1 t->t; t->t+1 t->t]
B = get_B(L, pDerivMatrix)
BHat = get_BHat(L, B)
symXi = get_symXi(L; substitute = true, xDef = t^2+2)
evaluate_xi(L, 1, symXi)
xi = [t->t^2+2; t->2t]
evaluate_xi(L, 1, xi)
U = VectorBoundaryForm([1 2; 3 4], [4 3; 2 1])
get_boundaryCondition(L, U, symXi)
get_boundaryCondition(L, U, xi)

# Complex entries
t = symbols("t")
(a,b) = (0,1)
symPFunctions = [t+im t*im t]
pFunctions = [t->t+im t->t*im t->t]
symL = SymLinearDifferentialOperator(symPFunctions, (a,b), t)
L = LinearDifferentialOperator(pFunctions, (a,b), symL)
n = 2
MCand = rand(Uniform(1.0,10.0), n, n)
NCand = rand(Uniform(1.0,10.0), n, n)
U = VectorBoundaryForm(MCand, NCand)
pDerivMatrix = [t->t+im t->t; t->t*im t->im]
construct_validAdjoint(L, U, pDerivMatrix)
B = get_B(L; pDerivMatrix = pDerivMatrix)
BHat = get_BHat(L, B)
Uc = get_Uc(U)
H = get_H(U, Uc)
J = get_J(BHat, H)
adjointU = get_adjoint(J)
n = convert(Int, size(J)[1]/2)
Pstar = J[(n+1):2n,1:n]
Qstar = J[(n+1):2n, (n+1):2n]
adjoint = VectorBoundaryForm(Pstar, Qstar)
rank(hcat(Pstar, Qstar))