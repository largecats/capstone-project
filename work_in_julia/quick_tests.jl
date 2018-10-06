#############################################################################
# Date created: 2018/10/06
# Description: Quick tests along the way.
#############################################################################
using SymPy

t = symbols("t")
symL = SymLinearDifferentialOperator([t+1 t+1 t+1], (0,1), t)
symDerivMatrix = get_symPDerivMatrix(symL, true)
symDerivMatrix = get_symPDerivMatrix(symL, false)
pStringMatrix = get_pStringMatrix(symL)
u, v = SymFunction("u")(t), SymFunction("v")(t)
symUVForm = get_symUvForm(symL, u, v, true)
symUVForm = get_symUvForm(symL, u, v, false)
get_symB(symL, true)
get_symB(symL, false)
get_symXi(symL)