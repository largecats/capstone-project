xiX = get_xi(L, x)
get_boundary_condition(L, U, xiX)
xiY = get_xi(L, y)
get_boundary_condition(L, U, xiY)

U = VectorBoundaryForm([1 2+im; 2 1+3im], [2+1im 3; 3 2])
U = VectorBoundaryForm([1 2; 2 4], [1 3; 2 6])
U = VectorBoundaryForm([1 2; 2 1], [2 3; 3 2])
Uc = find_Uc(U)
H = construct_H(U, Uc)

# function p(t)
#     return t + 1
# end
t = symbols("t")
p = t + 1
L = LinearDifferentialOperator([p,p], (0, 1))
function w(x)
    return x + 1
end
LinearDifferentialOperator([w,w], (0, 1))

u, v = SymFunction("u")(t), SymFunction("v")(t)
pStringMatrix = get_pString_matrix(L)
lambdify(deriv(p, t, 1))(1)
subs(deriv(p, t, 1), t, 3)
# Lambda(t, deriv(p, t, 1))
pFuncMatrix = get_pFunc_matrix(L, pStringMatrix, t)
subs(pFuncMatrix[1,2], t, 2)

uvForm = get_uv_form(L, t, u, v)
coeff1 = coeff(uvForm, deriv(u, t, 0)*deriv(conj(v), t, 0))
args(coeff1)
coeffMatrix = get_coefficient_matrix(L, uvForm, t, u, v)

bMatrix = get_B_matrix(L, uvForm, t, u, v, pStringMatrix, pFuncMatrix, coeffMatrix)
evaluate_matrix(bMatrix, t, 0)
bHatMatrix = get_B_hat(L, bMatrix, t)

# The kth derivative of a function u is encoded by the ordered pair (u, k)
struct Derivative
    u::Function
    k::Int

    Derivative(u::Function, k::Int) =
    try
        uDev = new(u, k)
        check_derivative_input(uDev)
        return uDev
    catch err
        return err
    end
end

# Checks whether the input degree k is valid
function check_derivative_input(uDev::Derivative)
    u, k = uDev.u, uDev.k
    if k < 0
        error("Degree of derivative must be >= 0")
    else
        return true
    end
end

# Need to find a way to encode Bjk in terms of the p functions so that the matrix B can be evaluated at a and b.
# I'm thinking of either deducing a formula for Bjk by checking n = 2, 3, 4, or implementing the Polish notation to figure out what Bjk are.
# Product rule of derivatives
function product_derivative(uDev::Derivative, vDev::Derivative, k)

end

using SymPy
x = symbols("x")
y = SymFunction("y")(x)
expr = 3 + x + x^2 + diff(y,x)*x*2, x
expr = 3 + x + x^2 + Derivative(y,x)*x*2, x
<<<<<<< HEAD
expr = 3 + x + x^2
=======
>>>>>>> master
subs(expr, diff(y,x),x)
subs(expr, Derivative(y,x), x)
coeff(expr,x)
print(args(expr))

t = Symbol("t")
u = SymFunction("u")(t)
Derivative(u, t)
diff(u)