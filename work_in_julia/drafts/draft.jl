# f(x) = e^(-im*x)
f(x)=x+im
f(x) = x
cF = Fun(f, 0..1) # Chebyshev approximation of f on [0,1]
# Fun(Chebyshev(【0.0,1.0】),[0.5, 0.5]) means that f can be approximated by a Chebyshev polynomial with coefficients 0.5, 0.5
# get coefficients of x as an array
coefficients(cF)
Fun(x->x, Chebyshev())
g = Fun(Chebyshev(),[1,2,3]) # represents 1*1 + 2*x + 3*cos2arccosx
l(x) = 1+2*x+3*cos(2*acos(x))
l(0)
g(0)
using Gadfly
plot(g, -1, 1)
h = im*x-x
roots(h)

integrate(t^2 * exp(t) * cos(t), t)
integrate(exp(-im*t)*(3t^2+cos(t)+1), t)
integrate(exp(-im*t)*(3t^2+cos(t)+1), (t,0,1))
SymPy.integrate(exp(-im*t)*(3t^2+cos(t)+1), (t,0,im))
6*e^(-im) - 2*im*e^(-im) + (e^(-im)*sin(1))/2 + (e^(-im)*cos(1))/2 + (im*e^(-im)*sin(1))/2 + 5*im

using NumericalMath
using HCubature
g(x) = exp(-im*x)*(3x^2+cos(x)+1)
g(x) = e^(-im*x)*(3x^2+cos(x)+1)
quadgk(g, 0, 1)
quadgk(g, 0, im)
# points = [0+0im, 1+0im]
# line_integral(fz, points)

t = symbols("t")
p = SymFunction("p")(t)
for counter = 2:6
    symL = SymLinearDifferentialOperator(repeat([p],outer=(1,counter)), (0, 1), t)
    u, v = SymFunction("u")(t), SymFunction("v")(t)
    symUvForm = symUv_form(symL, u, v)
    pStringMatrix = pString_matrix(symL)
    coeffMatrix = sym_coefficient_matrix(symL, symUvForm, u, v)
    Base.showarray(STDOUT, coeffMatrix, false)
end

t = symbols("t")
p = t + 1
symL = SymLinearDifferentialOperator([p p p], (0, 1), t)
pStringMatrix = pString_matrix(symL)
pSymDerivMatrix = pSymDeriv_matrix(symL)
u, v = SymFunction("u")(t), SymFunction("v")(t)
symUvForm = sym_uv_form(symL, u, v)
symCoeffMatrix = sym_coefficient_matrix(symL, symUvForm, u, v)

U = VectorBoundaryForm([1 0; 0 0], [0 0; 1 0])
rank_of_U(U)
Uc = get_Uc(U)
H = get_H(U, Uc)

function f0(x)
    return x+1
end

function f1(x)
    return 1
end
L = LinearDifferentialOperator([f0 f0 f0], (0,1))
pDerivMatrix = [f0 f1; f0 f1]
B = get_B(L, pStringMatrix, pDerivMatrix, symCoeffMatrix)
bHat = B_hat(L, B)
J = get_J(bHat, H)
adjointU = get_adjoint(J)

# x = t; xi = [x(t); x'(t); x(t); x'(t)]
function x0(t)
    return t + 2
end
function x1(t)
    return 1
end
xi = [x0; x1]
xiEval = evaluate_xi(L, xi)
get_boundary_condition(L, U, xi)

check_adjoint(L, U, adjointU, B)

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