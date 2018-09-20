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
subs(expr, diff(y,x),x)
subs(expr, Derivative(y,x), x)
coeff(expr,x)
print(args(expr))

t = Symbol("t")
u = SymFunction("u")(t)
Derivative(u, t)
diff(u)