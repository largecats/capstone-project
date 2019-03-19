
using SymPy
using PyCall
sympy = pyimport("sympy")
# using Roots
using Distributions
# using IntervalArithmetic
# using IntervalRootFinding
using ApproxFun

using Plots

using NBInclude
using QuadGK
import QuadGK.quadgk
# using HCubature
using ApproxFun
using Roots
using Gadfly
using PyPlot
# pygui(true)

TOL = 1e-05

DIGITS = 3

INFTY = 10

# t = symbols("t")

x = symbols("x")
sympyAddExpr = 1 + x

sympyMultExpr = 2*x

sympyPowerExpr = x^2

sympyExpExpr = e^x

function check_all(array, condition)
    for x in array
        if !condition(x)
            return false
        end
    end
    return true
end

array = [0,1,2,3]
condition = x -> x>2
check_all(array, condition)

function check_any(array, condition)
    for x in array
        if condition(x)
            return true
        end
    end
    return false
end

array = [0,1,2,3]
condition = x -> x>2
check_any(array, condition)

function set_tol(x::Union{Number, Array}, y::Union{Number, Array}; atol = TOL)
    if isa(x, Number) && isa(y, Number)
       return atol * mean([x y])
    elseif isa(x, Array) && isa(y, Array)
        if size(x) != size(y)
            throw(error("Array dimensions do not match"))
        end
        # Avoid InexactError() when taking norm()
        x = convert(Array{Complex}, x)
        y = convert(Array{Complex}, y)
        return atol * (norm(x,2) + norm(y,2))
    else
        throw(error("Invalid input"))
    end
end

x = 14
y = 21
println("set_tol($x, $y) = $(set_tol(x,y))")

A = [4 0.6; 3 2]
B = [5 1; 10 3]
set_tol(A, B)
println("set_tol($A, $B) = $(set_tol(A,B))")

function evaluate(func::Union{Function, SymPy.Sym, Number}, a::Number)
    if isa(func, Function)
        funcA = func(a)
    elseif isa(func, SymPy.Sym) # SymPy.Sym must come before Number because SymPy.Sym will be recognized as Number
        freeSymbols = free_symbols(func)
        if length(freeSymbols) > 1
            throw(error("func should be univariate"))
        elseif length(freeSymbols) == 1
            t = free_symbols(func)[1,1]
            if isa(a, SymPy.Sym) # if x is SymPy.Sym, do not convert result to Number to preserve pretty printing
                funcA = subs(func, t, a)
            else
                funcA = SymPy.N(subs(func, t, a))
            end
        else
            funcA = func
        end
    else # func is Number
        funcA = func
    end
    return funcA
end

func = x -> x+1
xVal = 2
println("func($xVal) = $(evaluate(func, xVal))")

x = symbols("x")
func = x+1
println("func($xVal) = $(evaluate(func, xVal))")

a = symbols("a")
println("func($a) = $(evaluate(func, a))")

x = symbols("x")
array = [2x 1; x^3 x]
a = symbols("a")
evaluate.(array, a)

function partition(n::Int)
    if n < 0
        throw(error("Non-negative n required"))
    end
    output = []
    for i = 0:n
        j = n - i
        push!(output, (i,j))
    end
    return output
end

n = 5
partition(5)

function get_deriv(u::Union{SymPy.Sym, Number}, k::Int)
    if k < 0
        throw(error("Non-negative k required"))
    end
    if isa(u, SymPy.Sym)
        freeSymbols = free_symbols(u)
        if length(freeSymbols) > 1
            throw(error("u should be univariate"))
        elseif length(freeSymbols) == 1
            t = freeSymbols[1]
            y = u
            for i = 1:k
                newY = diff(y, t)
                y = newY
            end
            return y
        else
            if k > 0
                return 0
            else
                return u
            end
        end
    else
        if k > 0
            return 0
        else
            return u
        end
    end
end

u = 3
k = 1
println("get_deriv($u, $k) = $(get_deriv(u, k))")

t = symbols("t")
u = t^2+2t+3
get_deriv(u, k)

function add_func(f::Union{Number, Function}, g::Union{Number, Function})
    function h(x)
        if isa(f, Number)
            if isa(g, Number)
                return f + g
            else
                return f + g(x)
            end
        elseif isa(f, Function)
            if isa(g, Number)
                return f(x) + g
            else
                return f(x) + g(x)
            end
        end
    end
    return h
end

f(x) = x^3+1
g(x) = 4x
x = 2
add_func(f, g)(x) == f(x) + g(x)

function mult_func(f::Union{Number, Function}, g::Union{Number, Function})
    function h(x)
        if isa(f, Number)
            if isa(g, Number)
                return f * g
            else
                return f * g(x)
            end
        elseif isa(f, Function)
            if isa(g, Number)
                return f(x) * g
            else
                return f(x) * g(x)
            end
        end
    end
    return h
end

f(x) = x^3+1
g(x) = 4x
x = 2
mult_func(f, g)(x) == f(x) * g(x)

function sym_to_func(sym::Union{SymPy.Sym, Number})
    try
        freeSymbols = free_symbols(sym)
        if length(freeSymbols) > 1
            throw(error("sym should be univariate"))
        else
            function func(x)
                if length(freeSymbols) == 0
                    result = SymPy.N(sym)
                else
                    result = SymPy.N(subs(sym, freeSymbols[1], x))
                end
                return result
            end
            return func
        end
    catch
        function func(x)
            return sym
        end
        return func
    end
end

t = symbols("t")
sym = t^2+t
tVal = 2
println("sym_to_func($sym)($tVal) = $(sym_to_func(sym)(tVal))")

sym = [1 t; t^2 t^3]
sym_to_func.(sym)[2,1](2)
println("sym_to_func($sym[2,1])($tVal) = $(sym_to_func.(sym[2,1])(tVal))")

sym = 2
println("sym_to_func($sym)($tVal) = $(sym_to_func(sym)(tVal))")

function prettyRound(x::Number; digits::Int = DIGITS)
    if isa(x, Int)
        return x
    elseif isa(x, Real)
        if isa(x, Rational) || isa(x, Irrational) # If x is rational or irrational numbers like e, pi
            return x
        elseif round(abs(x), digits) == floor(abs(x))
            return Int(round(x))
        else
            return round(x, digits)
            # return rationalize(x)
        end
    elseif isa(x, Complex)
        roundedReal = prettyRound(real(x), digits = digits)
        roundedComplex = prettyRound(imag(x), digits = digits)
        return roundedReal + im*roundedComplex
    else
        return round(x, digits)
    end
end

# x = 0.0000001*im
x = 4.854572304702339 - 49.023298458074294im
prettyRound(x)
println("prettyRound($x) = $(prettyRound(x))")

x = [1/3 1/4 0]
println("prettyRound($x) = $(prettyRound.(x))")

x = [1//3 1//4 0//1]
println("prettyRound($x) = $(prettyRound.(x))")

function prettyPrint(x::Union{Number, SymPy.Sym})
    expr = x
    if isa(expr, SymPy.Sym)
        prettyExpr = expr
        for a in sympy[:preorder_traversal](expr)
            if length(free_symbols(a)) == 0 && length(args(a)) == 0
                if !(a in [e, PI]) && length(intersect(args(a), [e, PI])) == 0 # keep the transcendental numbers as symbols
                    prettyA = prettyRound.(SymPy.N(a))
                    prettyExpr = subs(prettyExpr, (a, prettyA))
                end
            end
        end
    else
        prettyExpr = prettyRound.(expr)
        prettyExpr = convert(SymPy.Sym, prettyExpr)
    end
    return prettyExpr
end

t = symbols("t")
x = 1//3*t + e^(2t) + PI/2 + 1.00001*im + sympy[:sqrt](2) + 1//6
prettyPrint(x)

t = symbols("t")
x = [1//3*t PI; 1//6*t^2+t 1]
prettyPrint.(x)

x = [0.0+1.0*im 1.0+0.0*im -1.0+0.0*im 0.0-1.0*im]
prettyPrint.(x)

function is_approxLess(x::Number, y::Number; atol = TOL)
    return !isapprox(x,y; atol = atol) && x < y
end

x = 1
y = x + TOL/2
println("is_approxLess($x,$y) = $(is_approxLess(x,y))")

y = x + TOL*2
println("is_approxLess($x,$y) = $(is_approxLess(x,y))")

function is_approx(x::Number, y::Number; atol = TOL)
    return isapprox(x, y; atol = atol)
end

x = 1
y = x + TOL/2
println("is_approx($x,$y) = $(is_approx(x,y))")

y = x + TOL*2
println("is_approx($x,$y) = $(is_approx(x,y))")

function argument(z::Number)
    if angle(z) >= 0 # in [0,pi]
        return angle(z)
    else 
        # angle(z) is in (-pi, 0]
        # Shift it up to (pi,2pi]
        argument = 2pi + angle(z) # This is in (pi,2pi]
        if is_approx(argument, 2pi) # Map 2pi to 0
            return 0
        else
            return argument # This is now in [0,2pi)
        end
    end
end

z = 1
println("argument($z) = $(argument(z))")

z = -1-im
println("argument($z) = $(argument(z))")

function trace_contour(a::Number, n::Int, sampleSize::Int; infty = INFTY)
    lambdaVec = []
    for counter = 1:sampleSize
        x = rand(Uniform(-infty,infty), 1, 1)[1]
        y = rand(Uniform(-infty,infty), 1, 1)[1]
        lambda = x + y*im
        if real(a*lambda^n)>0
            append!(lambdaVec, lambda)
        end
    end
    Gadfly.plot(x=real(lambdaVec), y=imag(lambdaVec), Guide.xlabel("Re"), Guide.ylabel("Im"), Coord.Cartesian(ymin=-infty,ymax=infty, xmin=-infty, xmax=infty, fixed = true))
end

n = 2
a = 1
sampleSize = 1000
trace_contour(a, n, sampleSize)

n = 3
a = im
sampleSize = 1000
trace_contour(a, n, sampleSize)

n = 3
a = -im
sampleSize = 1000
trace_contour(a, n, sampleSize)

function get_distancePointLine(z::Number, theta::Number)
    if theta >= 2pi && theta < 0
        throw(error("Theta must be in [0,2pi)"))
    else
        if is_approx(argument(z), theta)
            return 0
        else
            x0, y0 = real(z), imag(z)
            if is_approx(theta, pi/2) || is_approx(theta, 3pi/2)
                return abs(x0)
            elseif is_approx(theta, 0) || is_approx(theta, 2pi)
                return abs(y0)
            else
                k = tan(theta)
                x = (y0+1/k*x0)/(k+1/k)
                y = k*x
                distance = norm(z-(x+im*y))
                return distance
            end
        end
    end
end

z = 1
theta = pi/3
println("get_distancePointLine($z, $theta) = $(get_distancePointLine(z, theta))") # sqrt(3)/2

struct StructDefinitionError <: Exception
    msg::String
end

struct SymLinearDifferentialOperator
    # Entries in the array should be SymPy.Sym or Number. SymPy.Sym seems to be a subtype of Number, i.e., Array{Union{Number,SymPy.Sym}} returns Array{Number}. But specifying symPFunctions as Array{Number,2} gives a MethodError when the entries are Sympy.Sym objects.
    symPFunctions::Array
    interval::Tuple{Number,Number}
    t::SymPy.Sym
    SymLinearDifferentialOperator(symPFunctions::Array, interval::Tuple{Number,Number}, t::SymPy.Sym) =
    try
        symL = new(symPFunctions, interval, t)
        check_symLinearDifferentialOperator_input(symL)
        return symL
    catch err
        throw(err)
    end
end

function check_symLinearDifferentialOperator_input(symL::SymLinearDifferentialOperator)
    symPFunctions, (a,b), t = symL.symPFunctions, symL.interval, symL.t
    for symPFunc in symPFunctions
        if isa(symPFunc, SymPy.Sym)
            if size(free_symbols(symPFunc)) != (1,) && size(free_symbols(symPFunc)) != (0,)
                throw(StructDefinitionError(:"Only one free symbol is allowed in symP_k"))
            end
        elseif !isa(symPFunc, Number)
            throw(StructDefinitionError(:"symP_k should be SymPy.Sym or Number"))
        end
    end
    return true
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)

# symL is an attribute of L that needs to be input by the user. There are checks to make sure symL is indeed the symbolic version of L.
# Principle: Functionalities of Julia Functions >= Functionalities of SymPy. If p_k has no SymPy representation, the only consequence should be that outputs by functions that take L as arugment has no symbolic expression. E.g., we allow L.pFunctions and L.symL.pFunctions to differ.
struct LinearDifferentialOperator
    pFunctions::Array # Array of julia functions or numbers representing constant functions
    interval::Tuple{Number,Number}
    symL::SymLinearDifferentialOperator
    LinearDifferentialOperator(pFunctions::Array, interval::Tuple{Number,Number}, symL::SymLinearDifferentialOperator) =
    try
        L = new(pFunctions, interval, symL)
        check_linearDifferentialOperator_input(L)
        return L
    catch err
        throw(err)
    end
end

# Assume symFunc has only one free symbol, as required by the definition of SymLinearDifferentialOperator. 
# That is, assume the input symFunc comes from SymLinearDifferentialOperator.
function check_func_sym_equal(func::Union{Function,Number}, symFunc, interval::Tuple{Number,Number}, t::SymPy.Sym) # symFunc should be Union{SymPy.Sym, Number}, but somehow SymPy.Sym gets ignored
    (a,b) = interval
    # Randomly sample 1000 points from (a,b) and check if func and symFunc agree on them
    for i = 1:1000
        # Check endpoints
        if i == 1
            x = a
        elseif i == 2
            x = b
        else
            x = rand(Uniform(a,b), 1)[1,1]
        end
        funcEvalX = evaluate.(func, x)
        if isa(symFunc, SymPy.Sym)
            symFuncEvalX = SymPy.N(subs(symFunc,t,x))
            # N() converts SymPy.Sym to Number
            # https://docs.sympy.org/latest/modules/evalf.html
            # subs() works no matter symFunc is Number or SymPy.Sym
        else
            symFuncEvalX = symFunc
        end
        tol = set_tol(funcEvalX, symFuncEvalX)
        if !isapprox(real(funcEvalX), real(symFuncEvalX); atol = real(tol)) ||
            !isapprox(imag(funcEvalX), imag(symFuncEvalX); atol = imag(tol))
            println("x = $x")
            println("symFunc = $symFunc")
            println("funcEvalX = $funcEvalX")
            println("symFuncEvalX = $symFuncEvalX")
            return false
        end
    end
    return true
end

# Check whether the inputs of L are valid.
function check_linearDifferentialOperator_input(L::LinearDifferentialOperator)
    pFunctions, (a,b), symL = L.pFunctions, L.interval, L.symL
    symPFunctions, t = symL.symPFunctions, symL.t
    # domainC = Complex(a..b, 0..0) # Domain [a,b] represented in the complex plane
    p0 = pFunctions[1]
    # p0Chebyshev = Fun(p0, a..b) # Chebysev polynomial approximation of p0 on [a,b]
    if !check_all(pFunctions, pFunc -> (isa(pFunc, Function) || isa(pFunc, Number)))
        throw(StructDefinitionError(:"p_k should be Function or Number"))
    elseif length(pFunctions) != length(symPFunctions)
        throw(StructDefinitionError(:"Number of p_k and symP_k do not match"))
    elseif (a,b) != symL.interval
        throw(StructDefinitionError(:"Intervals of L and symL do not match"))
    # # Assume p_k are in C^{n-k}. Check whether p0 vanishes on [a,b]. 
    # # roots() in IntervalRootFinding doesn't work if p0 is sth like t*im - 2*im. Neither does find_zero() in Roots.
    # # ApproxFun.roots() 
    # elseif (isa(p0, Function) && (!isempty(roots(p0Chebyshev)) || all(x->x>b, roots(p0Chebyshev)) || all(x->x<b, roots(p0Chebyshev)) || p0(a) == 0 || p0(b) == 0)) || p0 == 0 
    #     throw(StructDefinitionError(:"p0 vanishes on [a,b]"))
    elseif !all(i -> check_func_sym_equal(pFunctions[i], symPFunctions[i], (a,b), t), 1:length(pFunctions))
        # throw(StructDefinitionError(:"symP_k does not agree with p_k on [a,b]"))
        warn("symP_k does not agree with p_k on [a,b]") # Make this a warning instead of an error because the functionalities of Julia Functions may be more than those of SymPy objects; we do not want to compromise the functionalities of LinearDifferentialOperator because of the restrictions on SymPy.
    else
        return true
    end
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
# Direct construction
pFunctions = [t->1 t->t+1 t->t^2+t+1]
L = LinearDifferentialOperator(pFunctions, interval, symL)

struct VectorBoundaryForm
    M::Array # Why can't I specify Array{Number,2} without having a MethodError?
    N::Array
    VectorBoundaryForm(M::Array, N::Array) =
    try
        U = new(M, N)
        check_vectorBoundaryForm_input(U)
        return U
    catch err
        throw(err)
    end
end

# Check whether the input matrices that characterize U are valid
function check_vectorBoundaryForm_input(U::VectorBoundaryForm)
    # M, N = U.M, U.N
    # Avoid Inexact() error when taking rank()
    M = convert(Array{Complex}, U.M)
    N = convert(Array{Complex}, U.N)
    if !(check_all(U.M, x -> isa(x, Number)) && check_all(U.N, x -> isa(x, Number)))
        throw(StructDefinitionError(:"Entries of M, N should be Number"))
    elseif size(U.M) != size(U.N)
        throw(StructDefinitionError(:"M, N dimensions do not match"))
    elseif size(U.M)[1] != size(U.M)[2]
        throw(StructDefinitionError(:"M, N should be square matrices"))
    elseif rank(hcat(M, N)) != size(M)[1] # rank() throws weird "InexactError()" when taking some complex matrices
        throw(StructDefinitionError(:"Boundary operators not linearly independent"))
    else
        return true
    end
end

M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)

function get_L(symL::SymLinearDifferentialOperator)
    symPFunctions, (a,b), t = symL.symPFunctions, symL.interval, symL.t
    if check_all(symPFunctions, x->!isa(x, SymPy.Sym))
        pFunctions = symPFunctions
    else
        pFunctions = sym_to_func.(symPFunctions)
    end
    L = LinearDifferentialOperator(pFunctions, (a,b), symL)
    return L
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)

function get_URank(U::VectorBoundaryForm)
    # Avoid InexactError() when taking hcat() and rank()
    M = convert(Array{Complex}, U.M)
    N = convert(Array{Complex}, U.N)
    MHcatN = hcat(M, N)
    return rank(MHcatN)
end

M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)
get_URank(U)

function get_Uc(U::VectorBoundaryForm)
    try
        check_vectorBoundaryForm_input(U)
        n = get_URank(U)
        I = complex(eye(2*n))
        M, N = U.M, U.N
        MHcatN = hcat(M, N)
        # Avoid InexactError() when taking rank()
        mat = convert(Array{Complex}, MHcatN)
        for i = 1:(2*n)
            newMat = vcat(mat, I[i:i,:])
            newMat = convert(Array{Complex}, newMat)
            if rank(newMat) == rank(mat) + 1
                mat = newMat
            end
        end
        UcHcat = mat[(n+1):(2n),:]
        Uc = VectorBoundaryForm(UcHcat[:,1:n], UcHcat[:,(n+1):(2n)])
        return Uc
    catch err
        return err
    end
end

M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)
Uc = get_Uc(U)

function get_H(U::VectorBoundaryForm, Uc::VectorBoundaryForm)
    MHcatN = hcat(convert(Array{Complex}, U.M), convert(Array{Complex}, U.N))
    McHcatNc = hcat(convert(Array{Complex}, Uc.M), convert(Array{Complex}, Uc.N))
    H = vcat(MHcatN, McHcatNc)
    return H
end

M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)
Uc = get_Uc(U)
get_H(U, Uc)

function get_pDerivMatrix(L::LinearDifferentialOperator; symbolic = false)
    if symbolic
        symL = L.symL
        symPFunctions, t = symL.symPFunctions, symL.t
        n = length(symPFunctions)-1
        symPDerivMatrix = Array{SymPy.Sym}(n,n)
        pFunctionSymbols = symPFunctions
        for i in 0:(n-1)
            for j in 0:(n-1)
                index, degree = i, j
                symP = pFunctionSymbols[index+1]
                # If symP is not a Sympy.Sym object (e.g., is a Number instead), then cannot use get_deriv()
                if !isa(symP, SymPy.Sym)
                    if degree > 0
                        symPDeriv = 0
                    else
                        symPDeriv = symP
                    end
                else
                    symPDeriv = get_deriv(symP, degree)
                end
                symPDerivMatrix[i+1,j+1] = symPDeriv
            end
        end
        return symPDerivMatrix
    else
        symPDerivMatrix = get_pDerivMatrix(L; symbolic = true)
        n = length(L.pFunctions)-1
        pDerivMatrix = sym_to_func.(symPDerivMatrix)
    end
    return pDerivMatrix
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
# pFunctions = [t->1 t->t+1 t->t^2+t+1]
# L = LinearDifferentialOperator(pFunctions, interval, symL)
L = get_L(symL)

tVal = 2
pDerivMatrix = get_pDerivMatrix(L; symbolic = false)
println("pDerivMatrix($tVal) = $(evaluate.(pDerivMatrix,tVal))")

symPDerivMatrix = get_pDerivMatrix(L; symbolic = true)

function get_Bjk(L::LinearDifferentialOperator, j::Int, k::Int; symbolic = false, pDerivMatrix = get_pDerivMatrix(L; symbolic = symbolic))
    n = length(L.pFunctions)-1
    if j <= 0 || j > n || k <= 0 || k > n
        throw("j, k should be in {1, ..., n}")
    end
    sum = 0
    if symbolic
        symPDerivMatrix = get_pDerivMatrix(L; symbolic = true)
        for l = (j-1):(n-k)
            summand = binomial(l, j-1) * symPDerivMatrix[n-k-l+1, l-j+1+1] * (-1)^l
            sum += summand
        end
    else
        for l = (j-1):(n-k)
            summand = mult_func(binomial(l, j-1) * (-1)^l, pDerivMatrix[n-k-l+1, l-j+1+1])
            sum = add_func(sum, summand)
        end
    end
    return sum
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)

# pFunctions = [t->1 t->t+1 t->t^2+t+1]
# L = LinearDifferentialOperator(pFunctions, interval, symL)
L = get_L(symL)

# pDerivMatrix = [1 0; t->1+t t->1]
j, k = 1, 1
BjkSym = get_Bjk(L, j, k; symbolic = true)

tVal = 1
Bjk = get_Bjk(L, j, k; symbolic = false)
println("Bjk($tVal) = $(Bjk(tVal))")
println("BjkSym($tVal) = $(evaluate.(BjkSym, tVal))")

function get_B(L::LinearDifferentialOperator; symbolic = false, pDerivMatrix = get_pDerivMatrix(L; symbolic = symbolic))
    n = length(L.pFunctions)-1
    B = Array{Union{Function, Number, SymPy.Sym}}(n,n)
    for j = 1:n
        for k = 1:n
            B[j,k] = get_Bjk(L, j, k; symbolic = symbolic, pDerivMatrix = pDerivMatrix)
        end
    end
    return B
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)

BSym = get_B(L; symbolic = true)
# Only Array{SymPy.Sym} would be automatically pretty-printed
# Enfore pretty-printing
prettyPrint.(BSym)

B = get_B(L; symbolic = false)
tVal = 1
println("B($tVal) = $(evaluate.(B, tVal))")
println("BSym($tVal) = $(evaluate.(BSym, tVal))")

function get_BHat(L::LinearDifferentialOperator, B::Array)
#     if check_any(B, x->isa(x, SymPy.Sym))
#         throw("Entries of B should be Function or Number")
#     end
    pFunctions, (a,b) = L.pFunctions, L.interval
    n = length(pFunctions)-1
    BHat = Array{Complex}(2n,2n)
    BEvalA = evaluate.(B, a)
    BEvalB = evaluate.(B, b)
    BHat[1:n,1:n] = -BEvalA
    BHat[(n+1):(2n),(n+1):(2n)] = BEvalB
    BHat[1:n, (n+1):(2n)] = 0
    BHat[(n+1):(2n), 1:n] = 0
    return BHat
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)

B = get_B(L; symbolic = false)
get_BHat(L, B)

function get_J(BHat, H)
    n = size(H)[1]
    H = convert(Array{Complex}, H)
    J = (BHat * inv(H))'
    # J = convert(Array{Complex}, J)
    return J
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)

B = get_B(L; symbolic = false)
BHat = get_BHat(L, B)

M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)
Uc = get_Uc(U)
H = get_H(U, Uc)

get_J(BHat, H)

function get_adjointCand(J)
    n = convert(Int, size(J)[1]/2)
    J = convert(Array{Complex}, J)
    PStar = J[(n+1):2n,1:n]
    QStar = J[(n+1):2n, (n+1):2n]
    adjointU = VectorBoundaryForm(PStar, QStar)
    return adjointU
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)

B = get_B(L; symbolic = false)
BHat = get_BHat(L, B)

M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)
Uc = get_Uc(U)
H = get_H(U, Uc)

J = get_J(BHat, H)
adjoint = get_adjointCand(J)

function get_xi(L::LinearDifferentialOperator; symbolic = true, xSym= nothing)
    if symbolic
        t = L.symL.t
        n = length(L.pFunctions)-1
        symXi = Array{SymPy.Sym}(n,1)
        if isa(xSym, Void)
            throw(error("xSymrequired"))
        else
            for i = 1:n
                symXi[i] = get_deriv(xSym, i-1)
            end
            return symXi
        end
    else
        if isa(xSym, Void)
            throw(error("xSym required"))
        elseif !isa(xSym, SymPy.Sym) && !isa(xSym, Number)
            throw(error("xSym should be SymPy.Sym or Number"))
        else
            symXi = get_xi(L; symbolic = true, xSym = xSym)
            xi = sym_to_func.(symXi)
        end
    end
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)

xSym = t^2+2t
symXi = get_xi(L; symbolic = true, xSym = xSym)

xi = get_xi(L; symbolic=false, xSym = xSym)
tVal = 1
println("xi($tVal) = $(evaluate.(xi, tVal))")
println("symXi($tVal) = $(evaluate.(symXi, tVal))")

function get_Ux(L::LinearDifferentialOperator, U::VectorBoundaryForm, xSym)
    (a, b) = L.interval
    xi = get_xi(L; symbolic = false, xSym = xSym)
    xiEvalA = evaluate.(xi, a)
    xiEvalB = evaluate.(xi, b)
    M, N = U.M, U.N
    Ux = M*xiEvalA + N*xiEvalB
    return Ux
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)
M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)

xSym = t^2+2t
Ux = get_Ux(L, U, xSym)
println("U($symX) = $Ux")

function check_adjoint(L::LinearDifferentialOperator, U::VectorBoundaryForm, adjointU::VectorBoundaryForm, B::Array)
    (a, b) = L.interval
    M, N = U.M, U.N
    P, Q = (adjointU.M)', (adjointU.N)'
    # Avoid InexactError() when taking inv()
    BEvalA = convert(Array{Complex}, evaluate.(B, a))
    BEvalB = convert(Array{Complex}, evaluate.(B, b))
    left = M * inv(BEvalA) * P
    right = N * inv(BEvalB) * Q
#     println("left = $left")
#     println("right = $right")
    tol = set_tol(left, right)
    return all(i -> isapprox(left[i], right[i]; atol = tol), 1:length(left)) # Can't use == to deterimine equality because left and right are arrays of floats
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)
M = [1 0; 2 0]
N = [0 2; 0 1]
# M = [3.9 5.4; 1+2*im 2]
# N = [4.7 8.1 + im; 0.5*im 10]
U = VectorBoundaryForm(M, N)
Uc = get_Uc(U)

# Non-symbolic
B = get_B(L)
BHat = get_BHat(L, B)
H = get_H(U, Uc)
J = get_J(BHat, H)
adjointU = get_adjointCand(J)

check_adjoint(L, U, adjointU, B)

function get_adjointU(L::LinearDifferentialOperator, U::VectorBoundaryForm, pDerivMatrix=get_pDerivMatrix(L))
    B = get_B(L; pDerivMatrix = pDerivMatrix)
    BHat = get_BHat(L, B)
    Uc = get_Uc(U)
    H = get_H(U, Uc)
    J = get_J(BHat, H)
    adjointU = get_adjointCand(J)
    if check_adjoint(L, U, adjointU, B)
        return adjointU
    else
        throw(error("Adjoint found not valid"))
    end
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
# pFunctions = [t->1 t->t+1 t->t^2+t+1]
# L = LinearDifferentialOperator(pFunctions, interval, symL)
L = get_L(symL)
M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)

adjointU = get_adjointU(L, U)

# although the function body is the same as "power" and "others", this case is isolated because negative exponents, e.g., factor_list(e^(-im*x)), give PolynomialError('a polynomial expected, got exp(-I*x)',), while factor_list(cos(x)) runs normally
function separate_real_imaginary_exp(expr::SymPy.Sym)
    result = real(expr) + im*imag(expr)
    return result
end

x = symbols("x", real = true)
y = symbols("y", real = true)
expr = e^(x+im*y)
println("func($expr) = $(SymPy.func(expr))")
separate_real_imaginary_exp(expr)

# we won't be dealing with cases like x^(x^x)
function separate_real_imaginary_power(expr::SymPy.Sym)
    result = real(expr) + im*imag(expr)
    return result
end

x = symbols("x", real = true)
y = symbols("y", real = true)
expr = (x+im*y)^2
println("func($expr) = $(SymPy.func(expr))")
separate_real_imaginary_power(expr)

function separate_real_imaginary_mult(expr::SymPy.Sym)
    terms = args(expr)
    result = 1
    # if the expanded expression contains toplevel multiplication, the individual terms must all be exponentials or powers
    for term in terms
        # println("term = $term")
        # if term is exponential
        if SymPy.func(term) == SymPy.func(sympyExpExpr)
            termSeparated = separate_real_imaginary_exp(term)
        # if term is power (not sure if this case and the case below overlaps)
        elseif SymPy.func(term) == SymPy.func(sympyPowerExpr)
            termSeparated = separate_real_imaginary_power(term)
            # else, further split each product term into indivdual factors (this case also includes the case where term is a number, which would go into the "constant" below)
        else
            termSeparated = term # term is a number
#             (constant, factors) = factor_list(term)
#             termSeparated = constant
#             # separate each factor into real and imaginary parts and collect the product of separated factors
#             for (factor, power) in factors
#                 factor = factor^power
#                 termSeparated = termSeparated * (real(factor) + im*imag(factor))
#             end
        end
        # println("termSeparated = $termSeparated") 
        # collect the product of separated term, i.e., product of separated factors
        result = result * termSeparated
    end
    result = real(result) + im*imag(result)
    return result
end

x = symbols("x", real = true)
y = symbols("y", real = true)
expr = (x+im*y)*2x
println("func($expr) = $(SymPy.func(expr))")
separate_real_imaginary_mult(expr)

function separate_real_imaginary_add(expr::SymPy.Sym)
    x = symbols("x")
    # if the expanded expression contains toplevel addition, the individual terms must all be products or symbols
    terms = args(expr)
    result = 0
    # termSeparated = 0 # to avoid undefined error if there is no else (case incomplete)
    for term in terms
        # println("term = $term")
        # if term is a symbol
        if SymPy.func(term) == SymPy.func(x)
            termSeparated = term
        # if term is exponential
        elseif SymPy.func(term) == SymPy.func(sympyExpExpr)
            termSeparated = separate_real_imaginary_exp(term)
        # if term is a power
        elseif SymPy.func(term) == SymPy.func(sympyPowerExpr)
            termSeparated = separate_real_imaginary_power(term)
        # if term is a product
        elseif SymPy.func(term) == SymPy.func(sympyMultExpr)
            termSeparated = separate_real_imaginary_mult(term)
        # if term is a number
        else
            termSeparated = term
        end
        # println("termSeparated = $termSeparated")
        result = result + termSeparated
    end
    result = real(result) + im*imag(result)
    return result
end

x = symbols("x", real = true)
y = symbols("y", real = true)
expr = (x+im*y) + 2x
println("func($expr) = $(SymPy.func(expr))")
separate_real_imaginary_add(expr)

function separate_real_imaginary_power_mult_add(expr::SymPy.Sym)
    if SymPy.func(expr) == SymPy.func(sympyPowerExpr)
        result = separate_real_imaginary_power(expr)
    elseif SymPy.func(expr) == SymPy.func(sympyMultExpr)
        result = separate_real_imaginary_mult(expr)
        else #if SymPy.func(expr) == SymPy.func(sympyAddExpr)
        result = separate_real_imaginary_add(expr)
#     else
#         result = expr
    end
    return result
end

x = symbols("x", real = true)
y = symbols("y", real = true)
expr = (x+im*y) + 2x
println("func($expr) = $(SymPy.func(expr))")
separate_real_imaginary_power_mult_add(expr)

function separate_real_imaginary_others(expr::SymPy.Sym)
    # if the expanded expression is neither of the above, it must be a single term, e.g., x or cos(2x+1), which is a function wrapping around an expression; in this case, use the helper function to clean up the expression and feed it back to the function
    term = args(expr)[1]
    termCleaned = separate_real_imaginary_power_mult_add(term)
    result = subs(expr,args(expr)[1],termCleaned)
    result = real(result) + im*imag(result)
    return result
end

x = symbols("x", real = true)
y = symbols("y", real = true)
expr = cos(x+im*y)
println("func($expr) = $(SymPy.func(expr))")
separate_real_imaginary_others(expr)

function separate_real_imaginary(delta::SymPy.Sym)
    x = symbols("x", real = true)
    y = symbols("y", real = true)
    
    freeSymbols = free_symbols(delta)
    # check if delta has one and only one free symbol (e.g., global variable lambda)
    if length(freeSymbols) == 1
        lambda = freeSymbols[1]
        # substitute lambda with x+iy
        expr = subs(delta, lambda, x+im*y)
        # expand the new expression
        expr = expand(expr)
        
        if SymPy.func(expr) == SymPy.func(sympyPowerExpr)
#             println(expr)
#             println("power!")
            result = separate_real_imaginary_power(expr)
#             println("result = $result")
        elseif SymPy.func(expr) == SymPy.func(sympyAddExpr)
#             println(expr)
#             println("addition!")
            result = separate_real_imaginary_add(expr)
#             println("result = $result")
        elseif SymPy.func(expr) == SymPy.func(sympyMultExpr)
#             println(expr)
#             println("multiplication!")
            result = separate_real_imaginary_mult(expr)
#             println("result = $result")
        else
#             println(expr)
#             println("single term!")
            result = separate_real_imaginary_others(expr)
#             println("result = $result")
        end
        result = expand(result)
        return real(result) + im*imag(result)
        
    else
        throw("Delta has more than one variable")
    end
end

lambda = symbols("lambda")
delta = lambda + 1
separate_real_imaginary(delta)

delta = e^(lambda)
separate_real_imaginary(delta)

delta = lambda^2
separate_real_imaginary(delta)

delta = cos(lambda)
separate_real_imaginary(delta)

delta = cos(lambda)*e^(lambda)
separate_real_imaginary(delta)

delta = (lambda^3+lambda+2)*e^(lambda)
separatedDelta = separate_real_imaginary(delta)
prettyPrint(separatedDelta)
# Verified with WolframAlpha

function plot_levelCurves(bivariateDelta::SymPy.Sym; realFunc = real(bivariateDelta), imagFunc = imag(bivariateDelta), xRange = (-INFTY, INFTY), yRange = (-INFTY, INFTY), step = INFTY/1000, width = 1500, height = 1000)
    freeSymbols = free_symbols(bivariateDelta)
    x = symbols("x", real = true)
    y = symbols("y", real = true)
    
    xGridStep = (xRange[2] - xRange[1])/50
    yGridStep = (yRange[2] - yRange[1])/50
    
    if freeSymbols == [x, y]
        Plots.contour(xRange[1]:step:xRange[2], yRange[1]:step:yRange[2], realFunc, levels=[0], size = (width, height), tickfontsize = 20, seriescolor=:reds, transpose = false, linewidth = 4, linealpha = 1, xticks = xRange[1]:xGridStep:xRange[2], yticks = yRange[1]:yGridStep:yRange[2], grid = true, gridalpha = 0.5)
        Plots.contour!(xRange[1]:step:xRange[2], yRange[1]:step:yRange[2], imagFunc, levels=[0], size = (width, height), tickfontsize = 20, seriescolor=:blues, transpose = false, linewidth = 4, linealpha = 1, xticks = xRange[1]:xGridStep:xRange[2], yticks = yRange[1]:yGridStep:yRange[2], grid = true, gridalpha = 0.5)
    else
        Plots.contour(xRange[1]:step:xRange[2], yRange[1]:step:yRange[2], realFunc, levels=[0], size = (width, height), tickfontsize = 20, seriescolor=:reds, transpose = true, linewidth = 4, linealpha = 1, xticks = xRange[1]:xGridStep:xRange[2], yticks = yRange[1]:yGridStep:yRange[2], grid = true, gridalpha = 0.5)
        Plots.contour!(xRange[1]:step:xRange[2], yRange[1]:step:yRange[2], imagFunc, levels=[0], size = (width, height), tickfontsize = 20, seriescolor=:blues, transpose = true, linewidth = 4, linealpha = 1, xticks = xRange[1]:xGridStep:xRange[2], yticks = yRange[1]:yGridStep:yRange[2], grid = true, gridalpha = 0.5)
    end
end

lambda = symbols("lambda")
x = symbols("x", real = true)
y = symbols("y", real = true)
delta = lambda + 1
bivariateDelta = subs(delta, lambda, x+im*y)
# plot_levelCurves(separatedDelta) # somehow causes method error, probably because real(separatedDelta) is a function of x only and imag(separatedDelta) is a function of y only. In this case, we need to input the realFunc and imagFunc manually
# plot_levelCurves(bivariateDelta; realFunc = (x, y) -> x + 1, imagFunc = (x, y) -> y)

delta = (lambda^3+lambda+2)*e^(lambda)
bivariateDelta = subs(delta, lambda, x+im*y)
# plot_levelCurves(bivariateDelta)

function check_boundaryConditions(L::LinearDifferentialOperator, U::VectorBoundaryForm, fSym::Union{SymPy.Sym, Number})
    # Checks whether f satisfies the homogeneous boundary conditions
    Ux = get_Ux(L, U, fSym)
    return check_all(Ux, x->is_approx(x, 0))
end

t = symbols("t")
symPFunctions = [-1 0 0]
interval = (0,1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)
M = [1 0; 0 1]
N = [-1 0; 0 -1]
U = VectorBoundaryForm(M, N)
x = symbols("x")
fSym = sin(2*pi*x)
check_boundaryConditions(L, U, fSym)

function get_MPlusMinus(adjointU::VectorBoundaryForm; symbolic = false)
    # these are numeric matrices
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    if symbolic
        # return MPlus and MMinus as symbolic expressions with (the global variable) lambda as free variable
        lambda = symbols("lambda")
        alpha = e^(2*PI*im/n)
        MPlusMat = Array{SymPy.Sym}(n,n)
        for k = 1:n
            for j = 1:n
                sumPlus = 0
                for r = 0:(n-1)
                    summandPlus = (-im*alpha^(k-1)*lambda)^r * bStar[j,r+1]
                    summandPlus = expand(summandPlus)
                    sumPlus += summandPlus
                end
                # sumPlus = simplify(prettyPrint(sumPlus))
                MPlusMat[k,j] = sumPlus
            end
        end
        MPlusSym = MPlusMat
        MMinusMat = Array{SymPy.Sym}(n,n)
        for k = 1:n
            for j = 1:n
                sumMinus = 0
                for r = 0:(n-1)
                    summandMinus = (-im*alpha^(k-1)*lambda)^r * betaStar[j,r+1]
                    summandMinus = expand(summandMinus)
                    sumMinus += summandMinus
                end
                # sumMinus = simplify(prettyPrint(sumMinus))
                MMinusMat[k,j] = sumMinus
            end
        end
        MMinusSym = MMinusMat
        return (MPlusSym, MMinusSym)
    else
        alpha = e^(2pi*im/n)
        # if not symbolic, return MPlus and MMinus as functions of lambda
        function MPlus(lambda::Number)
            MPlusMat = Array{Number}(n,n)
            for k = 1:n
                for j = 1:n
                    sumPlus = 0
                    for r = 0:(n-1)
                        summandPlus = (-im*alpha^(k-1)*lambda)^r * bStar[j,r+1]
                        sumPlus += summandPlus
                    end
                    MPlusMat[k,j] = sumPlus
                end
            end
            return MPlusMat
        end
        function MMinus(lambda::Number)
            MMinusMat = Array{Number}(n,n)
            for k = 1:n
                for j = 1:n
                    sumMinus = 0
                    for r = 0:(n-1)
                        summandMinus = (-im*alpha^(k-1)*lambda)^r * betaStar[j,r+1]
                        sumMinus += summandMinus
                    end
                    MMinusMat[k,j] = sumMinus
                end
            end
            return MMinusMat
        end
    end
    return (MPlus, MMinus)
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)
M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)
adjointU = get_adjointU(L, U)

(MPlus, MMinus) = get_MPlusMinus(adjointU; symbolic = false)
(MPlusSym, MMinusSym) = get_MPlusMinus(adjointU; symbolic = true)
prettyPrint.(MPlusSym)

prettyPrint.(MMinusSym)

lambda = 1+im
println("MPlus($lambda) = $(MPlus(lambda))")
println("MMinus($lambda) = $(MMinus(lambda))")
println("MPlusSym($lambda) = $(prettyPrint.(evaluate.(MPlusSym, lambda)))")
println("MMinusSym($lambda) = $(prettyPrint.(evaluate.(MMinusSym, lambda)))")

function get_M(adjointU::VectorBoundaryForm; symbolic = false)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    if symbolic
        # return M as a symbolic expression with lambda as free variable
        lambda = symbols("lambda")
        alpha = e^(2*PI*im/n)
        (MPlusSym, MMinusSym) = get_MPlusMinus(adjointU; symbolic = true)
        MLambdaSym = Array{SymPy.Sym}(n,n)
        for k = 1:n
            for j = 1:n
                MLambdaSym[k,j] = simplify(MPlusSym[k,j] + MMinusSym[k,j] * e^(-im*alpha^(k-1)*lambda))
            end
        end
        MSym = MLambdaSym
        return MSym
    else
        alpha = e^(2pi*im/n)
        function M(lambda::Number)
            (MPlus, MMinus) = get_MPlusMinus(adjointU)
            MPlusLambda, MMinusLambda = MPlus(lambda), MMinus(lambda)
            MLambda = Array{Number}(n,n)
            for k = 1:n
                for j = 1:n
                    MLambda[k,j] = MPlusLambda[k,j] + MMinusLambda[k,j] * e^(-im*alpha^(k-1)*lambda)
                end
            end
            return MLambda
        end
        return M 
    end
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)

M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)

adjointU = get_adjointU(L, U)

MSym = prettyPrint.(get_M(adjointU; symbolic = true))

M = get_M(adjointU; symbolic = false)
lambda = 1+im
println("M($lambda) = $(M(lambda))")
println("MSym($lambda) = $(evaluate.(sym_to_func.(MSym), lambda))")

function get_delta(adjointU::VectorBoundaryForm; symbolic = false)
    if symbolic
        MSym = get_M(adjointU; symbolic = true)
        deltaSym = simplify(SymPy.det(MSym))
        return deltaSym
    else
       function delta(lambda::Number)
            M = get_M(adjointU)
            MLambda = convert(Array{Complex}, M(lambda))
            return det(MLambda)
        end
        return delta 
    end
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)

M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)

adjointU = get_adjointU(L, U)

deltaSym = prettyPrint(get_delta(adjointU; symbolic = true))

delta = get_delta(adjointU; symbolic = false)
lambda = 1+im
println("delta($lambda) = $(delta(lambda))")
println("deltaSym($lambda) = $(evaluate(deltaSym, lambda))")

function get_Xlj(adjointU::VectorBoundaryForm, l::Number, j::Number; symbolic = false)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    if symbolic
        MSym = get_M(adjointU; symbolic = true)
        MBlockSym = [MSym MSym; MSym MSym]
        XljSym = MBlockSym[(l+1):(l+1+n-2), (j+1):(j+1+n-2)]
        return XljSym
    else
        M = get_M(adjointU; symbolic = false)
        function Xlj(lambda::Number)
            MLambda = M(lambda)
            MLambdaBlock = [MLambda MLambda; MLambda MLambda]
            XljLambda = MLambdaBlock[(l+1):(l+1+n-2), (j+1):(j+1+n-2)]
            return XljLambda
        end
        return Xlj 
    end
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)

M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)

adjointU = get_adjointU(L, U)

l, j = 1, 1
XljSym = prettyPrint.(get_Xlj(adjointU, l, j; symbolic = true))

Xlj = get_Xlj(adjointU, l, j; symbolic = false)
lambda = 1+im
println("Xlj($lambda) = $(Xlj(lambda))")
println("XljSym($lambda) = $(evaluate.(XljSym, lambda))")

function get_alpha(m::Int, n::Int)
    if m == 1
        result = (-1)^n
    elseif m == 2
        result = (-1)^(n+1)*n^2
    else
        sum = 0
        for k = 1:(n-m+2)
            product = 1
            for j = k:(m+k-3)
                product *= (n-j)
            end
            sum += binomial(m+k-3, k-1) * product
        end
        result = (-1)^(n+m-1)*2^(m-2)*n*sum
    end
    return result
end

n = 3
println("alpha(1,$n) = $(get_alpha(1, n))")
println("alpha(2,$n) = $(get_alpha(2, n))")
println("alpha(3,$n) = $(get_alpha(3, n))")

function get_ChebyshevTermIntegral(n::Int; symbolic = true, c = symbols("c"))
    if symbolic
        if c == 0
            if n == 0
                expr = 2
            elseif n == 1
                expr = 0
            else
                expr = ((-1)^(n+1)-1)/(n^2-1)
            end
        else
            expr = 0
            for m = 1:(n+1)
                summand = get_alpha(m, n) * (e^(im*c)/(im*c)^m  + (-1)^(m+n)*e^(-im*c)/(im*c)^m)
                expr += summand
            end
        end
        expr = simplify(expr)
        return expr
    else
        function TTilde(c)
            if c == 0
                if n == 0
                    result = 2
                elseif n == 1
                    result = 0
                else
                    result = ((-1)^(n+1)-1)/(n^2-1)
                end
            else
                result = 0
                for i = 1:(n+1)
                    summand = get_alpha(i, n) * (e^(im*c)/(im*c)^i  + (-1)^(i+n)*e^(-im*c)/(im*c)^i)
                    result += summand
                end
            end
            return result
        end
        return TTilde
    end
end

# alpha = e^(2pi*im/1)
# lambda = symbols("lambda")
# c = alpha^(l-1)*lambda/2
chebTermIntSym = get_ChebyshevTermIntegral(2; symbolic = true)
prettyPrint(chebTermIntSym)

chebTermInt = get_ChebyshevTermIntegral(2; symbolic = false)
alpha = e^(2pi*im/1)
lambda = 1+im
l = 2
c = alpha^(l-1)*lambda/2
println("chebTermInt($(prettyPrint(c))) = $(chebTermInt(c))")
println("chebTermIntSym($(prettyPrint(c))) = $(evaluate(chebTermIntSym, c))")

function get_ChebyshevCoefficients(f::Union{Function,Number})
    fCheb = ApproxFun.Fun(f, 0..1) # Approximate f on [0,1] using chebyshev polynomials
    chebCoefficients = ApproxFun.coefficients(fCheb) # get coefficients of the Chebyshev polynomial
    return chebCoefficients
end

f(x) = x^2+1
fChebCoeffs = get_ChebyshevCoefficients(f)

function get_ChebyshevIntegral(l::Number, f::Union{Function, Number}; symbolic = false, lambda = nothing, alpha = nothing)
    fChebCoeffs = get_ChebyshevCoefficients(f)
    # Replace coefficients too close to 0 by 0
    # fChebCoeffs = [if is_approx(x, 0) 0 else x end for x in fChebCoeffs]
    if symbolic
        lambda = symbols("lambda")
        c = alpha^(l-1)*lambda/2
        integralSym = 0
        for m = 1:length(fChebCoeffs)
            fChebCoeff = fChebCoeffs[m]
            if is_approx(fChebCoeff, 0)
                continue
            else
            integralSym += fChebCoeffs[m] * get_ChebyshevTermIntegral(m-1; symbolic = true, c = c)
            end
        end
        integralSym = integralSym/(2*e^(im*c))
        integralSym = simplify(integralSym)
        return integralSym
    else
        if isa(lambda, Void) || isa(alpha, Void)
            throw("lambda, alpha required")
        else
            c = alpha^(l-1)*lambda/2
            integral = 0
            for m = 1:length(fChebCoeffs)
                fChebCoeff = fChebCoeffs[m]
                if is_approx(fChebCoeff, 0)
                    continue
                else
                    integral += fChebCoeff * get_ChebyshevTermIntegral(m-1; symbolic = false)(c)
                end
            end
            integral = integral/(2*e^(im*c))
            return integral
        end
    end
end

alpha = e^(2pi*im/1)
l = 2
# f(x) = x^2+1
f(x) = sin(2*pi*x)

chebIntSym = simplify(get_ChebyshevIntegral(l, f; symbolic = true, alpha = alpha))
prettyPrint(chebIntSym)

lambda = 1+im
chebInt = get_ChebyshevIntegral(l, f; symbolic = false, lambda = lambda, alpha = alpha)

# Directly compute the integral
g(x) = e^(-im*alpha^(l-1)*lambda*x)
directInt = quadgk(mult_func(g,f), 0, 1)[1]

println("chebIntegral($lambda) = $chebInt")
println("chebIntegralSym($lambda) = $(evaluate.(chebIntSym, lambda))")
println("directIntegral = $directInt")

function get_FPlusMinus(adjointU::VectorBoundaryForm; symbolic = false)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    if symbolic
        alpha = e^(2*PI*im/n)
        lambda = symbols("lambda")
        (MPlusSym, MMinusSym) = get_MPlusMinus(adjointU; symbolic = true)
        deltaSym = get_delta(adjointU; symbolic = true)
        function FPlusSym(f::Union{Function, Number})
            sumPlusSym = 0
            for l = 1:n
                summandPlusSym = 0
                for j = 1:n
                    XljSym = get_Xlj(adjointU, l, j; symbolic = true)
                    integralSym = get_ChebyshevIntegral(l, f; symbolic = true, lambda = lambda, alpha = alpha)
                    summandPlusSym = summandPlusSym += (-1)^((n-1)*(l+j)) * SymPy.det(XljSym) * MPlusSym[1,j] * integralSym
                end
                sumPlusSym += summandPlusSym
            end
            return simplify(1/(2*PI*deltaSym)*sumPlusSym)
        end
        function FMinusSym(f::Union{Function, Number})
            sumMinusSym = 0
            for l = 1:n
                summandMinusSym = 0
                for j = 1:n
                    XljSym = get_Xlj(adjointU, l, j; symbolic = true)
                    c = alpha^(l-1)*lambda/2
                    integralSym = get_ChebyshevIntegral(l, f; symbolic = true, lambda = lambda, alpha = alpha)
                    summandMinusSym = summandMinusSym += (-1)^((n-1)*(l+j)) * SymPy.det(XljSym) * MMinusSym[1,j] * integralSym
                end
                sumMinusSym += summandMinusSym
            end
            return simplify((-e^(-im*lambda))/(2*PI*deltaSym)*sumMinusSym)
        end
        return (FPlusSym, FMinusSym)
    else
        alpha = e^(2pi*im/n)
        (MPlus, MMinus) = get_MPlusMinus(adjointU; symbolic = false)
        function FPlus(lambda::Number, f::Union{Function, Number})
            MPlusLambda, MMinusLambda = MPlus(lambda), MMinus(lambda)
            M = get_M(adjointU; symbolic = false)
            MLambda = convert(Array{Complex}, M(lambda))
            deltaLambda = det(MLambda) # or deltaLambda = (get_delta(adjointU))(lambda)
            sumPlus = 0
            for l = 1:n
                summandPlus = 0
                for j = 1:n
                    Xlj = get_Xlj(adjointU, l, j; symbolic = false)
                    XljLambda = convert(Array{Complex}, Xlj(lambda))
                    integral = get_ChebyshevIntegral(l, f; symbolic = false, lambda = lambda, alpha = alpha)
                    summandPlus = summandPlus += (-1)^((n-1)*(l+j)) * det(XljLambda) * MPlusLambda[1,j] * integral
                end
                sumPlus += summandPlus
            end
            return 1/(2pi*deltaLambda)*sumPlus
        end
        function FMinus(lambda::Number, f::Union{Function, Number})
            MPlusLambda, MMinusLambda = MPlus(lambda), MMinus(lambda)
            M = get_M(adjointU; symbolic = false)
            MLambda = convert(Array{Complex}, M(lambda))
            deltaLambda = det(MLambda) # or deltaLambda = (get_delta(adjointU))(lambda)
            sumMinus = 0
            for l = 1:n
                summandMinus = 0
                for j = 1:n
                    Xlj = get_Xlj(adjointU, l, j)
                    XljLambda = convert(Array{Complex}, Xlj(lambda))
                    integral = get_ChebyshevIntegral(l, f; symbolic = false, lambda = lambda, alpha = alpha)
                    summandMinus = summandMinus += (-1)^((n-1)*(l+j)) * det(XljLambda) * MMinusLambda[1,j] * integral
                end
                sumMinus += summandMinus
            end
            return (-e^(-im*lambda))/(2pi*deltaLambda)*sumMinus
        end
        return (FPlus, FMinus)
    end
end

t = symbols("t")
symPFunctions = [1 t+1 t^2+t+1]
interval = (0, 1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)

M = [1 0; 2 0]
N = [0 2; 0 1]
U = VectorBoundaryForm(M, N)

adjointU = get_adjointU(L, U)

f(x) = x^2+1

# Keep lambda as a free symbol
(FPlusSym, FMinusSym) = get_FPlusMinus(adjointU; symbolic = true)
prettyPrint(FPlusSym(f))

lambda = 1+im
(FPlus, FMinus) = get_FPlusMinus(adjointU; symbolic = false)
println("FPlus($lambda, f) = $(FPlus(lambda, f))")
println("FMinus($lambda, f) = $(FMinus(lambda, f))")
println("FPlusSym($lambda, f) = $(evaluate.(FPlusSym(f), lambda))")
println("FMinusSym($lambda, f) = $(evaluate.(FMinusSym(f), lambda))")

function get_gammaAAngles(a::Number, n::Int; symbolic = false)
    # thetaA = argument(a)
    thetaA = angle(a)
    if symbolic
        thetaStartList = Array{SymPy.Sym}(n) # List of angles that characterize where domain sectors start
        thetaEndList = Array{SymPy.Sym}(n) # List of angles that characterize where domain sectors end
        k = symbols("k")
        counter = 0
        while (2pi*counter + pi/2 - thetaA)/n < 2pi
        # Substituting counter for k
        # while SymPy.N(subs((2PI*k + PI/2 - thetaA)/n, k, counter)) < 2pi
            thetaStart = (2*PI*counter - PI/2 - rationalize(thetaA/pi)*PI)/n
            thetaEnd = (2*PI*counter + PI/2 - rationalize(thetaA/pi)*PI)/n
            counter += 1
            thetaStartList[counter] = thetaStart
            thetaEndList[counter] = thetaEnd
        end
    else
        thetaStartList = Array{Number}(n)
        thetaEndList = Array{Number}(n)
        k = 0
        while (2pi*k + pi/2 - thetaA)/n < 2pi
            thetaStart = (2pi*k - pi/2 - thetaA)/n
            thetaEnd = (2pi*k + pi/2 - thetaA)/n
            k += 1
            thetaStartList[k] = thetaStart
            thetaEndList[k] = thetaEnd
        end
    end
    return (thetaStartList, thetaEndList)
end

n = 2
a = 1
gammaAAnglesSym = get_gammaAAngles(a, n; symbolic = true)
println("n=$n, a=$a, gammaAAnglesSym = $gammaAAnglesSym")

gammaAAngles = get_gammaAAngles(a, n; symbolic = false)
println("n=$n, a=$a, gammaAAngles = $gammaAAngles")

n = 3
a = im
gammaAAnglesSym = get_gammaAAngles(a, n; symbolic = true)
println("n=$n, a=$a, gammaAAnglesSym = $gammaAAnglesSym")

n = 3
a = -im
gammaAAnglesSym = get_gammaAAngles(a, n; symbolic = true)
println("n=$n, a=$a, gammaAAnglesSym = $gammaAAnglesSym")

function get_gammaAAnglesSplit(a::Number, n::Int; symbolic = false)
    (thetaStartList, thetaEndList) = get_gammaAAngles(a, n; symbolic = symbolic)
    # Split sectors that contain the positive half of the real line (angle = 0)
    zeroIndex = find(i -> ((is_approxLess(thetaStartList[i], 0) && is_approxLess(0, thetaEndList[i]))), 1:length(thetaStartList))
    if !isempty(zeroIndex)
        index = zeroIndex[1]
        # Insert 0 after thetaStart
        splice!(thetaStartList, (index+1):index, 0)
        # Insert 0 before thetaEnd
        splice!(thetaEndList, index:(index-1), 0)
    end
    # Split sectors that contain the negative half of the real line (angle = pi)
    piIndex = find(i -> ((is_approxLess(thetaStartList[i], pi) && is_approxLess(pi, thetaEndList[i]))), 1:length(thetaStartList))
    if !isempty(piIndex)
        index = piIndex[1]
        if symbolic
            # Insert pi after thetaStart
            splice!(thetaStartList, (index+1):index, pi)
            # Insert pi before thetaEnd
            splice!(thetaEndList, index:(index-1), pi)
        else
            # Use pi*1 instead of pi to avoid "<= not defined for Irrational{:pi}" error in get_gamma()
            splice!(thetaStartList, (index+1):index, pi*1)
            splice!(thetaEndList, index:(index-1), pi*1)
        end
    end
    return (thetaStartList, thetaEndList)
end

n = 2
a = 1
gammaAAnglesSplitSym = get_gammaAAnglesSplit(a, n; symbolic = true)
println("n=$n, a=$a, gammaAAnglesSplitSym = $gammaAAnglesSplitSym")
gammaAAnglesSplit = get_gammaAAnglesSplit(a, n; symbolic = false)
println("n=$n, a=$a, gammaAAnglesSplit = $gammaAAnglesSplit")

n = 3
a = im
gammaAAnglesSplitSym = get_gammaAAnglesSplit(a, n; symbolic = true)
println("n=$n, a=$a, gammaAAnglesSplitSym = $gammaAAnglesSplitSym")
gammaAAnglesSplit = get_gammaAAnglesSplit(a, n; symbolic = false)
println("n=$n, a=$a, gammaAAnglesSplit = $gammaAAnglesSplit")

function pointOnSector(z::Number, sectorAngles::Tuple{Number, Number})
    (startAngle, endAngle) = sectorAngles
    return is_approx(argument(z), startAngle) || is_approx(argument(z), endAngle) || is_approx(angle(z), startAngle) || is_approx(angle(z), endAngle)
end

z = 1+im
sectorAngles = (-pi/4, pi/4)
pointOnSector(z, sectorAngles)

function pointInSector(z::Number, sectorAngles::Tuple{Number, Number})
    (startAngle, endAngle) = sectorAngles
    # First check if z is on the sector boundary
    if pointOnSector(z, sectorAngles)
        return false
    else
        # angle(z) would work if it's in the sector with positive real parts and both positive and negative imaginary parts; argument(z) would work if it's in the sector with negative real parts and both positive and negative imaginary parts
        return (angle(z) > startAngle && angle(z) < endAngle) || (argument(z) > startAngle && argument(z) < endAngle) # no need to use is_approxLess because the case of approximatedly equal is already checked in pointOnSector
    end
end

z = 1+im
sectorAngles = (-pi/4, pi/4)
println("pointInSector($z, $sectorAngles) = $(pointInSector(z, sectorAngles))")

z = 1
println("pointInSector($z, $sectorAngles) = $(pointInSector(z, sectorAngles))")

z = -1
println("pointInSector($z, $sectorAngles) = $(pointInSector(z, sectorAngles))")

function pointExSector(z::Number, sectorAngles::Tuple{Number, Number})
    return !pointOnSector(z, sectorAngles) && !pointInSector(z, sectorAngles)
end

z = 1+im
sectorAngles = (-pi/4, pi/4)
println("pointExSector($z, $sectorAngles) = $(pointExSector(z, sectorAngles))")

z = 1
println("pointExSector($z, $sectorAngles) = $(pointExSector(z, sectorAngles))")

z = -1
println("pointExSector($z, $sectorAngles) = $(pointExSector(z, sectorAngles))")

function get_epsilon(zeroList::Array, a::Number, n::Int)
    (thetaStartList, thetaEndList) = get_gammaAAnglesSplit(a, n; symbolic = false)
    thetaStartEndList = collect(Iterators.flatten([thetaStartList, thetaEndList]))
    truncZeroList = []
    for zero in zeroList
        # If zero is interior to any sector, discard it
        if any(i -> pointInSector(zero, (thetaStartList[i], thetaEndList[i])), 1:n)
        else # If not, append it to truncZeroList
            append!(truncZeroList, zero)
        end
    end
    # List of distance between each zero and each line marking the boundary of some sector
    pointLineDistances = [get_distancePointLine(z, theta) for z in zeroList for theta in thetaStartEndList]
    if length(truncZeroList)>1
        # List of distance between every two zeroes
        pairwiseDistances = [norm(z1-z2) for z1 in zeroList for z2 in truncZeroList]
    else
        pairwiseDistances = []
    end
    distances = collect(Iterators.flatten([pairwiseDistances, pointLineDistances]))
    # Distances of nearly 0 could be instances where the zero is actually on some sector boundary
    distances = filter(x -> !is_approx(x, 0), distances)
    epsilon = minimum(distances)/3
    return epsilon
end

n = 2
a = 1
zeroList = [1+sqrt(3)*im, 2+2*sqrt(3)*im, 0+0*im, 0+5*im, 0-5*im]
get_epsilon(zeroList, a, n)

function get_nGonAroundZero(zero::Number, epsilon::Number, n::Int)
    z = zero
    theta = argument(zero)
    deltaAngle = 2pi/n
    vertices = []
    for i = 1:n
        newAngle = pi-deltaAngle*(i-1)
        vertex = z + epsilon*e^(im*(theta+newAngle))
        append!(vertices, vertex)
    end
    # vertices = vcat(vertices, vertices[1])
    return vertices
end

zero = im
epsilon = 1
n = 4
prettyPrint.(get_nGonAroundZero(zero, epsilon, n))

function get_gamma(a::Number, n::Int, zeroList::Array; infty = INFTY, nGon = 8)
    (thetaStartList, thetaEndList) = get_gammaAAnglesSplit(a, n; symbolic = false)
    nSplit = length(thetaStartList)
    gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus = [], [], [], []
    epsilon = get_epsilon(zeroList, a, n)
    for i in 1:nSplit
        thetaStart = thetaStartList[i]
        thetaEnd = thetaEndList[i]
        # Initialize the boundary of each sector with the ending boundary, the origin, and the starting boundary (start and end boundaries refer to the order in which the boundaries are passed if tracked counterclockwise)
        initialPath = [infty*e^(im*thetaEnd), 0+0*im, infty*e^(im*thetaStart)]
        if thetaStart >= 0 && thetaStart <= pi && thetaEnd >= 0 && thetaEnd <= pi # if in the upper half plane, push the boundary path to gamma_a+
            push!(gammaAPlus, initialPath) # list of lists
        else # if in the lower half plane, push the boundary path to gamma_a-
            push!(gammaAMinus, initialPath)
        end
    end
    # Sort the zeroList by norm, so that possible zero at the origin comes last. We need to leave the origin in the initial path unchanged until we have finished dealing with all non-origin zeros because we use the origin in the initial path as a reference point to decide where to insert the deformed path
    zeroList = sort(zeroList, lt=(x,y)->!isless(norm(x), norm(y)))
    for zero in zeroList
        # println(zero)
        # If zero is not at the origin
        if !is_approx(zero, 0+0*im)
            # Draw an n-gon around it
            vertices = get_nGonAroundZero(zero, epsilon, nGon)
            # If zero is on the boundary of some sector
            if any(i -> pointOnSector(zero, (thetaStartList[i], thetaEndList[i])), 1:nSplit)
                # Find which sector(s) zero is on
                indices = find(i -> pointOnSector(zero, (thetaStartList[i], thetaEndList[i])), 1:nSplit)
                # If zero is on the boundary of one sector
                if length(indices) == 1
                    # if vertices[2] is interior to any sector, include vertices on the other half of the n-gon in the contour approximation
                    z0 = vertices[2]
                    if any(i -> pointInSector(z0, (thetaStartList[i], thetaEndList[i])), 1:nSplit)
                        # Find which sector vertices[2] is in
                        index = find(i -> pointInSector(z0, (thetaStartList[i], thetaEndList[i])), 1:nSplit)[1]
                    else # if vertices[2] is exterior, include vertices on this half of the n-gon in the contour approximation
                        # Find which sector vertices[length(vertices)] is in
                        z1 = vertices[length(vertices)]
                        index = find(i -> pointInSector(z1, (thetaStartList[i], thetaEndList[i])), 1:nSplit)[1]
                    end
                    thetaStart = thetaStartList[index]
                    thetaEnd = thetaEndList[index]
                    # Find all vertices exterior to or on the boundary of this sector, which would form the nGonPath around the zero
                    nGonPath = vertices[find(vertex -> !pointInSector(vertex, (thetaStart, thetaEnd)), vertices)]
                    # If this sector is in the upper half plane, deform gamma_a+
                    if thetaStart >= 0 && thetaStart <= pi && thetaEnd >= 0 && thetaEnd <= pi
                        gammaAPlusIndex = find(path -> (is_approx(argument(zero), argument(path[1])) || is_approx(argument(zero), argument(path[length(path)]))), gammaAPlus)[1]
                        deformedPath = copy(gammaAPlus[gammaAPlusIndex])
                        if any(i -> is_approx(argument(zero), thetaStartList[i]) || is_approx(angle(zero), thetaStartList[i]), 1:nSplit) # if zero is on the starting boundary, insert the n-gon path after 0+0*im
                            splice!(deformedPath, length(deformedPath):(length(deformedPath)-1), nGonPath)
                        else # if zero is on the ending boundary, insert the n-gon path before 0+0*im
                            splice!(deformedPath, 2:1, nGonPath)
                        end
                        gammaAPlus[gammaAPlusIndex] = deformedPath
                    else # if sector is in the lower half plane, deform gamma_a-
                        # # Find all vertices interior to or on the boundary of this sector, which would form the nGonPath around the zero
                        # nGonPath = vertices[find(vertex -> !pointExSector(vertex, (thetaStart, thetaEnd)), vertices)]
                        gammaAMinusIndex = find(path -> (is_approx(argument(zero), argument(path[1])) || is_approx(argument(zero), argument(path[length(path)]))), gammaAMinus)[1]
                        deformedPath = copy(gammaAMinus[gammaAMinusIndex])
                        if any(i -> is_approx(argument(zero), thetaStartList[i]) || is_approx(angle(zero), thetaStartList[i]), 1:nSplit) # if zero is on the starting boundary, insert the n-gon path after 0+0*im
                            splice!(deformedPath, length(deformedPath):(length(deformedPath)-1), nGonPath)
                        else # if zero is on the ending boundary, insert the n-gon path before 0+0*im
                            splice!(deformedPath, 2:1, nGonPath) 
                        end
                        gammaAMinus[gammaAMinusIndex] = deformedPath
                    end
                else # If zero is on the boundary of two sectors, then it must be on the real line, and we need to deform two sectors
                    # Find out which vertices are in the lower half plane
                    nGonPath = vertices[find(vertex -> !pointInSector(vertex, (0, pi)), vertices)]
                    for index in indices
                        thetaStart = thetaStartList[index]
                        thetaEnd = thetaEndList[index]
                        # If this is the sector in the upper half plane, deform gamma_a+
                        if thetaStart >= 0 && thetaStart <= pi && thetaEnd >= 0 && thetaEnd <= pi
                            gammaAPlusIndex = find(path -> (is_approx(argument(zero), argument(path[1])) || is_approx(argument(zero), argument(path[length(path)]))), gammaAPlus)[1]
                            deformedPath = copy(gammaAPlus[gammaAPlusIndex])
                            if is_approx(argument(zero), argument(deformedPath[length(deformedPath)])) # if zero is on the starting boundary, insert the n-gon path after 0+0*im
                                splice!(deformedPath, length(deformedPath):(length(deformedPath)-1), nGonPath)
                            else # if zero is on the ending boundary, insert the n-gon path before 0+0*im
                                splice!(deformedPath, 2:1, nGonPath)
                            end
                            gammaAPlus[gammaAPlusIndex] = deformedPath
                        else # If this is the sector in the lower half plane, deform gamma_a-
                            gammaAMinusIndex = find(path -> (is_approx(argument(zero), argument(path[1])) || is_approx(argument(zero), argument(path[length(path)]))), gammaAMinus)[1]
                            deformedPath = copy(gammaAMinus[gammaAMinusIndex])
                            if is_approx(argument(zero), argument(deformedPath[length(deformedPath)])) # if zero is on the starting boundary, insert the n-gon path after 0+0*im
                                splice!(deformedPath, length(deformedPath):(length(deformedPath)-1), nGonPath)
                            else # if zero is on the ending boundary, insert the n-gon path before 0+0*im
                                splice!(deformedPath, 2:1, nGonPath)
                            end
                            gammaAMinus[gammaAMinusIndex] = deformedPath
                        end
                    end
                end
                # Sort each sector's path in the order in which they are integrated over
                gammaAs = [gammaAPlus, gammaAMinus]
                for j = 1:length(gammaAs)
                    gammaA = gammaAs[j]
                    for k = 1:length(gammaA)
                        inOutPath = gammaA[k]
                        originIndex = find(x->x==0+0*im, inOutPath)[1]
                        inwardPath = inOutPath[1:(originIndex-1)]
                        outwardPath = inOutPath[(originIndex+1):length(inOutPath)]
                        # Sort the inward path and outward path
                        if length(inwardPath) > 0
                            inwardPath = sort(inwardPath, lt=(x,y)->!isless(norm(x), norm(y)))
                        end
                        if length(outwardPath) > 0
                            outwardPath = sort(outwardPath, lt=(x,y)->isless(norm(x), norm(y)))
                        end
                        inOutPath = vcat(inwardPath, 0+0*im, outwardPath)
                        gammaA[k] = inOutPath
                    end
                    gammaAs[j] = gammaA 
                end
                gammaAPlus, gammaAMinus = gammaAs[1], gammaAs[2]
            # If zero is interior to any sector (after splitting by real line), ignore it
            # If zero is exterior to the sectors, avoid it
            elseif all(i -> pointExSector(zero, (thetaStartList[i], thetaEndList[i])), 1:nSplit)
                nGonPath = vcat(vertices, vertices[1]) # counterclockwise
                # If zero is in the upper half plane, add the n-gon path to gamma_0+
                if argument(zero) >= 0 && argument(zero) <= pi
                    push!(gamma0Plus, nGonPath)
                else # If zero is in the lower half plane, add the n-gon path to gamma_0-
                    push!(gamma0Minus, nGonPath)
                end
            end
        else # If zero is at the origin, we deform all sectors and draw an n-gon around the origin
            # deform each sector in gamma_a+
            for i = 1:length(gammaAPlus)
                deformedPath = gammaAPlus[i]
                # find the index of the zero at origin in the sector boundary path
                index = find(j -> is_approx(deformedPath[j], 0+0*im), 1:length(deformedPath))
                # If the origin is not in the path, then it has already been bypassed
                if isempty(index)
                else # If not, find its index
                    index = index[1]
                end
                # create a path around zero (origin); the origin will not be the first or the last point in any sector boundary because it was initialized to be in the middle, and only insertions are performed. Moreover, the boundary path has already been sorted into the order in which they will be integrated over, so squarePath defined below has deformedPath[index-1], deformedPath[index+1] in the correct order.
                squarePath = [epsilon*e^(im*argument(deformedPath[index-1])), epsilon*e^(im*argument(deformedPath[index+1]))]
                # replace the zero with the deformed path
                deleteat!(deformedPath, index) # delete the origin
                splice!(deformedPath, index:(index-1), squarePath) # insert squarePath into where the origin was at
                gammaAPlus[i] = deformedPath
            end
            # deform each sector in gamma_a-
            for i = 1:length(gammaAMinus)
                deformedPath = gammaAMinus[i]
                index = find(j -> is_approx(deformedPath[j], 0+0*im), 1:length(deformedPath))
                if isempty(index)
                else
                    index = index[1]
                end
                squarePath = [epsilon*e^(im*argument(deformedPath[index-1])), epsilon*e^(im*argument(deformedPath[index+1]))]
                deleteat!(deformedPath, index)
                splice!(deformedPath, index:(index-1), squarePath)
                gammaAMinus[i] = deformedPath
            end
            # Draw an n-gon around the origin and add to gamma_0+
            vertices = get_nGonAroundZero(zero, epsilon/2, nGon)
            nGonPath = vcat(vertices, vertices[1])
            push!(gamma0Plus, nGonPath)
        end
    end
    return (gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus)
end

n = 2
a = 1
zeroList = [3+3*sqrt(3)*im, 2+2*sqrt(3)*im, 0+0*im, 0+5*im, 0-5*im]
(gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus) = get_gamma(a, n, zeroList)

println("gammaAPlus = $gammaAPlus")
println("gammaAMinus = $gammaAMinus")
println("gamma0Plus = $gamma0Plus")
println("gamma0Minus = $gamma0Minus")

function plot_contour(gamma::Array; infty = INFTY)
    sectorPathList = Array{Any}(length(gamma),1)
    for i = 1:length(gamma)
        # For each sector path in the gamma contour, plot the points in the path and connect them in the order in which they appear in the path
        sectorPath = gamma[i]
        # labels = map(string, collect(1:1:length(sectorPath)))
        sectorPathList[i] = layer(x = real(sectorPath), y = imag(sectorPath), Geom.line(preserve_order=true))
    end
    coord = Coord.cartesian(xmin=-infty, xmax=infty, ymin=-infty, ymax=infty, fixed=true)
    Gadfly.plot(Guide.xlabel("Re"), Guide.ylabel("Im"), coord, sectorPathList...)
end

n = 3
a = im
zeroList = [3+3*sqrt(3)*im, 2+2*sqrt(3)*im, 0+0*im, 0+5*im, 0-5*im, 4]
(gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus) = get_gamma(a, n, zeroList)
gamma = collect(Iterators.flatten([gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus]))

plot_contour(gamma)

n = 2
a = 1
(gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus) = get_gamma(a, n, zeroList)
gamma = collect(Iterators.flatten([gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus]))

plot_contour(gamma)

n = 4
a = 1
(gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus) = get_gamma(a, n, zeroList)
gamma = collect(Iterators.flatten([gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus]))

plot_contour(gamma)

function solve_IBVP(L::LinearDifferentialOperator, U::VectorBoundaryForm, a::Number, zeroList::Array, f::Function; FPlusFunc = lambda->get_FPlusMinus(adjointU; symbolic = false)[1](lambda, f), FMinusFunc = lambda->get_FPlusMinus(adjointU; symbolic = false)[2](lambda, f), pDerivMatrix = get_pDerivMatrix(L), infty = INFTY)
    n = length(L.pFunctions)-1
    adjointU = get_adjointU(L, U, pDerivMatrix)
    (gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus) = get_gamma(a, n, zeroList)
    function q(x,t)
        integrandPlus(lambda) = e^(im*lambda*x)*e^(-a*lambda^n*t) * FPlusFunc(lambda)
        integrandMinus(lambda) = e^(im*lambda*x)*e^(-a*lambda^n*t) * FMinusFunc(lambda)
        println("integrandPlus = $(integrandPlus(1+im))")
        println("integrandMinus = $(integrandMinus(1+im))")
        # Integrate over individual paths in the Gamma contours
        # println("gamma0Plus = $gamma0Plus")
        integralGamma0Plus = 0
        for path in gamma0Plus
            if length(path) == 0
                path = [im,im]
            end
            integralGamma0Plus += quadgk(integrandPlus, path...)[1]
        end
        println("int_0_+ = $integralGamma0Plus")
        # println("gammaAPlus = $gammaAPlus")
        integralGammaAPlus = 0
        for path in gammaAPlus
            # println("path = $path")
            if length(path) == 0
                path = [im,im]
            end
#             pathIntegral = 0
#             for i in 1:(length(path)-1)
#                 pathIntegral += quadgk(integrandPlus, path[i], path[i+1])[1]
#             end
            integralGammaAPlus += quadgk(integrandPlus, path...)[1]
            # integralGammaAPlus += pathIntegral
        end
        println("int_a_+ = $integralGammaAPlus")
        integralGamma0Minus = 0
        for path in gamma0Minus
            if length(path) == 0
                path = [-im,-im]
            end
            integralGamma0Minus += quadgk(integrandMinus, path...)[1]
        end
        println("int_0_- = $integralGamma0Minus")
        integralGammaAMinus = 0
        for path in gammaAMinus
            if length(path) == 0
                path = [-im,-im]
            end
            integralGammaAMinus += quadgk(integrandMinus, path...)[1]
        end
        println("int_a_- = $integralGammaAMinus")
        return (integralGamma0Plus + integralGammaAPlus + integralGamma0Minus + integralGammaAMinus)
    end
    return q
end

n = 2
beta0 = -1
beta1 = -1
theta = 0
a = e^(im*theta)

t = symbols("t")
symPFunctions = [-1 0 0]
interval = (0,1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)
b = [1 0; 0 1]
beta = [beta0 0; 0 beta1]
U = VectorBoundaryForm(b, beta)
f(x) = sin(x*2*pi)
x = symbols("x")
fSym = sin(x*2*PI)
# fSym = x^2-1
# f(x) = x^2-1
check_boundaryConditions(L, U, fSym)

adjointU = get_adjointU(L, U)
delta = get_delta(adjointU; symbolic = true)
lambda = free_symbols(delta)[1]
x = symbols("x", real = true)
y = symbols("y", real = true)
bivariateDelta = subs(delta, lambda, x+im*y)
plot_levelCurves(bivariateDelta; width = 2500, height = 2000)

zeroList = [0, 2*pi, -2*pi]
(gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus) = get_gamma(a, n, zeroList)
gamma = collect(Iterators.flatten([gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus]))
plot_contour(gamma)

(FPlusSym, FMinusSym) = get_FPlusMinus(adjointU; symbolic = true)
FPlusSymF = prettyPrint(simplify(FPlusSym(f)))
# FPlusSymF = simplify(FPlusSym(f))

function FPlusSymFunc(lambda)
    result = (-pi*im*lambda^(21)*e^(2*im*lambda) + 2*pi*im*lambda^(21)*e^(im*lambda) - pi*im*lambda^(21) + 2*pi*lambda^(20)*e^(2*im*lambda) - 2*pi*lambda^(20) - 124.025*im*lambda^(19)*e^(2*im*lambda) + 248.05*im*lambda^19*e^(im*lambda) - 124.025*im*lambda^(19) + 248.05*lambda^(18)*e^(2*im*lambda) - 248.05*lambda^(18) - 4896.315*im*lambda^(17)*e^(2*im*lambda) + 9792.63*im*lambda^(17)*e^(im*lambda) - 4896.315*im*lambda^(17) + 9792.628*lambda^(16)*e^(2*im*lambda) - 9792.628*lambda^(16) - 193298.889*im*lambda^(15)*e^(2*im*lambda) + 386598.019*im*lambda^(15)*e^(im*lambda) - 193298.889*im*lambda^(15) + 386604.354*lambda^(14)*e^(2*im*lambda) - 386604.354*lambda^(14) - 7630806.896*im*lambda^(13)*e^(2*im*lambda) + 15261185.681*im*lambda^(13)e^(im*lambda) - 7630806.896*im*lambda^(13) + 15249211.297*lambda^(12)*e^(2*im*lambda) - 15249211.297*lambda^(12) - 301719045.236*im*lambda^(11)*e^(2*im*lambda) + 603799795.567*im*lambda^(11)*e^(im*lambda) - 301719045.236*im*lambda^(11) + 553484502.712*lambda^(10)*e^(2*im*lambda) - 553484502.712*lambda^(10) - 7150474691.827*im*lambda^9*e^(2*im*lambda) + 13888394926.889*im*lambda^9*e^(im*lambda) - 7150474691.827*im*lambda^9 - 94738383623.968*lambda^8*e^(2*im*lambda) + 94738383623.968*lambda^8 + 1428811180235.8*im*lambda^7*e^(2*im*lambda) - 3728406758015.12*im*lambda^7*e^(im*lambda) + 1428811180235.8*im*lambda^7 - 1612610005459.41*lambda^6*e^(2*im*lambda) + 1612610005459.41*lambda^6 + 6700533879857.01*im*lambda^5*e^(2*im*lambda) + 41890893311.889*im*lambda^5*e^(im*lambda) + 6700533879857.01*im*lambda^5 - 6739824792702.45*lambda^4*e^(2*im*lambda) + 6739824792702.45*lambda^4 - 238540562621.508*im*lambda^3*e^(2*im*lambda) - 3443965023.125*im*lambda^3*e^(im*lambda) - 238540562621.508*im*lambda^3 + 234852019362.493*lambda^2*e^(2*im*lambda) - 234852019362.493*lambda^2 - 4516731666.814*im*lambda*e^(2*im*lambda) + 1057517363.677*im*lambda*e^(im*lambda) - 4516731666.814*im*lambda + 3987972984.976*e^(2*im*lambda) - 3987972984.976)/(pi*lambda^22*(im*lambda*e^(2*im*lambda) + im*lambda - 2*e^(2*im*lambda)+2))
    return result
end

FMinusSymF = prettyPrint(FMinusSym(f))

function FMinusSymFunc(lambda)
    result = ((2*pi*im*lambda^(21)*e^(2*im*lambda) - 2*pi*im*lambda^(21) + 12.566*lambda^(20)*e^(2*im*lambda) - 12.566*lambda^(20) + 248.05*im*lambda^(19)*e^(2*im*lambda) - 248.05*im*lambda^(19) + 496.1*lambda^(18)*e^(2*im*lambda) - 496.1*lambda^(18) + 9792.63*im*lambda^17*e^(2*im*lambda) - 9792.63*im*lambda^(17) + 19585.268*lambda^(16)*e^(2*im*lambda) + 0.007*lambda^(16)*e^(im*lambda) - 19585.251*lambda^(16) + 386598.019*im*lambda^(15)*e^(2*im*lambda) - 0.014*im*lambda^(15)*e^(im*lambda) - 386598.052*im*lambda^(15) + 773170.281*lambda^(14)*e^(2*im*lambda) - 26.335*lambda^(14)*e^(im*lambda) - 773221.861*lambda^(14) + 15261185.681*im*lambda^(13)*e^(2*im*lambda) + 52.671*im*lambda^(13)*e^(im*lambda) - 15261082.521*im*lambda^(13) + 30564311.357*lambda^(12)*e^(2*im*lambda) + 49715.326*lambda^(12)*e^(im*lambda) - 30480225.046*lambda^(12) + 603799795.567*im*lambda^(11)*e^(2*im*lambda) - 99430.652*im*lambda^(11)*e^(im*lambda) - 603967968.19*im*lambda^(11) + 1313123613.035*lambda^(10)*e^(2*im*lambda) + 199615489.739*lambda^(10)*e^(im*lambda) - 1102411914.478*lambda^(10) + 13888394926.889*im*lambda^9*e^(2*im*lambda) - 399230979.478*im*lambda^9*e^(im*lambda) - 14309818324.003*im*lambda^9 + 245906746294.396*lambda^8*e^(2*im*lambda) + 435358870071.532*lambda^8*e^(im*lambda) + 189510319792.611*lambda^8 - 3728406758015.12*im*lambda^7*e^(2*im*lambda) - 870717740143.064*im*lambda^7*e^(im*lambda) + 2857572625841.11*im*lambda^7 - 9946584026385.15*lambda^6*e^(2*im*lambda) - 6721484900334.91*lambda^6*e^(im*lambda) + 3225374741327.32*lambda^6 + 41890893311.889*im*lambda^5*e^(2*im*lambda) + 13442969800669.8*im*lambda^5*e^(im*lambda) + 13400527676803.8*im*lambda^5 - 13238593725664.6*lambda^4*e^(2*im*lambda) + 240967733293.361*lambda^4*e^(im*lambda) + 13478679841319.2*lambda^4 - 3443965023.125*im*lambda^3*e^(2*im*lambda) - 481935466586.721*im*lambda^3*e^(im*lambda) - 476728266286.113*im*lambda^3 + 474749529073.638*lambda^2*e^(2*im*lambda) + 5045490348.653*lambda^2*e^(im*lambda) - 471819073452.339*lambda^2 + 1057517363.677*im*lambda*e^(2*im*lambda) - 10090980697.305*im*lambda*e^(im*lambda) - 6918428606.276*im*lambda + 7975945969.952*e^(2*im*lambda) - 7975945969.952) * e^(-im*lambda))/(pi*lambda^(22)*(2*im*lambda*e^(2*im*lambda) + 2*im*lambda - 4*e^(2*im*lambda)+4))
    return result
end

for lambda = 1:20
    tic()
    FPlusSymFunc(lambda)
    toc()
end

# f(x) = x^2-1
# function FPlusSymFunc(lambda)
#     result = 2 *(e^(im*lambda)-1) * (1/8*im*lambda^2*e^(im*lambda) + 1/4*lambda + 1/4*im*e^(im*lambda) - 0.25*im) * e^(-im*lambda)/(pi*lambda^3*(cos(lambda)-1))
#     return result
# end
# function FMinusSymFunc(lambda)
#     result = 2 *(e^(im*lambda)-1) * (1/8*im*lambda^2*e^(im*lambda) + 1/4*lambda + 1/4*im*e^(im*lambda) - 0.25*im) * e^(-2im*lambda)/(pi*lambda^3*(cos(lambda)-1))
#     return result
# end

q = solve_IBVP(L, U, a, zeroList, f; FPlusFunc = FPlusSymFunc, FMinusFunc = FMinusSymFunc)

q(1/2, 1/2)

t = 0.1
# Gadfly.plot(x -> real(q(x,t)), 0, 1)

# TBD

# TBD

# TBD
