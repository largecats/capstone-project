##########################################################################################################################################################
# Course: YSC4103 MCS Capstone
# Date created: 2018/10/07
# Name: Linfan XIAO
# Description: Implement the transform pair (2.15a), (2.15b) on page 10 of "Evolution PDEs and augmented eigenfunctions. Finite interval."
##########################################################################################################################################################
# Importing packages and modules
##########################################################################################################################################################
using SymPy
# using QuadGK
# using HCubature
using ApproxFun

include("C:\\Users\\LinFan Xiao\\Academics\\College\\Capstone\\work_in_julia\\construct_adjoint.jl")
##########################################################################################################################################################
# Helper functions
##########################################################################################################################################################
# Function addition (f + g)(x) := f(x) + g(x)
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

# Function multiplication (f * g)(x) := f(x) * g(x)
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
##########################################################################################################################################################
# Structs
##########################################################################################################################################################

##########################################################################################################################################################
# Functions
##########################################################################################################################################################
# Get M+, M- in (2.13a), (2.13b)
function get_MPlusMinus(adjointU::VectorBoundaryForm, lambda::Number)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    alpha = e^(2pi*im/n)
    MPlus = Array{Number}(n,n)
    MMinus = Array{Number}(n,n)
    for k = 1:n
        for j = 1:n
            sumPlus, sumMinus = 0, 0
            for r = 0:(n-1)
                summandPlus = (-im*alpha^(k-1)*lambda)^r * bStar[j,r+1]
                summandMinus = (-im*alpha^(k-1)*lambda)^r * betaStar[j,r+1]
                sumPlus += summandPlus
                sumMinus += summandMinus
            end
            MPlus[k,j] = sumPlus
            MMinus[k,j] = sumMinus
        end
    end
    return (MPlus, MMinus)
end

# Get M in (2.14)
function get_M(adjointU::VectorBoundaryForm, lambda::Number)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    alpha = e^(2pi*im/n)
    (MPlus, MMinus) = get_MPlusMinus(adjointU, lambda)
    M = Array{Number}(n,n)
    for k = 1:n
        for j = 1:n
            M[k,j] = MPlus[k,j] + MMinus[k,j] * e^(-im*alpha^(k-1)*lambda) # M+ is just M+(lambda)
        end
    end
    return M
end

# Get Xlj, which is the (n-1)*(n-1) submatrix of M with (1,1) entry the (l + 1, j + 1) entry of M
function get_Xlj(adjointU::VectorBoundaryForm, M::Array, lambda::Number, l::Number, j::Number)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    M = get_M(adjointU, lambda)
    Xlj = M[(l+1):(l+1+n-2), (j+1):(j+1+n-2)]
    return Xlj
end

# Get F+, F- in (2.16a), (2.16b)
function get_FPlusMinus(adjointU::VectorBoundaryForm, f::Union{Function, Number}, lambda::Number)
    # add constraints on lambda \in contour
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    alpha = e^(2pi*im/n)
    (MPlus, MMinus) = get_MPlusMinus(adjointU, lambda)
    M = get_M(adjointU, lambda)
    M = convert(Array{Complex}, M)
    delta = det(M)
    sumPlus, sumMinus = 0, 0
    for l = 1:n
        for j = 1:n
            Xlj = get_Xlj(adjointU, M, lambda, l, j)
            Xlj = convert(Array{Complex}, Xlj)
            g(x) = e^(-im*alpha^(l-1)*lambda*x)
            integrand = mult_func(g, f)
            summandPlus = (-1)^((n-1)*(l+j)) * det(Xlj) * MPlus[1,j] * quadgk(integrand, 0, 1)[1]
            summandMinus = (-1)^((n-1)*(l+j)) * det(Xlj) * MMinus[1,j] * quadgk(integrand, 0, 1)[1]
        end
        sumPlus += summandPlus
        sumMinus += summandMinus
    end
    FPlus = 1/(2pi*delta)*sumPlus
    FMinus = (-e^(-im*lambda))/(2pi*delta)*sumMinus
    return (FPlus, FMinus)
end

# Get F in (2.15a)
function get_F(adjointU::VectorBoundaryForm, f::Union{Function, Number}, lambda::Number)
    (FPlus, FMinus) = get_FPlusMinus(adjointU, f, lambda)
    # if lambda in GammaPlusDomain return FPlus else return FMinus
end

# Get f in (2.15b)
function get_f(Gamma)
    f(x) = e^(im*lambda*x) # * integral(Gamma, F(lambda))
    return f
end

# Get the nth term in the Chebyshev series
# The three values in xRange correspond to cases where x in [-1,1], x<-1, and x>1
# https://en.wikipedia.org/wiki/Chebyshev_polynomials
function get_ChebyshevSeriesTerm(n::Int64; xRange = "=", symbolic = true)
    if symbolic
        x = symbols("x")
        if xRange == "="
            return cos(n*acos(x))
        elseif xRange == "<"
            return cosh(n*acosh(x))
        elseif xRange == ">"
            return (-1)^n * cosh(n*acosh(-x))
        end
    else
        if xRange == "="
            return x -> cos(n*acos(x))
        elseif xRrange == "<"
            return x-> cosh(n*acosh(x))
        elseif xRange == ">"
            return x -> (-1)^n * cosh(n*acosh(-x))
        end
    end
end

# Get the Chebyshev approximation of a function in range [a,b] as a Julia Function object or its symbolic expression
function get_ChebyshevApproximation(f::Function, a::Number, b::Number; symbolic = true)
    fCheb = ApproxFun.Fun(f, a..b) # Approximate f on [a,b] using chebyshev polynomials
    chebCoefficients = ApproxFun.coefficients(fCheb) # get coefficients of the Chebyshev polynomial
    n = length(coefficients)
    if symbolic
        equal, smaller, greater = 0, 0, 0
        for i = 1:n
            chebCoefficient = chebCoefficients[i]
            equal += chebCoefficient * get_ChebyshevSeriesTerm(i; xRange = "=", symbolic = symbolic)
            smaller += chebCoefficient * get_ChebyshevSeriesTerm(i; xRange = "<", symbolic = symbolic)
            greater += chebCoefficient * get_ChebyshevSeriesTerm(i; xRange = ">", symbolic = symbolic)
        end
        # An array containing the symbolic expression of fChebApprox on [-1,1], (-\infty, -1), and (1, \infty)
        fChebApproxSymArray = [smaller, equal, greater]
        return fChebApproxSymArray
    else
        function fChebApprox(x)
            sum = 0
            for i = 1:n
                chebCoefficient = chebCoefficients[i]
                if x <= 1 || x >= -1
                    summand = mult_func(chebCoefficient, get_ChebyshevSeriesTerm(i; xRange = "=", symbolic = symbolic))
                elseif x < -1
                    summand = mult_func(chebCoefficient, get_ChebyshevSeriesTerm(i; xRange = "<", symbolic = symbolic))
                else # x > 1
                    summand = mult_func(chebCoefficient, get_ChebyshevSeriesTerm(i; xRange = ">", symbolic = symbolic))
                end
                sum = add_func(sum, summand)
            end
            println(typeof(sum(x)))
            return sum(x)
        end
        return fChebApprox
    end
end