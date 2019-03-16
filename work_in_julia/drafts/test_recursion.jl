using SymPy

function get_Kn(n; symbolic = true, z = symbols("z"))
    println("n = $n")
    if symbolic
        if n == 0
            expr = 0
        elseif n == 1
            expr = 1/z*(e^z-e^(-z))
        else
            expr = get_Kn(n-2; symbolic = true) + 2(e^z/z + (-1)^n*e^(-z)/z - (n-1)/z*get_Kn(n-1; symbolic = true, z = z))
        end
        # expr = prettyPrint(expr)
        return expr
    else
        function Kn(z)
            if z == 0
                throw("z cannot be 0")
            end
            if n == 0
                result = 0
            elseif n == 1
                result = 1/z*(e^z-e^(-z))
            else
                result = get_Kn(n-2; symbolic = false)(z) + 2*(e^z/z + (-1)^n*e^(-z)/z - (n-1)/z*get_Kn(n-1; symbolic = false)(z))
            end
            return result
        end
        return Kn
    end
end

for n = 1:20
    tic()
    get_Kn(n; symbolic = true)
    toc()
end
# n = 1:20: elapsed time: 23.162084804 seconds
# n = 20: elapsed time: 16.718179318 seconds

for n = 1:5
    tic()
    get_Kn(n; symbolic = false)(1+im)
    toc()
end
# n = 1:20: elapsed time: 43.844209807 seconds

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
            expr = 1/(im*c)*(n*get_Kn(n; symbolic = true, z = (-im*c)) + (-1)^n*e^(im*c) - e^(-im*c))
        end
        # expr = prettyPrint(expr)
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
                result = 1/(im*c)*(n*get_Kn(n; symbolic = false)(-im*c) + (-1)^n*e^(im*c) - e^(-im*c))
            end
            return result
        end
        return TTilde
    end
end

tic()
get_ChebyshevTermIntegral(20; symbolic = true)
toc()
# elapsed time: 28.342955768 seconds

tic()
get_ChebyshevTermIntegral(20; symbolic = false)(1+im)
toc()
# elapsed time: 30.182036329 seconds

using ApproxFun
function get_ChebyshevCoefficients(f::Union{Function,Number})
    fCheb = ApproxFun.Fun(f, 0..1) # Approximate f on [0,1] using chebyshev polynomials
    chebCoefficients = ApproxFun.coefficients(fCheb) # get coefficients of the Chebyshev polynomial
    return chebCoefficients
end

function get_ChebyshevIntegral(l::Number, f::Union{Function, Number}; symbolic = false, lambda = nothing, alpha = nothing)
    fChebCoeffs = get_ChebyshevCoefficients(f)
    if symbolic
        c = alpha^(l-1)*lambda/2
        integralSym = 0
        for m = 1:length(fChebCoeffs)
            integralSym += fChebCoeffs[m] * get_ChebyshevTermIntegral(m-1; symbolic = true, c = c)
        end
        integralSym = integralSym/(2*e^(im*c))
        integralSym = prettyPrint(integralSym)
        integralSym = simplify(integralSym)
        return integralSym
    else
        if isa(lambda, Void) || isa(alpha, Void)
            throw("lambda, alpha required")
        else
            c = alpha^(l-1)*lambda/2
            integral = 0
            for m = 1:length(fChebCoeffs)
                integral += fChebCoeffs[m] * get_ChebyshevTermIntegral(m-1; symbolic = false)(c)
            end
            integral = integral/(2*e^(im*c))
            return integral
        end
    end
end

alpha = e^(2pi*im/1)
l = 2
f(x) = sin(x*2*pi)
lambda = 1+im
tic()
get_ChebyshevIntegral(l, f; symbolic = false, lambda = lambda, alpha = alpha)
toc()
# elapsed time: 87.117664234 seconds

function get_alpha(i, n)
    if i == 1
        alpha = (-1)^n
    elseif i == 2
        alpha = (-1)^(n+1)*n^2
    else
        sum = 0
        for k = 1:(n-i+2)
            product = 1
            for j = k:(i+k-3)
                product *= (n-j)
            end
            sum += binomial(i+k-3, k-1) * product
        end
        alpha = (-1)^(n+i-1)*2^(i-2)*n*sum
    end
    return alpha
end

function get_ChebyshevTermIntegral2(n::Int; symbolic = true, c = symbols("c"))
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
            for i = 1:(n+1)
                summand = get_alpha(i, n) * (e^(im*c)/(im*c)^i  + (-1)^(i+n)*e^(-im*c)/(im*c)^i)
                expr += summand
            end
        end
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

get_ChebyshevTermIntegral(0)
get_ChebyshevTermIntegral2(0)

simplify(get_ChebyshevTermIntegral(1))
get_ChebyshevTermIntegral2(1)

tic()
get_ChebyshevTermIntegral2(20; symbolic = true)
toc()
# elapsed time: 0.102118166 seconds

tic()
get_ChebyshevTermIntegral2(20; symbolic = false)(1+im)
toc()
# elapsed time: 0.093723125 seconds



function get_ChebyshevIntegral2(l::Number, f::Union{Function, Number}; symbolic = false, lambda = nothing, alpha = nothing)
    fChebCoeffs = get_ChebyshevCoefficients(f)
    # fChebCoeffs = filter(x->!isapprox(x,0; atol = 1e-05), fChebCoeffs)
    if symbolic
        c = alpha^(l-1)*lambda/2
        integralSym = 0
        for m = 1:length(fChebCoeffs)
            integralSym += fChebCoeffs[m] * get_ChebyshevTermIntegral2(m-1; symbolic = true, c = c)
        end
        integralSym = integralSym/(2*e^(im*c))
        # integralSym = prettyPrint(integralSym)
        integralSym = simplify(integralSym)
        return integralSym
    else
        if isa(lambda, Void) || isa(alpha, Void)
            throw("lambda, alpha required")
        else
            c = alpha^(l-1)*lambda/2
            integral = 0
            for m = 1:length(fChebCoeffs)
                integral += fChebCoeffs[m] * get_ChebyshevTermIntegral2(m-1; symbolic = false)(c)
            end
            integral = integral/(2*e^(im*c))
            return integral
        end
    end
end

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

alpha = e^(2pi*im/1)
l = 2
lambda = 1 + im
using QuadGK
f(x) = sin(x*pi*2)
g(x) = e^(-im*alpha^(l-1)*lambda*x)
quadgk(mult_func(g,f), 0, 1)[1]

get_ChebyshevIntegral2(l, f; symbolic = false, lambda = lambda, alpha = alpha)