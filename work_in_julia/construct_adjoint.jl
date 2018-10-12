#############################################################################
# Course: YSC4103 MCS Capstone
# Date created: 2018/10/07
# Name: Linfan XIAO
# Description: Algorithm to construct a valid adjoint boundary condition from a given (homogeneous) boundary condition based on Chapter 11 in Theory of Ordinary Differential Equations (Coddington & Levinson). The implementation uses Julia functions as main objects but supports symbolic expressions in the form of Julia struct attributes. If a function can produce both symbolic and non-symbolic outputs, the output type is controlled by get_output(...; symbolic = true/false). If a function produces only symbolic outputs, it is usually called get_symOutput().
#############################################################################
# Importing packages
#############################################################################
using SymPy
# using Roots
using Distributions
using IntervalArithmetic
using IntervalRootFinding
#############################################################################
# Helper functions
#############################################################################
# Check whether all elements in a not necessarily homogeneous array satisfy a given condition.
function check_all(array, condition)
    for x in array
        if !condition(x)
            return false
        end
    end
    return true
end

# Set an appropriate tolerance when checking whether x \approx y
function set_tol(x::Number, y::Number)
    return 1e-05 * mean([x y])
end

# Set an appropriate tolerance when checking whether M \approx N
function set_tol_matrix(A::Array{Number}, B::Array{Number})
    if size(A) != size(B)
        throw(error("Matrix dimensions do not match"))
    end
    return 1e-05 * (norm(A,2) + norm(B,2))
end

# Evaluate function on x where the function is Function, SymPy.Sym, or Number.
function evaluate(func::Union{Function,Number}, x::Number, t=nothing)
    if isa(func, Function)
        return func(x)
    elseif isa(func, SymPy.Sym) # SymPy.Sym must come before Number because SymPy.Sym will be recognized as Number
        return subs(func, t, x)
    else
        return func
    end
end

# Assign string as variable name
function assign(s::AbstractString, v::Any)
    s=Symbol(s)
    @eval (($s) = ($v))
end

# Generate two-integer partitions of n
function partition(n::Int)
    output = []
    for i = 0:n
        j = n - i
        push!(output, (i,j))
    end
    return output
end

# Construct the symbolic expression for the kth derivative of u with respect to t
function symDeriv(u::SymPy.Sym, t::SymPy.Sym, k::Int)
    if k < 0
        throw(error("Only nonnegative degrees are allowed"))
    end
    y = u
    for i = 1:k
        newY = diff(y, t)
        y = newY
    end
    return y
end

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

# Evaluate a matrix at t=a.
# Entries of B may be Function, Number, or Sympy.Sym.
function evaluate_matrix(matrix::Array, a::Number, t=nothing)
    (m, n) = size(matrix)
    matrixA = Array{Number}(m,n)
    for i = 1:m
        for j = 1:n
            matrixA[i,j] = evaluate(matrix[i,j], a, t)
        end
    end
    return matrixA
end

# Implement a polynomial in the form of Julia function given an array containing coefficients of x^n, x^{n-1},..., x^2, x, 1.
function get_polynomial(coeffList::Array)
    polynomial = 0
    n = length(coeffList)-1
    for i in 0:n
        newTerm = t -> coeffList[i+1] * t^(n-i)
        polynomial = add_func(polynomial, newTerm)
    end
    return polynomial
end

# Get the kth derivative of a polynomial implemented above
function get_polynomialDeriv(coeffList::Array, k::Int)
    if k < 0
        throw(error("Only nonnegative degrees are allowed"))
    elseif k == 0
        newCoeffList = coeffList
    else
        for counter = 1:k
            n = length(coeffList)
            newCoeffList = hcat([0],[(n-i)*coeffList[i] for i in 1:(n-1)]')
            coeffList = newCoeffList
        end
    end
    return get_polynomial(newCoeffList)
end
#############################################################################
# Structs
#############################################################################
# A struct definition error type is the class of all errors in a struct definition
struct StructDefinitionError <: Exception
    msg::String
end

# A symbolic linear differential operator of order n is encoded by an 1 x (n+1) array of symbolic expressions and an interval [a,b].
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

# A linear differential operator of order n is encoded by an 1 x (n+1) array of functions, an interval [a,b], and its symbolic expression.
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
        funcEvalX = evaluate(func, x)
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
    domainC = Complex(a..b, 0..0) # Domain [a,b] represented in the complex plane
    p0 = pFunctions[1]
    if !check_all(pFunctions, pFunc -> (isa(pFunc, Function) || isa(pFunc, Number)))
        throw(StructDefinitionError(:"p_k should be Function or Number"))
    elseif length(pFunctions) != length(symPFunctions)
        throw(StructDefinitionError(:"Number of p_k and symP_k do not match"))
    elseif (a,b) != symL.interval
        throw(StructDefinitionError(:"Intervals of L and symL do not match"))
    # Assume p_k are in C^{n-k}. Check whether p0 vanishes on [a,b]. roots() doesn't work if p0 is sth like t+2*im
    elseif (isa(p0, Function) && (!isempty(roots(p0, domainC, Newton)) || p0(a) == 0 || p0(b) == 0)) || p0 == 0 
        throw(StructDefinitionError(:"p0 vanishes on [a,b]"))
    elseif !all(i -> check_func_sym_equal(pFunctions[i], symPFunctions[i], (a,b), t), 1:length(pFunctions))
        # throw(StructDefinitionError(:"symP_k does not agree with p_k on [a,b]"))
        warn("symP_k does not agree with p_k on [a,b]") # Make this a warning instead of an error because the functionalities of Julia Functions may be more than those of SymPy objects; we do not want to compromise the functionalities of LinearDifferentialOperator because of the restrictions on SymPy.
    else
        return true
    end
end

# A boundary condition Ux = 0 is encoded by an ordered pair of two matrices (M, N) whose entries are Numbers.
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
    M, N = U.M, U.N
    if !(check_all(M, x -> isa(x, Number)) && check_all(N, x -> isa(x, Number)))
        throw(StructDefinitionError(:"Entries of M, N should be Number"))
    elseif size(M) != size(N)
        throw(StructDefinitionError(:"M, N dimensions do not match"))
    elseif size(M)[1] != size(M)[2]
        throw(StructDefinitionError(:"M, N should be square matrices"))
    elseif rank(hcat(M, N)) != size(M)[1]
        throw(StructDefinitionError(:"Boundary operators not linearly independent"))
    else
        return true
    end
end

#############################################################################
# Functions
#############################################################################
# Calculate the rank of U, i.e., rank(M:N)
function rank_of_U(U::VectorBoundaryForm)
    M, N = U.M, U.N
    MHcatN = hcat(M, N)
    return rank(MHcatN)
end

# Find Uc, a complementary form of U
function get_Uc(U::VectorBoundaryForm)
    try
        check_vectorBoundaryForm_input(U)
        n = rank_of_U(U)
        I = eye(2*n)
        M, N = U.M, U.N
        MHcatN = hcat(M, N)
        mat = MHcatN
        for i = 1:(2*n)
            newMat = vcat(mat, I[i:i,:])
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

# Construct H from M, N, Mc, Nc
function get_H(U::VectorBoundaryForm, Uc::VectorBoundaryForm)
    MHcatN = hcat(U.M, U.N)
    McHcatNc = hcat(Uc.M, Uc.N)
    H = vcat(MHcatN, McHcatNc)
    return H
end

# Construct a matrix whose ij-entry is a string "pij" which denotes the jth derivative of p_i
function get_pStringMatrix(L::LinearDifferentialOperator)
    if isa(L, LinearDifferentialOperator)
        pFunctions = L.pFunctions
    else
        pFunctions = L.symPFunctions
    end
    n = length(pFunctions)-1
    pStringMatrix = Array{String}(n,n)
    for i in 0:(n-1)
        for j in 0:(n-1)
            pStringMatrix[i+1,j+1] = string("p", i,j)
        end
    end
    return pStringMatrix
end

# Construct a matrix whose ij-entry is the symbolic expression of the jth derivative of p_i.
# Functions with keyword arguments are defined using a semicolon in the signature.
function get_symPDerivMatrix(L::LinearDifferentialOperator; substitute = true)
    symL = L.symL
    symPFunctions, t = symL.symPFunctions, symL.t
    n = length(symPFunctions)-1
    symPDerivMatrix = Array{SymPy.Sym}(n,n)
    if substitute
        pFunctionSymbols = symPFunctions
    else
        pFunctionSymbols = [SymFunction(string("p", i))(t) for i in 0:(n-1)]
    end
    for i in 1:n
        for j in 1:n
            index, degree = i-1, j-1
            symPDeriv = pFunctionSymbols[index+1]
            symPDerivMatrix[i,j] = symDeriv(symPDeriv, t, degree)
        end
    end
    return symPDerivMatrix
end

# For L, the above matrix would need to be constructed by hand.
# pDerivMatrix = 

# Create the symbolic expression for [uv](t).
# If substitute is true: Substitute the p_k SymFunctions with SymPy.Sym definitions, e.g., substitute p0 by t + 1.
function get_symUvForm(L::LinearDifferentialOperator, u::SymPy.Sym, v::SymPy.Sym; substitute = true)
    symL = L.symL
    symPFunctions, t = symL.symPFunctions, symL.t
    n = length(symPFunctions)-1
    if substitute
        pFunctionSymbols = symPFunctions
    else
        pFunctionSymbols = [SymFunction(string("p", i))(t) for i in 0:(n-1)]
    end
    sum = 0
    for m = 1:n
        for (j,k) in partition(m-1)
            summand = (-1)^j * symDeriv(u, t, k) * symDeriv(pFunctionSymbols[n-m+1] * conj(v), t, j)
            sum += summand
        end
    end
    sum = expand(sum)
    return sum
end

# # Find symbolic expression for Bjk using explicit formula.
# # If substitute is true: Substitute the p_k SymFunctions with SymPy.Sym definitions, e.g., substitute p0 by t + 1.
# function get_symBjk(L::LinearDifferentialOperator, j::Int, k::Int; substitute = true)
#     n = length(L.pFunctions)-1
#     sum = 0
#     matrix = get_symPDerivMatrix(L; substitute = substitute)
#     for l = (j-1):(n-k)
#         summand = binomial(l, j-1) * matrix[n-k-l+1, l-j+1+1] * (-1)^l
#         sum += summand
#     end
#     return sum
# end

# # Find symbolic B using explicit formula.
# # If substitute is true: Substitute the p_k SymFunctions with SymPy.Sym definitions, e.g., substitute p0 by t + 1.
# function get_symB(L::LinearDifferentialOperator; substitute = true)
#     n = length(L.pFunctions)-1
#     B = Array{Any}(n,n)
#     for j = 1:n
#         for k = 1:n
#             B[j,k] = get_symBjk(L, j, k; substitute = substitute)
#         end
#     end
#     return B
# end

# Find Bjk using explicit formula
function get_Bjk(L::LinearDifferentialOperator, j::Int, k::Int; symbolic = false, substitute = true, pDerivMatrix::Array = [])
    n = length(L.pFunctions)-1
    sum = 0
    if symbolic
        symPDerivMatrix = get_symPDerivMatrix(L; substitute = substitute)
        for l = (j-1):(n-k)
            summand = binomial(l, j-1) * symPDerivMatrix[n-k-l+1, l-j+1+1] * (-1)^l
            sum += summand
        end
    else
        if isempty(pDerivMatrix)
            throw(error("pDerivMatrix required"))
        elseif size(pDerivMatrix) != (n,n)
            throw(error("Size of pDerivMatrix should be ($n,$n)"))
        end
        for l = (j-1):(n-k)
            summand = mult_func(binomial(l, j-1) * (-1)^l, pDerivMatrix[n-k-l+1, l-j+1+1])
            sum = add_func(sum, summand)
        end
    end
    return sum
end

# Construct the B matrix using explicit formula
function get_B(L::LinearDifferentialOperator; symbolic = false, substitute = true, pDerivMatrix::Array = [])
    n = length(L.pFunctions)-1
    B = Array{Union{Function, Number}}(n,n)
    for j = 1:n
        for k = 1:n
            B[j,k] = get_Bjk(L, j, k; symbolic = symbolic, substitute = substitute, pDerivMatrix = pDerivMatrix)
        end
    end
    return B
end

# Construct B_hat. Since all entries of B_hat are evaluated, BHat is a numeric matrix.
function get_BHat(L::LinearDifferentialOperator, B::Array)
    pFunctions, (a,b) = L.pFunctions, L.interval
    n = length(pFunctions)-1
    BHat = Array{Number}(2n,2n)
    BEvalA = evaluate_matrix(B, a)
    BEvalB = evaluate_matrix(B, b)
    BHat[1:n,1:n] = -BEvalA
    BHat[(n+1):(2n),(n+1):(2n)] = BEvalB
    BHat[1:n, (n+1):(2n)] = 0
    BHat[(n+1):(2n), 1:n] = 0
    return BHat
end

# Construct J = (B_hat * H^{(-1)})^*, where ^* denotes conjugate transpose
function get_J(BHat, H)
    n = size(H)[1]
    J = (BHat * inv(H))'
    return J
end

# Construct U+
function get_adjoint(J)
    n = convert(Int, size(J)[1]/2)
    Pstar = J[(n+1):2n,1:n]
    Qstar = J[(n+1):2n, (n+1):2n]
    adjoint = VectorBoundaryForm(Pstar, Qstar)
    return adjoint
end

# Construct the symbolic expression of \xi = [x; x'; x''; ...], an n x 1 vector of derivatives of x(t)
function get_symXi(L::LinearDifferentialOperator; substitute = false, xDef = nothing)
    n = length(L.pFunctions)-1
    t = symbols("t")
    symXi = Array{SymPy.Sym}(n,1)
    if !substitute
        xDef = SymFunction("x")(t)
    end
    for i = 1:n
        try
            symXi[i] = symDeriv(xDef,t,i-1)
        catch err
            if isa(err, MethodError)
                throw(error("Definition of x required"))
            end
        end
    end
    return symXi
end

# For L, \xi needs to be contructed by hand, e.g., xi = [t->t^2+1; t->2t; t->2; ...]

# Evaluate \xi at a.
# \xi could be Function, SymPy.Sym, or Number.
function evaluate_xi(L::LinearDifferentialOperator, a::Number, xi)
    if any(xDeriv->isa(xDeriv, SymPy.Sym), xi)
        t = L.symL.t
    else
        t = nothing
    end
    n = length(xi)
    xiEvalA = Array{Number}(n,1)
    for i = 1:n
        xiEvalA[i,1] = evaluate(xi[i,1],a,t)
    end
    return xiEvalA
end

# Get boundary condition Ux = M\xi(a) + N\xi(b)
function get_boundaryCondition(L::LinearDifferentialOperator, U::VectorBoundaryForm, xi)
    # xi cannot contain other types than Function or Number
    if !all(xDeriv->isa(xDeriv, Union{Function, Number}), xi)
        typeXi = typeof(xi)
        throw(error("xi of type $typeXi does not match L of type LinearDifferentialOperator"))
    end
    (a,b) = L.interval
    M, N = U.M, U.N
    xiEvalA = evaluate_xi(L, a, xi)
    xiEvalB = evaluate_xi(L, b, xi)
    Ux = M*xiEvalA + N*xiEvalB
    return Ux
end

# Check if U+ is valid (only works for homogeneous cases Ux=0)
function check_adjoint(L::LinearDifferentialOperator, U::VectorBoundaryForm, adjointU::VectorBoundaryForm, B::Array)
    (a,b) = L.interval
    M, N = U.M, U.N
    P, Q = (adjointU.M)', (adjointU.N)'
    BEvalA = evaluate_matrix(B, a)
    BEvalB = evaluate_matrix(B, b)
    left = M * inv(BEvalA) * P
    right = N * inv(BEvalB) * Q
    # println("M * inv(BEvalA) * P = $left")
    # println("N * inv(BEvalB) * Q = $right")
    tol = set_tol_matrix(left, right)
    return all(i -> isapprox(left[i], right[i]; atol = tol), length(left)) # Can't use == to deterimine equality because left and right are arrays of floats
end

# Find a valid adjoint
function construct_validAdjoint(L::LinearDifferentialOperator, U::VectorBoundaryForm, pDerivMatrix::Array)
    B = get_B(L; pDerivMatrix = pDerivMatrix)
    BHat = get_BHat(L, B)
    Uc = get_Uc(U)
    H = get_H(U, Uc)
    J = get_J(BHat, H)
    adjointU = get_adjoint(J)
    if check_adjoint(L, U, adjointU, B)
        return adjointU
    else
        throw(error("Adjoint found not valid"))
    end
end
