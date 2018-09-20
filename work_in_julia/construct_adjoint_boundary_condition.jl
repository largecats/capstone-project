=============================================================================
# Course: YSC4103 MCS Capstone
# Date created: 2018/09/19
# Name: Linfan XIAO
# Description: Algorithm to construct a valid adjoint boundary condition from a given boundary condition.
# Based on: Chapter 11, Theory of Ordinary Differential Equations (Coddington & Levinson)
=============================================================================
# Importing packages
=============================================================================
using SymPy
=============================================================================
# Defining types and functions
=============================================================================
# A linear differential operator of order n is encoded by an n x 1 array of symbolic expressions and an interval [a,b]
struct LinearDifferentialOperator
    pFunctions::Array
    interval::Tuple

    LinearDifferentialOperator(pFunctions::Array, interval::Tuple) =
    try
        L = new(pFunctions, interval)
        check_linearDifferentialOperator_input(L)
        return L
    catch err
        return err  
    end
end

function check_linearDifferentialOperator_input(L::LinearDifferentialOperator)
    pFunctions, (a,b) = L.pFunctions, L.interval
    if (any(x -> (typeof(x) != SymPy.Sym), pFunctions))
        error("SymPy.Sym type required")
    elseif subs(pFunctions[1], free_symbols(pFunctions[1])[1,1], a) == 0
        error("p0 vanishes on [a,b]")
    else
        return true
    end
end

# A boundary condition Ux = 0 is encoded by an ordered pair of two matrices (M, N)
struct VectorBoundaryForm
    M::Array
    N::Array

    VectorBoundaryForm(M::Array, N::Array) =
    try
        U = new(M, N)
        check_vectorBoundaryForm_input(U)
        return U
    catch err
        err
    end
end

# Checks whether the input matrices that characterize U are valid
function check_vectorBoundaryForm_input(U::VectorBoundaryForm)
    M, N = U.M, U.N
    MHcatN = hcat(M, N)
    if size(M) != size(N)
        error("Matrices' dimensions do not match")
    elseif size(M)[1] != size(M)[2]
        error("Square matrices required")
    elseif rank(MHcatN) != ndims(M)
        error("Boundary operators not linearly independent")
    else
        return true
    end
end

# Calculates the rank of U, i.e., rank(M:N)
function get_rank(U::VectorBoundaryForm)
    try
        check_vectorboundaryform_input(U)
        M, N = U.M, U.N
        MHcatN = hcat(M, N)
        return rank(MHcatN)
    catch err
        return err
    end
end

# Finds a complementary form, Uc, to U
function find_Uc(U::VectorBoundaryForm)
    try
        check_vectorboundaryform_input(U)
        n = get_rank(U)
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

# Constructs H from M, N, Mc, Nc
function construct_H(U::VectorBoundaryForm, Uc::VectorBoundaryForm)
    MHcatN = hcat(U.M, U.N)
    McHcatNc = hcat(Uc.M, Uc.N)
    H = vcat(MHcatN, McHcatNc)
    return H
end

# Using expressions
# Assigns string as variable name
function assign(s::AbstractString,v::Any)
    s=Symbol(s)
    @eval (($s) = ($v))
end

# Two-integer partitions of n
function partition(n::Int)
    output = []
    for i = 0:n
        j = n - i
        push!(output, (i,j))
    end
    return output
end

# Construct the symbol (a SymPy.Sym object) for the kth derivative of u with respect to t
function deriv(u::SymPy.Sym, t::SymPy.Sym, k::Int)
    if k < 0
        error("Only nonnegative degrees are allowed")
    end
    y = u
    for i = 1:k
        newY = diff(y, t)
        y = newY
    end
    return y
end

# Construct a matrix whose ij-entry is a string "pij" which denotes the jth derivative of p_i
function get_pString_matrix(L::LinearDifferentialOperator)
    pFunctions, (a,b) = L.pFunctions, L.interval
    n = length(pFunctions)

    pStringMatrix = Array{String}(n,n)
    for i in 0:(n-1)
        for j in 0:(n-1)
            pStringMatrix[i+1,j+1] = string("p", i,j)
        end
    end
    return pStringMatrix
end

# Construct a matrix whose ij-entry is the symbolic expression of the jth derivative of p_i
function get_pFunc_matrix(L::LinearDifferentialOperator, pStringMatrix::Array{String}, t::SymPy.Sym)
    pFunctions = L.pFunctions
    n = ndims(pStringMatrix)
    pFuncMatrix = Array{SymPy.Sym}(n,n)

    for i in 1:n
        for j in 1:n
            pString = pStringMatrix[i,j]
            pStringSplit = split(pString, "")
            index, degree = map(parse, pStringSplit[2:3])
            pFuncSymbol = pFunctions[index+1]
            # pFuncMatrix[i,j] = lambdify(deriv(pFuncSymbol, t, degree))
            pFuncMatrix[i,j] = deriv(pFuncSymbol, t, degree)
        end
    end
    return pFuncMatrix
end

# Create the symbolic expression for [uv](t)
function get_uv_form(L::LinearDifferentialOperator, t::SymPy.Sym, u::SymPy.Sym, v::SymPy.Sym)
    pFunctions, (a,b) = L.pFunctions, L.interval
    n = length(pFunctions)

    pFunctionSymbols = [SymFunction(string("p", i))(t) for i in 0:(n-1)]
    sum = 0
    for m = 1:n
        for (j,k) in partition(m-1)
            summand = (-1)^j * deriv(u, t, k) * deriv(pFunctionSymbols[n-m+1] * conj(v), t, j)
            sum += summand
        end
    end
    sum = expand(sum)
    pFunctionsStrings = []
    pMatrix = get_pString_matrix(L)
    for i = 0:(n-1) # index
        for d = reverse(0:(n-1)) # degree, must start substitution from the highest degree otherwise replacement would not behave as expected
            pTerm = deriv(pFunctionSymbols[i+1], t, d)
            pString = pMatrix[i+1,d+1]
            sum = subs(sum, pTerm, symbols(pString))
        end
    end
    return sum
end

# Construct the B matrix in general form
function get_coefficient_matrix(L::LinearDifferentialOperator, uvForm::SymPy.Sym, t::SymPy.Sym, u::SymPy.Sym, v::SymPy.Sym)
    n = length(L.pFunctions)
    coeffMatrix = Array{SymPy.Sym}(n,n)
    for j in 1:n
        for k in 1:n
            term = deriv(u, t, k-1) * deriv(conj(v), t, j-1)
            coefficient = coeff(uvForm, term)
            coeffMatrix[j,k] = coefficient
        end
    end
    return coeffMatrix
end

# Construct the matrix B from the coefficients of u^{(i)}v^{(j)} in [uv](t)
function get_B_matrix(L::LinearDifferentialOperator, uvForm::SymPy.Sym, t::SymPy.Sym, u::SymPy.Sym, v::SymPy.Sym, pStringMatrix::Array{String}, pFuncMatrix::Array{SymPy.Sym}, coeffMatrix::Array{SymPy.Sym})
    n = length(L.pFunctions)
    bMatrix = Array{SymPy.Sym}(n,n)

    for i in 1:n
        for j in 1:n
            # @printf("i = %g; j = %g\n", i,j)
            if i == n + 1 - j
                bMatrix[i,j] = (-1)^(i+1) * pFuncMatrix[1,1]
                # @printf("bMatrix[i,j] = %s\n", bMatrix[i,j])
            elseif i > j
                bMatrix[i,j] = 0
                # @printf("bMatrix[i,j] = %s\n", bMatrix[i,j])
            else
                bEntry = 0
                coefficient = coeffMatrix[i,j]
                coeffArgs = args(coefficient)
                if isempty(coeffArgs)
                    coeffList = [coefficient]
                else
                    coeffList = coeffArgs
                end
                for l in 1:length(coeffList)
                    term = coeffList[l]
                    indexMat = find(x -> (symbols(x) == term || symbols(x) == -term), pStringMatrix)
                    if isempty(indexMat)
                        continue
                    else
                        index = indexMat[1,1]
                        termWoCoeff = symbols(pStringMatrix[index])
                        pFunc = pFuncMatrix[index]
                        termCoeff = coeff(term, termWoCoeff)
                        bEntry += termCoeff * pFunc
                    end
                end
                bMatrix[i,j] = bEntry
                # @printf("bMatrix[i,j] = %s\n", bMatrix[i,j])
            end
        end
    end
    return bMatrix
end

# Evaluate matrix at t = a
function evaluate_matrix(matrix::Array{SymPy.Sym}, t::SymPy.Sym, a::Number)
    n = ndims(bMatrix)
    matrixA = Array{Any}(n,n)
    for i = 1:n
        for j = 1:n
            matrixA[i,j] = subs(matrix[i,j], t, a)
        end
    end
    return matrixA
end

# Construct B_hat
function get_B_hat(L::LinearDifferentialOperator, bMatrix::Array{SymPy.Sym}, t::SymPy.Sym)
    pFunctions, (a,b) = L.pFunctions, L.interval
    n = length(pFunctions)
    bHat = zeros(2n,2n)

    bMatrixA = evaluate_matrix(bMatrix, t, a)
    bMatrixB = evaluate_matrix(bMatrix, t, b)
    bHat[1:n,1:n] = -bMatrixA
    bHat[(n+1):(2n),(n+1):(2n)] = bMatrixB

    return bHat
end


=============================================================================
# Tests
=============================================================================
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
bHatEval = get_B_hat(L, bMatrix, t)

