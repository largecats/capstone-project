#############################################################################
# Course: YSC4103 MCS Capstone
# Date created: 2018/09/19
# Name: Linfan XIAO
# Description: Algorithm to construct a valid adjoint boundary condition from a given (homogeneous) boundary condition based on SymPy objects (symbolic math).
# Based on: Chapter 11, Theory of Ordinary Differential Equations (Coddington & Levinson)
#############################################################################
# Importing packages
#############################################################################
using SymPy
#############################################################################
# Defining types and functions
#############################################################################
# A linear differential operator of order n is encoded by an n x 1 array of symbolic expressions and an interval [a,b]
struct LinearDifferentialOperator
    pFunctions::Array{}
    interval::Tuple
    t::SymPy.Sym
    LinearDifferentialOperator(pFunctions::Array, interval::Tuple, t::SymPy.Sym) =
    try
        L = new(pFunctions, interval, t)
        check_linearDifferentialOperator_input(L)
        return L
    catch err
        return err  
    end
end

function check_linearDifferentialOperator_input(L::LinearDifferentialOperator)
    pFunctions, (a,b), t = L.pFunctions, L.interval, L.t
    p0 = pFunctions[1]
    zeroP0 = real_roots(p0)
    if any(x -> (typeof(x) != SymPy.Sym), pFunctions)
        error("SymPy.Sym type required")
    elseif length(find(x -> x == t, free_symbols(p0))) != 0 && length(zeroP0) != 0 && zeroP0[1,1] >= a && zeroP0[1,1] <= b
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
    # rank() does not work with SymPy.Sym matrices
    # elseif rank(MHcatN) != size(M)[1]
    #     error("Boundary operators not linearly independent")
    else
        return true
    end
end

# Calculates the rank of U, i.e., rank(M:N)
function rank_of_U(U::VectorBoundaryForm)
    try
        check_vectorBoundaryForm_input(U)
        M, N = U.M, U.N
        MHcatN = hcat(M, N)
        return rank(MHcatN)
    catch err
        return err
    end
end

# Finds a complementary form, Uc, to U
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

# Constructs H from M, N, Mc, Nc
function get_H(U::VectorBoundaryForm, Uc::VectorBoundaryForm)
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
function pString_matrix(L::LinearDifferentialOperator)
    pFunctions, (a,b) = L.pFunctions, L.interval
    n = length(pFunctions)-1

    pStringMatrix = Array{String}(n,n)
    for i in 0:(n-1)
        for j in 0:(n-1)
            pStringMatrix[i+1,j+1] = string("p", i,j)
        end
    end
    return pStringMatrix
end

# Construct a matrix whose ij-entry is the symbolic expression of the jth derivative of p_i
function pFunc_matrix(L::LinearDifferentialOperator, pStringMatrix::Array{String})
    pFunctions, t = L.pFunctions, L.t
    n = size(pStringMatrix)[1]
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
function uv_form(L::LinearDifferentialOperator, u::SymPy.Sym, v::SymPy.Sym)
    pFunctions, (a,b), t = L.pFunctions, L.interval, L.t
    n = length(pFunctions)-1

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
    pMatrix = pString_matrix(L)
    for i = 0:(n-1) # index
        for d = reverse(0:(n-1)) # degree, must start substitution from the highest degree otherwise replacement would not behave as expected
            pTerm = deriv(pFunctionSymbols[i+1], t, d)
            pString = pMatrix[i+1,d+1]
            # Must substitute p_k SymFunction terms with symbols, else cannot collect coefficients of u^{(i)}v^{(j)}
            sum = subs(sum, pTerm, symbols(pString))
        end
    end
    return sum
end

# Construct the B matrix from the coefficients of u^{(i)}v^{(j)} in [uv](t)
function get_B(L::LinearDifferentialOperator, uvForm::SymPy.Sym, u::SymPy.Sym, v::SymPy.Sym, pStringMatrix::Array{String}, pFuncMatrix::Array{SymPy.Sym}, coeffMatrix::Array{SymPy.Sym})
    t = L.t
    n = length(L.pFunctions)-1
    B = Array{SymPy.Sym}(n,n)

    for i in 1:n
        for j in 1:n
            # @printf("i = %g; j = %g\n", i,j)
            if i == n + 1 - j
                B[i,j] = (-1)^(i+1) * pFuncMatrix[1,1]
                # @printf("B[i,j] = %s\n", B[i,j])
            elseif i > j
                B[i,j] = 0
                # @printf("B[i,j] = %s\n", B[i,j])
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
                B[i,j] = bEntry
                # @printf("B[i,j] = %s\n", B[i,j])
            end
        end
    end
    return B
end

# Construct the B matrix in the general form using undefined SymFunctions
function coefficient_matrix(L::LinearDifferentialOperator, uvForm::SymPy.Sym, u::SymPy.Sym, v::SymPy.Sym)
    t = L.t
    n = length(L.pFunctions)-1
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

# Evaluate (the B) matrix at t = a
function evaluate_matrix(matrix::Array{SymPy.Sym}, t::SymPy.Sym, a::Number)
    n = size(B)[1]
    matrixA = Array{SymPy.Sym}(n,n)
    for i = 1:n
        for j = 1:n
            matrixA[i,j] = subs(matrix[i,j], t, a)
        end
    end
    return matrixA
end

# Construct B_hat
function B_hat(L::LinearDifferentialOperator, B::Array{SymPy.Sym})
    pFunctions, (a,b), t = L.pFunctions, L.interval, L.t
    n = length(pFunctions)-1
    bHat = Array{SymPy.Sym}(2n,2n)
    BEvalA = evaluate_matrix(B, t, a)
    BEvalB = evaluate_matrix(B, t, b)
    bHat[1:n,1:n] = -BEvalA
    bHat[(n+1):(2n),(n+1):(2n)] = BEvalB
    bHat[1:n, (n+1):(2n)] = 0
    bHat[(n+1):(2n), 1:n] = 0
    return bHat
end

# Construct J = (B_hat * H^{(-1)})^*
function get_J(bHat::Array{SymPy.Sym}, H)
    n = size(H)[1]
    J = (bHat * inv(H))'
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

# Get \xi
function get_xi(L::LinearDifferentialOperator, x::SymPy.Sym)
    pFunctions, t = L.pFunctions, L.t
    n = length(pFunctions)-1
    xi = Array{SymPy.Sym}(n,1)
    for i = 1:n
        xi[i,1] = deriv(x, t, i-1)
    end
    return xi
end

# Evaluate \xi at a and b
function evaluate_xi(L::LinearDifferentialOperator, xi::Array{SymPy.Sym})
    pFunctions, (a,b), t = L.pFunctions, L.interval, L.t
    n = length(pFunctions)-1
    xiEvalA = Array{SymPy.Sym}(n,1)
    xiEvalB = Array{SymPy.Sym}(n,1)
    for i = 1:n
        xiEvalA[i,1] = subs(xi[i,1], t, a)
        xiEvalB[i,1] = subs(xi[i,1], t, b)
    end
    return (xiEvalA, xiEvalB)
end

# Get boundary condition Ux = M\xi(a) + N\xi(b)
function get_boundary_condition(L::LinearDifferentialOperator, U::VectorBoundaryForm, xi::Array{SymPy.Sym})
    M, N = U.M, U.N
    (xiEvalA, xiEvalB) = evaluate_xi(L, xi)
    Ux = M*xiEvalA + N*xiEvalB
    return Ux
end

# Check if U+ is valid
function check_adjoint(L::LinearDifferentialOperator, U::VectorBoundaryForm, adjointU::VectorBoundaryForm, B::Array{SymPy.Sym})
    (a,b), t = L.interval, L.t
    M, N = U.M, U.N
    P, Q = (adjointU.M)', (adjointU.N)'
    BEvalA = evaluate_matrix(B, t, a)
    BEvalB = evaluate_matrix(B, t, b)
    return M * inv(BEvalA) * P == N * inv(BEvalB) * Q
end

#############################################################################
# Tests
#############################################################################
U = VectorBoundaryForm([1 0; 0 0], [0 0; 1 0])
Uc = get_Uc(U)
H = get_H(U, Uc)

t = symbols("t")
epsilon = symbols("e")
p0, p1, p2 = -epsilon, 1, 0
u, v = SymFunction("u")(t), SymFunction("v")(t)
x, y = t + 1, t + 2

L = LinearDifferentialOperator([p0, p1, p2], (0, 1), t)
pStringMatrix = pString_matrix(L)
pFuncMatrix = pFunc_matrix(L, pStringMatrix)

uvForm = uv_form(L, u, v)
coeffMatrix = coefficient_matrix(L, uvForm, u, v)
B = get_B(L, uvForm, u, v, pStringMatrix, pFuncMatrix, coeffMatrix)
evaluate_matrix(B, t, 0)
bHat = B_hat(L, B)
J = get_J(bHat, H)
adjointU = get_adjoint(J)
P, Q = (adjointU.M)', (adjointU.N)'

check_adjoint(L, U, adjointU, B)

#############################################################################
# Finding patterns
#############################################################################
# p = SymFunction("p")(t)
p = 1 + t
L = LinearDifferentialOperator([p, p, p, p, p, p], (0, 1), t)
pStringMatrix = pString_matrix(L)
pFuncMatrix = pFunc_matrix(L, pStringMatrix)
uvForm = uv_form(L, u, v)
coeffMatrix = coefficient_matrix(L, uvForm, u, v)

# Finding an explicit formula for B_{jk}
function find_Bjk(L::LinearDifferentialOperator, j::Int, k::Int, pStringMatrix::Array{String})
    n = length(L.pFunctions)-1
    sum = 0
    for l = (j-1):(n-k)
        summand = binomial(l, j-1) * symbols(pStringMatrix[n-k-l+1, l-j+1+1]) * (-1)^l
        sum += summand
    end
    return sum
end

for counter = 2:11
    println(counter)
    L = LinearDifferentialOperator(repmat([p],counter), (0, 1), t)
    uvForm = uv_form(L, u, v)
    pStringMatrix = pString_matrix(L)
    coeffMatrix = coefficient_matrix(L, uvForm, u, v)
    n = length(L.pFunctions)-1
    for j = 1:n
        for k = 1:n
            println(find_Bjk(L, j, k, pStringMatrix) == coeffMatrix[j, k])
        end
    end
end
