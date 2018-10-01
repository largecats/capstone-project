#############################################################################
# Course: YSC4103 MCS Capstone
# Date created: 2018/09/21
# Name: Linfan XIAO
# Description: Algorithm to construct a valid (homogeneous) adjoint boundary condition from a given (homogeneous) boundary condition based on SymPy and julia functions. The symbolic expressions are used to generate the coefficients Bjk.
# Based on: Chapter 11, Theory of Ordinary Differential Equations (Coddington & Levinson)
#############################################################################
# Importing packages
#############################################################################
using SymPy
using Roots
#############################################################################
# Defining types
#############################################################################
# A linear differential operator of order n is encoded by an n x 1 array of functions and an interval [a,b]
struct LinearDifferentialOperator
    pFunctions::Array # Array of julia functions
    interval::Tuple{Number,Number}
    LinearDifferentialOperator(pFunctions::Array, interval::Tuple{Number,Number}) =
    try
        L = new(pFunctions, interval)
        check_linearDifferentialOperator_input(L)
        return L
    catch err
        return err
    end
end

# Check whether the inputs of L are valid
function check_linearDifferentialOperator_input(L::LinearDifferentialOperator)
    pFunctions, (a,b) = L.pFunctions, L.interval
    p0 = pFunctions[1]
    if (isa(p0, Function) && length(find_zeros(p0, a, b)) != 0) || p0 == 0
        error("p0 vanishes on [a,b]") # find_zeros is computationally expensive and makes defining the LinearDifferentialOperator type slow
    else
        return true
    end
end

# A symbolic linear differential operator of order n is encoded by an n x 1 array of symbolic expressions and an interval [a,b]
struct SymLinearDifferentialOperator
    pFunctions::Array{SymPy.Sym,2}
    interval::Tuple{Number,Number} # Must be numbers for the check function to work
    t::SymPy.Sym
    SymLinearDifferentialOperator(pFunctions::Array{SymPy.Sym,2}, interval::Tuple{Number,Number}, t::SymPy.Sym) =
    try
        symL = new(pFunctions, interval, t)
        check_symLinearDifferentialOperator_input(symL)
        return symL
    catch err
        return err  
    end
end

# Check whether the inputs of symL are valid. No checks are needed if symL is only used to create coeffMatrix containing the expressions for Bjk.
function check_symLinearDifferentialOperator_input(symL::SymLinearDifferentialOperator)
    pFunctions, (a,b), t = symL.pFunctions, symL.interval, symL.t
    # p0 = pFunctions[1]
    # zeroP0 = real_roots(p0)
    # if length(find(x -> x == t, free_symbols(p0))) != 0 && length(zeroP0) != 0 && zeroP0[1,1] >= a && zeroP0[1,1] <= b
    #     error("p0 vanishes on [a,b]")
    # else
    #     return true
    # end
    return true
end

# A boundary condition Ux = 0 is encoded by an ordered pair of two matrices (M, N) with symbolic or numeric entries
struct VectorBoundaryForm
    M::Array # Why can't I specify Array{Number,2} without having a MethodError?
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

# Check whether the input matrices that characterize U are valid
function check_vectorBoundaryForm_input(U::VectorBoundaryForm)
    M, N = U.M, U.N
    MHcatN = hcat(M, N)
    if size(M) != size(N)
        error("Matrices' dimensions do not match")
    elseif size(M)[1] != size(M)[2]
        error("Square matrices required")
    elseif rank(MHcatN) != size(M)[1]
        error("Boundary operators not linearly independent")
    else
        return true
    end
end

#############################################################################
# Defining functions
#############################################################################
# Calculate the rank of U, i.e., rank(M:N)
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

# Find a complementary form, Uc, to U
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

# Assign string as variable name
function assign(s::AbstractString,v::Any)
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
function pString_matrix(L::Union{LinearDifferentialOperator, SymLinearDifferentialOperator})
    pFunctions = L.pFunctions
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
function pSymDeriv_matrix(symL::SymLinearDifferentialOperator)
    pFunctions, t = symL.pFunctions, symL.t
    n = length(pFunctions)-1
    pStringMatrix = pString_matrix(symL)
    pSymDerivMatrix = Array{SymPy.Sym}(n,n)

    for i in 1:n
        for j in 1:n
            pString = pStringMatrix[i,j]
            pStringSplit = split(pString, "")
            index, degree = map(parse, pStringSplit[2:3])
            pSymDeriv = pFunctions[index+1]
            pSymDerivMatrix[i,j] = deriv(pSymDeriv, t, degree)
        end
    end
    return pSymDerivMatrix
end

# For L, the above matrix would need to be constructed by hand.
# pDerivMatrix = 

# Create the symbolic expression for [uv](t)
function symUv_form(symL::SymLinearDifferentialOperator, u::SymPy.Sym, v::SymPy.Sym)
    t = symL.t
    n = length(symL.pFunctions)-1
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
    pMatrix = pString_matrix(symL)
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

# Construct the B matrix in the general form using undefined SymFunctions
function get_symB(symL::SymLinearDifferentialOperator, symUvForm::SymPy.Sym, u::SymPy.Sym, v::SymPy.Sym)
    t = symL.t
    n = length(symL.pFunctions)-1
    symB = Array{SymPy.Sym}(n,n)
    for j in 1:n
        for k in 1:n
            term = deriv(u, t, k-1) * deriv(conj(v), t, j-1)
            coefficient = coeff(symUvForm, term)
            symB[j,k] = coefficient
        end
    end
    return symB
end

# (f + g)(x) := f(x) + g(x)
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

# (f * g)(x) := f(x) * g(x)
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

# Find symbolic expression for Bjk using explicit formula
function get_symBjk_by_formula(symL::SymLinearDifferentialOperator, j::Int, k::Int, pStringMatrix::Array{String})
    n = length(symL.pFunctions)-1
    sum = 0
    for l = (j-1):(n-k)
        summand = binomial(l, j-1) * symbols(pStringMatrix[n-k-l+1, l-j+1+1]) * (-1)^l
        sum += summand
    end
    return sum
end

# Find symbolic B using explicit formula
function get_symB_by_formula(symL::SymLinearDifferentialOperator, pStringMatrix::Array{String})
    n = length(symL.pFunctions)-1
    B = Array{Any}(n,n)
    for j = 1:n
        for k = 1:n
            B[j,k] = get_symBjk_by_formula(symL, j, k, pStringMatrix)
        end
    end
    return B
end

# Construct the B matrix from the coefficients of u^{(i)}v^{(j)} in [uv](t)
function get_B(L::LinearDifferentialOperator, pStringMatrix::Array{String}, pDerivMatrix::Array, symB::Array{SymPy.Sym})
    n = length(L.pFunctions)-1
    B = Array{Any}(n,n)
    for i in 1:n
        for j in 1:n
            if i == n + 1 - j
                # B[i,j] = (-1)^(i+1) * pDerivMatrix[1,1]
                B[i,j] = mult_func((-1)^(i+1), pDerivMatrix[1,1])
            elseif i > n + 1 - j # i > j also works?
                B[i,j] = 0
            else
                function bEntry(x)
                    return 0
                end
                coefficient = symB[i,j]
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
                        pDeriv = pDerivMatrix[index]
                        termCoeff = coeff(term, termWoCoeff)
                        # bEntry += termCoeff * pDeriv
                        bEntry = add_func(mult_func(termCoeff, pDeriv), bEntry)
                    end
                end
                B[i,j] = bEntry
            end
        end
    end
    return B
end

# Find Bjk using explicit formula
function get_Bjk_by_formula(L::LinearDifferentialOperator, j::Int, k::Int, pDerivMatrix::Array)
    n = length(L.pFunctions)-1
    sum = 0
    for l = (j-1):(n-k)
        summand = binomial(l, j-1) * pDerivMatrix[n-k-l+1, l-j+1+1] * (-1)^l
        sum += summand
    end
    return sum
end

# Construct the B matrix using explicit formula
function get_B_by_formula(L::LinearDifferentialOperator, pDerivMatrix::Array)
    n = length(L.pFunctions)-1
    B = Array{Any}(n,n)
    for j = 1:n
        for k = 1:n
            B[j,k] = get_Bjk_by_formula(L, j, k, pDerivMatrix)
        end
    end
    return B
end

# Evaluate (the B) matrix at t = a
function evaluate_matrix(matrix::Array, a::Number)
    n = size(matrix)[1]
    matrixA = Array{Number}(n,n)
    for i = 1:n
        for j = 1:n
            if isa(matrix[i,j], Function)
                matrixA[i,j] = matrix[i,j](a)
            else # if isa(matrix[i,j], Number)
                matrixA[i,j] = matrix[i,j]
            end
        end
    end
    return matrixA
end

# Construct B_hat
function B_hat(L::LinearDifferentialOperator, B::Array)
    pFunctions, (a,b) = L.pFunctions, L.interval
    n = length(pFunctions)-1
    bHat = Array{Number}(2n,2n)
    BEvalA = evaluate_matrix(B, a)
    BEvalB = evaluate_matrix(B, b)
    bHat[1:n,1:n] = -BEvalA
    bHat[(n+1):(2n),(n+1):(2n)] = BEvalB
    bHat[1:n, (n+1):(2n)] = 0
    bHat[(n+1):(2n), 1:n] = 0
    return bHat
end

# Construct J = (B_hat * H^{(-1)})^*
function get_J(bHat, H)
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

# \xi needs to be contructed by hand
# xi = [x; x'; x''; ...]

# Evaluate \xi at a and b
function evaluate_xi(L::LinearDifferentialOperator, xi)
    pFunctions, (a,b) = L.pFunctions, L.interval
    n = length(pFunctions)-1
    xiEvalA = Array{Number}(n,1)
    xiEvalB = Array{Number}(n,1)
    for i = 1:n
        xiEvalA[i,1] = xi[i,1](a)
        xiEvalB[i,1] = xi[i,1](b)
    end
    return (xiEvalA, xiEvalB)
end

# Get boundary condition Ux = M\xi(a) + N\xi(b)
function get_boundary_condition(L::LinearDifferentialOperator, U::VectorBoundaryForm, xi)
    M, N = U.M, U.N
    (xiEvalA, xiEvalB) = evaluate_xi(L, xi)
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
    return M * inv(BEvalA) * P == N * inv(BEvalB) * Q
end

#############################################################################
# Processes
#############################################################################
function get_symL(L::LinearDifferentialOperator)
    pFunctions, (a,b) = L.pFunctions, L.interval
    n = length(pFunctions)-1
    t = symbols("t")
    p = SymFunction("p")(t)
    symL = SymLinearDifferentialOperator(repeat([p],outer=(1,n+1)), (a,b), t)
    return symL
end

function get_symB_from_L(L::LinearDifferentialOperator)
    symL = get_symL(L)
    t = symL.t
    u, v = SymFunction("u")(t), SymFunction("v")(t)
    symUvForm = symUv_form(symL, u, v)
    symBFromL = get_symB(symL, symUvForm, u, v)
    return symBFromL
end

function construct_valid_adjoint(L::LinearDifferentialOperator, U::VectorBoundaryForm, pDerivMatrix::Array)
    symL = get_symL(L)
    t = symL.t
    pStringMatrix = pString_matrix(symL)
    symBFromL = get_symB_from_L(L)
    B = get_B(L, pStringMatrix, pDerivMatrix, symBFromL)
    bHat = B_hat(L, B)
    Uc = get_Uc(U)
    H = get_H(U, Uc)
    J = get_J(bHat, H)
    adjointU = get_adjoint(J)
    if check_adjoint(L, U, adjointU, B)
        return adjointU
    else
        error("Adjoint found not valid")
    end
end

#############################################################################
# Tests
#############################################################################
# Constant p_k
y'' + 4y = 0, y(0)=y(2pi)=0
p0, p1, p2 = 1, 0, 4
p00, p10 = p0, p1
p01, p11 = 0, 0
pDerivMatrix = [p00 p01; p10 p11]
L = LinearDifferentialOperator([p0 p1 p2], (0,2pi))
U = VectorBoundaryForm([1 0; 0 0], [0 0; 1 0]) 
# TODO: 1. Need checks that connect L and U? e.g., dimensions should match? 2. Most boundary conditions give Us like this?
construct_valid_adjoint(L, U, pDerivMatrix)

# Variable p_k
# p0*y'' + p1*y' + p2*y = 0, y(0)=y(1)=0
function p0(t) return t + 1 end
function p1(t) return t^2 + 1 end
function p2(t) return 1 end
p00, p10 = p0, p1
function p01(t) return 1 end
function p11(t) return 2t end
pDerivMatrix = [p00 p01; p10 p11]
L = LinearDifferentialOperator([p0 p1 p2], (0,1))
U = VectorBoundaryForm([1 0; 0 0], [0 0; 1 0]) 
construct_valid_adjoint(L, U, pDerivMatrix)

# Check if the Bjk formula works using symbolic epxressions
t = symbols("t")
p = SymFunction("p")(t)
for counter = 2:11
    symL = SymLinearDifferentialOperator(repeat([p],outer=(1,counter)), (0, 1), t)
    u, v = SymFunction("u")(t), SymFunction("v")(t)
    symUvForm = symUv_form(symL, u, v)
    pStringMatrix = pString_matrix(symL)
    symB = get_symB(symL, symUvForm, u, v)
    symBByFormula = get_symB_by_formula(symL, pStringMatrix)

    # Base.showarray(STDOUT, symB, false)
    # Base.showarray(STDOUT, symBByFormula, false)
    println(symB == symBByFormula)
end

# Verify Bjk formula with julia functions
function p(t) return t + 1 end
L = LinearDifferentialOperator([p p p], (0, 1))
pStringMatrix = pString_matrix(L)
symB = get_symB_from_L(L)
p00, p10 = p, p
p01, p11 = 1, 1
pDerivMatrix = [p00 p01; p10 p11]
B = get_B(L, pStringMatrix, pDerivMatrix, symB)
B[1,1](2)
B[1,2](2)
B[2,1](2)
B[2,2]
