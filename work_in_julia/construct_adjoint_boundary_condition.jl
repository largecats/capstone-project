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
# A linear differential operator of order n is encoded by an n x 1 array of functions and an interval [a,b]
struct LinearDifferentialOperator
    pFunctions::Array{}
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
    pFunctions, (a, b) = L.pFunctions, L.interval
    p0 = pFunctions[1]
    if p0(a) == 0
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
        push!(output, (i, j))
    end
    return output
end

# Get the symbol (a SymPy.Sym object) for the kth derivative of u with respect to t
function deriv(u::SymPy.Sym, t::Symbol, k::Int)
    if k < 0
        error("Only nonnegative degrees are allowed")
    end
    y = u
    for i = 1:k
        newY = diff(y, t)
        # newY = Derivative(y, t)
        y = newY
    end
    return y
end

function create_p_matrix(L::LinearDifferentialOperator)
    pFunctions, (a, b) = L.pFunctions, L.interval
    n = length(pFunctions)

    mat = Array{String}(n,n)
    for i in 0:(n-1)
        for j in 0:(n-1)
            mat[i+1,j+1] = string("p", i, j)
        end
    end
    return mat
end

function get_uv_form(L::LinearDifferentialOperator)
    pFunctions, (a, b) = L.pFunctions, L.interval
    n = length(pFunctions)

    t = Symbol("t")
    u, v = SymFunction("u")(t), SymFunction("v")(t)
    pFunctionSymbols = [SymFunction(string("p", i))(t) for i in 0:(n-1)]
    sum = 0
    for m = 1:n
        for (j,k) in partition(m-1)
            summand = (-1)^j * deriv(u, t, k) * deriv(pFunctionSymbols[n-m+1] * conj(v), t, j)
            sum += summand
        end
    end
    println(sum)
    sum = expand(sum)
    println(sum)
    pFunctionsStrings = []
    pMatrix = create_p_matrix(L)
    for i = 0:(n-1) # index
        for d = reverse(0:(n-1)) # degree
            pTerm = deriv(pFunctionSymbols[i+1], t, d)
            println(pTerm)
            pString = pMatrix[i+1,d+1]
            println(pString)
            sum = subs(sum, pTerm, Symbol(pString))
        end
    end
    return sum
end

=============================================================================
# Tests
=============================================================================
U = VectorBoundaryForm([1 2+im; 2 1+3im], [2+1im 3; 3 2])
U = VectorBoundaryForm([1 2; 2 4], [1 3; 2 6])
U = VectorBoundaryForm([1 2; 2 1], [2 3; 3 2])
Uc = find_Uc(U)
H = construct_H(U, Uc)

function p(t)
    return t + 1
end
L = LinearDifferentialOperator([p,p], (0, 1))
create_p_matrix(L)
sum = get_uv_form(L)
sum = subs(sum, deriv(pFunctionSymbols[1], t, 0), Symbol("p00"))
subs(sum, deriv(pFunctionSymbols[1], t, 1), Symbol("p00"))
uvForm = get_uv_form(L)
print(args(uvForm))
coeff(uvForm, deriv(u, t, 0)*deriv(conj(v), t, 0))

for i in reverse(1:5)
    println(i)
end