=============================================================================
# Course: YSC4103 MCS Capstone
# Date created: 2018/09/19
# Name: Linfan XIAO
# Description: Algorithm to construct a valid adjoint boundary condition from a given boundary condition.
# Based on: Chapter 11, Theory of Ordinary Differential Equations (Coddington & Levinson)
=============================================================================
# Defining types and functions
=============================================================================
# A linear differential operator of order n is encoded by an n x 1 array of functions and an interval [a,b]
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
L = LinearDifferentialOperator([p], (-1, 1))

function u(x)
    x
end
Derivative(u, 1)
Derivative(u, -1)