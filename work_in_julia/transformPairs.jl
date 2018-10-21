##########################################################################################################################################################
# Course: YSC4103 MCS Capstone
# Date created: 2018/10/07
# Name: Linfan XIAO
# Description: Implement the transform pair (2.15a), (2.15b) on page 10 of "Evolution PDEs and augmented eigenfunctions. Finite interval."
##########################################################################################################################################################
# Importing packages
##########################################################################################################################################################
using SymPy
using QuadGK
using HCubature
##########################################################################################################################################################
# Helper functions
##########################################################################################################################################################

##########################################################################################################################################################
# Structs
##########################################################################################################################################################

##########################################################################################################################################################
# Functions
##########################################################################################################################################################
function get_MPlusMinus(adjointU::VectorBoundaryForm, lambda::Number)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    alpha = e^(2pi*im/n) # shouldn't this be 1.0?
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

function get_Xlj(adjointU::VectorBoundaryForm, M::Array, lambda::Number, l::Number, j::Number)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    M = get_M(adjointU, lambda)
    Xlj = M[(l+1):(l+1+n-2), (j+1):(j+1+n-2)]
    return Xlj
end

# A lot of redundant calculations because of the argument lambda...
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
            x(t) = t
            # complex integrand not working with current packages
            # summandPlus = (-1)^((n-1)*(l+j)) * det(Xlj) * MPlus[1,j] * quadgk(e^(-im*alpha^(l-1)*lambda*x)*f, 0, 1, rtol=1e-5)
            # summandMinus = (-1)^((n-1)*(l+j)) * det(Xlj) * MMinus[1,j] * quadgk(e^(-im*alpha^(l-1)*lambda*x)*f, 0, 1, rtol=1e-5)
            # summandPlus = (-1)^((n-1)*(l+j)) * det(Xlj) * MPlus[1,j] * hquadrature(e^(-im*alpha^(l-1)*lambda*x)*f, 0, 1)
        end
        sumPlus += summandPlus
        sumMinus += summandMinus
    end
    FPlus = 1/(2pi*delta)*sumPlus
    FMinus = (-e^(-im*lambda))/(2pi*delta)*sumMinus
    return (FPlus, FMinus)
end

function get_F(adjointU::VectorBoundaryForm, f::Union{Function, Number}, lambda::Number)
    (FPlus, FMinus) = get_FPlusMinus(adjointU, f, lambda)
    # if lambda in GammaPlusDomain return FPlus else return FMinus
end

function get_f(Gamma)
    f(x) = e^(im*lambda*x) # * integral(Gamma, F(lambda))
    return f
end