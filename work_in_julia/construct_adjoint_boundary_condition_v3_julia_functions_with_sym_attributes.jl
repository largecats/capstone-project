#############################################################################
# Course: YSC4103 MCS Capstone
# Date created: 2018/09/21
# Name: Linfan XIAO
# Description: Algorithm to construct a valid adjoint boundary condition from a given (homogeneous) boundary condition based on julia functions. The symbolic expressions are implemented as attributes of the julia structs.
# Based on: Chapter 11, Theory of Ordinary Differential Equations (Coddington & Levinson)
#############################################################################
# Importing packages
#############################################################################
using SymPy
using Roots
using Distributions
#############################################################################
# Defining types
#############################################################################
# A struct definition error type is the class of all errors in a struct definition
struct StructDefinitionError <: Exception
    msg::String
end

# A symbolic linear differential operator of order n is encoded by an 1 x (n+1) array of symbolic expressions and an interval [a,b].
struct SymLinearDifferentialOperator
    # Entries in the array should be SymPy.Sym or Number, but somehow SymPy.Sym is ignored, e.g., Array{Union{Number,SymPy.Sym}} returns Array{Number}.
    symPFunctions::Array
    interval::Tuple{Number,Number}
    t::SymPy.Sym
    SymLinearDifferentialOperator(symPFunctions::Array, interval::Tuple{Number,Number}, t::SymPy.Sym) =
    try
        symL = new(symPFunctions, interval, t)
        check_symLinearDifferentialOperator_input(symL)
        return symL
    catch err
        return err  
    end
end

# Check whether the inputs of symL are valid. 
# The symbolic expression has at most one free variable t (may have no free variable when the expression is constant and is represented by a Number).
function check_symLinearDifferentialOperator_input(symL::SymLinearDifferentialOperator)
    symPFunctions, (a,b), t = symL.symPFunctions, symL.interval, symL.t
    for symPFunc in symPFunctions
        if !(isa(symPFunc, SymPy.Sym) || isa(symPFunc, Number))
            throw(StructDefinitionError(:"Symbolic expression is either SymPy.Sym or Number"))
        elseif isa(symPFunc, SymPy.Sym) && size(free_symbols(symPFunc)) != (1,) && size(free_symbols(symPFunc)) != (0,)
            throw(StructDefinitionError(:"Only one free symbol is allowed"))
        end
    end
    return true
end

# A linear differential operator of order n is encoded by an 1 x (n+1) array of functions, an interval [a,b], and its symbolic expression.
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
        return err
    end
end

# Set an appropriate tolerance when checking whether x \approx y
function set_tol(x::Number, y::Number)
    return 1e-09 * mean([x y])
end

# Evaluate p_k on x when p_k is either Function or Number.
# Assume p_k comes from LinearDifferentialOperator.
function evaluate_p(p::Union{Function,Number}, x::Number)
    if isa(p, Function)
        return p(x)
    else
        return p
    end
end

# Assume symFunc has only one free symbol, as required by the definition of SymLinearDifferentialOperator. 
# That is, assume the input symFunc comes from SymLinearDifferentialOperator.
function check_func_sym_equal(func::Union{Function,Number}, symFunc, interval)
    # symFunc should be Union{SymPy.Sym, Number}, but then somehow SymPy.Sym gets ignored
    t = free_symbols(symFunc)[1,1]
    (a,b) = interval
    # Randomly sample 100 points from [a,b] and check if func and symFunc agree on them
    for i = 1:100
        x = rand(Uniform(a,b), 1)[1,1]
        funcEvalX = evaluate_p(func, x)
        # N() converts SymPy.Sym to Number
        # https://docs.sympy.org/latest/modules/evalf.html
        # subs() works no matter symFunc is Number or SymPy.Sym
        symFuncEvalX = N(subs(symFunc,t,x))
        tol = set_tol(funcEvalX, symFuncEvalX)
        if !isapprox(funcEvalX, symFuncEvalX; atol = tol)
            return false
        end
        i += 1
    end
    return true
end

function check_all(list, condition)
    for item in list
        
    end
end

# Check whether the inputs of L are valid.
function check_linearDifferentialOperator_input(L::LinearDifferentialOperator)
    pFunctions, (a,b), symL = L.pFunctions, L.interval, L.symL
    symPFunctions = symL.symPFunctions
    p0 = pFunctions[1]
    if !all(pFunc -> (isa(pFunc, Function) || isa(pFunc, Number)), pFunctions)
        throw(StructDefinitionError(:"Function is either Function or Number"))
    # Assume p_k are in C^{n-k}. Check whether p0 vanishes on [a,b].
    elseif (isa(p0, Function) && length(find_zeros(p0, a, b)) != 0 && (p0(a) == 0 || p0(b) == 0)) || p0 == 0 # find_zeros() finds zeros on (a,b)
        throw(StructDefinitionError(:"p0 vanishes on [a,b]"))
    elseif !all(i -> check_func_sym_equal(pFunctions[i], symPFunctions[i], (a,b)), 1:length(pFunctions))
        throw(StructDefinitionError(:"Symbolic expression does not agree with function"))
    else
        return true
    end
end

#############################################################################
# Tests
#############################################################################
# Test for SymLinearDifferentialOperator definition
function test_symLinearDifferentialOperator_def()
    t = symbols("t")
    symP = t + 1
    a = symbols("a")
    symP1 = a*t + 1
    symL = SymLinearDifferentialOperator([symP symP symP], (0,1) , t)
    try
        symL1 = SymLinearDifferentialOperator([symP symP symP1], (0,1), t)
    catch err
        return (err.msg == "Only one free symbol is allowed") &&
        (isa(err,StructDefinitionError))
    end
    return true
end

test_symLinearDifferentialOperator_def()

# Test for LinearDifferentialOperator definition
function test_linearDifferentialOperator_def()
    function p(t) return t + 1 end
    t = symbols("t")
    symP = t + 1
    symL = SymLinearDifferentialOperator([symP symP symP], (0,1) , t)
    LinearDifferentialOperator([p p p], (0,1), symL)

    symP1 = t + 2
    symL1 = SymLinearDifferentialOperator([symP symP symP1], (0,1) , t)
    try
        LinearDifferentialOperator([p p p], (0,1), symL1)
    catch err
        return (err.msg == "Symbolic expression does not agree with function") &&
        (isa(err,StructDefinitionError))
    end

    function p(t) return t end
    symL0 = SymLinearDifferentialOperator([t t t], (0,1) , t)
    try
        LinearDifferentialOperator([p p p], (0,1), symL0)
    catch err
        return (err.msg == "p0 vanishes on [a,b]") &&
        (isa(err,StructDefinitionError))
    end

    # Constant functions
    symL = SymLinearDifferentialOperator([1 1 1], (0,1) , t)
    LinearDifferentialOperator([1 1 1], (0,1), symL)

    function p(t) return 1 end
    symL = SymLinearDifferentialOperator([1 1 1], (0,1) , t)
    LinearDifferentialOperator([p p p], (0,1), symL)

    return true
end

test_linearDifferentialOperator_def()