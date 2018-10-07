#############################################################################
# Course: YSC4103 MCS Capstone
# Date created: 2018/10/07
# Name: Linfan XIAO
# Description: Unit tests for construct_adjoint_boundary_condition.jl.
#############################################################################
# Importing packages
#############################################################################
using SymPy
using Roots
using Distributions
#############################################################################
# Tests
#############################################################################
# Test the algorithm to generate valid adjoint U+
function test_generate_adjoint(n, k)
    results = [true]
    t = symbols("t")
    (a,b) = (0,1)

    println("Testing the algorithm to generate valid adjoint U+: Constant p_k")
    pFunctions = rand(Uniform(1.0,10.0), 1, (n+1))
    symPFunctions = pFunctions
    symL = SymLinearDifferentialOperator(symPFunctions, (a,b), t)
    L = LinearDifferentialOperator(pFunctions, (a,b), symL)
    MCand = rand(Uniform(1.0,10.0), n, n)
    NCand = rand(Uniform(1.0,10.0), n, n)
    U = VectorBoundaryForm(MCand, NCand)
    pDerivMatrix = eye(n)
    println("Testing: order of L = $n")
    for counter = 1:k
        for i = 1:n
            for j = 1:n
                if j == 1
                    pDerivMatrix[i,j] = pFunctions[i]
                else
                    pDerivMatrix[i,j] = 0
                end
            end
        end
        passed = false
        try
            adjoint = construct_validAdjoint(L, U, pDerivMatrix)
            passed = true
            append!(results, passed)
        catch err
            println("Failed with $err")
        end
        if passed
            println("Test $counter: Passed!")
        end
    end

    println("Testing the algorithm to generate valid adjoint U+: Variable p_k")
    # Generate variable p_k
    pFunctions = Array{Function}(1,n+1)
    symPFunctions = Array{Number}(1,n+1)
    pDerivMatrix = Array{Union{Function, Number}}(n,n)
    for i = 1:(n+1)
        # Each p_k is a polynomial function with random degree between 0 to 4 and random coefficients between 0 and 10
        pFunctionCoeffs = rand(Uniform(1.0,10.0), 1, (rand(1:5)))
        if i < n+1
            pDerivMatrix[i:i, 1:n] = [get_polynomialDeriv(pFunctionCoeffs, k) for k = 0:(n-1)]
        end
        pFunction = get_polynomial(pFunctionCoeffs)
        pFunctions[i] = pFunction
        symPFunction = sum([pFunctionCoeffs[i+1]*t^(length(pFunctionCoeffs)-1-i) for i in 0:(length(pFunctionCoeffs)-1)])
        symPFunctions[i] = symPFunction
    end
    symL = SymLinearDifferentialOperator(symPFunctions, (a,b), t)
    L = LinearDifferentialOperator(pFunctions, (a,b), symL)
    MCand = rand(Uniform(1.0,10.0), n, n)
    NCand = rand(Uniform(1.0,10.0), n, n)
    U = VectorBoundaryForm(MCand, NCand)
    println("Testing: order of L = $n")
    for counter = 1:k
        for i = 1:n
            for j = 1:n
                if j == 1
                    pDerivMatrix[i,j] = pFunctions[i]
                else
                    pDerivMatrix[i,j] = 0
                end
            end
        end
        passed = false
        try
            adjoint = construct_validAdjoint(L, U, pDerivMatrix)
            passed = true
            append!(results, passed)
        catch err
            println("Failed with $err")
        end
        if passed
            println("Test $counter: Passed!")
        end
    end

    return all(results)
end
for n = 1:10
    println(test_generate_adjoint(n, 10))
end

# Test the SymLinearDifferentialOperator definition
function test_symLinearDifferentialOperatorDef()
    results = [true]
    println("Testing definition of SymLinearDifferentialOperator: symP_k are SymPy.Sym")
    t = symbols("t")
    passed = false
    try
        SymLinearDifferentialOperator([t+1 t+1 t+1], (0,1), t)
        passed = true
    catch err
        println("Failed with $err")
    end
    if passed
        println("Passed!")
    end
    append!(results, passed)

    println("Testing definition of SymLinearDifferentialOperator: symP_k are Number")
    passed = false
    try
        SymLinearDifferentialOperator([1 1 1], (0,1), t)
        passed = true
    catch err
        println("Failed with $err")
    end
    if passed
        println("Passed!")
    end
    append!(results, passed)

    println("Testing definition of SymLinearDifferentialOperator: symP_k are SymPy.Sym and Number")
    passed = false
    try
        SymLinearDifferentialOperator([1 1 t+1], (0,1), t)
        passed = true
    catch err
        println("Failed with $err")
    end
    if passed
        println("Passed!")
    end
    append!(results, passed)

    println("Testing StructDefinitionError: symP_k should be SymPy.Sym or Number")
    passed = false
    try
        SymLinearDifferentialOperator(['s' 1 t+1], (0,1), t)
    catch err
        if isa(err,StructDefinitionError) && err.msg == "symP_k should be SymPy.Sym or Number"
            passed = true
            println("Passed!")
        else
            println("Failed with $err")
        end
    end
    if !passed
        println("Failed!")
    end
    append!(results, passed)

    println("Testing StructDefinitionError: Only one free symbol is allowed in symP_k")
    a = symbols("a")
    passed = false
    try
        SymLinearDifferentialOperator([t+1 t+1 a*t+1], (0,1), t)
    catch err
        if isa(err,StructDefinitionError) && err.msg == "Only one free symbol is allowed in symP_k"
            passed = true
            println("Passed!")
        else
            println("Failed with $err")
        end
    end
    if !passed
        println("Failed!")
    end
    append!(results, passed)

    return all(results)
end
test_symLinearDifferentialOperatorDef()

# Test the LinearDifferentialOperator definition
function test_linearDifferentialOperatorDef()
    results = [true]
    # Variable p_k
    println("Testing definition of LinearDifferentialOperator: p_k are variable Function")
    t = symbols("t")
    symL = SymLinearDifferentialOperator([t+1 t+1 t+1], (1,2), t)
    passed = false
    try
        L = LinearDifferentialOperator([t->t+1 t->t+1 t->t+1], (1,2), symL)
        passed = true
    catch err
        println("Failed with $err")
    end
    if passed
        println("Passed!")
    end
    append!(results, passed)

    # Constant coefficients
    println("Testing definition of LinearDifferentialOperator: p_k are Number")
    symL = SymLinearDifferentialOperator([1 1 1], (0,1), t)
    passed = false
    try
        LinearDifferentialOperator([1 1 1], (0,1), symL)
        passed = true
    catch err
        println("Failed with $err")
    end
    if passed
        println("Passed!")
    end
    append!(results, passed)

    println("Testing definition of LinearDifferentialOperator: p_k are constant Function")
    symL = SymLinearDifferentialOperator([1 1 1], (0,1), t)
    passed = false
    try
        LinearDifferentialOperator([t->1 t->1 t->1], (0,1), symL)
        passed = true
    catch err
        println("Failed with $err")
    end
    if passed
        println("Passed!")
    end
    append!(results, passed)

    # Mixed coefficients
    println("Testing definition of LinearDifferentialOperator: p_k are mixed")
    symL = SymLinearDifferentialOperator([1 1 t+1], (0,1), t)
    passed = false
    try
        LinearDifferentialOperator([1 t->1 t->t+1], (0,1), symL)
        passed = true
    catch err
        println("Failed with $err")
    end
    if passed
        println("Passed!")
    end
    append!(results, passed)

    println("Testing StructDefinitionError: p_k should be Function or Number")
    passed = false
    try
        LinearDifferentialOperator(['s' 1 1], (0,1), symL)
    catch err
        if err.msg == "p_k should be Function or Number" && (isa(err,StructDefinitionError))
            passed = true
            println("Passed!")
        else
            println("Failed with $err")
        end
    end
    if !passed
        println("Failed!")
    end
    append!(results, passed)

    println("Testing StructDefinitionError: Number of p_k and symP_k do not match")
    symL = SymLinearDifferentialOperator([1 1 t+1], (0,1), t)
    passed = false
    try
        LinearDifferentialOperator([1 t->1], (0,1), symL)
    catch err
        if err.msg == "Number of p_k and symP_k do not match" && (isa(err, StructDefinitionError))
            passed = true
            println("Passed!")
        else
            println("Failed with $err")
        end
    end
    if !passed
        println("Failed!")
    end
    append!(results, passed)

    println("Testing StructDefinitionError: p0 vanishes on [a,b]")
    function p2(t) return t end
    symL = SymLinearDifferentialOperator([t 1 2], (0,1), t)
    passed = false
    try
        LinearDifferentialOperator([t->t 1 2], (0,1), symL)
    catch err
        if err.msg == "p0 vanishes on [a,b]" && (isa(err, StructDefinitionError))
            passed = true
            println("Passed!")
        else
            println("Failed with $err")
        end
    end
    if !passed
        println("Failed!")
    end
    append!(results, passed)
    
    # # This is now a warning
    # println("Testing StructDefinitionError: symP_k does not agree with p_k on [a,b]")
    # symL = SymLinearDifferentialOperator([t+1 t+1 t+2], (0,1), t)
    # passed = false
    # try
    #     LinearDifferentialOperator([t->t+1 t->t+1 t->t+1], (0,1), symL)
    # catch err
    #     if err.msg == "symP_k does not agree with p_k on [a,b]" && (isa(err,StructDefinitionError))
    #         passed = true
    #         println("Passed!")
    #     else
    #         println("Failed with $err")
    #     end
    # end
    # if !passed
    #     println("Failed!")
    # end
    # append!(results, passed)

    return all(results)
end
test_linearDifferentialOperatorDef()

# Test the VectorBoundaryForm definition
function test_vectorBoundaryFormDef()
    results = [true]
    println("Testing the definition of VectorBoundaryForm")
    M = eye(3)
    N = M
    passed = false
    try
        VectorBoundaryForm(M, N)
        passed = true
    catch err
        println("Failed with $err")
    end
    if passed
        println("Passed!")
    end
    append!(results, passed)

    println("Testing StructDefinitionError: Entries of M, N should be Number")
    M = ['a' 2; 3 4]
    N = M
    passed = false
    try
        VectorBoundaryForm(M, N)
    catch err
        if err.msg == "Entries of M, N should be Number" && isa(err, StructDefinitionError)
            passed = true
            println("Passed!")
        else
            println("Failed with $err")
        end
    end
    if !passed
        println("Failed!")
    end
    append!(results, passed)

    println("Testing StructDefinitionError: M, N dimensions do not match")
    M = eye(2)
    N = eye(3)
    passed = false
    try
        VectorBoundaryForm(M, N)
    catch err
        if err.msg == "M, N dimensions do not match" && isa(err,StructDefinitionError)
            passed = true
            println("Passed!")
        else
            println("Failed with $err")
        end
    end
    if !passed
        println("Failed!")
    end
    append!(results, passed)

    println("Testing StructDefinitionError: M, N should be square matrices")
    M = [1 2]
    N = M
    passed = false
    try
        VectorBoundaryForm(M, N)
    catch err
        if err.msg == "M, N should be square matrices" && isa(err,StructDefinitionError)
            passed = true
            println("Passed!")
        else
            println("Failed with $err")
        end
    end
    if !passed
        println("Failed!")
    end
    append!(results, passed)

    println("Testing StructDefinitionError: Boundary operators not linearly independent")
    M = [1 2; 2 4]
    N = [3 4; 6 8]
    passed = false
    try
        VectorBoundaryForm(M, N)
    catch err
        if err.msg == "Boundary operators not linearly independent" && isa(err,StructDefinitionError)
            passed = true
            println("Passed!")
        end
    end
    if !passed
        println("Failed!")
    end
    append!(results, passed)

    return all(results)
end
test_vectorBoundaryFormDef()
