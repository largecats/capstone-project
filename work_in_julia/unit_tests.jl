#############################################################################
# Course: YSC4103 MCS Capstone
# Date created: 2018/10/07
# Name: Linfan XIAO
# Description: Unit tests for construct_adjoint_boundary_condition.jl.
#############################################################################
# Importing packages and modules
#############################################################################
include("C:\\Users\\LinFan Xiao\\Academics\\College\\Capstone\\work_in_julia\\construct_adjoint.jl")
using construct_adjoint
#############################################################################
# Helper functions
#############################################################################
function generate_symPFunctions(n; random = false, constant = false)
    global t = symbols("t")
    if random
        symPFunctions = Array{Number}(1,n+1)
        for i = 1:(n+1)
            seed = rand(0:1)
            if seed == 0 # constant
                symPFunction = rand(Uniform(1.0,10.0), 1, 1)[1]
            else # variable
                pFunctionCoeffs = rand(Uniform(1.0,10.0), 1, (rand(1:5)))
                symPFunction = sum([pFunctionCoeffs[i+1]*t^(length(pFunctionCoeffs)-1-i) for i in 0:(length(pFunctionCoeffs)-1)])
            end
            symPFunctions[i] = symPFunction
        end
    else
        if constant # constant
            symPFunctions = rand(Uniform(1.0,10.0), 1, (n+1))
        else # variable
            symPFunctions = Array{Number}(1,n+1)
            for i = 1:(n+1)
                # Each p_k is a polynomial function with random degree between 0 to 4 and random coefficients between 0 and 10
                pFunctionCoeffs = rand(Uniform(1.0,10.0), 1, (rand(1:5)))
                symPFunction = sum([pFunctionCoeffs[i+1]*t^(length(pFunctionCoeffs)-1-i) for i in 0:(length(pFunctionCoeffs)-1)])
                symPFunctions[i] = symPFunction
            end
        end
    end
    return symPFunctions
end

function generate_pFunctions(n; random = false, constant = false)
    if random
        pFunctions = Array{Union{Function, Number}}(1,n+1)
        pDerivMatrix = Array{Union{Function, Number}}(n,n)
        for i = 1:(n+1)
            seed = rand(0:1)
            if seed == 0 # constant
                pFunction = rand(Uniform(1.0,10.0), 1, 1)[1]
                if i < n+1
                    pDerivMatrix[i,1] = pFunction
                    pDerivMatrix[i:i, 2:n] = 0
                end
            else # variable
                pFunctionCoeffs = rand(Uniform(1.0,10.0), 1, (rand(1:5)))
                pFunction = get_polynomial(pFunctionCoeffs)
                if i < n+1
                    pDerivMatrix[i:i, 1:n] = [get_polynomialDeriv(pFunctionCoeffs, k) for k = 0:(n-1)]
                end
            end
            pFunctions[i] = pFunction
        end
    else
        if constant # constant
            pFunctions = rand(Uniform(1.0,10.0), 1, (n+1))
            pDerivMatrix = eye(n)
            for i = 1:n
                for j = 1:n
                    if j == 1
                        pDerivMatrix[i,j] = pFunctions[i]
                    else
                        pDerivMatrix[i,j] = 0
                    end
                end
            end
        else # variable
            pFunctions = Array{Union{Function, Number}}(1,n+1)
            pDerivMatrix = Array{Union{Function, Number}}(n,n)
            for i = 1:(n+1)
                # Each p_k is a polynomial function with random degree between 0 to 4 and random coefficients between 0 and 10
                pFunctionCoeffs = rand(Uniform(1.0,10.0), 1, (rand(1:5)))
                if i < n+1
                    pDerivMatrix[i:i, 1:n] = [get_polynomialDeriv(pFunctionCoeffs, k) for k = 0:(n-1)]
                end
                pFunction = get_polynomial(pFunctionCoeffs)
                pFunctions[i] = pFunction
            end
        end
    end
    return pFunctions, pDerivMatrix
end

function generate_pFunctionsAndSymPFunctions(n; random = false, constant = false)
    global t = symbols("t")
    if random
        pFunctions = Array{Union{Function, Number}}(1,n+1)
        symPFunctions = Array{Number}(1,n+1)
        pDerivMatrix = Array{Union{Function, Number}}(n,n)
        for i = 1:(n+1)
            seed = rand(0:1)
            if seed == 0 # constant
                pFunction = rand(Uniform(1.0,10.0), 1, 1)[1]
                symPFunction = pFunction
                if i < n+1
                    pDerivMatrix[i,1] = pFunction
                    pDerivMatrix[i:i, 2:n] = 0
                end
            else # variable
                pFunctionCoeffs = rand(Uniform(1.0,10.0), 1, (rand(1:5)))
                if i < n+1
                    pDerivMatrix[i:i, 1:n] = [get_polynomialDeriv(pFunctionCoeffs, k) for k = 0:(n-1)]
                end
                pFunction = get_polynomial(pFunctionCoeffs)
                symPFunction = sum([pFunctionCoeffs[i+1]*t^(length(pFunctionCoeffs)-1-i) for i in 0:(length(pFunctionCoeffs)-1)])
            end
            pFunctions[i] = pFunction
            symPFunctions[i] = symPFunction
        end
    else
        if constant # constant
            pFunctions = rand(Uniform(1.0,10.0), 1, (n+1))
            symPFunctions = pFunctions
            pDerivMatrix = Array{Union{Function, Number}}(n,n)
            for i = 1:n
                for j = 1:n
                    if j == 1
                        pDerivMatrix[i,j] = pFunctions[i]
                    else
                        pDerivMatrix[i,j] = 0
                    end
                end
            end
        else # variable
            t = symbols("t")
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
        end
    end
    
    return pFunctions, symPFunctions, pDerivMatrix
end
#############################################################################
# Tests
#############################################################################
# Test the algorithm to generate valid adjoint U+
function test_generate_adjoint(n, k)
    global results = [true]
    global t = symbols("t")
    global (a,b) = (0,1)

    for counter = 1:k
        println("Test $counter")
        println("Testing the algorithm to generate valid adjoint U+: Constant p_k")
        (pFunctions, symPFunctions, pDerivMatrix) = generate_pFunctionsAndSymPFunctions(n; random = false, constant = true)
        symL = SymLinearDifferentialOperator(symPFunctions, (a,b), t)
        L = LinearDifferentialOperator(pFunctions, (a,b), symL)
        MCand = rand(Uniform(1.0,10.0), n, n)
        NCand = rand(Uniform(1.0,10.0), n, n)
        U = VectorBoundaryForm(MCand, NCand)
        println("Testing: order of L = $n")
        passed = false
        try
            adjoint = construct_validAdjoint(L, U, pDerivMatrix)
            passed = true
            append!(results, passed)
        catch err
            println("Failed with $err")
        end
        if passed
            println("Passed!")
        end

        println("Testing the algorithm to generate valid adjoint U+: Variable p_k")
        # Generate variable p_k
        (pFunctions, symPFunctions, pDerivMatrix) = generate_pFunctionsAndSymPFunctions(n; random = false, constant = false)
        symL = SymLinearDifferentialOperator(symPFunctions, (a,b), t)
        L = LinearDifferentialOperator(pFunctions, (a,b), symL)
        MCand = rand(Uniform(1.0,10.0), n, n)
        NCand = rand(Uniform(1.0,10.0), n, n)
        U = VectorBoundaryForm(MCand, NCand)
        println("Testing: order of L = $n")
        try
            adjoint = construct_validAdjoint(L, U, pDerivMatrix)
            passed = true
            append!(results, passed)
        catch err
            println("Failed with $err")
        end
        if passed
            println("Passed!")
        end

        println("Testing the algorithm to generate valid adjoint U+: Constant or variable p_k")
        # Generate p_k
        (pFunctions, symPFunctions, pDerivMatrix) = generate_pFunctionsAndSymPFunctions(n; random = true)
        symL = SymLinearDifferentialOperator(symPFunctions, (a,b), t)
        L = LinearDifferentialOperator(pFunctions, (a,b), symL)
        MCand = rand(Uniform(1.0,10.0), n, n)
        NCand = rand(Uniform(1.0,10.0), n, n)
        U = VectorBoundaryForm(MCand, NCand)
        println("Testing: order of L = $n")
        try
            adjoint = construct_validAdjoint(L, U, pDerivMatrix)
            passed = true
            append!(results, passed)
        catch err
            println("Failed with $err")
        end
        if passed
            println("Passed!")
        end
    end

    return all(results)
end
# for n = 1:10
#     println(test_generate_adjoint(n, 10))
# end

# Test the SymLinearDifferentialOperator definition
function test_symLinearDifferentialOperatorDef(n, k)
    global results = [true]
    global t = symbols("t")
    global (a,b) = (0,1)

    for counter = 1:k
        println("Test $counter")
        println("Testing definition of SymLinearDifferentialOperator: symP_k are SymPy.Sym")
        symPFunctions = generate_symPFunctions(n; random = false, constant = false)
        passed = false
        try
            SymLinearDifferentialOperator(symPFunctions, (a,b), t)
            passed = true
        catch err
            println("Failed with $err")
        end
        if passed
            println("Passed!")
        end
        append!(results, passed)

        println("Testing definition of SymLinearDifferentialOperator: symP_k are Number")
        symPFunctions = generate_symPFunctions(n; random = false, constant = true)
        passed = false
        try
            SymLinearDifferentialOperator(symPFunctions, (a,b), t)
            passed = true
        catch err
            println("Failed with $err")
        end
        if passed
            println("Passed!")
        end
        append!(results, passed)

        println("Testing definition of SymLinearDifferentialOperator: symP_k are SymPy.Sym and Number")
        symPFunctions = generate_symPFunctions(n; random = true)
        passed = false
        try
            SymLinearDifferentialOperator(symPFunctions, (a,b), t)
            passed = true
        catch err
            println("Failed with $err")
        end
        if passed
            println("Passed!")
        end
        append!(results, passed)

        println("Testing StructDefinitionError: symP_k should be SymPy.Sym or Number")
        symPFunctions = hcat(generate_symPFunctions(n-1; random = true), ["str"])
        passed = false
        try
            SymLinearDifferentialOperator(symPFunctions, (a,b), t)
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
        r = symbols("r")
        passed = false
        try
            SymLinearDifferentialOperator([t+1 t+1 r*t+1], (a,b), t)
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
    end

    return all(results)
end
# for n = 1:10
#     println(test_symLinearDifferentialOperatorDef(n, 10))
# end

# Test the LinearDifferentialOperator definition
function test_linearDifferentialOperatorDef(n, k)
    global results = [true]
    global t = symbols("t")
    global (a,b) = (0,1)
    
    for counter = 1:k
        println("Test $k")

        # Variable p_k
        println("Testing definition of LinearDifferentialOperator: p_k are variable Function")
        (pFunctions, symPFunctions, pDerivMatrix) = generate_pFunctionsAndSymPFunctions(n; random = false, constant = false)
        symL = SymLinearDifferentialOperator(symPFunctions, (a,b), t)
        passed = false
        try
            L = LinearDifferentialOperator(pFunctions, (a,b), symL)
            passed = true
        catch err
            println("Failed with $err")
        end
        if passed
            println("Passed!")
        end
        append!(results, passed)

        # Constant coefficients
        println("Testing definition of LinearDifferentialOperator: p_k are Constants")
        (pFunctions, symPFunctions, pDerivMatrix) = generate_pFunctionsAndSymPFunctions(n; random = false, constant = true)
        symL = SymLinearDifferentialOperator(symPFunctions, (a,b), t)
        passed = false
        try
            LinearDifferentialOperator(pFunctions, (a,b), symL)
            passed = true
        catch err
            println("Failed with $err")
        end
        if passed
            println("Passed!")
        end
        append!(results, passed)

        # # Not necessary, since constant Functions are still Functions
        # println("Testing definition of LinearDifferentialOperator: p_k are constant Function")
        # symL = SymLinearDifferentialOperator([1 1 1], (0,1), t)
        # passed = false
        # try
        #     LinearDifferentialOperator([t->1 t->1 t->1], (0,1), symL)
        #     passed = true
        # catch err
        #     println("Failed with $err")
        # end
        # if passed
        #     println("Passed!")
        # end
        # append!(results, passed)

        # Mixed coefficients
        println("Testing definition of LinearDifferentialOperator: p_k are mixed")
        # (pFunctions, symPFunctions, pDerivMatrix) = generate_pFunctionsAndSymPFunctions(n; random = true)
        symL = SymLinearDifferentialOperator([1 1 t+1], (a,b), t)
        passed = false
        try
            LinearDifferentialOperator([1 t->1 t->t+1], (a,b), symL)
            passed = true
        catch err
            println("Failed with $err")
        end
        if passed
            println("Passed!")
        end
        append!(results, passed)

        println("Testing StructDefinitionError: p_k should be Function or Number")
        pFunctions = hcat(generate_pFunctions(n-1; random = true)[1], ["str"])
        passed = false
        try
            LinearDifferentialOperator(['s' 1 1], (a,b), symL)
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
        symL = SymLinearDifferentialOperator([1 1 t+1], (a,b), t)
        passed = false
        try
            LinearDifferentialOperator([1 t->1], (a,b), symL)
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

        println("Testing StructDefinitionError: Intervals of L and symL do not match")
        symL = SymLinearDifferentialOperator([1 1 t+1], (a,b), t)
        passed = false
        try
            LinearDifferentialOperator([1 t->1 t->t+1], (-b,a), symL)
        catch err
            if err.msg == "Intervals of L and symL do not match" && (isa(err, StructDefinitionError))
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
        symL = SymLinearDifferentialOperator([t 1 2], (a,b), t)
        passed = false
        try
            LinearDifferentialOperator([t->t 1 2], (a,b), symL)
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
    end

    return all(results)
end
# for n = 1:10
#     println(test_linearDifferentialOperatorDef(n, 10))
# end

# Test the VectorBoundaryForm definition
function test_vectorBoundaryFormDef(n, k)
    global results = [true]

    for counter = 1:k
        println("Test $counter")

        println("Testing the definition of VectorBoundaryForm")
        MReal = rand(Uniform(1.0,10.0), n, n)
        MComplex = rand(Uniform(1.0,10.0), n, n)
        M = MReal + im*MComplex
        NReal = rand(Uniform(1.0,10.0), n, n)
        NComplex = rand(Uniform(1.0,10.0), n, n)
        N = NReal + im*NComplex
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
        M = rand(Uniform(1.0,10.0), n, n-1)
        N = rand(Uniform(1.0,10.0), n, n)
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
        M = rand(Uniform(1.0,10.0), n, n-1)
        N = rand(Uniform(1.0,10.0), n, n-1)
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
        M = [1 2*im; 2 4*im]
        N = [3 4*im; 6 8*im]
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
    end

    return all(results)
end
# for n = 1:10
#     println(test_vectorBoundaryFormDef(n, 10))
# end