##########################################################################################################################################################
# Course: YSC4103 MCS Capstone
# Date created: 2019/02/01
# Name: Linfan XIAO
# Description: Approximate the roots of an exponential polynomial delta by visualizing them as the intersections of the level curves 
# real(delta) = 0 and imaginary(delta) = 0.
##########################################################################################################################################################
# Importing packages and modules
##########################################################################################################################################################
using SymPy
using Plots
##########################################################################################################################################################
# Global variables
##########################################################################################################################################################
lambda = symbols("lambda")

# declare x, y as real variables
x = symbols("x", real = true)
y = symbols("y", real = true)

sympyAddExpr = 1 + x
sympyMultExpr = 2*x
sympyPowerExpr = x^2
sympyExpExpr = e^x

infty = 50
##########################################################################################################################################################
# Helper functions
##########################################################################################################################################################

##########################################################################################################################################################
# Structs
##########################################################################################################################################################

##########################################################################################################################################################
# Main functions
##########################################################################################################################################################
# function to separate real and imaginary parts of an expression delta(lambda), which is an exponential polynomial in one variable

# helper function that deals with the case where func(expr) = func(sympyExpExpr)
# although the function body is the same as "power" and "others", this case is isolated because negative exponents, e.g., factor_list(e^(-im*x)), give PolynomialError('a polynomial expected, got exp(-I*x)',), while factor_list(cos(x)) runs normally
function separate_real_imaginary_exp(expr::SymPy.Sym)
    result = real(expr) + im*imag(expr)
    return result
end
# helper function that deals with the case where func(expr) = func(sympyPowerExpr)
# we won't be dealing with cases like x^(x^x)
function separate_real_imaginary_power(expr::SymPy.Sym)
    result = real(expr) + im*imag(expr)
    return result
end
# helper function that deals with the case where func(expr) = func(sympyMultExpr)
function separate_real_imaginary_mult(expr::SymPy.Sym)
    terms = args(expr)
    result = 1
    # if the expanded expression contains toplevel multiplication, the individual terms must all be exponentials or powers
    for term in terms
        println("term = $term")
        # if term is exponential
        if func(term) == func(sympyExpExpr)
            termSeparated = separate_real_imaginary_exp(term)
        # if term is power (not sure if this case and the case below overlaps)
        elseif func(term) == func(sympyPowerExpr)
            termSeparated = separate_real_imaginary_power(term)
            # else, further split each product term into indivdual factors (this case also includes the case where term is a number, which would go into the "constant" below)
        else
            termSeparated = term # term is a number
#             (constant, factors) = factor_list(term)
#             termSeparated = constant
#             # separate each factor into real and imaginary parts and collect the product of separated factors
#             for (factor, power) in factors
#                 factor = factor^power
#                 termSeparated = termSeparated * (real(factor) + im*imag(factor))
#             end
        end
        println("termSeparated = $termSeparated") 
        # collect the product of separated term, i.e., product of separated factors
        result = result * termSeparated
    end
    return result
end
# helper function that deals with the case where func(expr) = func(sympyAddExpr)
function separate_real_imaginary_add(expr::SymPy.Sym)
    # if the expanded expression contains toplevel addition, the individual terms must all be products or symbols
    terms = args(expr)
    result = 0
    # termSeparated = 0 # to avoid undefined error if there is no else (case incomplete)
    for term in terms
        println("term = $term")
        # if term is a symbol
        if func(term) == func(x)
            termSeparated = term
        # if term is exponential
        elseif func(term) == func(sympyExpExpr)
            termSeparated = separate_real_imaginary_exp(term)
        # if term is a power
        elseif func(term) == func(sympyPowerExpr)
            termSeparated = separate_real_imaginary_power(term)
        # if term is a product
        elseif func(term) == func(sympyMultExpr)
            termSeparated = separate_real_imaginary_mult(term)
        # if term is a number
        else
            termSeparated = term
        end
        println("termSeparated = $termSeparated")
        result = result + termSeparated
    end
    return result
end
# helper function that deals with the case where func(expr) != func(sympyPowerExpr), func(sympyAddExpr), func(sympyMultExpr)
function separate_real_imaginary_others(expr::SymPy.Sym)
    # if the expanded expression is neither of the above, it must be a single term, e.g., x or cos(2x+1), which is a function wrapping around an expression; in this case, use the helper function to clean up the expression and feed it back to the function
    term = args(expr)[1]
    termCleaned = separate_real_imaginary_power_add_mult(term)
    result = subs(expr,args(expr)[1],termCleaned)
    return result
end
# helper function that deals with the cases where func(expr) = func(sympyPowerExpr), func(sympyAddExpr), func(sympyMultExpr)
function separate_real_imaginary_power_add_mult(expr::SymPy.Sym)
    if func(expr) == func(sympyPowerExpr)
        result = separate_real_imaginary_power(expr)
    elseif func(expr) == func(sympyAddExpr)
        result = separate_real_imaginary_add(expr)
    elseif func(expr) == func(sympyMultExpr)
        result = separate_real_imaginary_mult(expr)
    end
    return result
end
# main function
function separate_real_imaginary(delta::SymPy.Sym)
    
    # check if delta has one and only one free symbol (e.g., global variable lambda)
    if length(free_symbols(delta)) == 1
        # substitute lambda as x+iy
        expr = subs(delta, lambda, x+im*y)
        # expand the new expression
        expr = expand(expr)
        
        if func(expr) == func(sympyPowerExpr)
#             println(expr)
#             println("power!")
            result = separate_real_imaginary_power(expr)
#             println("result = $result")
        elseif func(expr) == func(sympyAddExpr)
#             println(expr)
#             println("addition!")
            result = separate_real_imaginary_add(expr)
#             println("result = $result")
        elseif func(expr) == func(sympyMultExpr)
#             println(expr)
#             println("multiplication!")
            result = separate_real_imaginary_mult(expr)
#             println("result = $result")
        else
#             println(expr)
#             println("single term!")
            result = separate_real_imaginary_others(expr)
#             println("result = $result")
        end
        result = expand(result)
        return real(result) + im*imag(result)
        
    else
        throw("Delta has more than one variable")
    end
end

# function to plot the level curves real(Delta(lambda)) = 0 and imaginary(Delta(lambda)) = 0
function plot_level_curves(separatedDelta::SymPy.Sym; realFunc = real(separatedDelta), imagFunc = imag(separatedDelta), xRange = (-infty, infty), yRange = (-infty, infty), step = infty/1000)
    xGridStep = (xRange[2] - xRange[1])/50
    yGridStep = (yRange[2] - yRange[1])/50
    if free_symbols(separatedDelta) == [x, y]
       Plots.contour(xRange[1]:step:xRange[2], yRange[1]:step:yRange[2], realFunc, levels=[0], 
            size = (1500, 1000), tickfontsize = 20, seriescolor=:reds, transpose = false, 
            linewidth = 4, linealpha = 1, 
            xticks = xRange[1]:xGridStep:xRange[2], yticks = yRange[1]:yGridStep:yRange[2], 
            grid = true, gridalpha = 0.5)
        Plots.contour!(xRange[1]:step:xRange[2], yRange[1]:step:yRange[2], imagFunc, levels=[0], 
            size = (1500, 1000), tickfontsize = 20, seriescolor=:blues, transpose = false, 
            linewidth = 4, linealpha = 1,
            xticks = xRange[1]:xGridStep:xRange[2], yticks = yRange[1]:yGridStep:yRange[2], 
            grid = true, gridalpha = 0.5)
    else
        Plots.contour(xRange[1]:step:xRange[2], yRange[1]:step:yRange[2], realFunc, levels=[0], 
            size = (1500, 1000), tickfontsize = 20, seriescolor=:reds, transpose = true, 
            linewidth = 4, linealpha = 1, 
            xticks = xRange[1]:xGridStep:xRange[2], yticks = yRange[1]:yGridStep:yRange[2], 
            grid = true, gridalpha = 0.5)
        Plots.contour!(xRange[1]:step:xRange[2], yRange[1]:step:yRange[2], imagFunc, levels=[0], 
            size = (1500, 1000), tickfontsize = 20, seriescolor=:blues, transpose = true, 
            linewidth = 4, linealpha = 1, 
            xticks = xRange[1]:xGridStep:xRange[2], yticks = yRange[1]:yGridStep:yRange[2], 
            grid = true, gridalpha = 0.5)
    end
end

# function to approximate a zero of Delta given a small range (upon inspecting the level curves)
function approximate_zero(separatedDelta::SymPy.Sym, xRange::Tuple{Number,Number}, yRange::Tuple{Number,Number}; step::Number = (xRange[2]-xRange[1])/200)
    zeroCandidates = [(0.0,0.0)]
    xList = xRange[1]:step:xRange[2]
    yList = yRange[1]:step:yRange[2]
    for xVal in xList
        for yVal in yList
            realVal = subs(real(separatedDelta), (x, xVal), (y, yVal))
            imagVal = subs(imag(separatedDelta), (x, xVal), (y, yVal))
            if isapprox(realVal, 0; atol = step) && isapprox(imagVal, 0; atol = step)
                println((xVal, yVal))
                push!(zeroCandidates, (xVal, yVal))
            end
        end
    end
    deleteat!(zeroCandidates, 1)
    zeroCandidates
end