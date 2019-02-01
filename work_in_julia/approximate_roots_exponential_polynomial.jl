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

##########################################################################################################################################################
# Helper functions
##########################################################################################################################################################

##########################################################################################################################################################
# Structs
##########################################################################################################################################################

##########################################################################################################################################################
# Main functions
##########################################################################################################################################################

# function to separate real and imaginary parts of an expression delta(lambda), which is an exponential polynomial
function separate_real_imaginary(delta)
    # check if delta has one free symbol (that is lambda)
    if length(free_symbols(delta)) == 1
        # declare x, y as real variables
        x = symbols("x", real = true)
        y = symbols("y", real = true)
        # substitute lambda as x+iy
        expr = subs(delta, lambda, x+im*y)
        # Split the new expression into summands, if any
        summands = args(expr)
        if length(summands) == 1
            summands = [expr]
        end
        println(summands)
        # For each summand, further split into factors
        sumSeparated = 0
        for summand in summands
            (constant, factors) = factor_list(summand)
            println(factors)
            # Initialize separated summand as the constant coefficient
            summandSeparated = constant
            for (factor, power) in factors
                factor = factor^power
                summandSeparated = summandSeparated * (real(factor) + im*imag(factor))
            end
            sumSeparated = sumSeparated + summandSeparated
        end
        real(sumSeparated) + im*imag(sumSeparated)
    end
end

