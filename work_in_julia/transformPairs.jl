##########################################################################################################################################################
# Course: YSC4103 MCS Capstone
# Date created: 2018/10/07
# Name: Linfan XIAO
# Description: Implement the transform pair (2.15a), (2.15b) on page 10 of "Evolution PDEs and augmented eigenfunctions. Finite interval."
##########################################################################################################################################################
# Importing packages and modules
##########################################################################################################################################################
using SymPy
using QuadGK
# using HCubature
using ApproxFun
using Roots
using Gadfly
include("C:\\Users\\LinFan Xiao\\Academics\\College\\Capstone\\work_in_julia\\construct_adjoint.jl")
##########################################################################################################################################################
# Helper functions
# Assign string as variable name
function assign(s::AbstractString, v::Any)
    s=Symbol(s)
    @eval (($s) = ($v))
end
##########################################################################################################################################################
# Function addition (f + g)(x) := f(x) + g(x)
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

# Function multiplication (f * g)(x) := f(x) * g(x)
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
##########################################################################################################################################################
# Structs
##########################################################################################################################################################

##########################################################################################################################################################
# Functions
##########################################################################################################################################################
# Get M+, M- in (2.13a), (2.13b)
function get_MPlusMinus(adjointU::VectorBoundaryForm, lambda::Number)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    alpha = e^(2pi*im/n)
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

# Get M in (2.14)
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

# Get Xlj, which is the (n-1)*(n-1) submatrix of M with (1,1) entry the (l + 1, j + 1) entry of M
function get_Xlj(adjointU::VectorBoundaryForm, M::Array, lambda::Number, l::Number, j::Number)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    M = get_M(adjointU, lambda)
    MBlock = [M M; M M]
    Xlj = MBlock[(l+1):(l+1+n-2), (j+1):(j+1+n-2)]
    return Xlj
end

# Get F+, F- in (2.16a), (2.16b)
function get_FPlusMinus(adjointU::VectorBoundaryForm, f::Union{Function, Number}, lambda::Number)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    alpha = e^(2pi*im/n)
    (MPlus, MMinus) = get_MPlusMinus(adjointU, lambda)
    M = get_M(adjointU, lambda)
    M = convert(Array{Complex}, M)
    delta = det(M)
    sumPlus, sumMinus = 0, 0
    for l = 1:n
        summandPlus, summandMinus = 0, 0
        for j = 1:n
            Xlj = get_Xlj(adjointU, M, lambda, l, j)
            Xlj = convert(Array{Complex}, Xlj)
            g(x) = e^(-im*alpha^(l-1)*lambda*x)
            integrand = mult_func(g, f) # Assume f is a finite sum of Chebyshev polynomials
            summandPlus = (-1)^((n-1)*(l+j)) * det(Xlj) * MPlus[1,j] * (quadgk(integrand, 0, 1)[1])
            summandMinus = (-1)^((n-1)*(l+j)) * det(Xlj) * MMinus[1,j] * (quadgk(integrand, 0, 1)[1])
        end
        sumPlus += summandPlus
        sumMinus += summandMinus
    end
    FPlus = 1/(2pi*delta)*sumPlus
    FMinus = (-e^(-im*lambda))/(2pi*delta)*sumMinus
    return (FPlus, FMinus)
end

# Get F in (2.15a)
function get_F(adjointU::VectorBoundaryForm, f::Union{Function, Number}, lambda::Number)
    (FPlus, FMinus) = get_FPlusMinus(adjointU, f, lambda)
    # if lambda in GammaPlusDomain return FPlus else return FMinus
end

# Get f in (2.15b)
function get_f(Gamma)
    f(x) = e^(im*lambda*x) # * integral(Gamma, F(lambda))
    return f
end

# Returns a function that shifts a point in the interval [a,b] to a corresponding point in [-1,1]
function shift_interval(originalInterval::Tuple{Number,Number}; targetInterval = (-1,1), symbolic = true)
    (a,b) = originalInterval
    (c,d) = targetInterval
    if symbolic
        x = symbols("x")
        return (d-c)*(x-a)/(b-a)+c
    else
        return x -> (d-c)*(x-a)/(b-a)+c
    end
end

# Get the nth term in the Chebyshev series expansion on [-1,1] as a Julia Function object or its symbolic expression
# https://en.wikipedia.org/wiki/Chebyshev_polynomials
function get_ChebyshevSeriesTerm(n::Int64; symbolic = true)
    if symbolic
        x = symbols("x")
        return cos(n*acos(x))
    else
        return x -> cos(n*acos(x))
    end
end

# Get the Chebyshev approximation of a function in range [a,b] as a Julia Function object or its symbolic expression
function get_ChebyshevApproximation(f::Function, interval::Tuple{Number,Number}; symbolic = true)
    (a,b) = interval
    fCheb = ApproxFun.Fun(f, a..b) # Approximate f on [a,b] using chebyshev polynomials
    chebCoefficients = ApproxFun.coefficients(fCheb) # get coefficients of the Chebyshev polynomial
    n = length(chebCoefficients)
    if symbolic
        fChebApproxSym = 0
        for i = 1:n
            chebCoefficient = chebCoefficients[i]
            chebTermSym = get_ChebyshevSeriesTerm(i-1; symbolic = symbolic)
            x = free_symbols(chebTermSym)
            if !isempty(x)
                x = x[1,1] # Assume there is only one free symbol
            end
            # Chebyshev approximation is best on [-1,1]. Thus, we compute the Chebyshev approximation on [a,b] by shifting [a,b] to [-1,1] and evaluating the Chebyshev approximation there.
            # Here, we do that by replacing x in the symbolic expression with the symbolic expression after the interval shift.
            fChebApproxSym += chebCoefficient * subs(chebTermSym, x, shift_interval((a,b); symbolic = symbolic))
        end
        return fChebApproxSym
    else
        function fChebApprox(x)
            sum = 0
            for i = 1:n
                chebCoefficient = chebCoefficients[i]
                # Here, we compose the interval shift function with the Chebyshev term.
                chebTerm = get_ChebyshevSeriesTerm(i-1; symbolic = symbolic)
                summand = mult_func(chebCoefficient, x -> chebTerm(shift_interval((a,b); symbolic = symbolic)(x)))
                sum = add_func(sum, summand)
            end
            return sum(x)
        end
        function finalChebApprox(x)
            if x >= a && x <= b # x in [a,b]
                return fChebApprox(x)
            else # if x is outside [a,b], make the approximation 0
                return 0.0
            end
        end
        return finalChebApprox
    end
end

# Find the angles of the lines characterizing gammaA, boundary of the domain {\lambda\in \C: Re(a*\lambda^n)>0} in [0, 2pi)
function find_gammaAAngles(a::Number, n::Int; symbolic = true)
    # thetaA = argument(a)
    thetaA = angle(a)
    thetaStartList = Array{Number}(n,1) # List of where domain sectors start
    thetaEndList = Array{Number}(n,1) # List of where domain sectors end
    if symbolic
        k = symbols("k")
        counter = 0
        while N(subs((2pi*k + pi/2 - thetaA)/n, k, counter)) < 2pi
            theta1 = subs((2PI*k - PI/2 - rationalize(thetaA/pi)*PI)/n, k, counter)
            theta2 = subs((2PI*k + PI/2 - rationalize(thetaA/pi)*PI)/n, k, counter)
            counter += 1
            thetaStartList[counter] = theta1
            thetaEndList[counter] = theta2
        end
    else
        k = 0
        while (2pi*k + pi/2 - thetaA)/n < 2pi
            theta1 = (2pi*k - pi/2 - thetaA)/n
            theta2 = (2pi*k + pi/2 - thetaA)/n
            k += 1
            thetaStartList[k] = theta1
            thetaEndList[k] = theta2
        end
    end
    return (thetaStartList, thetaEndList)
end

# Plot the contour sectors by sampling points
function contour_tracing(a::Number, n::Int, sampleSize::Int)
    lambdaVec = []
    for counter = 1:sampleSize
        x = rand(Uniform(-10.0,10.0), 1, 1)[1]
        y = rand(Uniform(-10.0,10.0), 1, 1)[1]
        lambda = x + y*im
        if real(a*lambda^n)>0
            append!(lambdaVec, lambda)
        end
    end
    plot(x=real(lambdaVec), y=imag(lambdaVec), Coord.Cartesian(ymin=-10,ymax=10, xmin=-10, xmax=10))
end

# Find the argument of a complex number in [0,2pi)
function argument(z::Number)
    if angle(z) >= 0
        return angle(z)
    else # Shift from (-pi, pi] to [0,2pi)
        argument = 2pi + angle(z) # This is in (0,2pi]
        if isapprox(argument, 2pi)
            return 0
        else
            return argument
        end
    end
end

# Finds distance between a complex number and a line (given by an angle in [0,2pi)
function find_distancePointLine(z::Number, theta::Number)
    if theta >= 2pi && theta < 0
        throw(error("Theta must be in [0,2pi)"))
    else
        if isapprox(argument(z), theta)
            return 0
        else
            x0, y0 = real(z), imag(z)
            if isapprox(theta, pi/2) || isapprox(theta, 3pi/2)
                return abs(x0)
            elseif isapprox(theta, 0) || isapprox(theta, 2pi)
                return abs(y0)
            else
                k = tan(theta)
                x = (y0+1/k*x0)/(k+1/k)
                y = k*x
                distance = norm(z-(x+im*y))
                return distance
            end
        end
    end
end

# Returns the minimum of the pairwise distances between zeroes in zeroList
function find_epsilon(zeroList::Array, a::Number, n::Int)
    (thetaStartList, thetaEndList) = find_gammaAAngles(a, n, symbolic = false)
    # List of distances between each zero and each sector boundary
    pointLineDistances = [find_distancePointLine(z, theta) for z in zeroList for theta in collect(Iterators.flatten([thetaStartList, thetaEndList]))]
    pointLineDistances = pointLineDistances[pointLineDistances.>0]
    if length(zeroList)>1
        # List of distance between every two zeroes
        pairwiseDistances = [norm(z1-z2) for z1 in zeroList for z2 in zeroList]
        pairwiseDistances = pairwiseDistances[pairwiseDistances.>0]
    else
        pairwiseDistances = []
    end
    distances = collect(Iterators.flatten([pairwiseDistances, pointLineDistances]))
    epsilon = minimum(distances)/2
    return epsilon
end

# Returns an array of four complex numbers representing the vertices of a square around the zero; each vertex is of distance epsilon from the zero.
function draw_squareAroundZero(zero::Number, epsilon::Number, a::Number, n::Int)
    z = zero
    theta = argument(zero)
    z1 = z - epsilon*e^(im*theta)
    z2 = z + epsilon*e^(im*(theta+pi/2))
    z3 = z + epsilon*e^(im*theta)
    z4 = z + epsilon*e^(im*(theta-pi/2))
    return [z1, z2, z3, z4] # The four vertices of the square in clockwise order
end

# Find the contours gamma_a+, gamma_a-, gamma_0+, gamma_0-
function find_gamma(a::Number, n::Int, zeroList::Array, infinity::Number)
    (thetaStartList, thetaEndList) = find_gammaAAngles(a, n; symbolic = false)
    gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus = [], [], [], []
    epsilon = find_epsilon(zeroList, a, n)
    for i in 1:n
        thetaStart = thetaStartList[i]
        thetaEnd = thetaEndList[i]
        # Initialize the boundary of each sector with the ending boundary, the origin, and the starting boundary (start and end boundaries refer to the order in which the boundaries are passed if tracked counterclockwise)
        initialPath = [infinity*e^(im*thetaEnd), 0+0*im, infinity*e^(im*thetaStart)]
        if thetaStart >= 0 && thetaStart <= pi && thetaEnd >= 0 && thetaEnd <= pi # if in the upper half plane, push the boundary path to gamma_a+
            push!(gammaAPlus, initialPath) # list of lists
        else # if in the lower half plane, push the boundary path to gamma_a-
            push!(gammaAMinus, initialPath)
        end
    end
    # Sort the zeroList by norm, so that possible zero at the origin comes last. We need to leave the origin in the initial path unchanged until we have finished dealing with all non-origin zeros because we use the origin in the initial path as a reference point to decide where to insert the deformed path
    zeroList = sort(zeroList, lt=(x,y)->!isless(norm(x), norm(y)))
    for zero in zeroList
        # If zero is not at the origin
        if !isapprox(zero, 0+0*im)
            # Draw a square around it
            (z1, z2, z3, z4) = draw_squareAroundZero(zero, epsilon, a, n) # z1, z3 are on the sector boundary, z2, z4 are on the interior or exterior
            # If zero is on the boundary of some sector
            if any(i -> isapprox(argument(zero), thetaStartList[i]) || isapprox(argument(zero), thetaEndList[i]), 1:n)
                # if z2 is interior to any sector, include z4 in the contour approximation
                if any(i -> argument(z2)>thetaStartList[i] && argument(z2)<thetaEndList[i], 1:n)
                    # Find which sector z2 is in
                    index = find(i -> argument(z2)>thetaStartList[i] && argument(z2)<thetaEndList[i], 1:n)[1]
                    squarePath = [z1, z4, z3]
                    thetaStart = thetaStartList[index]
                    thetaEnd = thetaEndList[index]
                    # If this sector is in the upper half plane, deform gamma_a+
                    if thetaStart >= 0 && thetaStart <= pi && thetaEnd >= 0 && thetaEnd <= pi
                        deformedPath = gammaAPlus[index]
                        if any(i -> isapprox(argument(zero), thetaStartList[i]), 1:n) # if zero is on the starting boundary, insert the square path after 0+0*im
                            splice!(deformedPath, length(deformedPath):(length(deformedPath)-1), squarePath)
                        else # if zero is on the ending boundary, insert the square path before 0+0*im
                            splice!(deformedPath, 2:1, squarePath)
                        end
                        gammaAPlus[index] = deformedPath
                    else # if sector is in the lower half plane, deform gamma_a-
                        deformedPath = gammaAMinus[index-length(gammaAPlus)]
                        if any(i -> isapprox(argument(zero), thetaStartList[i]), 1:n) # if zero is on the starting boundary, insert the square path after 0+0*im
                            splice!(deformedPath, length(deformedPath):(length(deformedPath)-1), squarePath)
                        else # if zero is on the ending boundary, insert the square path before 0+0*im
                            splice!(deformedPath, 2:1, squarePath) 
                        end
                        gammaAMinus[index-length(gammaAPlus)] = deformedPath
                    end
                else # if z2 is exterior, include z2 in the contour approximation
                    # Find which sector z4 is in
                    index = find(i -> argument(z4)>thetaStartList[i] && argument(z4)<thetaEndList[i], 1:n)[1]
                    squarePath = [z1, z2, z3]
                    thetaStart = thetaStartList[index]
                    thetaEnd = thetaEndList[index]
                    # If this sector is in the upper half plane, deform gamma_a+
                    if thetaStart >= 0 && thetaStart <= pi && thetaEnd >= 0 && thetaEnd <= pi
                        deformedPath = gammaAPlus[index]
                        if any(i -> isapprox(argument(zero), thetaStartList[i]), 1:n) # if zero is on the starting boundary, insert the square path after 0+0*im
                            splice!(deformedPath, length(deformedPath):(length(deformedPath)-1), squarePath)
                        else # if zero is on the ending boundary, insert the square path before 0+0*im
                            splice!(deformedPath, 2:1, squarePath)
                        end
                        gammaAPlus[index] = deformedPath
                    else # if sector is in the lower half plane, deform gamma_a-
                        deformedPath = gammaAMinus[index-length(gammaAPlus)]
                        if any(i -> isapprox(zero, thetaStartList[i]), 1:n) # if zero is on the starting boundary, insert the square path after 0+0*im
                            splice!(deformedPath, length(deformedPath):(length(deformedPath)-1), squarePath)
                        else # if zero is on the ending boundary, insert the square path before 0+0*im
                            splice!(deformedPath, 2:1, squarePath)
                        end
                        gammaAMinus[index-length(gammaAPlus)] = deformedPath
                    end
                end
                # Sort each sector's path in the order in which they are integrated over
                gammaAs = [gammaAPlus, gammaAMinus]
                for j = 1:length(gammaAs)
                    gammaA = gammaAs[j]
                    for k = 1:length(gammaA)
                        inOutPath = gammaA[k]
                        originIndex = find(x->x==0+0*im, inOutPath)[1]
                        inwardPath = inOutPath[1:(originIndex-1)]
                        outwardPath = inOutPath[(originIndex+1):length(inOutPath)]
                        # Sort the inward path and outward path
                        if length(inwardPath) > 0
                            inwardPath = sort(inwardPath, lt=(x,y)->!isless(norm(x), norm(y)))
                        end
                        if length(outwardPath) > 0
                            outwardPath = sort(outwardPath, lt=(x,y)->isless(norm(x), norm(y)))
                        end
                        inOutPath = vcat(inwardPath, 0+0*im, outwardPath)
                        gammaA[k] = inOutPath
                    end
                    gammaAs[j] = gammaA 
                end
                gammaAPlus, gammaAMinus = gammaAs[1], gammaAs[2]
            # elseif any(i -> argument(zero)>thetaStartList[i] && argument(zero)<thetaEndList[i], 1:n) # If zero is in any of the sectors, ignore it
            elseif all(i -> argument(zero)<thetaStartList[i] || argument(zero)>thetaEndList[i], 1:n) # If zero is exterior to the sectors, avoid it
                squarePath = [z1, z2, z3, z4, z1] # counterclockwise
                # If zero is in the upper half plane, add squarePath to gamma_0+
                if argument(zero) >= 0 && argument(zero) <= pi
                    push!(gamma0Plus, squarePath)
                else # If zero is in the lower half plane, add squarePath to gamma_0-
                    push!(gamma0Minus, squarePath)
                end
            end
        else # If zero is at the origin, we deform all sectors
            # deform each sector in gamma_a+
            for i = 1:length(gammaAPlus)
                deformedPath = gammaAPlus[i]
                # find the index of the zero at origin in the sector boundary path
                index = find(j -> isapprox(deformedPath[j], 0+0*im), 1:length(deformedPath))
                # If the origin is not in the path, then it has already been bypassed
                if isempty(index)
                else # If not, find its index
                    index = index[1]
                end
                # create a path around zero (origin); the origin will not be the first or the last point in any sector boundary because it was initialized to be in the middle, and only insertions are performed. Moreover, the boundary path has already been sorted into the order in which they will be integrated over, so squarePath defined below has deformedPath[index-1], deformedPath[index+1] in the correct order.
                squarePath = [epsilon*e^(im*argument(deformedPath[index-1])), epsilon*e^(im*argument(deformedPath[index+1]))]
                # replace the zero with the deformed path
                deleteat!(deformedPath, index) # delete the origin
                splice!(deformedPath, index:(index-1), squarePath) # insert squarePath into where the origin was at
                gammaAPlus[i] = deformedPath
            end
            # deform each sector in gamma_a-
            for i = 1:length(gammaAMinus)
                deformedPath = gammaAMinus[i]
                if isempty(index)
                else
                    index = index[1]
                end
                squarePath = [epsilon*e^(im*argument(deformedPath[index-1])), epsilon*e^(im*argument(deformedPath[index+1]))]
                deleteat!(deformedPath, index)
                splice!(deformedPath, index:(index-1), squarePath)
                gammaAMinus[i] = deformedPath
            end
        end
    end
    return (gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus)
end

# Visualize gamma contour
function plot_contour(gamma)
    sectorPathList = Array{Any}(length(gamma),1)
    for i = 1:length(gamma)
        # For each sector path in the gamma contour, plot the points in the path and connect them in the order in which they appear in the path
        sectorPath = gamma[i]
        sectorPathList[i] = layer(x = real(sectorPath), y = imag(sectorPath), Geom.line(preserve_order=true))
    end
    coord = Coord.cartesian(xmin=-10.0, xmax=10.0, ymin=-10.0, ymax=10.0, fixed=true)
    plot(coord, sectorList...)
end