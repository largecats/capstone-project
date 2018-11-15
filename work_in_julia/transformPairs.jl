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
# Get M+, M- in (2.13a), (2.13b) as functions of lambda (for fixed adjointU)
function get_MPlusMinus(adjointU::VectorBoundaryForm)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    alpha = e^(2pi*im/n)
    function MPlus(lambda::Number)
        MPlusMat = Array{Number}(n,n)
        for k = 1:n
            for j = 1:n
                sumPlus = 0
                for r = 0:(n-1)
                    summandPlus = (-im*alpha^(k-1)*lambda)^r * bStar[j,r+1]
                    sumPlus += summandPlus
                end
                MPlusMat[k,j] = sumPlus
            end
        end
        return MPlusMat
    end
    function MMinus(lambda::Number)
        MMinusMat = Array{Number}(n,n)
        for k = 1:n
            for j = 1:n
                sumMinus = 0
                for r = 0:(n-1)
                    summandMinus = (-im*alpha^(k-1)*lambda)^r * betaStar[j,r+1]
                    sumMinus += summandMinus
                end
                MMinusMat[k,j] = sumMinus
            end
        end
        return MMinusMat
    end
    return (MPlus, MMinus)
end

# Get M in (2.14) as a function of lambda (for fixed adjointU)
function get_M(adjointU::VectorBoundaryForm)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    alpha = e^(2pi*im/n)
    function M(lambda::Number)
        (MPlus, MMinus) = get_MPlusMinus(adjointU)
        MPlusLambda, MMinusLambda = MPlus(lambda), MMinus(lambda)
        MLambda = Array{Number}(n,n)
        for k = 1:n
            for j = 1:n
                MLambda[k,j] = MPlusLambda[k,j] + MMinusLambda[k,j] * e^(-im*alpha^(k-1)*lambda)
            end
        end
        return MLambda
    end
    return M
end

# Get delta := det(M) as a function of M(lambda) (for fixed adjointU)
function get_delta(adjointU::VectorBoundaryForm)
    function delta(lambda::Number)
        M = get_M(adjointU)
        MLambda = convert(Array{Complex}, M(lambda))
        return det(MLambda)
    end
    return delta
end

# Get Xlj, which is the (n-1)*(n-1) submatrix of M(lambda) with (1,1) entry the (l + 1, j + 1) entry of M(lambda), as a function of lambda (for fixed adjointU, M, l, j)
function get_Xlj(adjointU::VectorBoundaryForm, M::Function, l::Number, j::Number)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    M = get_M(adjointU)
    function Xlj(lambda::Number)
        MLambda = M(lambda)
        MLambdaBlock = [MLambda MLambda; MLambda MLambda]
        XljLambda = MLambdaBlock[(l+1):(l+1+n-2), (j+1):(j+1+n-2)]
        return XljLambda
    end
    return Xlj
end

# Get F+_\lambda, F-_\lambda in (2.16a), (2.16b) as a function of f
function get_FPlusMinusLambda(adjointU::VectorBoundaryForm, lambda::Number)
    bStar, betaStar = adjointU.M, adjointU.N
    n = size(bStar)[1]
    alpha = e^(2pi*im/n)
    (MPlus, MMinus) = get_MPlusMinus(adjointU)
    MPlusLambda, MMinusLambda = MPlus(lambda), MMinus(lambda)
    M = get_M(adjointU)
    MLambda = convert(Array{Complex}, M(lambda))
    deltaLambda = det(MLambda) # or deltaLambda = (get_delta(adjointU))(lambda)
    function FPlusLambda(f)
        sumPlus = 0
        for l = 1:n
            summandPlus = 0
            for j = 1:n
                Xlj = get_Xlj(adjointU, M, l, j)
                XljLambda = convert(Array{Complex}, Xlj(lambda))
                g(x) = e^(-im*alpha^(l-1)*lambda*x)
                integrand = mult_func(g, f) # Assume f is a finite sum of Chebyshev polynomials
                summandPlus = (-1)^((n-1)*(l+j)) * det(XljLambda) * MPlusLambda[1,j] * (quadgk(integrand, 0, 1)[1])
            end
            sumPlus += summandPlus
        end
        return 1/(2pi*deltaLambda)*sumPlus
    end
    function FMinusLambda(f)
        sumMinus = 0
        for l = 1:n
            summandMinus = 0
            for j = 1:n
                Xlj = get_Xlj(adjointU, M, l, j)
                XljLambda = convert(Array{Complex}, Xlj(lambda))
                g(x) = e^(-im*alpha^(l-1)*lambda*x)
                integrand = mult_func(g, f) # Assume f is a finite sum of Chebyshev polynomials
                summandMinus = (-1)^((n-1)*(l+j)) * det(XljLambda) * MMinusLambda[1,j] * (quadgk(integrand, 0, 1)[1])
            end
            sumMinus += summandMinus
        end
        return (-e^(-im*lambda))/(2pi*deltaLambda)*sumMinus
    end
    return (FPlusLambda, FMinusLambda)
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

# Find the angles of the lines characterizing gammaA, boundary of the domain {\lambda\in \C: Re(a*\lambda^n)>0} in [-2pi, 2pi). We make the interval like this to ensure that [z is in a sector] <=> [angle(z) >= sectorStart && angle(z) <= sectorEnd]
function find_gammaAAngles(a::Number, n::Int; symbolic = false)
    # thetaA = argument(a)
    thetaA = angle(a)
    thetaStartList = Array{Number}(n) # List of where domain sectors start
    thetaEndList = Array{Number}(n) # List of where domain sectors end
    if symbolic
        k = symbols("k")
        counter = 0
        while N(subs((2pi*k + pi/2 - thetaA)/n, k, counter)) < 2pi
            thetaStart = subs((2PI*k - PI/2 - rationalize(thetaA/pi)*PI)/n, k, counter)
            thetaEnd = subs((2PI*k + PI/2 - rationalize(thetaA/pi)*PI)/n, k, counter)
            counter += 1
            thetaStartList[counter] = thetaStart
            thetaEndList[counter] = thetaEnd
            # append!(thetaStartList, thetaStart)
            # append!(thetaEndList, thetaEnd)
        end
    else
        k = 0
        while (2pi*k + pi/2 - thetaA)/n < 2pi
            thetaStart = (2pi*k - pi/2 - thetaA)/n
            thetaEnd = (2pi*k + pi/2 - thetaA)/n
            k += 1
            thetaStartList[k] = thetaStart
            thetaEndList[k] = thetaEnd
            # append!(thetaStartList, thetaStart)
            # append!(thetaEndList, thetaEnd)
        end
    end
    return (thetaStartList, thetaEndList)
end

function isApproxLess(x::Number, y::Number; atol = 1e-15)
    return !isapprox(x,y; atol = atol) && x < y
end

function isApprox(x::Number, y::Number; atol = 1e-15)
    return isapprox(x, y; atol = atol)
end

# Split sectors that contain the real line
function find_gammaAAnglesSplit(a::Number, n::Int; symbolic = false)
    (thetaStartList, thetaEndList) = find_gammaAAngles(a, n; symbolic = symbolic)
    zeroIndex = find(i -> ((isApproxLess(thetaStartList[i], 0) && isApproxLess(0, thetaEndList[i]))), 1:n)
    if !isempty(zeroIndex)
        index = zeroIndex[1]
        # Insert 0 after thetaStart
        splice!(thetaStartList, (index+1):index, 0)
        # Insert 0 before thetaEnd
        splice!(thetaEndList, index:(index-1), 0)
    end
    piIndex = find(i -> ((isApproxLess(thetaStartList[i], pi) && isApproxLess(pi, thetaEndList[i]))), 1:n)
    if !isempty(piIndex)
        index = piIndex[1]
        if symbolic
            # Insert pi after thetaStart
            splice!(thetaStartList, (index+1):index, PI)
            # Insert pi before thetaEnd
            splice!(thetaEndList, index:(index-1), PI)
        else
            # Insert pi after thetaStart
            splice!(thetaStartList, (index+1):index, pi*1)
            # Insert pi before thetaEnd
            splice!(thetaEndList, index:(index-1), pi*1)
        end
    end
    return (thetaStartList, thetaEndList)
end

# Plot the contour sectors by sampling points
function contour_tracing(a::Number, n::Int, infinity::Number, sampleSize::Int)
    lambdaVec = []
    for counter = 1:sampleSize
        x = rand(Uniform(-infinity,infinity), 1, 1)[1]
        y = rand(Uniform(-infinity,infinity), 1, 1)[1]
        lambda = x + y*im
        if real(a*lambda^n)>0
            append!(lambdaVec, lambda)
        end
    end
    plot(x=real(lambdaVec), y=imag(lambdaVec), Coord.Cartesian(ymin=-infinity,ymax=infinity, xmin=-infinity, xmax=infinity, fixed = true))
end

# Find the argument of a complex number in [0,2pi)
function argument(z::Number)
    if angle(z) >= 0 # in [0,pi]
        return angle(z)
    else # Shift from (-pi, 0] to [pi,2pi)
        argument = 2pi + angle(z) # This is in (pi,2pi]
        if isApprox(argument, 2pi) # Change to [pi,2pi)
            return 0
        else
            return argument
        end
    end
end

# Find distance between a complex number and a line (given by an angle in [0,2pi)
function find_distancePointLine(z::Number, theta::Number)
    if theta >= 2pi && theta < 0
        throw(error("Theta must be in [0,2pi)"))
    else
        if isApprox(argument(z), theta)
            return 0
        else
            x0, y0 = real(z), imag(z)
            if isApprox(theta, pi/2) || isApprox(theta, 3pi/2)
                return abs(x0)
            elseif isApprox(theta, 0) || isApprox(theta, 2pi)
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

# Returns the minimum of the pairwise distances between zeroes in zeroList that are not interior to any sector (since interior zeroes would not matter in any way)
function find_epsilon(zeroList::Array, a::Number, n::Int)
    (thetaStartList, thetaEndList) = find_gammaAAnglesSplit(a, n; symbolic = false)
    thetaStartEndList = collect(Iterators.flatten([thetaStartList, thetaEndList]))
    truncZeroList = []
    for zero in zeroList
        # If zero is interior to any sector, discard it
        if any(i -> pointInSector(zero, (thetaStartList[i], thetaEndList[i])), 1:n)
        else # If not, append it to truncZeroList
            append!(truncZeroList, zero)
        end
    end
    pointLineDistances = [find_distancePointLine(z, theta) for z in zeroList for theta in thetaStartEndList]
    if length(truncZeroList)>1
        # List of distance between every two zeroes
        pairwiseDistances = [norm(z1-z2) for z1 in zeroList for z2 in truncZeroList]
    else
        pairwiseDistances = []
    end
    distances = collect(Iterators.flatten([pairwiseDistances, pointLineDistances]))
    # Distances of nearly 0 could be instances where the zero is actually on some sector boundary
    distances = filter(x -> !isApprox(x, 0), distances)
    epsilon = minimum(distances)/3
    return epsilon
end

# Returns an array of four complex numbers representing the vertices of an n-gon around the zero; each vertex is of distance epsilon from the zero.
function draw_nGonAroundZero(zero::Number, epsilon::Number, n::Int)
    z = zero
    theta = argument(zero)
    deltaAngle = 2pi/n
    vertices = []
    for i = 1:n
        newAngle = pi-deltaAngle*(i-1)
        vertex = z + epsilon*e^(im*(theta+newAngle))
        append!(vertices, vertex)
    end
    # vertices = vcat(vertices, vertices[1])
    return vertices
end

# Determine whether a point z is on the boundary of a sector characterized by a start angle and an end angle
function pointOnSector(z::Number, sectorAngles::Tuple{Number, Number})
    (startAngle, endAngle) = sectorAngles
    return isApprox(argument(z), startAngle) || isApprox(argument(z), endAngle) || isApprox(angle(z), startAngle) || isApprox(angle(z), endAngle)
end

# Determine whether a point z is in the interior of a sector characterized by a start angle and an end angle
function pointInSector(z::Number, sectorAngles::Tuple{Number, Number})
    (startAngle, endAngle) = sectorAngles
    # First check if z is on the sector boundary
    if pointOnSector(z, sectorAngles)
        return false
    else
        # angle(z) would work if it's in the sector with positive real parts and both positive and negative imaginary parts; argument(z) would work if it's in the sector with negative real parts and both positive and negative imaginary parts
        return (angle(z) > startAngle && angle(z) < endAngle) || (argument(z) > startAngle && argument(z) < endAngle)
    end
end

# Determine whether a point z is in the exterior of a sector characterized by a start angle and an end angle
function pointExSector(z::Number, sectorAngles::Tuple{Number, Number})
    return !pointOnSector(z, sectorAngles) && !pointInSector(z, sectorAngles)
end

# Find the contours gamma_a+, gamma_a-, gamma_0+, gamma_0-
function find_gamma(a::Number, n::Int, zeroList::Array, infinity::Number; nGon = 8)
    (thetaStartList, thetaEndList) = find_gammaAAnglesSplit(a, n; symbolic = false)
    nSplit = length(thetaStartList)
    gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus = [], [], [], []
    epsilon = find_epsilon(zeroList, a, n)
    for i in 1:nSplit
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
        println(zero)
        # If zero is not at the origin
        if !isApprox(zero, 0+0*im)
            # Draw an n-gon around it
            vertices = draw_nGonAroundZero(zero, epsilon, nGon)
            # If zero is on the boundary of some sector
            if any(i -> pointOnSector(zero, (thetaStartList[i], thetaEndList[i])), 1:nSplit)
                # Find which sector(s) zero is on
                indices = find(i -> pointOnSector(zero, (thetaStartList[i], thetaEndList[i])), 1:nSplit)
                # If zero is on the boundary of one sector
                if length(indices) == 1
                    # if vertices[2] is interior to any sector, include vertices on the other half of the n-gon in the contour approximation
                    z0 = vertices[2]
                    if any(i -> pointInSector(z0, (thetaStartList[i], thetaEndList[i])), 1:nSplit)
                        # Find which sector vertices[2] is in
                        index = find(i -> pointInSector(z0, (thetaStartList[i], thetaEndList[i])), 1:nSplit)[1]
                    else # if vertices[2] is exterior, include vertices on this half of the n-gon in the contour approximation
                        # Find which sector vertices[length(vertices)] is in
                        z1 = vertices[length(vertices)]
                        index = find(i -> pointInSector(z1, (thetaStartList[i], thetaEndList[i])), 1:nSplit)[1]
                    end
                    thetaStart = thetaStartList[index]
                    thetaEnd = thetaEndList[index]
                    # Find all vertices exterior to or on the boundary of this sector, which would form the nGonPath around the zero
                    nGonPath = vertices[find(vertex -> !pointInSector(vertex, (thetaStart, thetaEnd)), vertices)]
                    # If this sector is in the upper half plane, deform gamma_a+
                    if thetaStart >= 0 && thetaStart <= pi && thetaEnd >= 0 && thetaEnd <= pi
                        deformedPath = gammaAPlus[index]
                        if any(i -> isApprox(argument(zero), thetaStartList[i]) || isApprox(angle(zero), thetaStartList[i]), 1:nSplit) # if zero is on the starting boundary, insert the n-gon path after 0+0*im
                            splice!(deformedPath, length(deformedPath):(length(deformedPath)-1), nGonPath)
                        else # if zero is on the ending boundary, insert the n-gon path before 0+0*im
                            splice!(deformedPath, 2:1, nGonPath)
                        end
                        gammaAPlus[index] = deformedPath
                    else # if sector is in the lower half plane, deform gamma_a-
                        deformedPath = gammaAMinus[index-length(gammaAPlus)]
                        if any(i -> isApprox(argument(zero), thetaStartList[i]) || isApprox(angle(zero), thetaStartList[i]), 1:nSplit) # if zero is on the starting boundary, insert the n-gon path after 0+0*im
                            splice!(deformedPath, length(deformedPath):(length(deformedPath)-1), nGonPath)
                        else # if zero is on the ending boundary, insert the n-gon path before 0+0*im
                            splice!(deformedPath, 2:1, nGonPath) 
                        end
                        gammaAMinus[index-length(gammaAPlus)] = deformedPath
                    end
                else # If zero is on the boundary of two sectors, then it must be on the real line, and we need to deform two sectors
                    # Find out which vertices are in the lower half plane
                    nGonPath = vertices[find(vertex -> !pointInSector(vertex, (0, pi)), vertices)]
                    for index in indices
                        thetaStart = thetaStartList[index]
                        thetaEnd = thetaEndList[index]
                        # If this is the sector in the upper half plane, deform gamma_a+
                        if thetaStart >= 0 && thetaStart <= pi && thetaEnd >= 0 && thetaEnd <= pi
                            gammaAPlusIndex = find(path -> (isApprox(argument(zero), argument(path[1])) || isApprox(argument(zero), argument(path[length(path)]))), gammaAPlus)[1]
                            deformedPath = copy(gammaAPlus[gammaAPlusIndex])
                            if isApprox(argument(zero), argument(deformedPath[length(deformedPath)])) # if zero is on the starting boundary, insert the n-gon path after 0+0*im
                                splice!(deformedPath, length(deformedPath):(length(deformedPath)-1), nGonPath)
                            else # if zero is on the ending boundary, insert the n-gon path before 0+0*im
                                splice!(deformedPath, 2:1, nGonPath)
                            end
                            gammaAPlus[gammaAPlusIndex] = deformedPath
                        else # If this is the sector in the lower half plane, deform gamma_a-
                            gammaAMinusIndex = find(path -> (isApprox(argument(zero), argument(path[1])) || isApprox(argument(zero), argument(path[length(path)]))), gammaAMinus)[1]
                            deformedPath = copy(gammaAMinus[gammaAMinusIndex])
                            if isApprox(argument(zero), argument(deformedPath[length(deformedPath)])) # if zero is on the starting boundary, insert the n-gon path after 0+0*im
                                splice!(deformedPath, length(deformedPath):(length(deformedPath)-1), nGonPath)
                            else # if zero is on the ending boundary, insert the n-gon path before 0+0*im
                                splice!(deformedPath, 2:1, nGonPath)
                            end
                            gammaAMinus[gammaAMinusIndex] = deformedPath
                        end
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
            # If zero is interior to any sector (after splitting by real line), ignore it
            # If zero is exterior to the sectors, avoid it
            elseif all(i -> pointExSector(zero, (thetaStartList[i], thetaEndList[i])), 1:nSplit)
                nGonPath = vcat(vertices, vertices[1]) # counterclockwise
                # If zero is in the upper half plane, add the n-gon path to gamma_0+
                if argument(zero) >= 0 && argument(zero) <= pi
                    push!(gamma0Plus, nGonPath)
                else # If zero is in the lower half plane, add the n-gon path to gamma_0-
                    push!(gamma0Minus, nGonPath)
                end
            end
        else # If zero is at the origin, we deform all sectors and draw an n-gon around the origin
            # deform each sector in gamma_a+
            for i = 1:length(gammaAPlus)
                deformedPath = gammaAPlus[i]
                # find the index of the zero at origin in the sector boundary path
                index = find(j -> isApprox(deformedPath[j], 0+0*im), 1:length(deformedPath))
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
                index = find(j -> isApprox(deformedPath[j], 0+0*im), 1:length(deformedPath))
                if isempty(index)
                else
                    index = index[1]
                end
                squarePath = [epsilon*e^(im*argument(deformedPath[index-1])), epsilon*e^(im*argument(deformedPath[index+1]))]
                deleteat!(deformedPath, index)
                splice!(deformedPath, index:(index-1), squarePath)
                gammaAMinus[i] = deformedPath
            end
            # Draw an n-gon around the origin and add to gamma_0+
            vertices = draw_nGonAroundZero(zero, epsilon/2, nGon)
            nGonPath = vcat(vertices, vertices[1])
            push!(gamma0Plus, nGonPath)
        end
    end
    return (gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus)
end

# Visualize gamma contour
function plot_contour(gamma::Array, infinity::Number)
    sectorPathList = Array{Any}(length(gamma),1)
    for i = 1:length(gamma)
        # For each sector path in the gamma contour, plot the points in the path and connect them in the order in which they appear in the path
        sectorPath = gamma[i]
        # labels = map(string, collect(1:1:length(sectorPath)))
        sectorPathList[i] = layer(x = real(sectorPath), y = imag(sectorPath), Geom.line(preserve_order=true))
    end
    coord = Coord.cartesian(xmin=-infinity, xmax=infinity, ymin=-infinity, ymax=infinity, fixed=true)
    plot(coord, sectorPathList...)
end


# Get F in (2.15a) as a function of lambda
function get_F(f::Function, adjointU::VectorBoundaryForm, a::Number, n::Number, zeroList::Array, infinity::Number)
    (gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus) = find_gamma(a, n, zeroList, infinity)
    function F(lambda::Number)
        (FPlusLambda, FMinusLambda) = get_FPlusMinusLambda(adjointU, lambda)
        # if lambda in gamma+ return FPlusLambda(f); if lambda in gamma- return FMinusLambda(f)
        if lambda in gammaAPlus || lambda in gamma0Plus
            return FPlusLambda(f)
        elseif lambda in gammaAMinus || lambda in gamma0Minus
            return FMinusLambda(f)
        end
    end
    return F
end

# Get f in (2.15b) as a function of x. The argument gamma may be cconstructed as follows.
# (gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus) = find_gamma(a, n, zeroList, infinity)
# gamma = collect(Iterators.flatten([gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus]))
function get_f(gamma::Array, F::Function)
    function f(x::Number)
        if x <= 1 && x >= 0
            return quadgk(e^(im*lambda*x)*F(lambda), gamma...)[1]
        else
            throw(error("f_x(F) is only defined for x in [0,1]"))
        end
    end
    return f
end

# Function to find (approximately) zeroes of delta(lambda)
function get_zeroList(delta::Function, infinity::Number)
    # deltaChebApprox = Fun(delta, ) # domain of lambda is the complex plane, approximated by the circle of radius infinity
end

# Function to solve a given IBVP with given L and U
function solve_IBVP(L::LinearDifferentialOperator, U::VectorBoundaryForm, pDerivMatrix::Array, lambda::Number, a::Number, zeroList::Array, infinity::Number, f::Function) # f is given
    n = length(L.pFunctions)-1
    adjointU = construct_validAdjoint(L, U, pDerivMatrix)
    F = get_F(f, adjointU, a, n, zeroList, infinity)
    function q(x,t)
        return f(e^(-a*lambda^n*F(lambda)))
    end
    return q
end