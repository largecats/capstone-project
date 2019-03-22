# Plotting

a = -im
n = 5
INFTY = 10
sampleSize = 10000
simulation = trace_contour(a, n, sampleSize; infty = INFTY)
Gadfly.draw(PDF("C:\\Users\\LinFan Xiao\\Academics\\College\\Capstone\\reports\\report 3 (initial thesis)\\simulation.pdf"), simulation)

# zeroList = [3+3*sqrt(3)*im, 2+2*sqrt(3)*im, 0+0*im, 0+5*im, 0-5*im, 4]
zeroList = [3+3*sqrt(3)*im, 2+2*sqrt(3)*im, 0+0*im, 0+5*im, 0-5*im, 3, -5, 4-4*sqrt(3)*im]
(gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus) = get_gamma(a, n, zeroList; infty = INFTY, nGon = 8)
gamma = collect(Iterators.flatten([gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus]))
myPlot = plot_contour(gamma; infty = INFTY)
Gadfly.draw(PDF("C:\\Users\\LinFan Xiao\\Academics\\College\\Capstone\\reports\\report 3 (initial thesis)\\contourPlot.pdf"), myPlot)

lambda = symbols("lambda")
x = symbols("x", real = true)
y = symbols("y", real = true)
delta = (lambda^3+lambda+2)*e^(lambda)
delta = cos(lambda)
# delta = cos(lambda)*e^(lambda)
delta = cos(lambda)*e^(lambda) + sin(lambda)
bivariateDelta = subs(delta, lambda, x+im*y)
p = plot_levelCurves(bivariateDelta; width = 2500, height = 2000)
Plots.pdf(p, "C:\\Users\\LinFan Xiao\\Academics\\College\\Capstone\\reports\\report 3 (initial thesis)\\levelCurves.pdf")

n = 2
beta0 = -1
beta1 = -1
theta = 0
a = e^(im*theta)

t = symbols("t")
symPFunctions = [-1 0 0]
interval = (0,1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)
b = [1 0; 0 1]
beta = [beta0 0; 0 beta1]
U = VectorBoundaryForm(b, beta)

adjointU = get_adjointU(L, U)
delta = get_delta(adjointU; symbolic = true)
x = symbols("x", real = true)
y = symbols("y", real = true)
bivariateDelta = subs(delta, lambda, x+im*y)
p = plot_levelCurves(bivariateDelta; width = 2500, height = 2000)
Plots.pdf(p, "C:\\Users\\LinFan Xiao\\Academics\\College\\Capstone\\reports\\report 3 (initial thesis)\\levelCurves.pdf")

n = 2
beta0 = -1
beta1 = -1
theta = 0
a = e^(im*theta)
t = symbols("t")
symPFunctions = [-1 0 0]
interval = (0,1)
symL = SymLinearDifferentialOperator(symPFunctions, interval, t)
L = get_L(symL)
b = [1 0; 0 1]
beta = [beta0 0; 0 beta1]
U = VectorBoundaryForm(b, beta)
f(x) = sin(x*2*pi)
x = symbols("x")
fSym = sin(x*2*PI)
check_boundaryConditions(L, U, fSym)

adjointU = get_adjointU(L, U)
delta = get_delta(adjointU; symbolic = true)
lambda = free_symbols(delta)[1]
x = symbols("x", real = true)
y = symbols("y", real = true)
bivariateDelta = subs(delta, lambda, x+im*y)
p = plot_levelCurves(bivariateDelta; width = 2500, height = 2000)
Plots.pdf(p, "C:\\Users\\LinFan Xiao\\Academics\\College\\Capstone\\reports\\report 3 (initial thesis)\\levelCurves.pdf")

zeroList = [0, 2*pi, -2*pi]
(gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus) = get_gamma(a, n, zeroList)
gamma = collect(Iterators.flatten([gammaAPlus, gammaAMinus, gamma0Plus, gamma0Minus]))
p = plot_contour(gamma)
Gadfly.draw(PDF("C:\\Users\\LinFan Xiao\\Academics\\College\\Capstone\\reports\\report 3 (initial thesis)\\contourPlot.pdf"), p)