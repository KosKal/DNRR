cd("INSERT_PATH")

using Distributions
# using Polynomials
using PolynomialRoots
using DataFrames

include("toolBox.jl")
include("polyDNRR.jl")

#########################################################
# Parameter Declaration
#########################################################

ss = 270

## used data sets
# w1, lam1 = 0.6, 1.8e04 # f_{2,1}, (ss200), dataseeds: 23,24,21,32 , etas = 3.5,4.9,6.3,7.5 (sampseeds=94,95,93,96)
# w1, lam1 = 0.6, 3e04 # f_{2,1}, (ss200), dataseed: 44 , eta = 2.7 (sampseed=99)

w1, lam1 = 0.6, 1e05 # f_{2,1}, (ss200), dataseeds: 23,24,21,32 , etas = 3.5,4.9,6.3,7.5 (sampseeds=94,95,93,96)

lam2 = 1e-02 * lam1

# map parameters and initial condition
theta, x0 = [0.05, 2.55, 0., -0.99], 1.
# theta, x0 = [1., 0., -1.71], 0.1

# corresponding deterministic orbit
xdet = zeros(ss)
xdet[1] = polyMap(theta, x0)
for i in 2:ss
  xdet[i] = polyMap(theta, xdet[i - 1])
end

# orbit contaminated with dynamical noise 
# first argument is  type of noise: 0 is parametric, 1-3 nonparametric
dataSeed = 2
data = genData(3, ss, theta, x0, w1, lam1, lam2, dataSeed)
eta = sqrt(w1*1/lam1+(1-w1)*1/lam2)/std(data)

#########################################################
# GSB Reconstruction - Denoising
######################################################### 

# Dynamical Noise Only

degree, k, gibbsIter, burnIn, dela, delb, wdist, gEps1, gEps2, thLow, thUp, zLow, zUp, papr, pbpr, thin, samplerSeed =
  5, ss, 150000, 50000, 1e05, 1e-03, 1e-03, 1e-03, 1e-03, -10.0, 10.0, -10.0, 10.0, 1., 1., 10, 999;

thetas, x0s, noise, clusters, ps, deltas, ys = 0, 0, 0, 0, 0, 0, 0;
gc()

filename = "/Results"
savelocation = string(pwd(), filename, "/seed$samplerSeed")
mkpath(savelocation)
writedlm(string(savelocation,"/xdet.txt"), xdet)

@time thetas = polyDNRR(data, degree, k, gibbsIter, burnIn, dela, delb, wdist, gEps1, gEps2, thLow, thUp, zLow, zUp, papr, pbpr, thin, samplerSeed);

thetas = 0;
gc()


