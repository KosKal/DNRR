cd("INSERT_PATH")

using Distributions
using PolynomialRoots
using DataFrames

include("toolBox.jl")
include("quadDNRR.jl")

#########################################################
# Parameter Declaration
#########################################################

ss = 1000

# same \eta, different \sigma
# w1, lam1 = 0.6, 4.8e04 # f_{2,1},/ seed: 13 
# w1, lam1 = 0.7, 3.5e04 # f_{2,2}, dseed: 116, sseed = 122, sp = 3e-06
# w1, lam1 = 0.8, 2.5e04 # f_{2,3}, dseed: 13, sseed = 1236, sp = 4.5e-06 
w1, lam1 = 0.9, 1.3e04 # f_{2,4}, dseed: 13, sseed = 1214, sp = 4.5e-06 

# different \eta, same \sigma
# w1, lam1 = 0.6, 2e04 # f_{2,1}, seed: 8152
# w1, lam1 = 0.7, 2e04 # f_{2,2}, seed: 1134
# w1, lam1 = 0.8, 2e04 # f_{2,3}, seed: 89
# w1, lam1 = 0.9, 2e04 # f_{2,4}, seed: 58


lam2 = 1e-02 * lam1 # noise mixture 3 {dynamical noise}
# w1, lam1, lam2 = 0.3, 2.5e03, 2.5e03 # noise mixture 3 {dynamical noise}


# # map parameters {for sampler}
# # theta, x0, x00 = [1.4, 0., 0.3, 0., -1., 0.], 0.5, 1.0
# theta = [1.31, 0., 0.23, 0., -1., 0.]
theta, x0, x00 = [1.38, 0., 0.27, 0., -1., 0.], 0.5, 0.5

# x0, x00 =  0.1, 0.1 # initial conditions

# map parameters {for data generation}
# thetamap = [1.4, -1., 0.3] # pars1 (standard)
# thetamap = [1., -1.31, 0.23] # pars2 (SJH)
# thetamap = [1.38, -1., 0.27] # pars3 (Kantz)

# corresponding deterministic orbit
xdet = zeros(ss)
xdet[1] = henon(theta, x0, x00)
xdet[2] = henon(theta, xdet[1], x0)
for i in 3:ss
  xdet[i] = henon(theta, xdet[i-1], xdet[i-2])
end

# orbit contaminated with dynamical noise (latent in State Space)
# first argument is  type of noise: 0 is parametric, 1-3 nonparametric
dataSeed = 13
data = genData(3, ss, theta, x0, x00, w1, lam1, lam2, dataSeed)
eta = sqrt(w1*1/lam1+(1-w1)*1/lam2)/std(data)

# read data set
path = "/home/kkaloudis/Documents/Mathematics/PhD/Julia/Prediction"
x1 = readdlm(string(path,"/lynx.txt"))
x=zeros(114)
for i in 1:114
	x[i] = x1[i]
end

#########################################################
# GSB Reconstruction - Denoising
#########################################################

# Dynamical Noise Only

# include("quadDPDeno.jl")
# path = "/home/kkaloudis/Documents/DNRR paper codes/full quadratic/Results"
# data = readdlm(string(path,"/data1.txt"))
# xdet = readdlm(string(path,"/xdet1.txt"))

k, gibbsIter, burnIn, dela, delb, wdist, gEps1, gEps2, thLow, thUp, zLow, zUp, papr, pbpr, thin, samplerSeed =
  0, 1000000, 1000000, 1e05, 1e-03, 1e02, 1e-03, 1e-03, -10.0, 10.0, -20.0, 20.0, 0.5, 0.5, 100, 76169;


filename = "/Results"
savelocation = string(pwd(), filename, "/seed$samplerSeed")
mkpath(savelocation)
writedlm(string(savelocation,"/xdet.txt"), xdet)

## DNRR application
@time thetas = quadDNRR(x, k, gibbsIter, burnIn, dela, delb, wdist, gEps1, gEps2, thLow,
  thUp, zLow, zUp, papr, pbpr, thin, samplerSeed);

thetas = 0;
gc()

# ################################################################
# ################################################################


