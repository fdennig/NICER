##############################################################################################################################################################
# NICER computational model used for risk, inequality and climate change
##############################################################################################################################################################
# File: Optimization.jl
# Content: template for main types of optimization of the NICER integrated assessment model
##############################################################################################################################################################
# Source: original program running in Julia 0.6.4 can be found at https://github.com/fdennig/NICER
# Changes introduced: all changes required to run with Julia 1.1.0
# File release: March 16th 2019 (BM)
##############################################################################################################################################################
# Execution requires following files:
#  ../Optimization.jl
#  ../Function_definitions.jl
#  ../createPrandom.jl
#  ../Data/certainPARAMETERS.jld
#  ../Data/dparam_i.jld
# Program execution via REPL command: include("Optimization.jl")
# Results available in folder:
#  ../Outputs
##############################################################################################################################################################


println("\n------------------------------\n   Program: optimization.jl\n------------------------------")

# Preliminaries
println("+ Loading of required libraries")
using HDF5, JLD, LinearAlgebra, CSV, DataFrames, NLopt, Distributions
# using Random  # recommanded by Julia 0.7, not required with Julia 1.1

# Select the regime for parameter randomization
  # 0 = no randomization (just uses means)
  # 1 = total randomization
  # 2 = High initial TFP growth vs. Low initial TFP growth
  # 3 = High initial decarbonization rate vs. Low initial decarbonization rate
  # 4 = High elasticity of income wrt damage vs. Low elasticity of income wrt damage
  # 5 = High climate sensitivity vs. Low climate sensitivity
  # 6 = High atmosphere to upper ocean transfer coefficient vs. Low atmosphere to upper ocean transfer coefficient
  # 7 = High initial world backstop price vs. Low initial world backstop price
  # 8 = High T7 coefficient vs. Low T7 coefficient
  # 9 = Deciles of TFP (all else fixed at means) (nsample = 10)
  # 95 same as 9, but includes medial (nsample = 11)
  # 10 = Deciles of decarbonization rates (all else fixed at means) (nsample = 10)
  # 105 same as 10, but includes medial (nsample = 11)
  # 11 = Deciles - High TFP and High decarb vs. Low TFP and Low decarb (technology spillovers?)
  # 12 = Deciles - High TFP and Low decarb vs. Low TFP and High decarb (substitutable tech?)
  # 13 = Deciles of elasticity of income wrt damage (ee)
  # 14 = Deciles - High TFP and High ee vs. Low TFP and Low ee
  # 15 = Deciles - High TFP and Low ee vs. Low TFP and High ee
  # 16 = DECILES of climate sensitivity
  # 165 same as 16, but includes medial (nsample = 11)
regime_select = 0

# Define exogenous parameters
println("+ Setup of parameters")
Tm = 32  # time period to consider, Tm <= 60
tm = 18  # (must be an integer) length of the tax vector to consider 
         # THIS WILL DIFFER FROM THE DIMENSION OF THE ACTUAL TAX VECTOR OBJECT, tax_length!
         # note that the tax vector contains 2 sets of taxes -
         # the first lm elements are common to both sets, the remainder of the vector is then split in two, one for each set of taxes (must of equal length, (tm - lm)/2)
backstop_same = "Y" # "N" is default - 
                    # choose "Y" if we want all the countries to face the same backstop prices over time
rho = 0.015 # PP[1].para[1] # discount rate
eta = 2 # PP[1].para[3] # inequality aversion/time smoothing parameter (applied to per capita consumption by region, quintile and time period)
nu = 2  # risk aversion parameter (applied over random draws)
exi = 1 # elasticity of damage with respect to income

# Now execute the whole code: select all and evaluate

###########################################################################################################################################################################################################

# Run Function Definitions File AND createPrandom first to build necessary parameters!
println("+ Loading of model function definitions")
folder =  pwd()
include("$folder/Function_definitions.jl")
include("$folder/createPrandom.jl")

PP = createP(regime_select; eeM = exi)
nsample = size(PP)[1]
pb = zeros(nsample,12,60)
for i=1:nsample
  pb[i,:,:] = PP[i].pb'
end

# Optimization of welfare function using NLopt package
idims = Int(max(round(nsample/2),1)) # bifurcates the random draws into two subsets

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#learning happens instantly. you know from period one which learning branch you are on
# 1. lm = 0
lm = 0
tax_length = 2*tm - lm
count = 0 # keep track of number of iterations

# Define function to be maximized (requires special format for NLopt package)
 model = "DICE"   # 3.5.19
# model = "NICE"   # 3.14.19
println("+ Model selected = ",model)

function welfaremax(x,grad) # x is the tax vector, grad is the gradient (unspecified here since we have no way of computing it)
  WW = tax2expectedwelfare(x,PP,rho,eta,nu,Tm,tm,lm,idims,model="$model")[1] #change to model="RICE" or "DICE" for different levels of aggregation
  global count
  count::Int += 1
  if count%100 == 0
    println("  iteration f_$count \n$x\n")
  end
  return WW
end

# Choose algorithm (gradient free method) and dimension of tax vector, tm+1 <= n <= Tm
println("+ Setup of optimization algorithm")
n = tax_length
opt = Opt(:LN_BOBYQA, n) # algorithm and dimension of tax vector, possible (derivative-free) algorithms: LN_COBYLA, LN_BOBYQA

#BM*Modified: 
#  squeeze(a,3)  => dropdims(a, dims=3)
#  maximum(A, 2) => maximum(A, dims=2)
#ub_lm = maximum(squeeze(maximum(pb,2),2),1)[2:lm+1]                                          # Julia 0.6
 ub_lm = maximum(dropdims(maximum(pb,dims=2),dims=2),dims=1)[2:lm+1]                          # Julia 1.x   
#ub_1 = maximum(squeeze(maximum(pb,2),2)[1:idims,:],1)[lm+2:tm+1]                             # Julia 0.6
 ub_1 = maximum(dropdims(maximum(pb,dims=2),dims=2)[1:idims,:],dims=1)[lm+2:tm+1]             # Julia 1.x
#ub_2 = maximum(squeeze(maximum(pb,2),2)[(idims+1):nsample,:],1)[lm+2:tm+1]                   # Julia 0.6
 ub_2 = maximum(dropdims(maximum(pb,dims=2),dims=2)[(idims+1):nsample,:],dims=1)[lm+2:tm+1]   # Julia 1.x

# lower bound is zero
lower_bounds!(opt, zeros(n))
upper_bounds!(opt, [ub_lm; ub_1; ub_2])

# Set maximization
max_objective!(opt, welfaremax)

# Set relative tolerance for the tax vector choice - vs. ftol for function value?
ftol_rel!(opt,0.00000000000005)

# Optimize! Initial guess defined above
println("+ Optimization loops\n")
(expected_welfare,tax_vector,ret) = optimize(opt, [ub_lm; ub_1; ub_2].*0.5)

# Extract the two tax vectors from the optimization object
taxes_a = tax_vector[1:tm]
taxes_b = [tax_vector[1:lm];tax_vector[tm+1:end]]

println("+ Calculation of variables from taxes")
c, K, T, E, M, mu, lam, D, AD, Y, Q = VarsFromTaxes(taxes_a, taxes_b, PP, nsample, model="$model")

TAX = maximum(PP[1].pb[1:Tm,:],dims=2)[:,1]
TAX[1] = 0
TAX[2:(tm+1)] = taxes_a
taxes_1=TAX
TAX[2:(tm+1)] = taxes_b
taxes_2 = TAX

println("+ Generation of output data file =")
# create .jld of S_lm
res = Results(regime_select,nsample,Tm,tm,lm,Regions,taxes_1,taxes_2,expected_welfare,c,K,T,E,M,mu,lam,D,AD,Y,Q,rho,eta,nu,PP)
filenm = string(regime_select)
# create dataframe of period by region by state data
dataP = FrameFromResults(res, Tm, nsample, Regions, idims)

# Creation CSV output file
   # writetable("$(pwd())/Outputs/Optima/meanOptimum$(model).csv", dataP)   # Julia 0.6
   # CSV.write("$(pwd())/Outputs/Optima/meanOptimum$(model).csv", dataP)   # Julia 1.x
CSV.write("$(pwd())/Outputs/meanOptimum$(model)Julia110.csv", dataP)
println("$(pwd())/Outputs/meanOptimum$(model)Julia110.csv")

# Creation JLD output file
#jldopen("$(pwd())/Outputs/meanOptimum$(model)Julia11.jld", "w") do file
#    write(file, "res", res)
#println("\n$(pwd())/Outputs/meanOptimum$(model)Julia11.jld")
#end


##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# # now you learn in period lm
# # 2. lm = 4
# lm = 4
# tax_length = 2*tm - lm
#
# count = 0 # keep track of number of iterations
#
# # Define function to be maximized (requires special format for NLopt package)
#
# function welfaremax(x,grad) # x is the tax vector, grad is the gradient (unspecified here since we have no way of computing it)
#   WW = tax2expectedwelfare(x,PP,rho,eta,nu,Tm,tm,lm,idims,model="$model")[1]
#   global count
#   count::Int += 1
#   println("  iteration f_$count \n$x\n")
#   return WW
# end
#
# # Choose algorithm (gradient free method) and dimension of tax vector, tm+1 <= n <= Tm
# n = tax_length
# opt = Opt(:LN_BOBYQA, n) # algorithm and dimension of tax vector, possible (derivative-free) algorithms: LN_COBYLA, LN_BOBYQA
#
# #BM-Julia-0.7: squeeze is deprecated in favor of dropdims
# #ub_lm = maximum(squeeze(maximum(pb,2),2),1)[2:lm+1]                                          # Julia 0.6
# ub_lm = maximum(dropdims(maximum(pb,dims=2),dims=2),dims=1)[2:lm+1]                           # Julia 1.x
# #ub_1 = maximum(squeeze(maximum(pb,2),2)[1:idims,:],1)[lm+2:tm+1]                             # Julia 0.6
# ub_1 = maximum(dropdims(maximum(pb,dims=2),dims=2)[1:idims,:],dims=1)[lm+2:tm+1]              # Julia 1.x
# #ub_2 = maximum(squeeze(maximum(pb,2),2)[(idims+1):nsample,:],1)[lm+2:tm+1]                   # Julia 0.6
# ub_2 = maximum(dropdims(maximum(pb,dims=2),dims=2)[(idims+1):nsample,:],dims=1)[lm+2:tm+1]    # Julia 1.x
#
# # lower bound is zero
# lower_bounds!(opt, zeros(n))
# upper_bounds!(opt, [ub_lm; ub_1; ub_2])
#
# # Set maximization
# max_objective!(opt, welfaremax)
#
# # Set relative tolerance for the tax vector choice - vs. ftol for function value?
# ftol_rel!(opt,0.00000000000005)
#
# # Optimize! Initial guess defined above
# (expected_welfare,tax_vector,ret) = optimize(opt, [ub_lm; ub_1; ub_2].*0.5)
#
# # Extract the two tax vectors from the optimization object
# taxes_a = tax_vector[1:tm]
# taxes_b = [tax_vector[1:lm];tax_vector[tm+1:end]]
#
# c, K, T, E, M, mu, lam, D, AD, Y, Q = VarsFromTaxes(taxes_a, taxes_b, PP, nsample, model="$model")
#
# TAX = maximum(PP[1].pb[1:Tm,:],dims=2)[:,1]
# TAX[1] = 0
# TAX[2:(tm+1)] = taxes_a
# taxes_1=TAX
# TAX[2:(tm+1)] = taxes_b
# taxes_2 = TAX
#
# #create Results data type as well as dataframe
# reslm = Results(regime_select,nsample,Tm,tm,lm,Regions,taxes_1,taxes_2,expected_welfare,c,K,T,E,M,mu,lam,D,AD,Y,Q,rho,eta,nu,PP)
#
# # create dataframe of period by region by state data
# using DataFrames
# dataPlm = FrameFromResults(reslm, Tm, nsample, Regions, idims)
#
# Creation CSV output file
# CSV.write("$(pwd())/Outputs/meanOptimum$(model)-Julia11Plm.csv", dataPlm)
#

##!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# # you never learn which branch you are on
# # 3. lm = tm
# lm= tm
# tax_length = 2*tm - lm
#
# count = 0 # keep track of number of iterations
#
# # Define function to be maximized (requires special format for NLopt package)
#
# function welfaremax(x,grad) # x is the tax vector, grad is the gradient (unspecified here since we have no way of computing it)
#   WW = tax2expectedwelfare(x,PP,rho,eta,nu,Tm,tm,lm,idims)[1]
#   global count
#   count::Int += 1
#   println("  iteration f_$count \n$x\n")
#   return WW
# end
##
# # Choose algorithm (gradient free method) and dimension of tax vector, tm+1 <= n <= Tm
# n = tax_length
# opt = Opt(:LN_BOBYQA, n) # algorithm and dimension of tax vector, possible (derivative-free) algorithms: LN_COBYLA, LN_BOBYQA
#
# #BM-Julia-0.7: squeeze is deprecated in favor of dropdims
# #ub_lm = maximum(squeeze(maximum(pb,2),2),1)[2:lm+1]                                          # Julia 0.6
# ub_lm = maximum(dropdims(maximum(pb,dims=2),dims=2),dims=1)[2:lm+1]                           # Julia 1.x
# #ub_1 = maximum(squeeze(maximum(pb,2),2)[1:idims,:],1)[lm+2:tm+1]                             # Julia 0.6
# ub_1 = maximum(dropdims(maximum(pb,dims=2),dims=2)[1:idims,:],dims=1)[lm+2:tm+1]              # Julia 1.x
# #ub_2 = maximum(squeeze(maximum(pb,2),2)[(idims+1):nsample,:],1)[lm+2:tm+1]                   # Julia 0.6
# ub_2 = maximum(dropdims(maximum(pb,dims=2),dims=2)[(idims+1):nsample,:],dims=1)[lm+2:tm+1]    # Julia 1.x
#
# # lower bound is zero
# lower_bounds!(opt, zeros(n))
# upper_bounds!(opt, [ub_lm; ub_1; ub_2])
#
# # Set maximization
# max_objective!(opt, welfaremax)
#
# # Set relative tolerance for the tax vector choice - vs. ftol for function value?
# ftol_rel!(opt,0.00000000000005)
#
# # Optimize! Initial guess defined above
# (expected_welfare,tax_vector,ret) = optimize(opt, [ub_lm; ub_1; ub_2].*0.5)
#
# # Extract the two tax vectors from the optimization object
# taxes_a = tax_vector[1:tm]
# taxes_b = [tax_vector[1:lm];tax_vector[tm+1:end]]
#
# c, K, T, E, M, mu, lam, D, AD, Y, Q = VarsFromTaxes(taxes_a, taxes_b, PP, nsample, model="$model")
#
# TAX = maximum(PP[1].pb[1:Tm,:],dims=2)[:,1]
# TAX[1] = 0
# TAX[2:(tm+1)] = taxes_a
# taxes_1=TAX
# TAX[2:(tm+1)] = taxes_b
# taxes_2 = TAX
#
# #create Results data type as well as dataframe
# restm = Results(regime_select,nsample,Tm,tm,lm,Regions,taxes_1,taxes_2,expected_welfare,c,K,T,E,M,mu,lam,D,AD,Y,Q,rho,eta,nu,PP)
# dataPtm = FrameFromResults(restm, Tm, nsample, Regions, idims)
#
# Creation CSV output file
# CSV.write("$(pwd())/Outputs/meanOptimum$(model)-Julia11Ptm.csv", dataPtm)
#

###########################################################################################

#SS = Array(Results,3)    >>??     SS = Array{Float64,2}(undef, Results, 3)
#SS = [res reslm restm]
