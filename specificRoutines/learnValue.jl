###########################################################################################################################################################################################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Preliminaries
using HDF5, JLD, DataFrames, DataArrays, StatsBase, NLopt, Roots

folder = pwd()
include("$folder/Function_definitions.jl")
include("$folder/createPrandom.jl")
#Regimes
  # 9 = Deciles of TFP (all else fixed at means) (nsample = 10)
  # 10 = Deciles of decarbonization rates (all else fixed at means) (nsample = 10)
  # 16 = DECILES of climate sensitivity
  # 17 = uncertainty (100 random draws) on the quadratic damage coefficient (full correlation)
# Define exogenous parameters
Tm = 32
tm = 31
rho = 0.015 # discount rate
eta = 2 # inequality aversion/time smoothing parameter (applied to per capita consumption by region, quintile and time period)
nu = 2 # risk aversion parameter (applied over random draws)
eeM = 1 # elasticity of damage with respect to income
models = "NICE"
regimes = [9, 10, 16, 17, 18, 19]
lmv = [0 3 tm]
exi = [-1 0 1]
#compute optima
# pre allocate the dataframes and results array
consN = DataFrame(State = Int[], Region = String[], Year = Int[], tax = Float64[], cq1 = Float64[], cq2 = Float64[], cq3 = Float64[], cq4 = Float64[], cq5 = Float64[], K = Float64[], E = Float64[], mu = Float64[], lam = Float64[], D = Float64[], Y = Float64[], L = Float64[], Regime = Float64[], LearnP = Float64[])
consD = DataFrame(State = Int[], Region = String[], Year = Int[], tax = Float64[], c = Float64[], K = Float64[], E = Float64[], mu = Float64[], lam = Float64[], D = Float64[], Y = Float64[], L = Float64[], Regime = Float64[], LearnP = Float64[])
resArray = Array{Results}(3,6,3)
for k=1:3 # loop over value of XI, NO over models NICE and DICE
  for j = 1:6 # loop over regimes 9, 10, 16, 17
  # optimisation preliminaries
    PP = createP(regimes[j];eeM=exi[k]) # createP(regime_select; backstop_same = "Y", gy0M = dparam_i["gy0"][2]', gy0sd = ones(12)*0.0059, sighisTM = dparam_i["sighisT"][2]', sighisTsd = ones(12)*0.0064, xi1M = 3, xi1sd = 1.4, eeM = 1, eesd = 0.3)
    nsample=length(PP)
    pb = zeros(nsample,12,60)
    for ii=1:nsample
      pb[ii,:,:] = PP[ii].pb'
    end
    idims = Int(max(round(nsample/2),1)) # bifurcates the random draws into two subsets
    for i = 1:3
      lm=lmv[i] # loop over learning period
      tax_length = 2*tm - lm
      global count = 0 # keep track of number of iterations
  # Define function to be maximized (requires special format for NLopt package)
      function welfaremax(x,grad) # x is the tax vector, grad is the gradient (unspecified here since we have no way of computing it)
        WW = tax2expectedwelfare(x,PP,rho,eta,nu,Tm,tm,lm,idims,model="$(models)")[1] #change to model="RICE" or "DICE" for different levels of aggregation
        count += 1
        if count%100==0
          println("f_$count($x)")
        end
        return WW
      end
      n = tax_length
      opt = Opt(:LN_BOBYQA, n) # algorithm and dimension of tax vector, possible (derivative-free) algorithms: LN_COBYLA, LN_BOBYQA
      ub_lm = maximum(squeeze(maximum(pb,2),2),1)[2:lm+1]
      ub_1 = maximum(squeeze(maximum(pb,2),2)[1:idims,:],1)[lm+2:tm+1]
      ub_2 = maximum(squeeze(maximum(pb,2),2)[(idims+1):nsample,:],1)[lm+2:tm+1]
      # lower bound is zero
      lower_bounds!(opt, zeros(n))
      upper_bounds!(opt, [ub_lm; ub_1; ub_2])
      # Set maximization
      max_objective!(opt, welfaremax)
      # Set relative tolerance for the tax vector choice - vs. ftol for function value?
      ftol_rel!(opt,0.000000000005)
      # Optimize! Initial guess defined above
      (expected_welfare,tax_vector,ret) = optimize(opt, [ub_lm; ub_1; ub_2].*0.5)

      # Extract the two tax vectors from the optimization object
      taxes_a = tax_vector[1:tm]
      taxes_b = [tax_vector[1:lm];tax_vector[tm+1:end]]
      # get all endogenous variables
      c, K, T, E, M, mu, lam, D, AD, Y, Q = VarsFromTaxes(taxes_a, taxes_b, PP, nsample, model="$(models)")
      # get taxes for both learning branches
      taxes_1 = maximum(PP[1].pb[1:Tm,:],2)[:,1]
      taxes_1[1] = 0
      taxes_1[2:(tm+1)] = taxes_a
      taxes_2 = maximum(PP[1].pb[1:Tm,:],2)[:,1]
      taxes_2[1] = 0
      taxes_2[2:(tm+1)] = taxes_b

      # create Results variable
      res = Results(regimes[j],nsample,Tm,tm,lm,Regions,taxes_1,taxes_2,expected_welfare,c,K,T,E,M,mu,lam,D,AD,Y,Q,rho,eta,nu,PP)
      resArray[k,j,i] = res
      # create dataframe of period by region by state data
      dataP = FrameFromResults(res, Tm, nsample, Regions, idims)
      dataP[:Regime] = ones(size(dataP)[1])*regimes[j]
      dataP[:LearnP] = ones(size(dataP)[1])*(2015+10*lm)
      delete!(dataP, [:ID, :A, :sigma, :th1, :pb, :EL])
      if models == "NICE"
        consN = append!(consN, dataP)
      elseif models == "DICE"
        consD = append!(consD, dataP)
      end
    end
  end
end

learnN = zeros(6,3,2)
for j=1:6
  for k=1:3
    for i=1:2
      x = resArray[k,j,i].EWelfare - resArray[k,j,i+1].EWelfare
      learnN[j,k,i] = ((1+(1-eta)*x/sum(repmat(resArray[1,j,k].PP[1].L[1,:]/5,1,5).*resArray[1,j,k].c[1,:,:,1].^(1-eta)))^(1/(1-eta))-1)*100
    end
  end
end

learn = DataFrame(Model = repmat([-1; 0; 1]',6,1)[:], Regime = repmat(regimes,3), l2015lto2045=learnN[:,:,1][:], l2045tolend =learnN[:,:,2][:])
writetable("$(pwd())/Outputs/valueOfLearning/valueOfLearningNICE2.csv",learn)
save("$(pwd())/Outputs/valueOfLearning/resArrayNICE2.jld", "resArray", resArray)
# learnD = zeros(5,2)
# for j=1:5
#   for k=1:2
#     x = resArray[2,j,k].EWelfare - resArray[2,j,k+1].EWelfare
#     learnD[j,k] = (1+(1-eta)*x/sum(repmat(resArray[1,j,k].PP[1].L[1,:]/5,1,5).*resArray[1,j,k].c[1,:,:,1].^(1-eta)))^(1/(1-eta))-1
#   end
# end
# learn = DataFrame(Model = repmat(["NICE"; "DICE"]',5,1)[:], Regime = repmat(regimes,2), l2015l2025=[learnN[:,1]; learnD[:,1]], l2024lend = l2015l2025=[learnN[:,2]; learnD[:,2]])
# save("$(pwd())/Outputs/valueOfLearning/resArray.jld", "resArray", resArray)
# writetable("$(pwd())/Outputs/valueOfLearning/consD.csv", consD)
# writetable("$(pwd())/Outputs/valueOfLearning/consN.csv", consN)
# writetable("$(pwd())/Outputs/valueOfLearning/valueOfLearning.csv",learn)
