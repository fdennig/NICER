###########################################################################################################################################################################################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Preliminaries
using HDF5, JLD, DataFrames, DataArrays, StatsBase, NLopt, Roots

folder = pwd()
include("$folder/Function_definitions.jl")
include("$folder/createPrandom.jl")
#Regimes
  #20: Normally distributed damage elasticity with mean=eeM and sd=eesd
Tm = 32
tm = 31
rho = 0.015 # discount rate
eta = 2 # inequality aversion/time smoothing parameter (applied to per capita consumption by region, quintile and time period)
nu = 2 # risk aversion parameter (applied over random draws)
eeM = 1 # elasticity of damage with respect to income
model="NICE"
regime = 20
lmv = [0 3 tm]
eeMv = [0.3 0.4 0.5]
eesdv = [0.2 0.3 0.4]
#compute optima
# pre allocate the dataframes and results array
consN = DataFrame(State = Int[], Region = String[], Year = Int[], tax = Float64[], cq1 = Float64[], cq2 = Float64[], cq3 = Float64[], cq4 = Float64[], cq5 = Float64[], K = Float64[], E = Float64[], mu = Float64[], lam = Float64[], D = Float64[], Y = Float64[], L = Float64[], Regime = Float64[], LearnP = Float64[])
resArray = Array{Results}(3,3,3)

# optimisation preliminaries
for (j,eeM) in enumerate(eeMv)
  for (k,eesd) in enumerate(eesdv)
    PP = createP(regime;eeM=eeM, eesd=eesd) # createP(regime_select; backstop_same = "Y", gy0M = dparam_i["gy0"][2]', gy0sd = ones(12)*0.0059, sighisTM = dparam_i["sighisT"][2]', sighisTsd = ones(12)*0.0064, xi1M = 3, xi1sd = 1.4, eeM = 1, eesd = 0.3)
    nsample=length(PP)
    pb = zeros(nsample,12,60)
    for ii=1:nsample
      pb[ii,:,:] = PP[ii].pb'
    end
    idims = Int(max(round(nsample/2),1)) # bifurcates the random draws into two subsets
    for (i, lm) in enumerate(lmv)
      tax_length = 2*tm - lm
      global count = 0 # keep track of number of iterations
    # Define function to be maximized (requires special format for NLopt package)
      function welfaremax(x,grad) # x is the tax vector, grad is the gradient (unspecified here since we have no way of computing it)
        WW = tax2expectedwelfare(x,PP,rho,eta,nu,Tm,tm,lm,idims,model="$(model)")[1] #change to model="RICE" or "DICE" for different levels of aggregation
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
      c, K, T, E, M, mu, lam, D, AD, Y, Q = VarsFromTaxes(taxes_a, taxes_b, PP, nsample, model="$(model)")
      # get taxes for both learning branches
      taxes_1 = maximum(PP[1].pb[1:Tm,:],2)[:,1]
      taxes_1[1] = 0
      taxes_1[2:(tm+1)] = taxes_a
      taxes_2 = maximum(PP[1].pb[1:Tm,:],2)[:,1]
      taxes_2[1] = 0
      taxes_2[2:(tm+1)] = taxes_b

      # create Results variable
      res = Results(regime,nsample,Tm,tm,lm,Regions,taxes_1,taxes_2,expected_welfare,c,K,T,E,M,mu,lam,D,AD,Y,Q,rho,eta,nu,PP)
      resArray[i,j,k] = res
      # create dataframe of period by region by state data
      dataP = FrameFromResults(res, Tm, nsample, Regions, idims)
      dataP[:Regime] = ones(size(dataP)[1])*regime
      dataP[:LearnP] = ones(size(dataP)[1])*(2015+10*lm)
      delete!(dataP, [:ID, :A, :sigma, :th1, :pb, :EL])
      consN = append!(consN, dataP)
    end
  end
end


learnN = zeros(2,3,3)
for i=1:2 #(i, lm) in enumerate(lmv)
  for (j, eeM) in enumerate(eeMv)
    for (k, eesd) in enumerate(eesdv)
      x = resArray[i,j,k].EWelfare - resArray[i+1,j,k].EWelfare
      learnN[i,j,k] = ((1+(1-eta)*x/sum(repmat(resArray[i,j,k].PP[1].L[1,:]/5,1,5).*resArray[i,j,k].c[1,:,:,1].^(1-eta)))^(1/(1-eta))-1)*100
    end
  end
end

learn = DataFrame(Mean = repmat(eeMv',3)[:], SD = repmat(eesdv,3)[:], l2015lto2045=learnN[1,:,:][:], l2045tolend =learnN[2,:,:][:])
writetable("$(pwd())/Outputs/valueOfLearning/valueOfLearningElasticity.csv",learn)
# save("$(pwd())/Outputs/valueOfLearning/resArrayNICE2.jld", "resArray", resArray)
# learnD = zeros(5,2)
# for j=1:5
#   for k=1:2
#     x = resArray[2,j,k].EWelfare - resArray[2,j,k+1].EWelfare
#     learnD[j,k] = (1+(1-eta)*x/sum(repmat(resArray[1,j,k].PP[1].L[1,:]/5,1,5).*resArray[1,j,k].c[1,:,:,1].^(1-eta)))^(1/(1-eta))-1
#   end
# end
# save("$(pwd())/Outputs/valueOfLearning/resArray.jld", "resArray", resArray)
# writetable("$(pwd())/Outputs/valueOfLearning/consD.csv", consD)
# writetable("$(pwd())/Outputs/valueOfLearning/consN.csv", consN)
# writetable("$(pwd())/Outputs/valueOfLearning/valueOfLearning.csv",learn)
