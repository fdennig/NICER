using HDF5, JLD, DataFrames, DataArrays, StatsBase, NLopt, Gadfly, Distributions
f(A,names) = DataFrame(Any[A[:,i] for i = 1:size(A,2)], map(Symbol,names))
saveFolder = "$(pwd())/Outputs/March2018" ############################################ SET THIS TO YOUR SAVING DESTINATION
folderDat = "$(pwd())/Outputs/valueOfLearning/tenbranches" ###### set this location of the directory that contains resultseeM1.jld
folder = pwd()
include("$(pwd())/Function_definitions.jl")
include("$(pwd())/createPrandom.jl")
include("$(pwd())/optimisations.jl")
results = load("$folderDat/resultseeM1.jld")
resA = results["resArray"]# resA = resAL["resArray"]
eeMv = [0 1]
regimes = [19 9 161 18]
names = ["damage" "TFP" "climateSensitivity" "convergenceRate"]
Tm=32
tm=22
lmv = [0 3 tm]
nsample = 10
resArray = Array{Results10}(length(eeMv)+1,length(regimes),length(lmv))


for (k,lm) in enumerate(lmv) #loop over learning dates
for (j,regime) in enumerate(regimes) #loop over regimes
  #run the DICE model
  i=3
  PP = createP(regime)
  pred1 = resA[3,j,k,1].taxes[2:lm+1,1]*0.95
  pred2 = ones(tm-lm,nsample)
  for jj = 1:nsample
    pred2[:,jj] = resA[3,j,k,1].taxes[lm+2:tm+1,jj]*0.9+0.01
  end
  pre = [pred1; pred2[:]]
    resArray[i,j,k] = optimiseNICER10(PP,regime,lm=lm,tm=tm,inite=pre,model="DICE")
    save("$(saveFolder)resultsX.jld","resArray",resArray)
  #run NICE for e=0, and e=1
  for (i,eeM) in enumerate(eeMv)
    PP = createP(regime, eeM = eeM)
    pred1 = resA[i+1,j,k,1].taxes[2:lm+1,1]*0.95
    pred2 = ones(tm-lm,nsample)
    for jj = 1:nsample
      pred2[:,jj] = resA[i+1,j,k,1].taxes[lm+2:tm+1,jj]*0.9+0.01
    end
    pre = [pred1; pred2[:]]
    resArray[i,j,k] = optimiseNICER10(PP,regime,lm=lm,tm=tm,inite=pre,model="NICE")
    save("$(saveFolder)resultsX.jld","resArray",resArray)
  end
end
end