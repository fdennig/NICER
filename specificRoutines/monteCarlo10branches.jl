using HDF5, JLD, DataFrames, DataArrays, StatsBase, NLopt, Roots, Gadfly, Distributions
f(A,names) = DataFrame(Any[A[:,i] for i = 1:size(A,2)], map(Symbol,names))
folder = pwd()
folderLearning = "$(pwd())/Outputs/valueOfLearning"
include("$(pwd())/Function_definitions.jl")
include("$(pwd())/createPrandom.jl")
include("$(pwd())/optimisations.jl")
results = load("$folderLearning/tenbranches/resultseeM1.jld")# resAL = load("$folderLearning/resArray10branchesmarch15.jld")
resA = results["resArray"]# resA = resAL["resArray"]
# initial = load("$folderLearning/tenBranches/initials.jld")
# initials = initial["initials"]
# WELF = readtable("$folderLearning/tenBranches/welfares10branches.csv")
eeMv = [-1 0 1]
regimes = [19 9 16 18]
Tm=32
tm=22
lmv = [0 3 tm]
runs = 8 #15 #2
ND = Normal()
            # resArray = Array{Results10}(length(eeMv),length(regimes),length(lmv),runs)
            # initials = Array{Vector{Float64}}(length(eeMv),length(regimes),length(lmv),runs)
            # WELF = DataFrame(Elasticity = Int[], Regime = Int[], Learnperiod = Int[], Welfare = Float64[])
for j=1:8 #1:runs
for g=2:2 #2:3 #(g,eeM) in enumerate(eeMv)
  eeM=eeMv[g]
for (i,lm) in enumerate(lmv)
  lm = lmv[i]
for (h,regime) in enumerate(regimes)
  # regime = regimes[h]
PP = createP(regime, eeM = eeM)
nsample = length(PP)
  if j==1
    pred1 = resA[1,h,i,1].taxes[2:lm+1,1]*0.9
    pred2 = zeros(tm-lm,nsample)
    for jj = 1:nsample
      pred2[:,jj] = resA[1,h,i,1].taxes[lm+2:tm+1,jj]*0.9+0.01
    end
  else
    pred1 = resArray[g,h,i,j-1].taxes[2:lm+1,1]*0.9+rand(ND,lm)
    pred2 = zeros(tm-lm,nsample)
    for jj = 1:nsample
      pred2[:,jj] = (resArray[g,h,i,j-1].taxes[lm+2:tm+1,jj]+0.01)*0.99 + max(0,rand(ND,tm-lm))
    end
  end
  pre = [pred1; pred2[:]]
  resArray[g,h,i,j] = optimiseNICER10(PP,regime,lm=lm,tm=tm,inite=pre)
  lp = 2015+lm*10
  WELFA = DataFrame(Elasticity = eeM, Regime = regime, Learnperiod = lp, Welfare = resArray[g,h,i,j].EWelfare)
  WELF = append!(WELF,WELFA)
  writetable("$(pwd())/Outputs/valueOfLearning/tenbranches/welfares10branches.csv",WELF)
  initials[g,h,i,j] = pre
  using FileIO; save(File(format"JLD","$(pwd())/Outputs/valueOfLearning/tenbranches/initialsB"),"initials",initials)
  using FileIO; save(File(format"JLD","$(pwd())/Outputs/valueOfLearning/tenbranches/resultsB"),"resArray",resArray)
  # save("$(pwd())/Outputs/valueOfLearning/tenbranches/results.jld","resArray",resArray)
  println("e=$(eeM)lp=$(lp)regime=$(regime)run=$j")
end
end
end
end


resArray = load("$(pwd())/Outputs/valueOfLearning/tenbranches/resultsB", "resArray")
legend=["x", "y1", "y2", "y3", "y4", "y5", "y6", "y7", "y8", "y9", "y10"]
tm = 22
p = Array{Any}(4,3)
for i=1:4
  for j=1:3
    taxs = resArray[2,i,j,1].taxes[1:tm+1,:]
    df = f([convert(Array,10*(0:tm)+2005) taxs], legend)
    p[i,j] = plot(melt(df, :x), x = :x, y = :value, color = :variable, Geom.line, Geom.point,Guide.title("Regime $(regimes[i])"))
  end
end
draw(PDF("$(pwd())/Outputs/valueOfLearning/workinprogress/tenBranches5.pdf", 20inch, 12inch),
        vstack(hstack(p[1,1],p[2,1],p[3,1],p[4,1]),hstack(p[1,2],p[2,2],p[3,2],p[4,2]),hstack(p[1,3],p[2,3],p[3,3],p[4,3])))

JLD.save("$(pwd())/Outputs/valueOfLearning/resArray10branchesmarch15.jld", "resArray", resArray)
