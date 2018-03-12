using DataFrames, DataArrays, HDF5, JLD

df = readtable("$(pwd())/Outputs/valueOfLearning/tenbranches/welfarestenbranches.csv")  #"$(pwd())/Outputs/valueOfLearning/welfaresxi10m1.csv"
rNames = ["TFP", "CS", "CR", "Damage"]
rCodes = [9, 16, 18, 19]
dfd=DataFrame(Elasticity = Int[], Regime = String[], WLearn2015 = Float64[], WLearn2045 = Float64[], WLearn2325 = Float64[])
for i=[1,0,-1]
  for j=1:4
    x = maximum(df[(df[:Elasticity].==i)&(df[:Regime].==rCodes[j])&(df[:Learnperiod].==2015),:Welfare])
    y = maximum(df[(df[:Elasticity].==i)&(df[:Regime].==rCodes[j])&(df[:Learnperiod].==2045),:Welfare])
    z = maximum(df[(df[:Elasticity].==i)&(df[:Regime].==rCodes[j])&(df[:Learnperiod].==2325),:Welfare])
    push!(dfd, [i, rNames[j], x, y, z])
  end
end
function learnValue(dW, res, eta)
vl = ((1+(1-eta)*dW/sum(repmat(res.PP[1].L[1,:]/5,1,5).*res.c[1,:,:,1].^(1-eta))).^(1/(1-eta))-1)*100
end
folder = pwd()
include("$folder/Function_definitions.jl")
resA = load("$(pwd())/Outputs/valueOfLearning/tenbranches/results.jld","results") #load("$(pwd())/Outputs/valueOfLearning/resArrayNICE.jld","resArray")
res = resA["resArray"][1,1,1,1]
dfd[:Value2015to2045] = learnValue(dfd[:WLearn2015] - dfd[:WLearn2045], res,2)
dfd[:Value2045to2325] = learnValue(dfd[:WLearn2045] - dfd[:WLearn2325], res,2)
output = DataFrame(Elasticity=dfd[:Elasticity], Uncertainty = dfd[:Regime], V2015to2045=dfd[:Value2015to2045], V2045to2325=dfd[:Value2045to2325])

writetable("$(pwd())/Outputs/valueOfLearning/tenbranchesVoLtable.csv",output)

using YTables
println(latex(output))
