using HDF5, JLD, DataFrames, DataArrays, Gadfly, Distributions
folder = pwd()
folderLearning = "$(pwd())/Outputs/valueOfLearning"
include("$(pwd())/Function_definitions.jl")
include("$(pwd())/createPrandom.jl")
include("$(pwd())/optimisations.jl")
resArray = load("$(pwd())/Outputs/valueOfLearning/tenbranches/resultsB", "resArray")
eev = [-1 0 1];
regimes = [19 9 16 18]
lmv = [0 2045];
f(A,names) = DataFrame(Any[A[:,i] for i = 1:size(A,2)], map(Symbol,names))
nemas = ["$i" for i=1:10];
year = DataFrame(year=10*(0:31)+2005);
eWelfare = Array{Float64}(2,4,2,8);
for i=2:2
    for j=1:4
        for k = 1:2
            for l=1:8
                eWelfare[i-1,j,k,l] = resArray[i,j,k,l].EWelfare
            end
            l = indmax(eWelfare[i-1,j,k,:])
            df = [year f(resArray[i,j,k,l].taxes,nemas)]
            writetable("$(pwd())/Outputs/valueOfLearning/tenbranches/feb2018/tenTaxes_e_$(eev[i])_D_$(regimes[j])_l_$(lmv[k]).csv", df)
        end
    end
end

