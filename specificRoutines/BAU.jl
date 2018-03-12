###########################################################################################################################################################################################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

using HDF5, JLD
# Set user name to run auxiliary files
user = "francis" # or "francis" or "marc" as the case may be

regime_select = 0
idims = 1
Tm = 32
tm = 18 #Tm-1
lm=tm
backstop_same = "Y" # "N"
model = "RICE"
rho=0.015
eta = 2
# Now execute the whole code: select all and evaluate
###########################################################################################################################################################################################################
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

folder = pwd()
include("$folder/Function_definitions.jl")
include("$folder/createPrandom.jl")



# Extract the two tax vectors from the optimization object
taxes_a = zeros(tm)
taxes_b = taxes_a
expected_welfare = tax2welfare(taxes_1, PP[1], rho, eta, Tm, model="RICE")

c, K, T, E, M, mu, lam, D, AD, Y, Q = VarsFromTaxes(taxes_a, taxes_b, PP, nsample)

TAX = maximum(PP[1].pb[1:Tm,:],2)[:,1]
TAX[1] = 0
TAX[2:(tm+1)] = taxes_a
taxes_1=TAX
TAX[2:(tm+1)] = taxes_b
taxes_2 = TAX

# create .jld of S_lm
res = Results(regime_select,nsample,Tm,tm,lm,Regions,taxes_1,taxes_2,expected_welfare,c,K,T,E,M,mu,lam,D,AD,Y,Q,rho,eta,nu,PP)
filenm = string(regime_select)
using DataFrames
# create dataframe of period by region by state data
dataP = FrameFromResults(res, Tm, nsample, Regions, idims)

jldopen("$(pwd())/Outputs/BAU/BAUparameterstm15.jld", "w") do file
    write(file, "PP", PP[1])
end
jldopen("$(pwd())/Outputs/BAU/BAUtm15.jld", "w") do file
    write(file, "BAU", res)
end

writetable("$(pwd())/Outputs/BAU/BAUtm15.csv", dataP)
