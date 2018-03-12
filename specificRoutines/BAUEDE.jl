# Preliminaries
using HDF5, JLD, DataFrames, DataArrays, Gadfly, Cairo

folder = pwd()
# Run Function Definitions File
include("$(pwd())/Function_definitions.jl")
include("$(pwd())/createPrandom.jl")

eeMv = [0 1]
regimes = [175 95 165 185] #damage TFP climateSensitivity convergenceRate
names = ["damage" "TFP" "Csensi" "convergenceRate"]
decs = [0.95 0.85 0.75 0.65 0.55 0.5 0.45 0.35 0.25 0.15 0.05] 
decsA = [0.95,0.85,0.75,0.65,0.55,0.5,0.45,0.35,0.25,0.15,0.05] 
Tm=32
tm=Tm
lm = tm
taxes_a = zeros(tm)
taxes_b = taxes_a
nsample = 11
model = "RICE"
rho=0.015
eta = 2
nu = eta
idims = 5
l=2095

output = Array{DataFrames.DataFrame,2}(2,4)
output2=Array{Results,2}(2,4)
for (j, regime) in enumerate(regimes)
for (i, ee) in enumerate(eeMv)
PP = createP(regime, eeM = ee)
expected_welfare = tax2expectedwelfare(taxes_a, PP, rho, eta, nu, Tm, tm, lm, idims; model="NICE")
c, K, T, E, M, mu, lam, D, AD, Y, Q = VarsFromTaxes(taxes_a, taxes_b, PP, nsample)
res = Results(regime,nsample,Tm,tm,lm,Regions,taxes_a,taxes_b,expected_welfare,c,K,T,E,M,mu,lam,D,AD,Y,Q,rho,eta,nu,PP)
dataDist = FrameFromResults(res, Tm, nsample, Regions, idims)
byState = groupby(dataDist,:State) #groups by state
fin = hcat([by(df, :Year, dg -> (dot(dg[:cq1].^(1-eta)+dg[:cq2].^(1-eta)+dg[:cq3].^(1-eta)+dg[:cq4].^(1-eta)+dg[:cq5].^(1-eta),dg[:L]/5)./sum(dg[:L])).^(1/(1-eta)) )[:x1] for df in byState]...)
#yields a Tm*nsample dataArray with the EDEs for each time period
EDE = reshape(repmat(convert(Array, fin),12,1),Tm*nsample*12,1)[:,1]  #creates a Tm*nsample*12 long vector
dataDist[:EDE] = EDE #puts into the dataframe
output[i,j] = dataDist
end
end

saveFolder = "$(pwd())/Outputs/BAU"
for i=1:2
    for j=1:4
df = output[i,j]
df2095 = df[df[:Year].==2095,:]
dfplot = DataFrame(c = [df2095[:cq1];df2095[:cq2];df2095[:cq2];df2095[:cq4];df2095[:cq5]])
one = plot(dfplot, x="c",Geom.histogram)
draw(PDF("$(saveFolder)/$(names[j])E$(eeMv[i]).pdf", 6inch, 4inch), one)
    end
end
i=1; j=1;
df0 = output[1,j]
df1 = output[2,j]
dfEDE0 = df0[(df0[:Year].==2095)&(df0[:Region].=="USA"),:]
dfEDE1 = df1[(df1[:Year].==2095)&(df1[:Region].=="USA"),:]
dfplot = DataFrame(decs = decsA[11:-1:1],EDE0 = dfEDE0[:EDE],EDE1=dfEDE1[:EDE])
layer(x=rand(10), y=rand(10), Geom.point)
one = plot(layer(x=dfplot[:EDE0],y=dfplot[:decs],Geom.line),
 layer(x=dfplot[:EDE1],y=dfplot[:decs],Geom.line))  
 draw(PDF("$(saveFolder)/$(names[j])E$(eeMv[i])EDE.pdf", 6inch, 4inch), one)      