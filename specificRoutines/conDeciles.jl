# Preliminaries
using HDF5, JLD, DataFrames, DataArrays, Gadfly, Cairo

folder = pwd()
# Run Function Definitions File
include("$folder/Function_definitions.jl")

# 95 same as 9, but includes median (nsample = 11)
# 105 same as 10, but includes medial (nsample = 11)
# 165 same as 16, but includes medial (nsample = 11)
# 175  deciles of quadratic damage coefficient
Tm = 32
tm = 18
lm = tm
tax_length = 2*tm - lm
backstop_same = "Y" # "N" is default - choose "Y" if we want all the countries to face the same backstop prices over time

rho = 0.015 # PP[1].para[1] # discount rate
eta = 2 # PP[1].para[3] # inequality aversion/time smoothing parameter (applied to per capita consumption by region, quintile and time period)
nu = 2 # risk aversion parameter (applied over random draws)
eeM = 1 # elasticity of damage with respect to income

regimes = [95 105 165 175]
names = ["TFP" "Decarb" "Csensi" "Damage"]
decs = [0.95 0.85 0.75 0.65 0.55 0.5 0.45 0.35 0.25 0.15 0.05] 

l=2095

model = "NICE"
res = load("$(folder)/Outputs/Optima/meanOptimum$(model).jld", "res")
dataP = readtable("$(folder)/Outputs/Optima/meanOptimum$(model).csv")
tax_vector = res.taxes_1[2:(tm+1)]
include("$folder/createPrandom.jl")

for i=1:4
  regime_select = regimes[i]
  PP = createP(regime_select)
  nsample = size(PP)[1]
  idims = Int(max(round(nsample/2),1))
  expected_welfare = tax2expectedwelfare(tax_vector,PP,rho,eta,nu,Tm,tm,lm,idims,model="$model")[1]

  c, K, T, E, M, mu, lam, D, AD, Y, Q = VarsFromTaxes(tax_vector, tax_vector, PP, nsample, model="$model")

  # create Results structure
  resDist = Results(regime_select,nsample,Tm,tm,lm,Regions,res.taxes_1,res.taxes_2,expected_welfare,c,K,T,E,M,mu,lam,D,AD,Y,Q,rho,eta,nu,PP)
  # create dataframe of period by region by state data
  global dataDist = FrameFromResults(resDist, Tm, nsample, Regions, idims)
  byState = groupby(dataDist,:State) #groups by state
  fin = hcat([by(df, :Year, dg -> (dot(dg[:cq1].^(1-eta)+dg[:cq2].^(1-eta)+dg[:cq3].^(1-eta)+dg[:cq4].^(1-eta)+dg[:cq5].^(1-eta),dg[:L]/5)./sum(dg[:L])).^(1/(1-eta)) )[:x1] for df in byState]...)
  #yields a Tm*nsample dataArray with the EDEs for each time period
  EDE = reshape(repmat(convert(Array, fin),12,1),Tm*nsample*12,1)[:,1]  #creates a Tm*nsample*12 long vector
  dataDist[:EDE] = EDE #puts into the dataframe
  assign1 = parse("global data$(names[i])N = dataDist")
  eval(assign1)
  global xCDF = dataDist[(dataDist[:Region].=="USA")&(dataDist[:Year].==l),:EDE]
  assign2 = parse("global x$(names[i]) = xCDF")
  eval(assign2)
end
#plots
pTFP = plot(x=xTFP, y=decs, Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on TFP growth"))
pDecarb = plot(x=xDecarb, y=(1-decs), Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on decarbonisation rate"))
pCsensi = plot(x=xCsensi, y=(1-decs), Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on climate sensitivity"))
pDamage = plot(x=xDamage, y=(1-decs), Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on damage"))
pCombinedexTFP = plot(layer(x=xDecarb, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "red"))),
                layer(x=xCsensi, y=(1-decs), Geom.line),
                layer(x=xDamage, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "purple"))),
                Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on all"),
                Guide.manual_color_key("", ["TFP", "Decarbonisation", "C sensitivity", "Damage"], ["green", "red", "deepskyblue", "purple"])
            )
pCombined = plot(layer(x=xTFP, y=decs, Geom.line, Theme(default_color=parse(Compose.Colorant, "green"))),
                layer(x=xDecarb, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "red"))),
                layer(x=xCsensi, y=(1-decs), Geom.line),
                layer(x=xDamage, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "purple"))),
                Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on all"),
                Guide.manual_color_key("", ["TFP", "Decarbonisation", "C sensitivity", "Damage"], ["green", "red", "deepskyblue", "purple"])
            )
#save
saveFolder = "$(pwd())/Outputs/CDFs2"
draw(PDF("$(saveFolder)/EDE$(l)TFPRisk$(model).pdf", 6inch, 4inch), pTFP)
draw(PDF("$(saveFolder)/EDE$(l)DecarbRisk$(model).pdf", 6inch, 4inch), pDecarb)
draw(PDF("$(saveFolder)/EDE$(l)CsensiRisk$(model).pdf", 6inch, 4inch), pCsensi)
draw(PDF("$(saveFolder)/EDE$(l)DamageRisk$(model).pdf", 6inch, 4inch), pDamage)
draw(PDF("$(saveFolder)/EDE$(l)CombinedRisk$(model).pdf", 6inch, 4inch), pCombined)
draw(PDF("$(saveFolder)/EDE$(l)CombinedexTFPRisk$(model).pdf", 6inch, 4inch), pCombinedexTFP)

model = "DICE"

#keep the NICE policy, but calculate the DICE (average income) distributions
for i=1:4
  regime_select = regimes[i]
  PP = createP(regime_select)
  nsample = size(PP)[1]
  idims = Int(max(round(nsample/2),1))
  expected_welfare = tax2expectedwelfare(tax_vector,PP,rho,eta,nu,Tm,tm,lm,idims,model="$model")[1]

  c, K, T, E, M, mu, lam, D, AD, Y, Q = VarsFromTaxes(tax_vector, tax_vector, PP, nsample, model="DICE")

  # create Results structure
  resDist = Results(regime_select,nsample,Tm,tm,lm,Regions,res.taxes_1,res.taxes_2,expected_welfare,c,K,T,E,M,mu,lam,D,AD,Y,Q,rho,eta,nu,PP)
  # create dataframe of period by region by state data
  global dataDist = FrameFromResults(resDist, Tm, nsample, Regions, idims)
  byState = groupby(dataDist,:State) #groups by state
  fin = hcat([by(df, :Year, dg -> dot(dg[:c],dg[:L])./sum(dg[:L]))[:x1] for df in byState]...)
  #yields a Tm*nsample dataArray with the EDEs for each time period
  EDE = reshape(repmat(convert(Array, fin),12,1),Tm*nsample*12,1)[:,1]  #creates a Tm*nsample*12 long vector
  dataDist[:EDE] = EDE #puts into the dataframe
  assign1 = parse("global data$(names[i])D = dataDist")
  eval(assign1)
  global xCDF = dataDist[(dataDist[:Region].=="USA")&(dataDist[:Year].==l),:EDE]
  assign2 = parse("global x$(names[i]) = xCDF")
  eval(assign2)
end
#plots
pTFP = plot(x=xTFP, y=decs, Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on TFP growth"))
pDecarb = plot(x=xDecarb, y=(1-decs), Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on decarbonisation rate"))
pCsensi = plot(x=xCsensi, y=(1-decs), Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on climate sensitivity"))
pDamage = plot(x=xDamage, y=(1-decs), Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on damage"))
pCombinedexTFP = plot(layer(x=xDecarb, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "red"))),
                layer(x=xCsensi, y=(1-decs), Geom.line),
                layer(x=xDamage, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "purple"))),
                Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on EDE"),
                Guide.manual_color_key("", ["TFP", "Decarbonisation", "C sensitivity", "Damage"], ["green", "red", "deepskyblue", "purple"])
            )
pCombined = plot(layer(x=xTFP, y=decs, Geom.line, Theme(default_color=parse(Compose.Colorant, "green"))),
                layer(x=xDecarb, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "red"))),
                layer(x=xCsensi, y=(1-decs), Geom.line),
                layer(x=xDamage, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "purple"))),
                Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on EDE"),
                Guide.manual_color_key("", ["TFP", "Decarbonisation", "C sensitivity", "Damage"], ["green", "red", "deepskyblue", "purple"])
            )
#save
saveFolder = "$(pwd())/Outputs/CDFs2"
draw(PDF("$(saveFolder)/EDE$(l)TFPRisk$(model)_NICEpolicy.pdf", 6inch, 4inch), pTFP)
draw(PDF("$(saveFolder)/EDE$(l)DecarbRisk$(model)_NICEpolicy.pdf", 6inch, 4inch), pDecarb)
draw(PDF("$(saveFolder)/EDE$(l)CsensiRisk$(model)_NICEpolicy.pdf", 6inch, 4inch), pCsensi)
draw(PDF("$(saveFolder)/EDE$(l)DamageRisk$(model)_NICEpolicy.pdf", 6inch, 4inch), pDamage)
draw(PDF("$(saveFolder)/EDE$(l)CombinedRisk$(model)_NICEpolicy.pdf", 6inch, 4inch), pCombined)
draw(PDF("$(saveFolder)/EDE$(l)CombinedRiskexTFP$(model)_NICEpolicy.pdf", 6inch, 4inch), pCombinedexTFP)


IITFP =dataTFPD[(dataTFPD[:Region].=="USA")&(dataTFPD[:Year].==l),:EDE]-dataTFPN[(dataTFPN[:Region].=="USA")&(dataTFPN[:Year].==l),:EDE]
IICsensi =dataCsensiD[(dataCsensiD[:Region].=="USA")&(dataCsensiD[:Year].==l),:EDE]-dataCsensiN[(dataCsensiN[:Region].=="USA")&(dataCsensiN[:Year].==l),:EDE]
IIDecarb =dataDecarbD[(dataDecarbD[:Region].=="USA")&(dataDecarbD[:Year].==l),:EDE]-dataDecarbN[(dataDecarbN[:Region].=="USA")&(dataDecarbN[:Year].==l),:EDE]
IIDamage =dataDamageD[(dataDamageD[:Region].=="USA")&(dataDamageD[:Year].==l),:EDE]-dataDamageN[(dataDamageN[:Region].=="USA")&(dataDamageN[:Year].==l),:EDE]

pCombinedII = plot(layer(x=IITFP, y=decs, Geom.line, Theme(default_color=parse(Compose.Colorant, "green"))),
                layer(x=IIDecarb, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "red"))),
                layer(x=IICsensi, y=(1-decs), Geom.line),
                layer(x=IIDamage, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "purple"))),
                Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Absolute Inequality"),
                Guide.manual_color_key("", ["TFP", "Decarbonisation", "C sensitivity", "Damage"], ["green", "red", "deepskyblue", "purple"])
            )

draw(PDF("$(saveFolder)/II$(l)CombinedInequality_NICEpolicy.pdf", 6inch, 4inch), pCombinedII)

RITFP =1-dataTFPD[(dataTFPD[:Region].=="USA")&(dataTFPD[:Year].==l),:EDE].\dataTFPN[(dataTFPN[:Region].=="USA")&(dataTFPN[:Year].==l),:EDE]
RICsensi = 1-dataCsensiD[(dataCsensiD[:Region].=="USA")&(dataCsensiD[:Year].==l),:EDE].\dataCsensiN[(dataCsensiN[:Region].=="USA")&(dataCsensiN[:Year].==l),:EDE]
RIDecarb =1-dataDecarbD[(dataDecarbD[:Region].=="USA")&(dataDecarbD[:Year].==l),:EDE].\dataDecarbN[(dataDecarbN[:Region].=="USA")&(dataDecarbN[:Year].==l),:EDE]
RIDamage =1- dataDamageD[(dataDamageD[:Region].=="USA")&(dataDamageD[:Year].==l),:EDE].\dataDamageN[(dataDamageN[:Region].=="USA")&(dataDamageN[:Year].==l),:EDE]

pCombinedRI = plot(layer(x=RITFP, y=decs, Geom.line, Theme(default_color=parse(Compose.Colorant, "green"))),
                layer(x=RIDecarb, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "red"))),
                layer(x=RICsensi, y=(1-decs), Geom.line),
                layer(x=RIDamage, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "purple"))),
                Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Relative Inequality"),
                Guide.manual_color_key("", ["TFP", "Decarbonisation", "C sensitivity", "Damage"], ["green", "red", "deepskyblue", "purple"])
            )

draw(PDF("$(saveFolder)/RI$(l)CombinedInequality_NICEpolicy.pdf", 6inch, 4inch), pCombinedRI)

#now load the DICE optimum to use the DICE policy with the DICE (average income) distributions
res = load("$(folder)/Outputs/Optima/meanOptimum$(model).jld", "res")
dataP = readtable("$(folder)/Outputs/Optima/meanOptimum$(model).csv")
tax_vector = res.taxes_1[2:(tm+1)]

for i=1:4
  regime_select = regimes[i]
  PP = createP(regime_select)
  nsample = size(PP)[1]
  idims = Int(max(round(nsample/2),1))
  expected_welfare = tax2expectedwelfare(tax_vector,PP,rho,eta,nu,Tm,tm,lm,idims,model="$model")[1]

  c, K, T, E, M, mu, lam, D, AD, Y, Q = VarsFromTaxes(tax_vector, tax_vector, PP, nsample, model="DICE")

  # create Results structure
  resDist = Results(regime_select,nsample,Tm,tm,lm,Regions,res.taxes_1,res.taxes_2,expected_welfare,c,K,T,E,M,mu,lam,D,AD,Y,Q,rho,eta,nu,PP)
  # create dataframe of period by region by state data
  global dataDist = FrameFromResults(resDist, Tm, nsample, Regions, idims)
  byState = groupby(dataDist,:State) #groups by state
  fin = hcat([by(df, :Year, dg -> dot(dg[:c],dg[:L])./sum(dg[:L]))[:x1] for df in byState]...)
  #yields a Tm*nsample dataArray with the EDEs for each time period
  EDE = reshape(repmat(convert(Array, fin),12,1),Tm*nsample*12,1)[:,1]  #creates a Tm*nsample*12 long vector
  dataDist[:EDE] = EDE #puts into the dataframe
  assign1 = parse("global data$(names[i]) = dataDist")
  eval(assign1)
  global xCDF = dataDist[(dataDist[:Region].=="USA")&(dataDist[:Year].==l),:EDE]
  assign2 = parse("global x$(names[i]) = xCDF")
  eval(assign2)
end
#plots
pTFP = plot(x=xTFP, y=decs, Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on TFP growth"))
pDecarb = plot(x=xDecarb, y=(1-decs), Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on decarbonisation rate"))
pCsensi = plot(x=xCsensi, y=(1-decs), Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on climate sensitivity"))
pDamage = plot(x=xDamage, y=(1-decs), Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on damage"))
pCombinedexTFP = plot(layer(x=xDecarb, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "red"))),
                layer(x=xCsensi, y=(1-decs), Geom.line),
                layer(x=xDamage, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "purple"))),
                Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on all"),
                Guide.manual_color_key("", ["TFP", "Decarbonisation", "C sensitivity", "Damage"], ["green", "red", "deepskyblue", "purple"])
            )
pCombined = plot(layer(x=xTFP, y=decs, Geom.line, Theme(default_color=parse(Compose.Colorant, "green"))),
                layer(x=xDecarb, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "red"))),
                layer(x=xCsensi, y=(1-decs), Geom.line),
                layer(x=xDamage, y=(1-decs), Geom.line, Theme(default_color=parse(Compose.Colorant, "purple"))),
                Guide.xlabel("EDE consumption in $l"), Guide.ylabel("CDF"), Guide.title("Risk on all"),
                Guide.manual_color_key("", ["TFP", "Decarbonisation", "C sensitivity", "Damage"], ["green", "red", "deepskyblue", "purple"])
            )
#save
saveFolder = "$(pwd())/Outputs/CDFs2"
draw(PDF("$(saveFolder)/EDE$(l)TFPRisk$(model).pdf", 6inch, 4inch), pTFP)
draw(PDF("$(saveFolder)/EDE$(l)DecarbRisk$(model).pdf", 6inch, 4inch), pDecarb)
draw(PDF("$(saveFolder)/EDE$(l)CsensiRisk$(model).pdf", 6inch, 4inch), pCsensi)
draw(PDF("$(saveFolder)/EDE$(l)DamageRisk$(model).pdf", 6inch, 4inch), pDamage)
draw(PDF("$(saveFolder)/EDE$(l)CombinedRisk$(model).pdf", 6inch, 4inch), pCombined)
draw(PDF("$(saveFolder)/EDE$(l)CombinedexTFPRisk$(model).pdf", 6inch, 4inch), pCombinedexTFP)