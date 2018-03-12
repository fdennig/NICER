function optimiseNICER(PP, regime; rho=0.015, eta=2, nu=2, Tm=32, tm=31, lm=0, model="NICE", Regions=["USA", "OECD Europe", "Japan", "Russia", "Non-Russia Eurasia", "China", "India", "Middle East", "Africa", "Latin America", "OHI", "Other non-OECD Asia"]')
  nsample=length(PP)
  pb = zeros(nsample,12,60)
  for ii=1:nsample
    pb[ii,:,:] = PP[ii].pb'
  end
  idims = Int(max(round(nsample/2),1))
  n = 2*tm - lm
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
  res = Results(regime,nsample,Tm,tm,lm,Regions,taxes_1,taxes_2,expected_welfare,c,K,T,E,M,mu,lam,D,AD,Y,Q,rho,eta,nu,PP)
end

function optimiseNICER10(PP, regime; rho=0.015, eta=2, nu=2, Tm=32, tm=31, lm=0, model="NICE", inite = 0, Regions=["USA", "OECD Europe", "Japan", "Russia", "Non-Russia Eurasia", "China", "India", "Middle East", "Africa", "Latin America", "OHI", "Other non-OECD Asia"]')
  nsample=length(PP)
  #length of argument
  n = nsample*tm - (nsample-1)*lm
  global count = 0 # keep track of number of iterations
  # Define function to be maximized (requires special format for NLopt package)
  function welfaremax(x,grad) # x is the tax vector, grad is the gradient (unspecified here since we have no way of computing it)
    global taxes = zeros(tm,nsample) #produces the nsample different tax paths from x
    if lm<tm
      cuts = tm:(tm-lm):n
      taxes[:,1] = x[1:tm]
      for jj=1:(nsample-1)
        taxes[:,jj+1] = [x[1:lm];x[cuts[jj]+1:cuts[jj+1]]]
      end
    else
      taxes = repmat(x,1,nsample)
    end
    WW = tax2expectedwelfare10(taxes,PP,rho,eta,nu,Tm,tm,lm,model="$(model)")[1] #change to model="RICE" or "DICE" for different levels of aggregation
    count += 1
    if count%100==0
      println("f_$count($taxes)")
    end
    return WW
  end
  opt = Opt(:LN_BOBYQA, n) # algorithm and dimension of tax vector, possible (derivative-free) algorithms: LN_COBYLA, LN_BOBYQA
  pb = zeros(nsample,12,60)
  for ii=1:nsample
    pb[ii,:,:] = PP[ii].pb'
  end
  ub_lm = maximum(squeeze(maximum(pb,2),2),1)[2:lm+1]
  ub_n = zeros(tm-lm,nsample)
  for jj = 1:nsample
    ub_n[:,jj] = squeeze(maximum(pb,2),2)[1,lm+2:tm+1]
  end
  # lower bound is zero
  lower_bounds!(opt, zeros(n))
  upper_bounds!(opt, [ub_lm;ub_n[:]])
  # Set maximization
  max_objective!(opt, welfaremax)
  # Set relative tolerance for the tax vector choice - vs. ftol for function value?
  ftol_rel!(opt,0.00000000000000005)
  xtol_rel!(opt,1e-13)
  #Initial guess: RANDOMIZED
  # init = Array{Float64}(n)
  # for k=1:n
  #   init[k]=rand(Uniform(0,[ub_lm;ub_n[:]][k]),1)[1]
  # end
  if inite == 0
    init = [ub_lm;ub_n[:]]*0.8
  else
    init = inite
  end
  print(typeof(init))
  # Optimize!
  (expected_welfare,tax_vector,ret) = optimize(opt, init)

  # get all endogenous variables
  c, K, T, E, M, mu, lam, D, AD, Y, Q = VarsFromTaxes10(taxes, PP, nsample, model="$(model)", Tm=Tm)
  # get taxes for both learning branches
  TAXES = zeros(Tm, nsample)
  for jj = 1:nsample
    TAXES[2:Tm,jj] = maximum(PP[1].pb[2:Tm,:],2)[:,1]
    TAXES[2:(tm+1),jj] = taxes[:,jj]
  end
  # create Results variable
  res = Results10(regime,nsample,Tm,tm,lm,Regions,TAXES,expected_welfare,c,K,T,E,M,mu,lam,D,AD,Y,Q,rho,eta,nu,PP,ret)
end
