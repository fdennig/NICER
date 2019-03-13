##############################################################################################################################################################
# NICER computational model used for risk, inequality and climate change
##############################################################################################################################################################
# File: Function_definitions.jl
# Content: part of the NICER integrated assessment model
##############################################################################################################################################################
# Source: original program running in Julia 0.6.4 can be found at https://github.com/fdennig/NICER
# Changes introduced: all changes required to run with Julia 1.1.0
# File release: March 13th 2019 (BM)
##############################################################################################################################################################
# Program required by:
#  ../Optimization.jl
#  ../createPrandom.jl
##############################################################################################################################################################


function backstop(Th,RL,pw,du,dd,tau,nsample)
  #creates the nsample,12,60 array of the backstop price
  #inputs:
  #   Th (scalar): initial to final backstop price ratio
  #   RL Array(1,12):Region to world backstop price ratio
  #   pw (scalar): price of world backstop
  #   du (scalar): rate of decline of backstop price before tau
  #   dd (scalar): rate of decline of backstop price after tau
  #   tau (scalar): period at which rate of decline changes
  T = 60
  I = 12
  taut = convert(Int,(tau-1995)/10)
  p0 = pw.*RL
  pb = zeros(nsample,I,T) # note order of dimensions!
  pb[:,:,1] = p0
  for t = 2:taut
    pb[:,:,t] = Th*p0 + (1-du)*(pb[:,:,t-1]-Th*p0)
  end
  for t = (taut+1):T
    pb[:,:,t] = pb[:,:,t-1]*dd
  end
  return pb
end

function sig(gT,delsig,sighisT,adj15,Y0,E0,nsample)
  #creates the nsample,12,60 array of (unmitigated/BAU) emissions to output ratio
  #inputs:
  #   gT
  T = 60
  I = 12
  sigma = zeros(nsample,I,T) # note order of dimensions
  E000 = E0./1000
  sigma[:,:,1] = repeat(E000./Y0,nsample)
  sigma[:,:,2] = broadcast(*,broadcast(*,sigma[:,:,1],exp.(sighisT.*10)),adj15)
  compdelsig = ((1 .- broadcast(^,1-delsig,2:T)')./delsig) .-1 # creates all the compounding of delsig broadcast(.^,1-delsig,[2:T])'
      #Julia_0_6: compdelsig = ((1.-broadcast(.^,1-delsig,2:T)')./delsig).-1
      
  for t = 3:T
    G_ = exp.((ones(nsample,I)*(t-2)*gT + (sighisT.-gT)*compdelsig[1,t-2])*10) #sum(compdelsig[1:t-1]))*10)
      #Julia_0_6: G_ = exp((ones(nsample,I)*(t-2)*gT + (sighisT.-gT)*compdelsig[1,t-2])*10)
    sigma[:,:,t] = sigma[:,:,2].*G_
  end
  return sigma
end

function population(Pop0,poprates,nsample)
  T = 60
  I = 12
  L = zeros(nsample,I,T) # note order of dimensions
  L[:,:,1] = repeat(Pop0,nsample)
  for t = 2:31
    L[:,:,t] = L[:,:,t-1].*exp.(repeat(poprates[t-1,:]'*10,nsample))
      #Julia_0_6: L[:,:,t] = L[:,:,t-1].*exp(repmat(poprates[t-1,:]'*10,nsample))
  end
  for t = 32:T
    L[:,:,t] = L[:,:,31]
  end
  return L
end

function forcingEx(Fex2000,Fex2100)
  T = 60
  Fex = zeros(1,T)
  Fex = Fex2000*ones(1,T) + 0.1*(0:T-1)'*(Fex2100-Fex2000)
  for t = 12:T
    Fex[t] = 0.3
  end
  return Fex
end

function tfactorp(A0,gy0,tgl,delA,gamma,Crate,Cratio,y0,nsample)
  T = 60
  I = 12
  tfp = zeros(nsample,I,T) # note order of dimensions
  tfp[:,:,1] = repeat(A0,nsample)
  tfp[:,:,2] = broadcast(*,tfp[:,:,1],exp.(10*(1-gamma)*gy0))
      #Julia_0_6: tfp[:,:,2] = broadcast(*,tfp[:,:,1],exp(10*(1-gamma)*gy0))
  compdelA = repeat(exp.(-delA*(1:(T-2))'),nsample) # creates all the compounded delA values
      #Julia_0_6: compdelA = repmat(exp(-delA*(1:(T-2))'),nsample)
  gtUS = (1-gamma)*(tgl*ones(nsample,T-2) + (repeat(gy0[:,1],1,(T-2)) - tgl*ones(nsample,T-2)).*compdelA) # growth rates in US, periods 3 to 60
  cgtus = cumsum(gtUS, dims=2)
      #Julia_0_6: cgtus = cumsum(gtUS,2)
  tfp[:,1,3:T] = tfp[:,1,2].*exp.(cgtus.*10) # USA is correct
      #Julia_0_6: tfp[:,1,3:T] = tfp[:,1,2].*exp(cgtus.*10)
  fac = zeros(nsample,I-1)
  for i = 1:nsample
    fac[i,:] = log.(y0[1]./y0[2:I]') + log.(Cratio) + 10*(gy0[i,1].-gy0[i,2:I])'
      #Julia_0_6: fac[i,:] = log(y0[1]./y0[2:I]') + log(Cratio) + 10*(gy0[i,1].-gy0[i,2:I])'
  end
  k = (1 .-Crate).^(0:T-3)'
      #Julia_0_6: k = (1-Crate).^(0:T-3)'
  kR = zeros(nsample,I-1,T-2)
  for i = 1:nsample
    kR[i,:,:] = Crate[i].*(1-gamma)*0.1*(fac[i,:]*k[i,:]')
  end
  gtUS_ = permutedims(cat(gtUS,gtUS,gtUS,gtUS,gtUS,gtUS,gtUS,gtUS,gtUS,gtUS,gtUS;dims=3),[1 3 2]) # adds third dimension to gtUS (manual at the moment since I-1 = 11 is fixed)
      #Julia_0_6: gtUS_ = permutedims(cat(3,gtUS,gtUS,gtUS,gtUS,gtUS,gtUS,gtUS,gtUS,gtUS,gtUS,gtUS),[1 3 2])
  gtR = gtUS_ + kR
  cgtR = cumsum(gtR, dims=3)
      #Julia_0_6: cgtR = cumsum(gtR,3)
  tfp[:,2:I,3:T] = tfp[:,2:I,2].*exp.(10*cgtR)
  return tfp
end

function landuse(EL0,delL)
  T = 60
  EL = EL0'.*(1-delL).^(0:T-1)'
  return EL
end

function elasticity2attribution(e, shares)
  # damage attribution vector for given income elasticity of damage and income distribution
	da = shares.^e
  d = zeros(5,12)
	for i = 1:12
		d[:,i] = da[:,i]/sum(da[:,i])
	end
	return d
end

function damage(temp, psi)
  # maps atmospheric temperature to damage
	# calculate all regions' damage term
  # input:
  #   temp (scalar): atmospheric temperature (relative to preindustrial)
  #   psi Array(3,12): damage function coefficients
	D = (psi[1, :].*temp + psi[2, :].*temp^2 + (psi[3,:].*temp).^7).*0.01
end

function tempforcing(mat, fex, xi, transition, stock)
  # temperature cycle forced by carbon mass
	#inputs:
  #   mat (scalar): atmospheric carbon mass
  #   fex (scalar): exogenous forcing
  #   xi Array(1,7): forcing eqn parameters
  #   transition Array(2,2): stochastic temperature flow matrix
  #   stock Array(1,2): temperature stock in atmosphere and oceans
	forcing = xi[2]*(xi[1]*log2((mat+0.000001)/xi[6])+fex)
	T = stock*transition + [1 0]*forcing
end

function Mflow(stock, flow, transition)
  # carbon cycle forced by emissions
	#inputs:
  #   stock Array(1,3): carbon stock in three reservoirs
  #   flow (scalar): atmospheric carbon emissions
  #   transition Array(3,3): stochasting carbon cycle matrix
	M = stock*transition + 10*[1 0 0]*flow
end

function fromtax(tax,P,Tm)
  # this is the NICE model
  # maps the carbon tax (length(tax) < Tm) to consumption (Tmx12x5) using the parameters in the parameter draw P
  
  #consumption
  c = Array{Float64}(undef, Tm, 12, 5)
        #Julia_0_6: c = Array(Float64, Tm, 12, 5)
  cbar = Array{Float64}(undef, Tm, 12)
        #Julia_0_6: cbar = Array(Float64, Tm, 12)
  
  #capital
  K = Array{Float64}(undef, Tm, 12)
        #Julia_0_6:   K = Array(Float64, Tm, 12)
  K[1,:] = P.K0
  
  #temperature
  T = Array{Float64}(undef, Tm, 2)
       #Julia_0_6:   T = Array(Float64, Tm, 2)
	T[1, :] = P.T0
  T[2, :] = P.T1
  
  #emissions
  E = Array{Float64}(undef, Tm, 12)
       #Julia_0_6:   E = Array(Float64, Tm, 12)
  E[1,:] = P.E0/1000

  #carbon mass
  M = Array{Float64}(undef, Tm, 3)
       #Julia_0_6:   M = Array(Float64, Tm, 3)
	M[1, :] = P.M0
  M[2, :] = P.M1
  
	#savings
	s1 = P.para[4]/(1+P.para[1])^10
	S = ones(Tm,12).*s1

  #tax
  # TAX = [0; tax; maximum(P.pb,2)[(length(tax)+2):end]]
  TAX = maximum(P.pb,dims=2)
      #Julia_0_6:   TAX = maximum(P.pb,2)
	TAX[1] = 0
  TAX[2:length(tax)+1] = tax
	#mitition rate, abatement cost, damage, deflator
  mu = Array{Float64}(undef, Tm, 12)
       #Julia_0_6: mu = Array(Float64, Tm, 12)
	lam = Array{Float64}(undef, Tm, 12)
       #Julia_0_6: lam = Array(Float64, Tm, 12)
	D = Array{Float64}(undef, Tm, 12)
       #Julia_0_6: D = Array(Float64, Tm, 12)
  AD = Array{Float64}(undef, Tm, 12)
       #Julia_0_6: AD = Array(Float64, Tm, 12)

  #output
	Y = Array{Float64}(undef, Tm, 12)
       #Julia_0_6: Y = Array(Float64, Tm, 12)
  Q = Array{Float64}(undef, Tm, 12)
       #Julia_0_6: Q = Array(Float64, Tm, 12)

  #Period 1
  mu[1,:] = max.(min.((TAX[1]./P.pb[1,:]).^(1/(P.th2-1)),1),0) # mu between 0 and 1 element by element
      #Julia_0_6: mu[1,:] = max(min((TAX[1]./P.pb[1,:]).^(1/(P.th2-1)),1),0)
	lam[1,:] = max.(min.(P.th1[1,:].*mu[1,:].^P.th2,1),0) # lam between 0 and 1 element by element
      #Julia_0_6:   lam[1,:] = max(min(P.th1[1,:].*mu[1,:].^P.th2,1),0)
  D[1,:] = damage(T[1,1],P.psi)
  AD[1,:] = (1 .-lam[1,:])./(1 .+D[1,:])
      #Julia_0_6: AD[1,:] = (1-lam[1,:])./(1+D[1,:])
      #Julia_1_1_case
  Y[1,:] = P.A[1,:].*P.L[1,:].^(1-P.para[4]).*K[1,:].^P.para[4]
	Q[1,:] = AD[1,:].*Y[1,:]
  cbar[1,:] = (1 .-S[1,:]).*Q[1,:]./P.L[1,:]
      #Julia_0_6: cbar[1,:] = (1-S[1,:]).*Q[1,:]./P.L[1, :]
      #Julia_1_1_case

	#quintile consumptions period 1
	for i = 1:5
    c[1,:,i] = 5*cbar[1,:].*((1 .+D[1, :]).*P.q[i, :] - D[1, :].*P.d[i, :])
        #Julia_0_6: c[1,:,i] = 5*cbar[1,:].*((1+D[1,:]).*P.q[i,:] - D[1,:].*P.d[i,:])
        #Julia_1_1_case
  end

  # Period 2
  K[2, :] = max.(S[1,:].*Q[1, :]*10,0) # prevent negative capital (note, this will not bind at the optimum, but prevents the optmization from crashing)
        #Julia_0_6: K[2,:] = max(S[1,:].*Q[1,:]*10,0)
  Y[2, :] = P.A[2,:].*P.L[2,:].^(1-P.para[4]).*K[2,:].^P.para[4]
  mu[2, :] = max.(min.((TAX[2]./P.pb[2,:]).^(1/(P.th2-1)),1),0)
        #Julia_0_6: mu[2,:] = max(min((TAX[2]./P.pb[2,:]).^(1/(P.th2-1)),1),0)
  E[2, :] = (1 .-mu[2, :]).*P.sigma[2,:].*Y[2, :]
        #Julia_0_6: E[2,:] = (1 - mu[2, :]).*P.sigma[2, :].*Y[2, :]
      #Julia_1_1_case
  M[3, :] = Mflow(M[2,:]', sum(E[2, :] + P.EL[2, :]), P.TrM)
  lam[2, :] = max.(min.(P.th1[2, :].*mu[2, :].^P.th2,1),0)
        #Julia_0_6: lam[2, :] = max(min(P.th1[2, :].*mu[2, :].^P.th2,1),0)
	D[2, :] = damage(T[2, 1], P.psi)
  AD[2, :] = (1 .-lam[2, :])./(1 .+D[2, :])  
        #Julia_0_6: AD[2, :] = (1-lam[2, :])./(1+D[2, :])
        #Julia_1_1_case
	Q[2,:] = AD[2,:].*Y[2,:]
  cbar[2,:] = (1 .-S[2, :]).*Q[2, :]./P.L[2, :]
        #Julia_0_6: cbar[2,:] = (1-S[2, :]).*Q[2, :]./P.L[2, :]
        #Julia_1_1_case
   
	#quintile consumptions period 2
	for i = 1:5
    c[2,:,i] = max.(5*cbar[2,:].*((1 .+D[2, :]).*P.q[i, :] - D[2, :].*P.d[i, :]), P.tol)
        #Julia_0_6: c[2,:,i] = max(5*cbar[2,:].*((1+D[2, :]).*P.q[i, :] - D[2, :].*P.d[i, :]), P.tol)
      #Julia_1_1_case
	end
  K[3, :] = max.(S[2, :].*Q[2, :].*10,0) # prevent negative capital (note, this will not bind at the optimum, but prevents the optmization from crashing)
        #Julia_0_6: K[3, :] = max(S[2, :].*Q[2, :].*10,0)

  #periods 3 to Tm-1
	for t = 3:(Tm - 1)
    Y[t,:] = P.A[t,:].*P.L[t,:].^(1-P.para[4]).*K[t,:].^P.para[4]
    mu[t, :] =  max.(min.((TAX[t]./P.pb[t,:]).^(1/(P.th2-1)),1),0)
        #Julia_0_6: mu[t, :] =  max(min((TAX[t]./P.pb[t,:]).^(1/(P.th2-1)),1),0)
    lam[t, :] = max.(min.(P.th1[t, :].*mu[t, :].^P.th2,1),0) 
        #Julia_0_6: lam[t, :] = max(min(P.th1[t, :].*mu[t, :].^P.th2,1),0)
    E[t, :] = (1 .-mu[t, :]).*P.sigma[t, :].*Y[t, :]
        #Julia_0_6: E[t, :] = (1 - mu[t, :]).*P.sigma[t, :].*Y[t, :]
        #Julia_1_1_case
    M[t+1, :] = Mflow(M[t,:]', sum(E[t, :] + P.EL[t, :]), P.TrM)
		Mbar = (M[t+1, 1] +M[t, 1])/2
		T[t, :] = tempforcing(Mbar, P.Fex[t], P.xi, P.TrT, T[t-1, :]')
		D[t, :] = damage(T[t, 1], P.psi)
    AD[t, :] = (1 .-lam[t, :])./(1 .+D[t, :])
        #Julia_0_6: AD[t, :] = (1-lam[t, :])./(1+D[t, :])
        #Julia_1_1_case
		Q[t, :] = AD[t, :].*Y[t, :]                        
    cbar[t,:] = (1 .-S[t, :]).*Q[t, :]./P.L[t, :] 
        #Julia_0_6: cbar[t,:] = (1-S[t, :]).*Q[t, :]./P.L[t, :]
        #Julia_1_1_case

		for i = 1:5
      c[t, :, i] = max.(5*cbar[t,:].*((1 .+D[t, :]).*P.q[i, :] - D[t, :].*P.d[i, :]), P.tol)
        #Julia_0_6: c[t, :, i] = max(5*cbar[t,:].*((1+D[t, :]).*P.q[i, :] - D[t, :].*P.d[i, :]), P.tol)
        #Julia_1_1_case
		end
    K[t+1, :] = max.(S[t, :].*Q[t, :]*10,0) # prevent negative capital (note, this will not bind at the optimum, but prevents the optimization from crashing)
        #Julia_0_6: K[t+1, :] = max(S[t, :].*Q[t, :]*10,0)
    end
  
  # Period Tm
  Y[Tm, :] = P.A[Tm,:].*P.L[Tm,:].^(1-P.para[4]).*K[Tm,:].^P.para[4]
  mu[Tm, :] =  max.(min.((TAX[Tm]./P.pb[Tm,:]).^(1/(P.th2-1)),1),0)
        #Julia_0_6: mu[Tm, :] =  max(min((TAX[Tm]./P.pb[Tm,:]).^(1/(P.th2-1)),1),0)
  lam[Tm, :] = max.(min.(P.th1[Tm, :].*mu[Tm, :].^P.th2,1),0)
        #Julia_0_6: lam[Tm, :] = max(min(P.th1[Tm, :].*mu[Tm, :].^P.th2,1),0)
	T[Tm, :] = tempforcing(M[Tm, 1], P.Fex[Tm], P.xi, P.TrT, T[Tm-1, :]')
	D[Tm, :] = damage(T[Tm, 1], P.psi)
  AD[Tm, :] = (1 .-lam[Tm, :])./(1 .+D[Tm, :]) 
        #Julia_0_6: AD[Tm, :] = (1-lam[Tm, :])./(1+D[Tm, :])
        #Julia_1_1_case
	Q[Tm, :] = AD[Tm, :].*Y[Tm, :]
  cbar[Tm,:] = (1 .-S[Tm, :]).*Q[Tm, :]./P.L[Tm, :]
        #Julia_0_6:  cbar[Tm,:] = (1-S[Tm, :]).*Q[Tm, :]./P.L[Tm, :]
        #Julia_1_1_case
	for i = 1:5
    c[Tm, :, i] =  max.(5*cbar[Tm,:].*((1 .+D[Tm, :]).*P.q[i, :] - D[Tm, :].*P.d[i, :]), P.tol)
        #Julia_0_6: c[Tm, :, i] =  max(5*cbar[Tm,:].*((1+D[Tm, :]).*P.q[i, :] - D[Tm, :].*P.d[i, :]), P.tol)
        #Julia_1_1_case
	end

  return c,K,T,E,M,mu,lam,D,AD,Y,Q,cbar
end

function welfareN(c, L, rho, eta, Tm)
  R = 1 ./(1+rho).^(10 .*(0:(Tm-1)))
      #Julia_0_6: R = 1./(1+rho).^(10.*(0:(Tm-1)))
  A = Array{Float64}(undef, Tm, 12, 5)
      #Julia_0_6: A = Array(Float64, Tm, 12, 5)
	for i = 1:5
		A[:,:,i] = 0.2*L[1:Tm,:].*c[:,:,i].^(1-eta)
	end
  B = squeeze(sum(sum(A, dims=3), dims=2), dims=3)'
      #Julia_0_6: B = squeeze(sum(sum(A,3),2),3)'
	W = (B*R/(1-eta))[1]
  return W
end

function welfareR(c, L, rho, eta, Tm)
  R = 1 ./(1+rho).^(10 .*(0:(Tm-1)))
      #Julia_0_6: R = 1./(1+rho).^(10.*(0:(Tm-1)))
      #Julia_1_1_case
	A = L[1:Tm,:].*c[:,:].^(1-eta)
  B = sum(A, dims=2)'
      #Julia_0_6: B = sum(A,2)'
	W = ((B*R)/(1-eta))[1]
  return W
end

function welfareD(c,L,rho,eta,Tm)
  R = 1 ./(1+rho).^(10 .*(0:(Tm-1)))
      #Julia_0_6: R = 1./(1+rho).^(10.*(0:(Tm-1)))
      #Julia_1_1_case
  A = sum(L[1:Tm,:].*c[:,:], dims=2)'
      #Julia_0_6: A = sum(L[1:Tm,:].*c[:,:],2)'
  B = ((A.^(1-eta)*R)/(1-eta))
  W = B[1]
end

function W2EW(x,nu,eta)
  mean(((x.(1-eta)).^(1/(1-eta))).^(1-nu)./(1-nu))
end

function tax2welfare(tax, P, rho, eta, Tm; model="NICE")
  if model == "NICE"
    c = fromtax(tax, P, Tm)[1]
    W = welfareN(c, P.L, rho, eta, Tm)
  elseif model == "RICE"
    c = fromtax(tax, P, Tm)[12]
    W = welfareR(c,P.L,rho,eta,Tm)
  elseif model == "DICE"
    c = fromtax(tax, P, Tm)[12]
    W = welfareD(c,P.L,rho,eta,Tm)
  end
  return W
end

function tax2expectedwelfare(tax, P, rho, eta, nu, Tm, tm, lm, idims; model="NICE")
  nsample=length(P)
  
  if model == "NICE"
    # println("   + NICE in tax2expectedwelfare")
    c = zeros(Tm,12,5,nsample) # will contain per capita consumption at time t, in region I, in quintile q, for random draw n
    for i = 1:idims
        c[:,:,:,i] = fromtax(tax[1:tm],P[i],Tm)[1] # only consider tm length since we want to create a tax vector of particular length
    end
    for i = idims+1:length(P) # NB length(P) = nsample
        c[:,:,:,i] = fromtax([tax[1:lm];tax[tm+1:end]],P[i],Tm)[1]
    end

    R = 1 ./(1+rho).^(10 .*(0:(Tm-1))) # discount factors for each time period
      #Julia_0_6: R = 1./(1+rho).^(10.*(0:(Tm-1)))
      #Julia_1_1_case
    D = zeros(Tm,12,5,nsample)
    # Convert consumption to per capita discounted utility at time t, in region I (weighted by population), in quintile q, for random draw n
    for t = 1:Tm
      D[t,:,:,:] = ((c[t,:,:,:].^(1-eta)).*R[t])./(1-eta)
    end

    D_ = zeros(Tm,12,5,nsample)
    for i = 1:nsample
      D_[:,:,:,i] = D[:,:,:,i].*cat(P[i].L[1:Tm,:],P[i].L[1:Tm,:],P[i].L[1:Tm,:],P[i].L[1:Tm,:],P[i].L[1:Tm,:];dims=3)/5
        #Julia_0_6: D_[:,:,:,i] = D[:,:,:,i].*cat(3,P[i].L[1:Tm,:],P[i].L[1:Tm,:],P[i].L[1:Tm,:],P[i].L[1:Tm,:],P[i].L[1:Tm,:])/5
    end
    D = D_
    # Now sum over quintiles to get per capita discounted utility at time t, in region I, in random draw n
    B1 = sum(D, dims=3)
      #Julia_0_6: B1 = sum(D,3)
    # Now sum over regions to get per capita discounted utility at time t, in random draw n
    B2 = sum(B1, dims=2)
      #Julia_0_6: B2 = sum(B1,2)
    # Now sum over time to get per capita lifetime discounted utility in random draw n, and undo the concavity to get a "certainty equivalent" consumption measure
    B3 = (sum(B2, dims=1).*(1-eta)).^(1/(1-eta))
      #Julia_0_6: B3 = (sum(B2,1).*(1-eta)).^(1/(1-eta))
    # Now sum over random draws with the risk adjustment (nu) to get total world welfare (normalizing by nsample)
    W = sum(B3.^(1-nu))*(1/(1-nu))./nsample
    # W = (0.33*B3[1,1,1,1].^(1-nu) + 0.67*B3[1,1,1,2].^(1-nu))*(1/(1-nu)) # Test for unequal probabilities effect on learning...
    return W,c
  
  elseif model == "RICE"
    # println("   + NICE in tax2expectedwelfare")
    c = zeros(Tm,12,nsample) # will contain per capita consumption at time t, in region I, in quintile q, for random draw n
    for i = 1:idims
        c[:,:,i] = fromtax(tax[1:tm],P[i],Tm)[12] # only consider tm length since we want to create a tax vector of particular length
    end
    for i = idims+1:length(P) # NB length(P) = nsample
        c[:,:,i] = fromtax([tax[1:lm];tax[tm+1:end]],P[i],Tm)[12]
    end

    R = 1 ./(1+rho).^(10 .*(0:(Tm-1))) # discount factors for each time period
      #Julia_0_6:   R = 1./(1+rho).^(10.*(0:(Tm-1)))
      #Julia_1_1_case
    D = zeros(Tm,12,nsample)
    # Convert consumption to per capita discounted utility at time t, in region I (weighted by population), in quintile q, for random draw n
    for t = 1:Tm #convert per capita consumption to discounted utility
      D[t,:,:] = ((c[t,:,:].^(1-eta)).*R[t])./(1-eta)
    end
    for i = 1:nsample # weight discounted utility by regional population
      D[:,:,i] = D[:,:,i].*P[i].L[1:Tm,:]
    end
    # Now sum over regions to get per capita discounted utility at time t, in random draw n
    B2 = sum(D, dims=2)
      #Julia_0_6: B2 = sum(D,2) 
    # Now sum over time to get per capita lifetime discounted utility in random draw n, and undo the concavity to get a "certainty equivalent" consumption measure
    B3 = (sum(B2,dims=1).*(1-eta)).^(1/(1-eta))
      #Julia_0_6: B3 = (sum(B2,1).*(1-eta)).^(1/(1-eta))
    # Now sum over random draws with the risk adjustment (nu) to get total world welfare (normalizing by nsample)
    W = sum(B3.^(1-nu))*(1/(1-nu))./nsample
    return W,c
  
  elseif model == "DICE"
    # println("   + DICE in tax2expectedwelfare")
    c = zeros(Tm,nsample) # will contain per capita consumption at time t, in region I, in quintile q, for random draw n
    for i = 1:idims
        c[:,i] = sum(fromtax(tax[1:tm],P[i],Tm)[12].*P[i].L[1:Tm,:],dims=2)./sum(P[i].L[1:Tm,:],dims=2) # only consider tm length since we want to create a tax vector of particular length
            #Julia_0_6: c[:,i] = sum(fromtax(tax[1:tm],P[i],Tm)[12].*P[i].L[1:Tm,:],2)./sum(P[i].L[1:Tm,:],2)
    end
    for i = idims+1:length(P) # NB length(P) = nsample
        c[:,i] = sum(fromtax([tax[1:lm];tax[tm+1:end]],P[i],Tm)[12].*P[i].L[1:Tm,:],dims=2)./sum(P[i].L[1:Tm,:],dims=2)
            #Julia_0_6: c[:,i] = sum(fromtax([tax[1:lm];tax[tm+1:end]],P[i],Tm)[12].*P[i].L[1:Tm,:],2)./sum(P[i].L[1:Tm,:],2) 
    end

    R = 1 ./(1+rho).^(10 .*(0:(Tm-1))) # discount factors for each time period
      #Julia_0_6: R = 1./(1+rho).^(10.*(0:(Tm-1)))
      #Julia_1_1_case
    D = zeros(Tm,nsample)
    for t = 1:Tm #convert per capita consumption to discounted utility
      D[t,:] = ((c[t,:].^(1-eta)).*R[t])./(1-eta)
    end
    for i=1:nsample # weight discounted utility by global population
      D[:,i] = D[:,i].*sum(P[i].L[1:Tm,:],dims=2)
          #Julia_0_6: D[:,i] = D[:,i].*sum(P[i].L[1:Tm,:],2)
    end
    # Now sum over time to get per capita lifetime discounted utility in random draw n
    B3 = (sum(D,dims=1).*(1-eta)).^(1/(1-eta))
      #Julia_0_6: B3 = (sum(D,1).*(1-eta)).^(1/(1-eta))
    # Now sum over random draws with the risk adjustment (nu) to get total world welfare (normalizing by nsample)
    W = sum(B3.^(1-nu))*(1/(1-nu))./nsample
    return W,c
  end
  
end

function tax2expectedwelfare10(tax, P, rho, eta, nu, Tm, tm, lm; model="NICE")
  nsample=length(P)
  if model == "NICE"
    c = zeros(Tm,12,5,nsample) # will contain per capita consumption at time t, in region I, in quintile q, for random draw n
    for i = 1:nsample
        c[:,:,:,i] = fromtax(tax[:,i],P[i],Tm)[1] # only consider tm length since we want to create a tax vector of particular length
    end

    R = 1 ./(1+rho).^(10 .*(0:(Tm-1))) # discount factors for each time period
      #Julia_0_6: R = 1./(1+rho).^(10.*(0:(Tm-1)))
      #Julia_1_1_case
    D = zeros(Tm,12,5,nsample)
    # Convert consumption to per capita discounted utility at time t, in region I (weighted by population), in quintile q, for random draw n
    for t = 1:Tm
      D[t,:,:,:] = ((c[t,:,:,:].^(1-eta)).*R[t])./(1-eta)
    end
    D_ = zeros(Tm,12,5,nsample)
    for i = 1:nsample
      D_[:,:,:,i] = D[:,:,:,i].*cat(P[i].L[1:Tm,:],P[i].L[1:Tm,:],P[i].L[1:Tm,:],P[i].L[1:Tm,:],P[i].L[1:Tm,:];dims=3)/5
      #Julia_0_6: D_[:,:,:,i] = D[:,:,:,i].*cat(3,P[i].L[1:Tm,:],P[i].L[1:Tm,:],P[i].L[1:Tm,:],P[i].L[1:Tm,:],P[i].L[1:Tm,:])/5
      #Julia_1_1_case
    end
    D = D_
    # Now sum over quintiles to get per capita discounted utility at time t, in region I, in random draw n
    B1 = sum(D,dims=3)
        #Julia_0_6: B1 = sum(D,3)
    # Now sum over regions to get per capita discounted utility at time t, in random draw n
    B2 = sum(B1,dims=2)
        #Julia_0_6: B2 = sum(B1,2)
    # Now sum over time to get per capita lifetime discounted utility in random draw n, and undo the concavity to get a "certainty equivalent" consumption measure
    B3 = (sum(B2,dims=1).*(1-eta)).^(1/(1-eta))
        #Julia_0_6: B3 = (sum(B2,1).*(1-eta)).^(1/(1-eta))
    # Now sum over random draws with the risk adjustment (nu) to get total world welfare (normalizing by nsample)
    W = sum(B3.^(1-nu))*(1/(1-nu))./nsample
    # W = (0.33*B3[1,1,1,1].^(1-nu) + 0.67*B3[1,1,1,2].^(1-nu))*(1/(1-nu)) # Test for unequal probabilities effect on learning...
    return W,c

  elseif model == "RICE"
    c = zeros(Tm,12,nsample) # will contain per capita consumption at time t, in region I, for random draw n
    for i = 1:nsample
        c[:,:,i] = fromtax(tax[:,i],P[i],Tm)[12] # only consider tm length since we want to create a tax vector of particular length
    end
    R = 1 ./(1+rho).^(10 .*(0:(Tm-1))) # discount factors for each time period
      #Julia_0_6: R = 1./(1+rho).^(10.*(0:(Tm-1)))
      #Julia_1_1_case
    D = zeros(Tm,12,nsample)
    # Convert consumption to per capita discounted utility at time t, in region I (weighted by population), in quintile q, for random draw n
    for t = 1:Tm #convert per capita consumption to discounted utility
      D[t,:,:] = ((c[t,:,:].^(1-eta)).*R[t])./(1-eta)
    end
    for i = 1:nsample # weight discounted utility by regional population
      D[:,:,i] = D[:,:,i].*P[i].L[1:Tm,:]
    end
    # Now sum over regions to get per capita discounted utility at time t, in random draw n
    B2 = sum(D,dims=2)
      #Julia_0_6: B2 = sum(D,2)
    # Now sum over time to get per capita lifetime discounted utility in random draw n, and undo the concavity to get a "certainty equivalent" consumption measure
    B3 = (sum(B2,dims=1).*(1-eta)).^(1/(1-eta))
      #Julia_0_6: B3 = (sum(B2,1).*(1-eta)).^(1/(1-eta))
    # Now sum over random draws with the risk adjustment (nu) to get total world welfare (normalizing by nsample)
    W = sum(B3.^(1-nu))*(1/(1-nu))./nsample
    return W,c

  elseif model == "DICE"
    c = zeros(Tm,nsample) # will contain per capita consumption at time t, in region I, in quintile q, for random draw n
    for i = 1:nsample
        c[:,i] = sum(fromtax(tax[:,i],P[i],Tm)[12].*P[i].L[1:Tm,:],dims=2)./sum(P[i].L[1:Tm,:],dims=2) # only consider tm length since we want to create a tax vector of particular length
            #Julia_0_6: c[:,i] = sum(fromtax(tax[:,i],P[i],Tm)[12].*P[i].L[1:Tm,:],2)./sum(P[i].L[1:Tm,:],2)
    end
    R = 1 ./(1+rho).^(10 .*(0:(Tm-1))) # discount factors for each time period
      #Julia_0_6: R = 1./(1+rho).^(10.*(0:(Tm-1)))
      #Julia_1_1_case
    D = zeros(Tm,nsample)
    for t = 1:Tm #convert per capita consumption to discounted utility
      D[t,:] = ((c[t,:].^(1-eta)).*R[t])./(1-eta)
    end
    for i=1:nsample # weight discounted utility by global population
      D[:,i] = D[:,i].*sum(P[i].L[1:Tm,:],dims=2)
         #Julia_0_6: D[:,i] = D[:,i].*sum(P[i].L[1:Tm,:],2)
    end
    # Now sum over time to get per capita lifetime discounted utility in random draw n
    B3 = (sum(D,dims=1).*(1-eta)).^(1/(1-eta))
      #Julia_0_6: B3 = (sum(D,1).*(1-eta)).^(1/(1-eta))
    # Now sum over random draws with the risk adjustment (nu) to get total world welfare (normalizing by nsample)
    W = sum(B3.^(1-nu))*(1/(1-nu))./nsample
    return W,c
  end
end

function welfare2c_bar(W, L, rho, eta, nu, Tm)
  R = 1 ./(1+rho).^(10 .*(0:(Tm-1))) # discount factors for each time period
      #Julia_0_6: R = 1./(1+rho).^(10.*(0:(Tm-1)))
      #Julia_1_1_case
  D = sum(R.*L)
  cbar = (((1-nu)*W)^(1/(1-nu)))/((D)^(1/(1-eta)))
  return cbar
end

function VarsFromTaxes(taxes_1, taxes_2, PP, nsample; model="NICE")
  # Create storage objects
  if (model == "RICE") | (model == "DICE")
    c = Array{Float64}(undef, Tm, 12, nsample)
      #Julia_0_6: c = Array(Float64, Tm, 12, nsample)
  else
    c = Array{Float64}(undef, Tm, 12, 5, nsample)
      #Julia_0_6: c = Array(Float64, Tm, 12, 5, nsample)
  end
  K = Array{Float64}(undef, Tm, 12, nsample)
      #Julia_0_6: K = Array(Float64, Tm, 12, nsample)
  T = Array{Float64}(undef, Tm, 2, nsample)
      #Julia_0_6: T = Array(Float64, Tm, 2, nsample)
  E = Array{Float64}(undef, Tm, 12, nsample)
      #Julia_0_6: E = Array(Float64, Tm, 12, nsample)
  M = Array{Float64}(undef, Tm, 3, nsample)
      #Julia_0_6: M = Array(Float64, Tm, 3, nsample)
  mu = Array{Float64}(undef, Tm, 12, nsample)
      #Julia_0_6: mu = Array(Float64, Tm, 12, nsample)
  lam = Array{Float64}(undef, Tm, 12, nsample)
      #Julia_0_6: lam = Array(Float64, Tm, 12, nsample)
  D = Array{Float64}(undef, Tm, 12, nsample)
      #Julia_0_6: D = Array(Float64, Tm, 12, nsample)
  AD = Array{Float64}(undef, Tm, 12, nsample)
      #Julia_0_6: AD = Array(Float64, Tm, 12, nsample)
  Y = Array{Float64}(undef, Tm, 12, nsample)
      #Julia_0_6:  Y = Array(Float64, Tm, 12, nsample)
  Q = Array{Float64}(undef, Tm, 12, nsample)
      #Julia_0_6: Q = Array(Float64, Tm, 12, nsample)


  # Store data
  for i = 1:Int(max(round(nsample/2),1))
    if (model == "RICE") | (model == "DICE")
      c[:,:,i] = fromtax(taxes_1,PP[i],Tm)[12]
    else
      c[:,:,:,i] = fromtax(taxes_1,PP[i],Tm)[1]
    end
    K[:,:,i] = fromtax(taxes_1,PP[i],Tm)[2]
    T[:,:,i] = fromtax(taxes_1,PP[i],Tm)[3]
    E[:,:,i] = fromtax(taxes_1,PP[i],Tm)[4]
    M[:,:,i] = fromtax(taxes_1,PP[i],Tm)[5]
    mu[:,:,i] = fromtax(taxes_1,PP[i],Tm)[6]
    lam[:,:,i] = fromtax(taxes_1,PP[i],Tm)[7]
    D[:,:,i] = fromtax(taxes_1,PP[i],Tm)[8]
    AD[:,:,i] = fromtax(taxes_1,PP[i],Tm)[9]
    Y[:,:,i] = fromtax(taxes_1,PP[i],Tm)[10]
    Q[:,:,i] = fromtax(taxes_1,PP[i],Tm)[11]
  end

  for i = (Int(max(round(nsample/2),1))+1):nsample
    if (model == "RICE") | (model == "DICE")
      c[:,:,i] = fromtax(taxes_2,PP[i],Tm)[12]
    else
      c[:,:,:,i] = fromtax(taxes_2,PP[i],Tm)[1]
    end
    K[:,:,i] = fromtax(taxes_2,PP[i],Tm)[2]
    T[:,:,i] = fromtax(taxes_2,PP[i],Tm)[3]
    E[:,:,i] = fromtax(taxes_2,PP[i],Tm)[4]
    M[:,:,i] = fromtax(taxes_2,PP[i],Tm)[5]
    mu[:,:,i] = fromtax(taxes_2,PP[i],Tm)[6]
    lam[:,:,i] = fromtax(taxes_2,PP[i],Tm)[7]
    D[:,:,i] = fromtax(taxes_2,PP[i],Tm)[8]
    AD[:,:,i] = fromtax(taxes_2,PP[i],Tm)[9]
    Y[:,:,i] = fromtax(taxes_2,PP[i],Tm)[10]
    Q[:,:,i] = fromtax(taxes_2,PP[i],Tm)[11]
  end
  return c, K, T, E, M, mu, lam, D, AD, Y, Q
end

function VarsFromTaxes10(taxes, PP, nsample; model = "NICE", Tm=32)
  # Create storage objects
  if (model == "RICE") | (model == "DICE")
    c = Array{Float64}(undef, Tm, 12, nsample)
  else
    c = Array{Float64}(undef, Tm, 12, 5, nsample)
  end
  K = Array{Float64}(undef, Tm, 12, nsample)
  T = Array{Float64}(undef, Tm, 2, nsample)
  E = Array{Float64}(undef, Tm, 12, nsample)
  M = Array{Float64}(undef, Tm, 3, nsample)
  mu = Array{Float64}(undef, Tm, 12, nsample)
  lam = Array{Float64}(undef, Tm, 12, nsample)
  D = Array{Float64}(undef, Tm, 12, nsample)
  AD = Array{Float64}(undef, Tm, 12, nsample)
  Y = Array{Float64}(undef, Tm, 12, nsample)
  Q = Array{Float64}(undef, Tm, 12, nsample)

  # Store data
  for i = 1:nsample
    if (model == "RICE") | (model == "DICE")
      c[:,:,i] = fromtax(taxes[:,i],PP[i],Tm)[12]
    else
      c[:,:,:,i] = fromtax(taxes[:,i],PP[i],Tm)[1]
    end
    K[:,:,i], T[:,:,i], E[:,:,i], M[:,:,i], mu[:,:,i], lam[:,:,i], D[:,:,i], AD[:,:,i], Y[:,:,i], Q[:,:,i] = fromtax(taxes[:,i],PP[i],Tm)[2:11]
  end
  return c, K, T, E, M, mu, lam, D, AD, Y, Q
end

# Create storage object
mutable struct Results      #Julia_0_6: type Results
  regime
  nsample
  Tm
  tm
  lm
  Regions
  taxes_1
  taxes_2
  EWelfare
  c
  K
  T
  E
  M
  mu
  lam
  D
  AD
  Y
  Q
  rho
  eta
  nu
  PP
end

mutable struct Results10      #Julia_0_6: type Results
  regime
  nsample
  Tm
  tm
  lm
  Regions
  taxes
  EWelfare
  c
  K
  T
  E
  M
  mu
  lam
  D
  AD
  Y
  Q
  rho
  eta
  nu
  PP
  optiRet
end

function FrameFromResults(res, Tm, nsample, Regions, idims)
  # set up dataframe with periods, regions, State
  if length(size(res.c)) > 2
    dataP = DataFrame(ID = 1:(Tm*12*nsample), 
                      State = reshape(repeat(collect(1:nsample)',Tm*12),Tm*12*nsample,1)[:,1], 
                      Region = repeat(reshape(repeat(Regions,Tm),Tm*12,1),nsample)[:,1], 
                      Year = repeat(repeat(10 .*(0:Tm-1) .+2005,12),nsample))
      #Julia_0_6:     Year = repmat(repmat(10*(0:Tm-1)+2005,12),nsample))
      #Julia_1_1_case
    # add taxes (to the correct states)
    dataP[:tax] = [repeat(repeat(res.taxes_1,12),idims);repeat(repeat(res.taxes_2,12),nsample-idims)]
    dataP[:T] = reshape(repeat(res.T[:,1,:],12),Tm*12*nsample)
    # add consumption quintiles
    if length(size(res.c)) == 4
      confield = [:cq1, :cq2, :cq3, :cq4, :cq5]
      cquintiles=reshape(permutedims(res.c,[1 2 4 3]),Tm*12*nsample,5)
      m=1
      for field in confield
        dataP[Symbol(field)] = cquintiles[:,m]
        m+=1
      end
    elseif length(size(res.c)) == 3
      cons = reshape(res.c,Tm*12*nsample,1)[:,1]
      dataP[:c] = cons
    end

    # add remaining endogenous variables
    for field in [:K,:E,:mu,:lam,:D,:Y]
      dataP[Symbol(field)] = reshape(getfield(res,field),Tm*12*nsample)
    end
    # add exogenous variables
    y = Array{Float64}(undef,Tm*12*nsample,6)
    x = Array{Float64}(undef,Tm*12,nsample)
    k=1
    for field in [:L,:A,:sigma,:th1,:pb,:EL]
      for m in 1:nsample
        x[:,m] = reshape(getfield(res.PP[m],field)[1:Tm,:],Tm*12)
      end
      y[:,k] = reshape(x,Tm*12*nsample)
      dataP[Symbol(field)] = y[:,k]
      k+=1
    end
  end
  return dataP
end

# Define Deep as the type object that will hold all the random parameter draws in createP
mutable struct Deep     #Julia_0_6: type Deep 
  gy0
  sighisT
  TrM12
  xi1
  psi7
  pw
  ee
  psi2
  Crate
end

#Region Labels
Regions = ["USA" "OECD Europe" "Japan" "Russia" "Non-Russia Eurasia" "China" "India" "Middle East" "Africa" "Latin America" "OHI" "Other non-OECD Asia"]

# Define PP_ as the type that will hold the parameters returned by creatP
struct PP_
  para::Array{Float64,2} # 1x4 vector, constant across nsample, regions, time
  L::Array{Float64,2} # TxI array
  A::Array{Float64} # TxI array
  sigma::Array{Float64} # TxI array
  th1::Array{Float64} # TxI array
  th2::Float64 # scalar (constant)
  pb::Array{Float64}# TxI array
  EL::Array{Float64} # TxI array
  Fex::Array{Float64} # 1xT array
  TrM::Array{Float64} # 3x3 array
  xi::Array{Float64} # 1x7 array
  TrT::Array{Float64} # 2x2 array
  psi::Array{Float64} # 3xI array
  T0::Array{Float64} # 1x2 array (constant)
  T1::Array{Float64} # 1x2 array (constant)
  M0::Array{Float64} # 1x3 array (constant)
  M1::Array{Float64} #1x3 array (constant)
  K0::Array{Float64} # 1xI array
  E0::Array{Float64} # 1x12 vector with 2005 emissions
  R::Array{Float64} # 1xT array
  q::Array{Float64} # 5x12 array (constant)
  d::Array{Float64} # 5x12 array
  tol::Float64 # scalar (constant)
end
