include("RD.jl")
include("lattice.jl")
#logarithmically ranged values
function logrange(x1, x2, n)
	return [10^y for y in range(log10(x1), log10(x2), length=n)]
end

function run_script()
	L=50
	mid::Int=L*(L/2)+L/2
	nn,centers=gen_hex_lattice(L)
	memb0=zeros(Int,L*L,4)	
	totmol=L*L*10
	totA=0.5*totmol
	totB=0.5*totmol
	unit_length=1.0
	Ac=unit_length^2#area associated to each cell, set unit lengthscale
	Stot=L*L*Ac #area of each cell is assumed to be 1, setting unit lengthscale
	#notice that setting volume/surface ratio somehow sets geometry of the problem (membrane+cytosolic reservoir)
	#theta=100000.0
	#V=Stot/theta
	#V=4*(pi/3)*(Stot/(4*pi))^(3/2) #considering a spherical reservoir V=4 *pi R^3 /3 where R=sqrt(S/4*pi) 
	theta=0.2126944621086619#Stot/V
	println(theta)
	V=Stot/theta
	T=20000.0
	Tmeas=10.0
	Nsave=10
	d_timescale=1.0#this sets unit timescale
	dA=d_timescale
	dB=d_timescale
	dEA=d_timescale
	dEB=d_timescale
	#rates in the theory (do not correspond excatly to those used in the simulations), some slight dimensional changes are needed
	kAa_th=1.0
	kAd_th=1.0*kAa_th
	kAc_th=1.0
	kBa_th=1.0
	kBd_th=1.0*kBa_th
	kBc_th=logrange(logrange(0.1,10,17)[4],logrange(0.1,10.0,17)[14],23)[13]*kAc_th
	KMM_th=1.0
	#rates to implement
	kAc=kAc_th
	kBc=kBc_th
	kAa=kAa_th/V
	kAd=kAd_th
	kBa=kBa_th/V
	kBd=kBd_th
	KMM=KMM_th*Ac
	seed=22
	for i in 1:totA
		memb0[rand(1:(L*L)),1]+=1
	end
	for i in 1:totB
		memb0[rand(1:(L*L)),2]+=1	
	end
	EA_tot_n=Int(0.1*L*L)
	EB_tot_n=Int(0.1*L*L)
	EA_mem=Int(floor(EA_tot_n*(theta/(kAd_th/kAa_th))*(totA/Stot)/(1+((theta/(kAd_th/kAa_th))*(totA/Stot)))))
	EB_mem=Int(floor(EB_tot_n*(theta/(kBd_th/kBa_th))*(totB/Stot)/(1+((theta/(kBd_th/kBa_th))*(totB/Stot)))))
	for i in 1:EA_mem
		memb0[rand(1:(L*L)),3]+=1
	end
	for i in 1:EB_mem
		memb0[rand(1:(L*L)),4]+=1
	end
	#set all enzymes in the memebrane initially
	cyto0::Vector{Int64}=[EA_tot_n-EA_mem,EB_tot_n-EB_mem]
	EB_tot=EB_mem/V+cyto0[2]/V #these are in theory volume concentration
	EA_tot=EA_mem/V+cyto0[1]/V
	rho_0=kBc_th*EB_tot/(kAc_th*EA_tot)
	k=KMM_th*Stot/(totmol)
	rho_p=((kBd_th/kBa_th)+((totmol)/V))/(kAd_th/kAa_th)
	rho_m=(kBd_th/kBa_th)/((kAd_th/kAa_th)+((totmol)/V))
	if rho_0 < rho_m
		phi_eq=-1
	elseif rho_0 > rho_p
		phi_eq=1
	else
		#phi_eq=(rho_0*((totA+totB)/(Stot)+2*kAd*V/(kAa*Stot))-((totA+totB)/(Stot)+2*kBd*V/(kBa*Stot)))/(rho_0+1)
		#phi_eq=(rho_0*(1+2*kAd/(kAa*(totA+totB)))-(1+2*kBd/(kBa*(totA+totB))))/(rho_0+1)
		#phi_eq=(rho_0*(1+2*kAd*V/(kAa*(totA+totB)))-(1+2*kBd*V/(kBa*(totA+totB))))/(rho_0+1)
		phi_eq=(Stot/totmol)*(rho_0*(((totmol)/V)+2*kAd_th/kAa_th)-(((totmol)/V)+2*kBd_th/kBa_th))/(Stot*(rho_0+1)/V)
	end
	function phieq(x)
		return (x*(1+2*kAd*V/(kAa*(totA+totB)))-(1+2*kBd*V/(kBa*(totA+totB))))/(x+1)
	end
	xvals=(0.01:0.01:10.0)
	#display(plot(xvals,phieq.(xvals),xscale=:log10))
	#display(vline!([rho_m,rho_p, rho_0]))
	#sleep(5)
	println("estimated rho0 =", rho_0)
	println("estimated rho+ =", rho_p)
	println("estimated rho- =", rho_m)
	println("estimated k =", k)
	println("bistability region between ", k/(1+k) ," and ", (1+k/k))
	println("physical region between ", rho_0/rho_p, " and ", rho_0/rho_m)
	println("estimatted phi_eq/c= ", phi_eq)
	#plot_diagram(rho_p,rho_m,rho_0,k)
	println("now run reaction diffusion")
	@time memb0,cyto0 = run_RD(memb0,cyto0,nn,centers,T,Tmeas,dA,dB,dEA,dEB,kAc,kBc,kAa,kAd,kBa,kBd,KMM,rho_0,seed,Nsave,"test/")
end

run_script()


