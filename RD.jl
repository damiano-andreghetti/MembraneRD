using CavityTools, JLD, Random
#THE FOLLOWING IS TO MAKE PLOT, COMMENT FOR CLUSTER USAGE
#using ColorSchemes,Plots,LaTeXStrings 


function run_RD(cfg_memb::Matrix{Int}, cfg_cyto::Vector{Int}, neig::Matrix{Int}, centers::Matrix{Float64}, T::Float64, Tmeas::Float64, dA::Float64, dB::Float64,dEA::Float64,dEB::Float64,kAc::Float64, kBc::Float64, kAa::Float64,kAd::Float64,kBa::Float64,kBd::Float64, KMM::Float64,rho_0::Float64, seed::Int,Nsave::Int,fold_files::String)
	ran_ng = Random.Xoshiro(seed)
	#cfg_memb[i,1]=number of A moelcules in cell i
	#cfg_memb[i,2]=number of B moelcules in cell i
	#cfg_memb[i,3]=number of EA moelcules in cell i
	#cfg_memb[i,4]=number of EB moelcules in cell i
	#cfg_cyto[1]=number of EA molecules in cytosol
	#cfg_cyto[2]=number of EB molecules in cytosol
	#here define quantities used for emasures
	measures=zeros(Int(T/Tmeas)+2,5)
	Nc::Int=length(cfg_memb)/4 
	cnt_meas::Int=1
	phi_field=zeros(Float64,Nc)
	phi_av::Float64=0.0
	binder_cumulant::Float64=0.0
	#here define other quantities used for simulation
	t::Float64=0.0
	last_meas::Float64=0.0
	#setup probabilities of events occuring in a specific site (note that THESE ARE NOT RATES ACTUALLY, this will
	#be used conditioned to the fact that that kind of event will occur, just to know at which site
	R_diffA=ExponentialQueue(Nc) #diffusivity of specie A
	R_diffB=ExponentialQueue(Nc) #diffusivity of specie B
	R_catA=ExponentialQueue(Nc)	#catalysis B+EA->A+EA
	R_catB=ExponentialQueue(Nc) #catalysis A+EB->B+EB
	#attachment of EA: use the same of diffusion of A, they are both proprtional to [A]
	R_detEA=ExponentialQueue(Nc) #detachment of EA
	#attachment of EB: use the same of diffusion of B, they are both proprtional to [B]
	R_detEB=ExponentialQueue(Nc) #detachment of EB
	#for diffusion of enzymes EA  (or EB) use the same of detachment of enymes, they're both proportional to [EA] (or [EB])
	for i in 1:Nc
		R_diffA[i]=cfg_memb[i,1]		
		R_diffB[i]=cfg_memb[i,2]		
		R_catA[i]=cfg_memb[i,3].*(cfg_memb[i,2]./(cfg_memb[i,2].+KMM))		
		R_catB[i]=cfg_memb[i,4].*(cfg_memb[i,1]./(cfg_memb[i,1].+KMM))
		R_detEA[i]=cfg_memb[i,3]
		R_detEB[i]=cfg_memb[i,4]
	end
	#now setup actual rates, for "global" events (e.g. "rate for any diffusion of specie A taking place in whole membrane")
	R=ExponentialQueue(10)
	R[1]=dA*R_diffA.sum[end]  #calling R_diffA.sum returns cumsum of the accumulator, the last element is thus the total sum
	R[2]=dB*R_diffB.sum[end]
	R[3]=kAc*R_catA.sum[end]
	R[4]=kBc*R_catB.sum[end]
	R[5]=kAa*cfg_cyto[1]*R_diffA.sum[end]
	R[6]=kAd*R_detEA.sum[end]
	R[7]=kBa*cfg_cyto[2]*R_diffB.sum[end]
	R[8]=kBd*R_detEB.sum[end]
	R[9]=dEA*R_detEA.sum[end]
	R[10]=dEB*R_detEB.sum[end]
	#plot()
	cnt_save::Int=1
	while t <T
		#extract an event
		i,dt=peek(R, rng=ran_ng)
		if t+dt>=T
			break
		end
		#make measures here
		while last_meas<t+dt
			#make measure  here
			#in some denominators I addee + 10^-10 in roder to avoid NaNs
			totA=sum(cfg_memb[:,1])
			totB=sum(cfg_memb[:,2])
			#calcualte phi at each cell (normalized, goes from -1, to 1)
			phi_field.=((cfg_memb[:,2].-cfg_memb[:,1])./(cfg_memb[:,1].+cfg_memb[:,2] .+ 10^-10))
			phi_av=(totB-totA)/(totA+totB)
			println("T=",last_meas, " and <phi>/c=", phi_av)
			#measure phi
			measures[cnt_meas,1]=phi_av
			#measure Binder cumulant
			binder_cumulant=1-(1/3)*((sum((phi_field .-phi_av).^4)/Nc)/((sum((phi_field .-phi_av).^2)/Nc)^2 + 10^-10))
			measures[cnt_meas,2]=binder_cumulant
			#measure enzymes in cytosolic reservoir
			measures[cnt_meas,3]=cfg_cyto[1]
			measures[cnt_meas,4]=cfg_cyto[2]
			println("Binder cumulant: ", measures[cnt_meas,2])			
			#this saves configuration every Nsave measures
			if cnt_save==Nsave
				save(fold_files*"config_T="*string(round(last_meas,digits=3))*".jld",compress=true, "membrane", cfg_memb, "cytosol", cfg_cyto)
				cnt_save=0
			end
			last_meas+=Tmeas
			cnt_save+=1
			#the following is to check for convergence of rho to 1
			rho=rho_0*((kAd/kAa)+totA)/((kBd/kBa)+totB)
			measures[cnt_meas,5]=rho
			cnt_meas+=1
			#this saves measures
			save(fold_files*"measures.jld","phi",measures[1:cnt_meas,1],"binder_cumulant", measures[1:cnt_meas,2],"EA_cyto",measures[1:cnt_meas,3],"EB_cyto",measures[1:cnt_meas,4],"rho",measures[1:cnt_meas,5])
			#FOLLOWING IS TO PLOT, COMMENT FOR CLUSTER USAGE, REQUIRES LIBRARIES Plots, ColorSchemes, LaTeXStrings
			#display(scatter(centers[:,1],centers[:,2],clim=(-1,1),zcolor=(cfg_memb[:,2].-cfg_memb[:,1])./(cfg_memb[:,1].+cfg_memb[:,2]),color=:vik,title="T="*string(round(last_meas,digits=3)),markerstrokewidths=0.0,markershape=:hexagon,markersize=4.0,label=L"\varphi",aspect_ratio=:equal))
			#to plot rho
			#scatter(collect(1:cnt_meas).*Tmeas,measures[1:cnt_meas,5])
			#savefig(string(last_meas)*".png")
			#lattice 100x100 -> markersize=1.9
			#lattice 50x50 -> markersize=4.0
			#lattice 20x20 ->markersize=10.6
		end		
		#make event happen (i.e. sample the site where event will occur and change number of molecules and update rates)
		if i==1
			#diffusion of specie A
			dep = peekevent(R_diffA, rng=ran_ng)
			arr=rand(ran_ng,neig[dep,:]) #arrival is chosen uniformly between its neighbours
			cfg_memb[dep,1]-=1
			cfg_memb[arr,1]+=1
			#update rates for diffusion of specie A, and for catalysis and for enzyme detachment (both  global and local)
			R_diffA[dep]=cfg_memb[dep,1]
			R_diffA[arr]=cfg_memb[arr,1]
			R_catB[dep]=cfg_memb[dep,4]*cfg_memb[dep,1]/(cfg_memb[dep,1]+KMM)
			R_catB[arr]=cfg_memb[arr,4]*cfg_memb[arr,1]/(cfg_memb[arr,1]+KMM)
			R[1]=dA*R_diffA.sum[end]
			R[4]=kBc*R_catB.sum[end]
			R[5]=kAa*cfg_cyto[1]*R_diffA.sum[end]
		elseif i==2
			#diffusion of specie B
			dep = peekevent(R_diffB, rng=ran_ng)
			arr=rand(ran_ng,neig[dep,:]) #arrival is chosen uniformly between its neighbours
			cfg_memb[dep,2]-=1
			cfg_memb[arr,2]+=1
			#update rates for diffusion of specie B, and for catalysis and for enzyme detachment (both  global and local)
			R_diffB[dep]=cfg_memb[dep,2]
			R_diffB[arr]=cfg_memb[arr,2]
			R_catA[dep]=cfg_memb[dep,3]*cfg_memb[dep,2]/(cfg_memb[dep,2]+KMM)
			R_catA[arr]=cfg_memb[arr,3]*cfg_memb[arr,2]/(cfg_memb[arr,2]+KMM)
			R[2]=dB*R_diffB.sum[end]
			R[3]=kAc*R_catA.sum[end]
			R[7]=kBa*cfg_cyto[2]*R_diffB.sum[end]
		elseif i==3
			#B+EA->A+EA
			site = peekevent(R_catA, rng=ran_ng)
			cfg_memb[site,1]+=1
			cfg_memb[site,2]-=1
			#update rates for diffusivity of specie A and B, catalysis of both species and attachment of enzymes
			R_diffA[site]=cfg_memb[site,1]
			R_diffB[site]=cfg_memb[site,2]
			R_catA[site]=cfg_memb[site,3]*cfg_memb[site,2]/(cfg_memb[site,2]+KMM)
			R_catB[site]=cfg_memb[site,4]*cfg_memb[site,1]/(cfg_memb[site,1]+KMM)
			R[1]=dA*R_diffA.sum[end]  
			R[2]=dB*R_diffB.sum[end]
			R[3]=kAc*R_catA.sum[end]
			R[4]=kBc*R_catB.sum[end]
			R[5]=kAa*cfg_cyto[1]*R_diffA.sum[end]
			R[7]=kBa*cfg_cyto[2]*R_diffB.sum[end]
		elseif i==4
			#A+EB->B+EB
			site = peekevent(R_catB, rng=ran_ng)
			cfg_memb[site,2]+=1
			cfg_memb[site,1]-=1
			#update rates for diffusivity of specie A and B, catalysis of both species and attachment of enzymes
			R_diffA[site]=cfg_memb[site,1]
			R_diffB[site]=cfg_memb[site,2]
			R_catA[site]=cfg_memb[site,3]*cfg_memb[site,2]/(cfg_memb[site,2]+KMM)
			R_catB[site]=cfg_memb[site,4]*cfg_memb[site,1]/(cfg_memb[site,1]+KMM)
			R[1]=dA*R_diffA.sum[end]  
			R[2]=dB*R_diffB.sum[end]
			R[3]=kAc*R_catA.sum[end]
			R[4]=kBc*R_catB.sum[end]
			R[5]=kAa*cfg_cyto[1]*R_diffA.sum[end]
			R[7]=kBa*cfg_cyto[2]*R_diffB.sum[end]
		elseif i==5
			#attachment of EA from cytosol
			site = peekevent(R_diffA, rng=ran_ng)
			cfg_memb[site,3]+=1
			cfg_cyto[1]-=1
			#set again rates for catalysis EA, detachment of EA, attachment of EA, also diffusion of EA
			R_catA[site]=cfg_memb[site,3]*cfg_memb[site,2]/(cfg_memb[site,2]+KMM)
			R_detEA[site]=cfg_memb[site,3]
			R[1]=dA*R_diffA.sum[end]
			R[3]=kAc*R_catA.sum[end]
			R[5]=kAa*cfg_cyto[1]*R_diffA.sum[end]
			R[6]=kAd*R_detEA.sum[end]
			R[9]=dEA*R_detEA.sum[end]
		elseif i==6
			#detachment of EA
			site = peekevent(R_detEA, rng=ran_ng)
			cfg_memb[site,3]-=1
			cfg_cyto[1]+=1
			#set again rates for catalysis EA, detachment of EA, attachment of EA, also diffusion of EA
			R_catA[site]=cfg_memb[site,3]*cfg_memb[site,2]/(cfg_memb[site,2]+KMM)
			R_detEA[site]=cfg_memb[site,3]
			R[1]=dA*R_diffA.sum[end]
			R[3]=kAc*R_catA.sum[end]
			R[5]=kAa*cfg_cyto[1]*R_diffA.sum[end]
			R[6]=kAd*R_detEA.sum[end]
			R[9]=dEA*R_detEA.sum[end]
		elseif i==7
			#attachment of EB
			site = peekevent(R_diffB, rng=ran_ng)
			cfg_memb[site,4]+=1
			cfg_cyto[2]-=1
			#set again rates for catalysis EB, detachment of EB, attachment of EB, also diffusion of EB
			R_catB[site]=cfg_memb[site,4]*cfg_memb[site,1]/(cfg_memb[site,1]+KMM)
			R_detEB[site]=cfg_memb[site,4]
			R[2]=dB*R_diffB.sum[end]
			R[4]=kBc*R_catB.sum[end]
			R[7]=kBa*cfg_cyto[2]*R_diffB.sum[end]
			R[8]=kBd*R_detEB.sum[end]
			R[10]=dEB*R_detEB.sum[end]
		elseif i==8
			#detachment of EB
			site = peekevent(R_detEB, rng=ran_ng)
			cfg_memb[site,4]-=1
			cfg_cyto[2]+=1
			#set again rates for catalysis EB, detachment of EB, attachment of EB, also diffusion of EB
			R_catB[site]=cfg_memb[site,4]*cfg_memb[site,1]/(cfg_memb[site,1]+KMM)
			R_detEB[site]=cfg_memb[site,4]
			R[2]=dB*R_diffB.sum[end]
			R[4]=kBc*R_catB.sum[end]
			R[7]=kBa*cfg_cyto[2]*R_diffB.sum[end]
			R[8]=kBd*R_detEB.sum[end]
			R[10]=dEB*R_detEB.sum[end]
		elseif i==9
			#diffusion of EA
			dep = peekevent(R_detEA, rng=ran_ng)
			arr=rand(ran_ng,neig[dep,:]) #arrival is chosen uniformly between its neighbours
			cfg_memb[dep,3]-=1
			cfg_memb[arr,3]+=1
			#update rates for detachment of EA locally, and for catalysis and for enzyme diffusion (both  global and local)
			R_detEA[dep]=cfg_memb[dep,3]
			R_detEA[arr]=cfg_memb[arr,3]
			R_catA[dep]=cfg_memb[dep,3]*cfg_memb[dep,2]/(cfg_memb[dep,2]+KMM)
			R_catA[arr]=cfg_memb[arr,3]*cfg_memb[arr,2]/(cfg_memb[arr,2]+KMM)
			R[9]=dEA*R_detEA.sum[end]
			R[3]=kAc*R_catA.sum[end]
			R[6]=kAd*R_detEA.sum[end]
		elseif i==10
			#diffusion of EB
			dep = peekevent(R_detEB, rng=ran_ng)
			arr=rand(ran_ng,neig[dep,:]) #arrival is chosen uniformly between its neighbours
			cfg_memb[dep,4]-=1
			cfg_memb[arr,4]+=1
			#update rates for detachment of EA locally, and for catalysis and for enzyme diffusion (both  global and local)
			R_detEB[dep]=cfg_memb[dep,4]
			R_detEB[arr]=cfg_memb[arr,4]
			R_catB[dep]=cfg_memb[dep,4]*cfg_memb[dep,1]/(cfg_memb[dep,1]+KMM)
			R_catB[arr]=cfg_memb[arr,4]*cfg_memb[arr,1]/(cfg_memb[arr,1]+KMM)
			R[10]=dEB*R_detEB.sum[end]
			R[4]=kBc*R_catB.sum[end]
			R[8]=kBd*R_detEB.sum[end]
		end
		#advance time of dt
		t+=dt
	end
	return cfg_memb, cfg_cyto
end
