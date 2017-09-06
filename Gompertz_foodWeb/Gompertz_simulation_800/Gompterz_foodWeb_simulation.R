###### ------------------------------------------------------------------------------------------------------------------ ########
###### --- FBarraquand 29/08/2017 // Compertz food web dyn, taking the Jacobian of the previous LVR models as base ------ ########
###### --- FBarraquand 18/04/2017 // Lotka-Volterra-Ricker models (food web here) with G. Certain and A. Gaardmark ------ ########
###### --- Adapted from a previous code of mine for competitive systems with various sites ------------------------------ ########
###### --- Revised and adapted 22/08/2017 to better match the simulation experiment of Certain et al. 2017's main text -- ########
###### ------------------------------------------------------------------------------------------------------------------ ########

rm(list=ls())
graphics.off()
set.seed(42)
setwd('/home/frederic/Documents/MAR_modelling/resubmission/MAR_foodWeb/Gompertz_foodWeb/Gompertz_simulation_800')
getwd()

tmax=2100 ## Number of time units. 2000 timesteps total for the simulation, and we select 1000 time steps for model fitting

nsites = 100 ## Could be sites or repeats, depending on whether parameters can differ between sites
### or are identical between sites (the default here, which makes nsites repeats of the same food web)

###### ----------------- Definition of an interaction (adjacency) matrix ------------------------------------------------- ########

### This will define number of species as a by-product
interaction_matrix = rbind(c(1,0,1,0,0,0,0,0,0,0,0,0),
                            c(0,1,0,0,0,0,0,1,0,0,0,0),
                            c(1,0,1,0,1,1,0,0,0,0,0,0),
                            c(0,0,0,1,1,1,1,0,0,0,0,0),
                            c(0,0,1,1,1,0,0,0,1,1,1,0),
                           c(0,0,1,1,0,1,0,0,0,0,0,1),
                           c(0,0,0,1,0,0,1,0,0,0,0,1),
                           c(0,1,0,0,0,0,0,1,0,0,0,1),
                           c(0,0,0,0,1,0,0,0,1,1,0,0),
                           c(0,0,0,0,1,0,0,0,1,1,0,0),
                           c(0,0,0,0,1,0,0,0,0,0,1,0),
                           c(0,0,0,0,0,1,1,1,0,0,0,1))
### Gatun lake interaction network
### Source: Aufderheide, H., Rudolf, L., Gross, T. & Lafferty, K.D. (2013) How to predict community responses to perturbations in the face of imperfect knowledge and network complexity. Proceedings of the Royal Society B: Biological Sciences 280, 20132355.

## Plot this
interaction_matrix
nspecies=nrow(interaction_matrix)
write.csv(interaction_matrix,file="DataFoodWeb_qualitativeMatrix.csv")

### Ecological interaction strength parameters
# defined later - we use a Jacobian of an LVR model

###### ----------------- Environmental covariables simulation -------------------------------------------------------------------- ########

### Variable 1, temperature like - variation between sites around a joint trend
# y1noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=tmax,sd=sqrt(0.05) )
## Remark FB // Adapted from code for simulated plankton dynamics 
y1site <-matrix(0, nrow=tmax,ncol=nsites,byrow=TRUE)
# sigma_y1<-0.05
# for (ksite in 1:nsites){ y1site[,ksite] <- y1noise + rnorm(tmax, mean = 0, sd = sigma_y1)}
## Output this data
# y1site
# matplot(y1site)
# matlines(y1site)
### Here's what is added to r 
colMeans(y1site)

### Variable 2, rainfall-like - variation in time differing between the k sites
# Based on previous model 
y2site<-matrix(0, nrow=tmax,ncol=nsites,byrow=TRUE)
for (ksite in 1:nsites){ 
  y2noise<-arima.sim(model=list(ar=c(0.1, 0.2, 0.1,0.5,-0.1)), n=tmax,sd=sqrt(0.1) ) # more noisy than first because no subsequent sampling
  plot(1:tmax,y2noise,type="o")
  y2site[,ksite] <- y2noise # changes for every site
}
y2site
matplot(y2site)
matlines(y2site)
colMeans(y2site)

### Variable 3, site "quality" - fixed variation between sites
## sigma_y3= 0.1 ## Low because we are using log-normal noise here. 
y3site <- rep(0,nsites) #rnorm(nsites, mean = 0, sd = sigma_y3)
### hist(y3site)
### Currently set to zero so that all sites are repeats. 

###### ------------------ Definition of species affected by environmental effects ------------------------------------------------ ########
# We decide that 1/3 of the species are affected by the forcing variable(s) 
# (to be similar to the 2-species experiment without exagerating complexity)
# Species 1, 2, 4 and 11 (some at the base, middle, and top of the food web)
environmental_forcing = c(1,1,0,1,0,0,0,0,0,0,1,0)

####### Variable covariate effects across sites or repeats
## Scaling of covariate effects (multiplier)
# scaling_covar = 0.2 #0.1 #0 
scaling_covar=1

# Initialize
q=matrix(0,nsites,nspecies) # the effect of the env variable(s)
local_cov=matrix(0,nsites,nspecies) # what is added to the growth rate considering static env. variables
N=matrix(0,nsites,nspecies) # matrix of site-specific abundances (for later)

for (ksite in 1:nsites)
{
  environmental_forcing_strength = runif(12,0.01,0.5)
  q[ksite,]=environmental_forcing_strength*environmental_forcing ## element-by-element multiplication so there's only the forced species that are affected
  q[ksite,11] = 0.001 # species 11 is a top predator and its mortality is only mildly modulated. 
  local_cov[ksite,] = scaling_covar*q[ksite,]*(mean(y1site[,ksite])+mean(y2site[,ksite])+y3site[ksite]) 
  ## Mind that q is varying across sites
  }
write.csv(q,file="DataFoodWeb_envEffects.csv")


###### ----------------- Loop over sites/repeats -------------------------------------------------------------------- ########
ksite=1
while (ksite <= nsites){

    print(paste("site/repeat =",ksite))
  
    ## FB 24/08/2017 Re-arranged the logical order there
    ## Initializing the conditions for feasibility
    PressFeasible<-F
    IsFeasible<-F
    TryPress=0
    
    while(PressFeasible==F){## PressFeasible checks that the eq. can be perturbed and still be feasible
      TryPress=TryPress+1
      print(paste("TryPress ",TryPress))
      if (TryPress>5){
        IsFeasible<-F ## Re-drawing a parameter set at random.  
        TryPress=1 ## Counter to do the press back to 1
        print(paste("TryPress ",TryPress))
        } 
      ## If the current feasible point cannot be perturbed than it is not feasible in a practical sense.. 
      ## IsFeasible checks that the multispecies equilibrium point is feasible (N*>0). 
      TryFeasible=0
      while(IsFeasible==F){
      TryFeasible=TryFeasible + 1
      print(paste("TryFeasible",TryFeasible))
      
      ### Definition of growth rates - we still need those here even though they don't have the same meaning
      r_prey<-c(2,1,1) #we can't put this at random // species 1 needs very high growth rate and species 2 cannot be too productive
      # Setting the first intrinsic GR to 3 works well but produces boom-bust dynamics with sometimes ln(density) = -40 ...
      # works better if the two other prey have equal GRs
      r_middle_predator<-runif(6,0.01,0.05) #let's say they reproduce poorly on smthg else
      r_top_predator<-c(-0.001,-0.0001,-0.001) #runif(3,-0.01,0) ### much better
      r<-c(r_prey[1],r_prey[2],r_middle_predator[1],r_prey[3],r_middle_predator[2:6],r_top_predator[1:3])
      write.csv(r,file=paste("intrinsicGR/DataFoodWeb_intrinsicGrowthRates",ksite,".csv",sep=""))

      ### Interaction parameters
      B_path = "/home/frederic/Documents/MAR_modelling/resubmission/MAR_foodWeb/LVR_foodWeb_partiallyTunedParameters/MAR_evaluation/Jacobian/"
      # B matrix
      B=as.matrix(read.csv(paste(B_path,"Jacobian",ksite,".csv",sep="")))
      B=B[,-1]
      write.csv(B,file=paste("alpha/DataFoodWeb_alphaMatrix",ksite,".csv",sep="")) #Not the same alpha but we call it like that in the text
      if (max(abs(eigen(B)$values))<1){IsFeasible<-T} 

      } #end of condition on IsFeasible
 

### Noise level for the dynamics. 
sigma=0.3 # So that Sigma^2= 0.1

#### ------------------------------------------------ Computing the dynamics ---------------------------------------------------- #########

index_time=1:tmax
  
  ### Initialize the abundances matrix
  Y=matrix(1,nrow=tmax,ncol=nspecies)
  y=matrix(1,nrow=tmax,ncol=nspecies)
  Z=matrix(1,nrow=tmax,ncol=nspecies)
  epsilon<-matrix(0, nrow=tmax-1,ncol=nspecies,byrow=TRUE) # Noise on growth rates
  
  # Initial abundances
  Y[1,]=abs(rnorm(nspecies,1,1)) # abundance vector 
  y[1,]=log(Y[1,]) # natural log-abundance vector
  Z[1,]=Y[1,]/sum(Y[1,]) # relative abundance vector
  
  # Loop over time
  for (t in 1:(tmax-1)){
    epsilon[t,]<-rnorm(nspecies,mean = 0, sd = sigma) # Residual noise on log-population growth rates
    ycovar =   q[ksite,] *(y1site[t+1,ksite] + y2site[t+1,ksite] + y3site[ksite]) # q is a vector of environmental effects of the covariates that is multiplied by the scalar value of covariates at time t
    y[t+1,] = r + scaling_covar*ycovar+epsilon[t,] + B %*%y[t,] # Gompertz dynamics - mind that r does not have the same meaning
    Y[t+1,] = exp(y[t+1,])
    Z[t+1,]=Y[t+1,]/sum(Y[t+1,]) # relative abundance vector
  } # end of loop in time 
 
  matplot(y)
  matlines(y)
  matplot(Z)
  matlines(Z)
  
  ### Local growth rates
  local_GR=r+local_cov[ksite,]
  N[ksite,]=exp(solve(B) %*% (-local_GR))
  ## warning - it might be that small differences there render the abundances negative

	### Theoretical computation of the PRESS perturbation
	PRESS=1 ## More generates non-feasible equilibria
	if (det(diag(nspecies)-B)!=0){PressFeasible<-T} #Always doable if B-I invertible
    } #end of loop on Press Feasible
	LDiff= solve(diag(nspecies)-B) %*% (q[ksite,]*PRESS) 
	N_star_bis = exp(LDiff+log(N[ksite,]))
	Diff=N_star_bis-N[ksite,] #N[ksite,] is a little more precise than N_star but should not change much

  ### Fill temporary data structures for each site
  DataFoodWebTemp_wTime_abs = data.frame(ksite*rep(1,tmax),index_time,Y,as.numeric(y1site[,ksite]),as.numeric(y2site[,ksite]),y3site[ksite]*rep(1,tmax))
  # Compute empirical means over time (i.e., temporal average of abundances on a log-scale, back-transformed to regular scale)
  #Y_mean=exp(colMeans(y))
  #Z_mean=Y_mean/sum(Y_mean)
  
  # Compute equilibrium model (theoretic) derived values
  Y_theor=N[ksite,]
  #Z_theor=Y_theor/sum(Y_theor)
  
  # Fill data structures
  DataFoodWebTemp_theoreticMeans_abs = data.frame(ksite,Y_theor[1],Y_theor[2],Y_theor[3],Y_theor[4],Y_theor[5],Y_theor[6],Y_theor[7],Y_theor[8],Y_theor[9],Y_theor[10],Y_theor[11],Y_theor[12],as.numeric(mean(y1site[,ksite])),as.numeric(y3site[ksite]))
  DataFoodWebTemp_theorPRESS = data.frame(1:nspecies,ksite*rep(1,nspecies),N[ksite,],N_star_bis,log(N[ksite,]),log(N_star_bis),Diff,LDiff)
  
   ### Write this in total dataframes
  if (ksite==1){ 	# Create data structures
    DataFoodWeb_wTime_abs=DataFoodWebTemp_wTime_abs
    DataFoodWeb_theoreticMeans_abs=DataFoodWebTemp_theoreticMeans_abs
    DataFoodWeb_theorPRESS = DataFoodWebTemp_theorPRESS

    } else {
      # Add to data structures
      DataFoodWeb_wTime_abs=rbind(DataFoodWeb_wTime_abs,DataFoodWebTemp_wTime_abs)
      DataFoodWeb_theoreticMeans_abs=rbind(DataFoodWeb_theoreticMeans_abs,DataFoodWebTemp_theoreticMeans_abs)
      DataFoodWeb_theorPRESS =rbind(DataFoodWeb_theorPRESS, DataFoodWebTemp_theorPRESS)
    }
    ksite=ksite+1
  
}# end of loop on sites
      
   
### Name variables
names(DataFoodWeb_wTime_abs)=c("Site","Time_index","Species1","Species2","Species3","Species4","Species5","Species6","Species7","Species8","Species9","Species10","Species11","Species12","Abiotic_var1","Abiotic_var2","Abiotic_var3")
names(DataFoodWeb_theoreticMeans_abs)=c("Site","Species1","Species2","Species3","Species4","Species5","Species6","Species7","Species8","Species9","Species10","Species11","Species12","Abiotic_var1","Abiotic_var3")
names(DataFoodWeb_theorPRESS)=c("Species","Site","Nstar","Nstar_prime","LNstar","LNstar_prime","DeltaN","DeltaLN")

### Output to file
write.csv(DataFoodWeb_wTime_abs,file="DataFoodWeb_wTime_abs.csv")
write.csv(DataFoodWeb_theoreticMeans_abs,file="DataFoodWebTemp_theoreticMeans_abs.csv")
write.csv(DataFoodWeb_theorPRESS,file="DataFoodWeb_theorPRESS.csv")

par(mfrow=c(1,1))

### Look at last site
pdf(file="CommunityDynamics.pdf",width=15,height=7)
par(cex=1.5)
matplot(Y[800:900,],ylab="True abundance",xlab="Time")
matlines(Y[800:900,])
dev.off()

pdf(file="CommunityDynamics_logScale.pdf",width=15,height=7)
par(cex=1.5)
matplot(y[800:900,],ylab="log(abundance)",xlab="Time")
matlines(y[800:900,])
dev.off()



