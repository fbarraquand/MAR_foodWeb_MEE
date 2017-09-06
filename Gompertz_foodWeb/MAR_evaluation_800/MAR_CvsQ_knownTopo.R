#########################################################################################################
##### FBarraquand 22/08/2017 - Extraction of fitted MAR parameters and comparison to simulated models ###
##### Lotka-Volterra-Ricker food web models with G. Certain.  #####################
#########################################################################################################

#### Check how Q is recovered by matrix C

rm(list=ls())
graphics.off()
set.seed(42)


nsites=100 #100


## Real, simulated parameters
simul_path = "/home/frederic/Documents/MAR_modelling/resubmission/MAR_foodWeb/Gompertz_foodWeb/Gompertz_simulation_800/"
# Q matrix
q=as.matrix(read.csv(paste(simul_path,"DataFoodWeb_envEffects.csv",sep="")))
q=q[,-1]

C_matrix=q #Initializing

for (ksite in 1:nsites){ ### Loop on sites or repeats 
  
  ### Extraction of estimated and real (simulated) parameters
  
  ## Estimated parameters
  estim_path = "/home/frederic/Documents/MAR_modelling/resubmission/MAR_foodWeb/Gompertz_foodWeb/MAR_estimation_knownTopology_800/"
  # C matrix
  C=as.matrix(read.csv(paste(estim_path,"estimatedC/C_point",ksite,".csv",sep="")))
  C=C[,-1] #remove first column
  C # point estimate
  nspecies=nrow(C)
  # Confidence intervals for C
  C_upper=as.matrix(read.csv(paste(estim_path,"estimatedC/C_upper",ksite,".csv",sep="")))
  C_upper=C_upper[,-1] #remove first column
  C_upper # upper bound C CI
  C_lower=as.matrix(read.csv(paste(estim_path,"estimatedC/C_lower",ksite,".csv",sep="")))
  C_lower=C_lower[,-1] #remove first column
  C_lower # lower bound C CI
  ### Computations on C
  C_signif = C
  C_signif[(C_lower<0)&(C_upper>0)] = 0
  C_signif ### Matrix of statistically significant coefficients at a 5% level
  
  C_matrix[ksite,] = C #C_signif
  
  ## Real values are lines of the Q matrix
  q[ksite,]
 
}

pdf(file="EnvEffect_knownTopo_LVR.pdf",width=8,height=8)
par(mfrow=c(2,2))
plot(q[,1],C_matrix[,1],xlab ="q_1",ylab = "c_1",xlim=c(-0.1,0.6),ylim=c(-0.1,0.6))
abline(0,1,lwd=2,col="black")
abline(lm(C_matrix[,1] ~q[,1]),lwd=2,col="blue",lty="dashed")
plot(q[,2],C_matrix[,2],xlab ="q_2",ylab = "c_2",xlim=c(-0.1,0.6),ylim=c(-0.1,0.6))
abline(0,1,lwd=2,col="black")
abline(lm(C_matrix[,1] ~q[,1]),lwd=2,col="blue",lty="dashed")
plot(q[,4],C_matrix[,4],xlab ="q_4",ylab = "c_4",xlim=c(-0.1,0.6),ylim=c(-0.1,0.6))
abline(0,1,lwd=2,col="black")
abline(lm(C_matrix[,1] ~q[,1]),lwd=2,col="blue",lty="dashed")
plot(q[,11],C_matrix[,11],xlab ="q_11",ylab = "c_11")
abline(0,1,lwd=2,col="black")
abline(lm(C_matrix[,1] ~q[,1]),lwd=2,col="blue",lty="dashed")
dev.off()
