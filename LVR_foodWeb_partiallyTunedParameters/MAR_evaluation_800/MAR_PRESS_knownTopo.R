#########################################################################################################
##### FBarraquand 22/08/2017 - Extraction of fitted MAR parameters and comparison to simulated models ###
##### Lotka-Volterra-Ricker food web models with G. Certain.  #####################
#########################################################################################################

#### Compare PRESS predictions

rm(list=ls())
graphics.off()
set.seed(42)

nsites=100 #100
nspecies=12

PRESS=1

simul_path = "/home/frederic/Documents/MAR_modelling/resubmission/MAR_foodWeb/LVR_foodWeb_partiallyTunedParameters/LVR_simulation_800/"
# alpha matrix

Press_data=as.data.frame(read.csv(paste(simul_path,"DataFoodWeb_theorPRESS.csv",sep="")))
head(Press_data)

PRESS_pred=PRESS_true=matrix(NA,nrow=nsites,ncol=nspecies)

for (ksite in 1:nsites){ ### Loop on sites or repeats 
  
  ### Extraction of estimated and real (simulated) parameters
  
  ## Estimated parameters B
  estim_path = "/home/frederic/Documents/MAR_modelling/resubmission/MAR_foodWeb/LVR_foodWeb_partiallyTunedParameters/MAR_estimation_knownTopology_800/"
  # B matrix
  B=as.matrix(read.csv(paste(estim_path,"estimatedB/B_point",ksite,".csv",sep="")))
  B=B[,-1] #remove first column
  B # point estimate
  
  ### All that stuff may serve no purpose #########################################################
  # Confidence intervals for B
  B_upper=as.matrix(read.csv(paste(estim_path,"estimatedB/B_upper",ksite,".csv",sep="")))
  B_upper=B_upper[,-1] #remove first column
  B_upper # upper bound B CI
  B_lower=as.matrix(read.csv(paste(estim_path,"estimatedB/B_lower",ksite,".csv",sep="")))
  B_lower=B_lower[,-1] #remove first column
  B_lower # lower bound B CI
  ### Computations on B
  B_signif = B
  B_signif[(B_lower<0)&(B_upper>0)] = 0
  B_signif ### Matrix of statistically significant coefficients at a 5% level
  ### All that stuff may serve no purpose #########################################################
  
  ## Estimated parameters C 
  estim_path = "/home/frederic/Documents/MAR_modelling/resubmission/MAR_foodWeb/LVR_foodWeb_partiallyTunedParameters/MAR_estimation_knownTopology_800/"
  # C matrix
  C=as.matrix(read.csv(paste(estim_path,"estimatedC/C_point",ksite,".csv",sep="")))
  C=C[,-1] #remove first column
  C # point estimate
 
  ### All that stuff may serve no purpose #########################################################
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
  ### All that stuff may serve no purpose #########################################################
  
  ## Press prediction and comparison
  PRESS_pred[ksite,]<-c(solve(diag(12)-B)%*%(C*PRESS))     
  delta<-Press_data$DeltaLN[Press_data$Site==ksite] 
  PRESS_true[ksite,]<-delta  
}

pdf(file="PRESS_prediction_LVR_knownTopo.pdf",width=8,height=12)
par(mfrow=c(3,4),pty="s")
for (ks in 1:nspecies){
  plot(PRESS_true[,ks],PRESS_pred[,ks],xlab="True Delta ln(N)",ylab="MAR predicted Delta ln(N)",main=paste("Species = ",ks,sep=""))
  abline(0,1,lwd=2,col="black")
  }
dev.off()
  
