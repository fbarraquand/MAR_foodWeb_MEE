################ MAR analyses of a simulated food web dynamics ##################
##### FBarraquand 18/07/2017, revised 22/08/2017 - Uses MARSS  ##################
##### Fits the model with known topology --------------------- ##################
##### Lotka-Volterra-Ricker food web models with G. Certain -- ##################
#################################################################################

rm(list=ls())
graphics.off()
set.seed(42)

### Estimation of an interaction matrix (biotic drivers) and environmental effects on PGR (abiotic)

### Loading useful packages
library('MARSS')
library('stringr')

## Utilitary function for extracting parameter values from a MAR model fit
getparam.fct<-function(m.ci,Bzero,Czero,ObsErr=T){ 
  #from GCertain -- works only if B and C are *not* specified as "uncontrained", which is the case here
  nsp<-ncol(Bzero)
  nenv<-ncol(Czero)
  
  Bfit<-array(NA,dim=c(nsp,nsp,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
  Bfit[,,1]<-Bzero
  Bfit[,,1][is.na(m.ci$call$model$B==0)]<-m.ci$par$B
  Bfit[,,2]<-Bzero
  Bfit[,,2][is.na(m.ci$call$model$B==0)]<-m.ci$par.lowCI$B
  Bfit[,,3]<-Bzero
  Bfit[,,3][is.na(m.ci$call$model$B==0)]<-m.ci$par.upCI$B
  Bfit[,,4]<-Bzero
  Bfit[,,4][is.na(m.ci$call$model$B==0)]<-m.ci$par.se$B
  
  Cfit<-array(NA,dim=c(nsp,nenv,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
  Cfit[,,1]<-Czero
  Cfit[,,1][is.na(m.ci$call$model$C==0)]<-m.ci$par$U
  Cfit[,,2]<-Czero
  Cfit[,,2][is.na(m.ci$call$model$C==0)]<-m.ci$par.lowCI$U
  Cfit[,,3]<-Czero
  Cfit[,,3][is.na(m.ci$call$model$C==0)]<-m.ci$par.upCI$U
  Cfit[,,4]<-Czero
  Cfit[,,4][is.na(m.ci$call$model$C==0)]<-m.ci$par.se$U
  
  Sigmafit<-array(NA,dim=c(nsp,nsp,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
  Sigmafit[,,1]<-diag(nsp)*c(m.ci$par$Q)
  Sigmafit[,,2]<-diag(nsp)*c(m.ci$par.lowCI$Q)
  Sigmafit[,,3]<-diag(nsp)*c(m.ci$par.upCI$Q)
  Sigmafit[,,4]<-diag(nsp)*c(m.ci$par.se$Q)
  
  fittedparam<-list(Bfit,Cfit,Sigmafit)
  
  if (ObsErr==T){
    ObsErrfit<-array(NA,dim=c(nsp,nsp,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
    ObsErrfit[,,1]<-diag(nsp)*c(m.ci$par$R)
    ObsErrfit[,,2]<-diag(nsp)*c(m.ci$par.lowCI$R)
    ObsErrfit[,,3]<-diag(nsp)*c(m.ci$par.upCI$R)
    ObsErrfit[,,4]<-diag(nsp)*c(m.ci$par.se$R)
    fittedparam<-list(Bfit,Cfit,Sigmafit,ObsErrfit)}
  
  return(fittedparam)
}


### Estimation of MAR/VAR matrix model, using MARSS

DBall=read.csv("/home/frederic/Documents/MAR_modelling/resubmission/MAR_foodWeb/Gompertz_foodWeb/Gompertz_simulation_800/DataFoodWeb_wTime_abs.csv") 
head(DBall)
nrow(DBall)

nsites=length(unique(DBall$Site)) #number of sites or repeats

################## Looping on the sites ####################
for (ksite in 1:nsites){

DB=DBall[DBall$Site==ksite,] ## Select a site
head(DB)
DB=DB[DB$Time_index %in% 1001:1800,]
head(DB) ## select 800 points between 1000 and 2000

abundance_mat=as.matrix(DB[,4:15]) ### Create matrix with time series of abundances
### From code with CP
sp=c("Species1","Species2","Species3","Species4","Species5","Species6","Species7","Species8","Species9","Species10","Species11","Species12")
nspecies=length(sp)
iter_min=100 #200 for Griffiths ; but it seems that about 30 are enough in our case
iter_estimate=10 #not used for now
iter_scenario=10
iter_boot=1000 #Ok for Griffiths
BFGS=FALSE

cov3_tot=c("Abiotic_var2")#,"Abiotic_var1","Abiotic_var3")

### Could use this with site instead of repeat
dates=DB$Time_index
dates_bis=dates #already interpolated
tab_sp=DB[,sp] # species abundance table
tab_cov=DB[,cov3_tot] # covariate table
tab_cov_bis=tab_cov   # new table keeping track of the original covariate, because we scale tab_cov afterwards

#Setting MARSS model -- known topology of the adjacency matrix
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

B1=matrix(list(0),nrow=length(sp),ncol=length(sp),dimnames=list(1:12,1:12))
for (i in 1:length(sp)){
  s=sp[i]
  for (j in 1:length(sp)){
    s2=sp[j]
    if(interaction_matrix[i,j]!=0){
      #B1[i,j]=paste("alpha_",i,",",j,sep="")
      B1[i,j]=paste(s2,s,sep="") ### beware, paste(s2,s,sep="") means we have the impacting species first  
    }else{
      if(s==s2){
        B1[i,j]=paste(s2,s,sep="")
      }
    }
  }
}

### Defining the environmental effect C
C1=matrix(list("c1","c2",0,"c4",0,0,0,0,0,0,"c11",0),nspecies,1) # vector because 1 env variable only here

### Defining other parameters
U1="zero"
Q1="diagonal and unequal"
Z1=diag(1,length(sp),length(sp))
A1="zero"
R1="zero"
V1=diag(1,length(sp))
pi1="zero"

aalpha=0.05
#cntl.list=list(conv.test.slope.tol=0.001,minit=iter_min,maxit=500,abstol=0.001)#,MCInit=TRUE,silent=2)#,numInits=iter_estimate,numInitSteps=10) #Took values from Griffiths for minit,maxit,abstol. conv.test.slope is recommended in MARSS User Guide itself, numInits and numInitStep are taken as rules of thumbs (quick search on the Internet, numInits=500 elsewhere but it seems really big to me
cntl.list=list(conv.test.slope.tol=0.5,minit=iter_min,maxit=500,abstol=0.05)
#Without season and with EM algorithm

tab_sp=log(tab_sp)
tab_sp=t(scale(tab_sp, scale = FALSE)) ### log and center the abundance data
tab_cov=t(scale(tab_cov_bis, scale = FALSE)) ### same for the covariates // centered but not scaled here for comparison to Gregoire's code
### Covariate means are already very close to zero though. 
### For the PRESS prediction, we take into account those very small differences in covariate value.
### For the short-term forecasts, it is best to compare the prediction of the MAR with centered covar
rownames(tab_sp)=sp
rownames(tab_cov)=cov3_tot
c1=tab_cov
model.list=list(B=B1,U=U1,Q=Q1,Z=Z1,A=A1,R=R1,x0=pi1,V0=V1,tinitx=1,C=C1,c=c1) ### Specifying the parameters to MARSS
fit_log=MARSS(tab_sp, method="kem",model=model.list,control=cntl.list) ### Fit of MAR model. 
fit_log$par$B ## 

### Need to add the confidence intervals absolutely... 
cis=MARSSparamCIs(fit_log,method="hessian",alpha=aalpha)
### This is freaking slow even with the Hessian, damn. 

### Gregoire's code for extraction
Bzero=matrix(0,nrow=12,ncol=12)
Czero=matrix(0,12,1)
fittedparam<-getparam.fct(cis,Bzero=Bzero,Czero=Czero,ObsErr = F)

################# The community matrix B ###########################################
Bfit.stock<- fittedparam[[1]][,,1]
Bfit.stock.low <- fittedparam[[1]][,,2]
Bfit.stock.up  <- fittedparam[[1]][,,3]
Bfit.stock.se  <- fittedparam[[1]][,,4]
### Send to file
write.csv(Bfit.stock,file=paste("estimatedB/B_point",ksite,".csv",sep=""))
write.csv(Bfit.stock.low,file=paste("estimatedB/B_lower",ksite,".csv",sep=""))
write.csv(Bfit.stock.up,file=paste("estimatedB/B_upper",ksite,".csv",sep=""))
write.csv(Bfit.stock.se,file=paste("estimatedB/B_SE",ksite,".csv",sep=""))
####################################################################################

############### The environmental effects on growth rates C ##############################
Cfit.stock    <- fittedparam[[2]][,,1]
Cfit.stock.low <- fittedparam[[2]][,,2]
Cfit.stock.up  <- fittedparam[[2]][,,3]
Cfit.stock.se  <- fittedparam[[2]][,,4]

# Code to export the C matrix 
write.csv(Cfit.stock,file=paste("estimatedC/C_point",ksite,".csv",sep=""))
write.csv(Cfit.stock.low,file=paste("estimatedC/C_lower",ksite,".csv",sep=""))
write.csv(Cfit.stock.up,file=paste("estimatedC/C_upper",ksite,".csv",sep=""))
write.csv(Cfit.stock.se,file=paste("estimatedC/C_SE",ksite,".csv",sep=""))
#################################################################################################


#################### Process error variances ####################################################
Omegafit.stock   <- fittedparam[[3]][,,1]
Omegafit.stock.low <- fittedparam[[3]][,,2]
Omegafit.stock.up <- fittedparam[[3]][,,3]
Omegafit.stock.se <- fittedparam[[3]][,,4]

# Code to export Omega
write.csv(Omegafit.stock,file=paste("Omega/Omega_point",ksite,".csv",sep=""))
write.csv(Omegafit.stock.low,file=paste("Omega/Omega_lower",ksite,".csv",sep=""))
write.csv(Omegafit.stock.up,file=paste("Omega/Omega_upper",ksite,".csv",sep=""))
write.csv(Omegafit.stock.se,file=paste("Omega/Omega_SE",ksite,".csv",sep=""))
####################################################################################################

} #end of loop on sites

#################################################################################################





