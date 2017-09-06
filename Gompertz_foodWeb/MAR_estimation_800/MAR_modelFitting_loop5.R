################ MAR analyses of a simulated food web dynamics ##################
##### FBarraquand 18/07/2017, revised 22/08/2017 - Uses MARSS  ##################
##### Fits the full, unconstrained model --------------------- ##################
##### Lotka-Volterra-Ricker food web models with G. Certain -- ##################
#################################################################################

rm(list=ls())
graphics.off()
set.seed(42)

### Estimation of an interaction matrix (biotic drivers) and environmental effects on PGR (abiotic)

### Loading useful packages
library('MARSS')
library('stringr')

### Estimation of MAR/VAR matrix model, using MARSS 

### Modified from plankton analyses with CP ###

DBall=read.csv("/home/frederic/Documents/MAR_modelling/resubmission/MAR_foodWeb/Gompertz_foodWeb/Gompertz_simulation_800/DataFoodWeb_wTime_abs.csv") 
head(DBall)
nrow(DBall)

nsites=length(unique(DBall$Site)) #number of sites or repeats

################## Looping on the sites ####################
for (ksite in 81:100){

DB=DBall[DBall$Site==ksite,] ## Select a site
head(DB)
DB=DB[DB$Time_index %in% 1001:1800,]
head(DB) ## 

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


#Setting MARSS model -- unconstrained matrix
B1="unconstrained"
U1="zero"
Q1="diagonal and unequal"
Z1=diag(1,length(sp),length(sp))
A1="zero"
R1="zero"
V1=diag(1,length(sp))
pi1="zero"
C1="unconstrained"
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

###################################################################################################################################
#### Code to extract the B matrix -- adapted from doi: https://doi.org/10.1101/171264 (a little bit complex for what it does here)
####################################################################################################################################
Bnew=Blower=Bupper=B_SE=matrix(0,length(sp),length(sp))# 
nom=dimnames(fit_log$par$B)[[1]] #B matrix coefficients

for (n in 1:length(nom)){
  if(grepl("^\\(",nom[n])){ #default when using unconstrained or diagonal #In this case, it begins with a parenthesis
    i=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][2]) #Coefficients are called (i,j)
    j=as.numeric(strsplit(nom[n],split="[\\(,\\)]")[[1]][3])
    if(is.na(i)){ #i and j might not be numeric
      i=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][2]))))
      j=which(unlist(lapply(sp,function(x) grepl(x,strsplit(nom[n],split="[\\(,\\)]")[[1]][3]))))
    }
  }else{ #When I defined my own matrix, I called the effect of species X on species Y "XY" ; but this means X in column, Y in row
    a=str_split(nom[n],sp) #remove the sp, you have a list of length S, with 1,2 or 3 values ("XY" when nothing, "" "Y" when you have index Y, "" "" "" when the interaction is "XX" and you have index X)
    for (abis in 1:length(a)){
      if(length(a[[abis]])==2){
        if((a[[abis]][1]=="")&&(a[[abis]][2]!="")){
          j=abis
        }
        if((a[[abis]][2]=="")&&(a[[abis]][1]!="")){
          i=abis
        }
      }else if(length(a[[abis]])==3){ 
        if((a[[abis]][2]=="")&&(a[[abis]][1]=="")&&(a[[abis]][3]=="")){
          j=abis
          i=abis
        }
      }
    }
  }
  Bnew[i,j]=fit_log$par$B[n]
  Blower[i,j]=cis$par.lowCI$B[n]
  Bupper[i,j]=cis$par.upCI$B[n]
  B_SE[i,j]=cis$par.se$B[n]
}
############################################## end of B extraction loop #######################

### Plotting now
image(Bnew)
image(Blower)
image(Bupper)

### Send to file
write.csv(Bnew,file=paste("estimatedB/B_point",ksite,".csv",sep=""))
write.csv(Blower,file=paste("estimatedB/B_lower",ksite,".csv",sep=""))
write.csv(Bupper,file=paste("estimatedB/B_upper",ksite,".csv",sep=""))
write.csv(B_SE,file=paste("estimatedB/B_SE",ksite,".csv",sep=""))

################### Code to extract the C matrix ################################################
fit_log$par$C

#Covariate effects -- from Coralie's code (but it plots the whole thing too...)
Cnew<-Cupper<-Clower<-C_SE<-rep(NA,nspecies)
  for (i in 1:length(sp)){
      Cnew[i]=fit_log$par$U[i] ### Shouldn't this be fit_log$par$C[i]? 
      ### I know we specify C=C1 and U="zero", and we get back the params in U... the results are correct but weirdly presented
      Cupper[i]=cis$par.upCI$U[i]
      Clower[i]=cis$par.lowCI$U[i]
      C_SE[i]=cis$par.se$U[i]
    }
write.csv(Cnew,file=paste("estimatedC/C_point",ksite,".csv",sep=""))
write.csv(Clower,file=paste("estimatedC/C_lower",ksite,".csv",sep=""))
write.csv(Cupper,file=paste("estimatedC/C_upper",ksite,".csv",sep=""))
write.csv(C_SE,file=paste("estimatedC/C_SE",ksite,".csv",sep=""))
#################################################################################################

################## Code to extract Sigma (also called Omega) #################################
# In that case no problem with Greg's code. 
nsp=nspecies
m.ci=cis
Sigmafit<-array(NA,dim=c(nsp,nsp,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
Sigmafit[,,1]<-diag(nsp)*c(m.ci$par$Q)
Sigmafit[,,2]<-diag(nsp)*c(m.ci$par.lowCI$Q)
Sigmafit[,,3]<-diag(nsp)*c(m.ci$par.upCI$Q)
Sigmafit[,,4]<-diag(nsp)*c(m.ci$par.se$Q)
write.csv(Sigmafit[,,1],file=paste("Omega/Omega_point",ksite,".csv",sep=""))
write.csv(Sigmafit[,,2],file=paste("Omega/Omega_lower",ksite,".csv",sep=""))
write.csv(Sigmafit[,,3],file=paste("Omega/Omega_upper",ksite,".csv",sep=""))
write.csv(Sigmafit[,,4],file=paste("Omega/Omega_SE",ksite,".csv",sep=""))
} #end of loop on sites

#################################################################################################

### Utilitary function for extracting parameter values from a MAR model fit
# getparam.fct<-function(m.ci,Bzero,Czero,ObsErr=T){
#   nsp<-ncol(Bzero)
#   nenv<-ncol(Czero)
#   
#   Bfit<-array(NA,dim=c(nsp,nsp,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
#   Bfit[,,1]<-Bzero
#   Bfit[,,1][is.na(m.ci$call$model$B==0)]<-m.ci$par$B
#   Bfit[,,2]<-Bzero
#   Bfit[,,2][is.na(m.ci$call$model$B==0)]<-m.ci$par.lowCI$B
#   Bfit[,,3]<-Bzero
#   Bfit[,,3][is.na(m.ci$call$model$B==0)]<-m.ci$par.upCI$B
#   Bfit[,,4]<-Bzero
#   Bfit[,,4][is.na(m.ci$call$model$B==0)]<-m.ci$par.se$B
#   
#   Cfit<-array(NA,dim=c(nsp,nenv,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
#   Cfit[,,1]<-Czero
#   Cfit[,,1][is.na(m.ci$call$model$C==0)]<-m.ci$par$U
#   Cfit[,,2]<-Czero
#   Cfit[,,2][is.na(m.ci$call$model$C==0)]<-m.ci$par.lowCI$U
#   Cfit[,,3]<-Czero
#   Cfit[,,3][is.na(m.ci$call$model$C==0)]<-m.ci$par.upCI$U
#   Cfit[,,4]<-Czero
#   Cfit[,,4][is.na(m.ci$call$model$C==0)]<-m.ci$par.se$U
#   
#   Sigmafit<-array(NA,dim=c(nsp,nsp,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
#   Sigmafit[,,1]<-diag(nsp)*c(m.ci$par$Q)
#   Sigmafit[,,2]<-diag(nsp)*c(m.ci$par.lowCI$Q)
#   Sigmafit[,,3]<-diag(nsp)*c(m.ci$par.upCI$Q)
#   Sigmafit[,,4]<-diag(nsp)*c(m.ci$par.se$Q)
#   
#   fittedparam<-list(Bfit,Cfit,Sigmafit)
#   
#   if (ObsErr==T){
#     ObsErrfit<-array(NA,dim=c(nsp,nsp,4),dimnames=list(NULL,NULL,c("fit","low","up","se")))
#     ObsErrfit[,,1]<-diag(nsp)*c(m.ci$par$R)
#     ObsErrfit[,,2]<-diag(nsp)*c(m.ci$par.lowCI$R)
#     ObsErrfit[,,3]<-diag(nsp)*c(m.ci$par.upCI$R)
#     ObsErrfit[,,4]<-diag(nsp)*c(m.ci$par.se$R)
#     fittedparam<-list(Bfit,Cfit,Sigmafit,ObsErrfit)}
#   
#   return(fittedparam)
# }




