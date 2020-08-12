# Do Monte Carlo experiments for 4 sites that are rougly bivariate LN (high PPCC, low M.Skew, M.Jurt around 8)
# I have already calculated the parameters for each site based on B bootstrap replicates 
# Parameters are averaged over all B bootstrap replicatse
# Parameters are LN3 estimators (except mean of O, mu_O. use product moment estimator)

# Barber, Lamontagne, Vogel
# A Monte Carlo experiment evaluating estimators of Efficiency
# Nash-Sutcliffe Efficiency (NSE), Kling-Gupta Efficiency (KGE), Pool et al. (2018) Estimator (POE)
# Introducing new estimator of efficiency: Barber-Lamontagne Efficiency (LBE)
# 2.11.2019

rm(list = ls())  # Remove everything in global enviroclcnment
graphics.off()   # Close all open plots
setwd("C:\\Users\\jlamon02\\Dropbox\\Caitline_Code\\ForJon\\")
source(file="C:\\Users\\jlamon02\\Dropbox\\Caitline_Code\\ForJon\\functions3BT.R", local=FALSE, echo=FALSE, print.eval=FALSE)
S <- read.csv(file='c:/Users/jlamon02/Dropbox/Caitline_Code/ForJon/S1.csv',header=TRUE,sep=",")
O <- read.csv(file='c:/Users/jlamon02/Dropbox/Caitline_Code/ForJon/O1.csv',header=TRUE,sep=",")
period <- read.csv(file='c:/Users/jlamon02/Dropbox/Caitline_Code/ForJon/period1.csv',header=TRUE,sep=",")

mu_O_month <- matrix(nrow=1, ncol=12)
Co_month <- matrix(nrow=1, ncol=12)
delta_month <- matrix(nrow=1, ncol=12)
theta_month <- matrix(nrow=1, ncol=12)
rho_month <- matrix(nrow=1, ncol=12)
tau_O_month<- matrix(nrow=1, ncol=12)
tau_S_month<- matrix(nrow=1, ncol=12)
mu_u_month<- matrix(nrow=1, ncol=12)
mu_v_month<- matrix(nrow=1, ncol=12)
sd_u_month<- matrix(nrow=1, ncol=12)
sd_v_month<- matrix(nrow=1, ncol=12)
rho_log_month<- matrix(nrow=1, ncol=12)
# mixture moment estimators
mu_mix_O_month<- matrix(nrow=1, ncol=12)
var_mix_O_month<- matrix(nrow=1, ncol=12)
mu_mix_S_month<- matrix(nrow=1, ncol=12)
var_mix_S_month<- matrix(nrow=1, ncol=12)
mu_mix_SO_month<- matrix(nrow=1, ncol=12)
mu_u_mix<- matrix(nrow=1, ncol=12)
mu_v_mix<- matrix(nrow=1, ncol=12)
tau_O_mix<- matrix(nrow=1, ncol=12)
tau_S_mix<- matrix(nrow=1, ncol=12)
sd_u_mix<- matrix(nrow=1, ncol=12)
sd_v_mix<- matrix(nrow=1, ncol=12)
mu_mix_O_month_MC<- matrix(nrow=1, ncol=12)
var_mix_O_month_MC_PART<- matrix(nrow=1, ncol=12)
mu_mix_S_month_MC<- matrix(nrow=1, ncol=12)
var_mix_S_month_MC_PART<- matrix(nrow=1, ncol=12)
mu_mix_SO_month_MC<- matrix(nrow=1, ncol=12)
rho_mix<- matrix(nrow=1, ncol=12)

mu_mix_O_MC <- matrix(nrow=1, ncol=1)
var_mix_O_MC <- matrix(nrow=1, ncol=1)
mu_mix_S_MC <- matrix(nrow=1, ncol=1)
var_mix_S_MC <- matrix(nrow=1, ncol=1)
mu_mix_SO_MC <- matrix(nrow=1, ncol=1)
theta_mix_MC <- matrix(nrow=1, ncol=1)
delta_mix_MC <- matrix(nrow=1, ncol=1)
Co_mix_MC <- matrix(nrow=1, ncol=1)
r_mix_MC <- matrix(nrow=1, ncol=1)
LBE_mix_MC <- matrix(nrow=1, ncol=1)
LBEprime_mix_MC <- matrix(nrow=1, ncol=1)
mu_mix_O <- matrix(nrow=1, ncol=1)
var_mix_O <- matrix(nrow=1, ncol=1)
mu_mix_S <- matrix(nrow=1, ncol=1)
var_mix_S <- matrix(nrow=1, ncol=1)

mu_mix_SO <- matrix(nrow=1, ncol=1)

theta_mix <- matrix(nrow=1, ncol=1)
delta_mix <- matrix(nrow=1, ncol=1)
Co_mix <- matrix(nrow=1, ncol=1)
Cs_mix <- matrix(nrow=1, ncol=1)
r_mix <- matrix(nrow=1, ncol=1)
LBE_mix <- matrix(nrow=1, ncol=1)
LBEprime_mix <- matrix(nrow=1, ncol=1)
rho <- matrix(nrow=1, ncol=1)


allData <- cbind(O,S,period)
colnames(allData, do.NULL = FALSE)
colnames(allData) <- c("O","S","month")
for (m in 1:12){ # For calculations based on monthly data
  oneSite_month <- subset(allData, allData$month== m) # Pull data for individual site's month
  LN3params_month <- parms(oneSite_month[,1:2])
  mu_O_month[m] <- LN3params_month[1] # real space mean
  rho_month[m] <- LN3params_month[2]
  Co_month[m] <- LN3params_month[3]
  delta_month[m] <- LN3params_month[4]
  theta_month[m] <- LN3params_month[5]
  tau_O_month[m] <- LN3params_month[6]
  tau_S_month[m] <- LN3params_month[7]
  mu_u_month[m] <- LN3params_month[8]
  mu_v_month[m] <- LN3params_month[9]
  sd_u_month[m] <- LN3params_month[10]
  sd_v_month[m] <- LN3params_month[11]
  rho_log_month[m] <- LN3params_month[12]
      
  # Per Vogel's suggestion (11/1/2019) if tau values are negative set tau to zero. This means fitting LN2 at these sites.
  if (tau_O_month[m]<0 | tau_S_month[m]<0){
    tau_O_month[m] <- 0
    tau_S_month[m] <- 0
    LN2count_month <- LN2count_month + 1
    }
      
  # mixture moment estimators
  mu_mix_O_month[m] <- tau_O_month[m]+exp(mu_u_month[m]+sd_u_month[m]^2/2)
  var_mix_O_month[m] <- (exp(2*mu_u_month[m]+sd_u_month[m]^2)*(exp(sd_u_month[m]^2)-1))
  mu_mix_S_month[m] <- tau_S_month[m]+exp(mu_v_month[m]+sd_v_month[m]^2/2)
  var_mix_S_month[m] <- (exp(2*mu_v_month[m]+sd_v_month[m]^2)*(exp(sd_v_month[m]^2)-1))
  mu_mix_SO_month[m] <- (mu_mix_S_month[m]*mu_mix_O_month[m]+rho_month[m]*sqrt(var_mix_S_month[m])*sqrt(var_mix_O_month[m]))
}

# mixture moments from RAW data 
mu_mix_O <- 1/12*sum(mu_mix_O_month)
var_mix_O <- 1/12*sum(var_mix_O_month+mu_mix_O_month^2)-mu_mix_O^2
mu_mix_S <- 1/12*sum(mu_mix_S_month)
var_mix_S <- 1/12*sum(var_mix_S_month+mu_mix_S_month^2)-mu_mix_S^2
    
mu_mix_SO <- 1/12*sum(mu_mix_SO_month)
    
theta_mix <- sqrt(var_mix_S)/sqrt(var_mix_O)
delta_mix <- 1-mu_mix_S/mu_mix_O
Co_mix <- sqrt(var_mix_O)/mu_mix_O
Cs_mix <- sqrt(var_mix_S)/mu_mix_S
#r1_mix[i] <- (1/nrow(oneSite)*sum(oneSite[,3]*oneSite[,4])-mu_mix_O[i]*mu_mix_S[i])/(sqrt(var_mix_O[i]*var_mix_S[i]))
r_mix <- (mu_mix_SO-mu_mix_O*mu_mix_S)/(sqrt(var_mix_O*var_mix_S))
    
LBE_mix <- 2*theta_mix*r_mix-theta_mix^2-delta_mix^2/Co_mix^2
LBEprime_mix <- 1-sqrt(delta_mix^2+(theta_mix-1)^2+(r_mix-1)^2)