
############## Stedinger's Lower Bound Estimator used for BOOTSTRAP DATA ##############
# don't think this is needed anymore
# DON'T ACTUALLY NEED THIS FUNCTION FOR BOOTSTRAP!

# Function to calculate Stedinger's lower bound estimator

# INPUT 
# data should by matrix with 3 columns
# column 2 should be observations, column 3 should be simulations
# Should be RAW data (not logged)

# OUTPUT 
# A vector with 2 values
# Value 1: lower bound estimator for observations
# Value 2: lower bound estimator for simulations

lwrBoundBT <- function(data){
  
  if ((min(data)+max(data)-2*median(data))>0){
    tau <- (min(data)*max(data)-median(data)^2)/(min(data)+max(data)-2*median(data))     # Compute tau_O for each MC trial
  } else {
    tau <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }
  return(tau)
}


############## Stedinger's Lower Bound Estimator used for REAL DATA ##############
# not needed anymore
# DON'T ACTUALLY NEED THIS FUNCTION FOR BOOTSTRAP!

# Function to calculate Stedinger's lower bound estimator

# INPUT 
# data should by matrix with 3 columns
# column 2 should be observations, column 3 should be simulations
# Should be RAW data (not logged)

# OUTPUT 
# A vector with 2 values
# Value 1: lower bound estimator for observations
# Value 2: lower bound estimator for simulations

lwrBound <- function(data){

  if ((min(data[,2])+max(data[,2])-2*median(data[,2]))>0){
    tau_O <- (min(data[,2])*max(data[,2])-median(data[,2])^2)/(min(data[,2])+max(data[,2])-2*median(data[,2]))     # Compute tau_O for each MC trial
  } else {
    tau_O <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }
  
  
  if ((min(data[,3])+max(data[,3])-2*median(data[,3]))>0){
    tau_S <- (min(data[,3])*max(data[,3])-median(data[,3])^2)/(min(data[,3])+max(data[,3])-2*median(data[,3]))     # Compute tau_O for each MC trial
  } else {
    tau_S <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }

  tau <- cbind(tau_O, tau_S)
  return(tau)
}


############### Mardia Skew and Kurtosis ###################
mardia <- function(data, cov = TRUE, tol = 1e-25){
  n <- nrow(data)
  data <- scale(data, scale = FALSE)
  S <- ((nrow(data) - 1)/nrow(data)) * cov(data)
  D <- data %*% solve(S, tol = tol) %*% t(data)
  skew <-  nrow(data) * (sum(D^3)/nrow(data)^2)/6
  kurt <- ((sum(diag((D^2)))/nrow(data)) - ncol(data) * (ncol(data) + 2)) * sqrt(nrow(data)/(8 * ncol(data) * (ncol(data) + 2)))
  result = rbind(skew,kurt)
  result
}



################ Efficiency Estimators ##################
# INPUT
# data should be a matrix with 2 columns, should be RAW data (not logged)
# column 1:observations, column 2: simulations

# OUTPUT
# returns matrix with efficiency estimators for a single site
# LBE, NSE, LNSE, LBE', KGE

eff <- function(data){
  S <- data[,2]   # Specify similuations
  O <- data[,1]   # Specify observations
  # Get sample size of station data (how many observed and PRMS simulated streamflow data)
  n <- nrow(data)
  
  # Calculate correlation coefficient (calculate both r and r*=r1)
  r <- cor(S,O)   # Calculate correlation coefficient (rho) between simulations and observations
  
  mu_O <- mean(O)  # Calculate mean of observations,mu_O (mean of column 1)
  mu_S <- mean(S)  # Calculate mean of simulations, mu_S (mean of column 2)
  
  # Product moment estimators
  delta <- (mu_O-mu_S)/mu_O   # Calculate delta from mean obs and mean sims
  theta <- sd(S)/sd(O)   # Calculate theta = stdev(S)/stdev(O)
  Co <- sd(O)/mu_O  # Caluclate Co as stdev/mean (different than how it's estimated for LBE)
  
  
  
  ############## NSE ##############
  # Nash-Sutcliffe Efficiency (NSE)
  NSE <- 1-((1/n)*sum((S-O)^2))/((1/(n-1))*sum((O-mean(O))^2))   # Barber et al. 2019 eqn. 8b
  
  # Compute NSE using natural logs. Call this LogNSE
  LogNSE <- 1-((1/n)*sum((log(S)-log(O))^2))/((1/(n-1))*sum((log(O)-mean(log(O)))^2))
  
  
  
  ############## KGE ##############
  alpha <- sd(S)/sd(O)        # Estimator of alpha = estimator of theta = std. dev. simulations / std. dev. observations
  beta <- mean(S)/mean(O)
  KGE <- 1-sqrt((beta-1)^2+(alpha-1)^2+(r-1)^2)   # Barber et al. 2019 eqn 14
  
  
  
  ############## RNP (PVSE) ##############
  # Pool Efficiency (RNP) (Taken from Pool et al. 2018 Supplemental Material R script)
  fdc.sim = sort(S / (mean(S) * n))     # Calculate normalized flow duration curves
  fdc.obs = sort(O / (mean(O) * n))
  alpha_RNP = 1 - 0.5 * sum(abs(fdc.sim - fdc.obs))     # Calculate alpha component
  # Beta = Same as KGE beta
  r_RNP = cor(S, O, method="spearman")     # Calculate r component
  RNP <- 1 - sqrt((alpha_RNP - 1)^2 + (beta - 1)^2 + (r_RNP - 1)^2)     # Return Non-Parametric Efficiency value
  
  
  
  ############## LBE ##############
  # Eq. 16
  # Calculate tau for observations and simulations
  if ((min(O)+max(O)-2*median(O))>0){
    tau_O <- (min(O)*max(O)-median(O)^2)/(min(O)+max(O)-2*median(O))     # Compute tau_O for each MC trial
  } else {
    tau_O <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }
  
  if ((min(S)+max(S)-2*median(S))>0){
    tau_S <- (min(S)*max(S)-median(S)^2)/(min(S)+max(S)-2*median(S))     # Compute tau_S for each MC trial
  } else {
    tau_S <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }
  
  #### ADDED 11.23.2019 to be consistent with Vogel's suggestion (11/1/2019) if tau values are negative set tau to zero. This means fitting LN2 at these sites.
  # If can't reproduce past results comment this out
  if (tau_O<0 | tau_S<0){
    tau_O <- 0
    tau_S <- 0
  }
  
  u <- log(O-tau_O)      # Compute u based on observations of LN2 to use for LBE equation
  # changed to sample variance on 11.23.2019
  var_u <- (1/(n-1))*sum((u-mean(u))^2)   # sample variance of 
  v <- log(S-tau_S)      # Compute v based on observations of LN2 to use for LBE equation
  var_v <- (1/(n-1))*sum((v-mean(v))^2)   # Population variance of v
  
  # Estimator of correlation coefficient, r* (r1),  using Stedinger (1981) eqn 2
  # Stedinger (1981) eqn 3: estimator of variance of log of observations and simulations (CHECK!!!)
  s2_yOyS <- 1/n*sum((u-mean(u))*(v-mean(v))) 
  s2_yO <- 1/n*sum((u-mean(u))^2)   # Stedinger (1981) eqn 3: estimator of variance of log of observations
  s2_yS <- 1/n*sum((v-mean(v))^2)    # Stedinger (1981) eqn 3: estimator of variance of log of simulations
  r1 <- (exp(s2_yOyS)-1)/sqrt((exp(s2_yO)-1)*(exp(s2_yS)-1)) # Stedinger (1981) eqn 2: estimator of correlation coefficient, r* (r1)
  
  # LN3 estimators
  Co_LBE <- (sqrt(exp(2*mean(u)+var_u)*(exp(var_u)-1)))/(tau_O+exp(mean(u)+var_u/2)) # LBE estimate of coefficient of variation of observations
  theta_LBE <- sqrt((exp(2*mean(v)+var_v)*(exp(var_v)-1))/(exp(2*mean(u)+var_u)*(exp(var_u)-1))) # LBE estimate of theta
  delta_LBE <- 1-(tau_S+exp(mean(v)+var_v/2))/(tau_O+exp(mean(u)+var_u/2)) # LBE estimate of delta
  
  LBE <- 2*theta_LBE*r1-theta_LBE^2-delta_LBE^2/Co_LBE^2 # Barber et al. 2019 eq. 16
  LBEprime <- 1-sqrt((-1*delta_LBE)^2+(theta_LBE-1)^2+(r1-1)^2)
  
  effResults <- cbind(LBE, NSE, LogNSE, LBEprime, KGE, RNP)

  return(effResults)
}



################ Efficiency Estimators - MODIFIED input ##################
# May not need anymore
# INPUT
# data should be in two different vectors, should be RAW data (not logged)
# x:observations, y: simulations

# OUTPUT
# returns matrix with efficiency estimators for a single site
# LBE, NSE, LNSE, LBE', KGE

eff2 <- function(x,y){
  S <- y   # Specify similuations
  O <- x   # Specify observations
  # Get sample size of station data (how many observed and PRMS simulated streamflow data)
  n <- nrow(x)
  
  # Calculate correlation coefficient (calculate both r and r*=r1)
  r <- cor(S,O)   # Calculate correlation coefficient (rho) between simulations and observations
  
  mu_O <- mean(O)  # Calculate mean of observations,mu_O (mean of column 1)
  mu_S <- mean(S)  # Calculate mean of simulations, mu_S (mean of column 2)
  
  # Product moment estimators
  delta <- (mu_O-mu_S)/mu_O   # Calculate delta from mean obs and mean sims
  theta <- sd(S)/sd(O)   # Calculate theta = stdev(S)/stdev(O)
  Co <- sd(O)/mu_O  # Caluclate Co as stdev/mean (different than how it's estimated for LBE)
  
  
  
  ############## NSE ##############
  # Nash-Sutcliffe Efficiency (NSE)
  NSE <- 1-((1/n)*sum((S-O)^2))/((1/(n-1))*sum((O-mean(O))^2))   # Barber et al. 2019 eqn. 8b
  
  # Compute NSE using natural logs. Call this LogNSE
  LogNSE <- 1-((1/n)*sum((log(S)-log(O))^2))/((1/(n-1))*sum((log(O)-mean(log(O)))^2))
  
  
  
  ############## KGE ##############
  alpha <- sd(S)/sd(O)        # Estimator of alpha = estimator of theta = std. dev. simulations / std. dev. observations
  beta <- mean(S)/mean(O)
  KGE <- 1-sqrt((beta-1)^2+(alpha-1)^2+(r-1)^2)   # Barber et al. 2019 eqn 14
  
  
  
  ############## RNP ##############
  # Pool Efficiency (RNP) (Taken from Pool et al. 2018 Supplemental Material R script)
  fdc.sim = sort(S / (mean(S) * n))     # Calculate normalized flow duration curves
  fdc.obs = sort(O / (mean(O) * n))
  alpha_RNP = 1 - 0.5 * sum(abs(fdc.sim - fdc.obs))     # Calculate alpha component
  # Beta = Same as KGE beta
  r_RNP = cor(S, O, method="spearman")     # Calculate r component
  RNP <- 1 - sqrt((alpha_RNP - 1)^2 + (beta - 1)^2 + (r_RNP - 1)^2)     # Return Non-Parametric Efficiency value
  
  
  
  ############## LBE ##############
  # Eq. 16
  # Calculate tau for observations and simulations
  if ((min(O)+max(O)-2*median(O))>0){
    tau_O <- (min(O)*max(O)-median(O)^2)/(min(O)+max(O)-2*median(O))     # Compute tau_O for each MC trial
  } else {
    tau_O <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }
  
  if ((min(S)+max(S)-2*median(S))>0){
    tau_S <- (min(S)*max(S)-median(S)^2)/(min(S)+max(S)-2*median(S))     # Compute tau_S for each MC trial
  } else {
    tau_S <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }
  
  u <- log(O-tau_O)      # Compute u based on observations of LN2 to use for LBE equation
  var_u <- (1/n)*sum((u-mean(u))^2)   # Population variance of 
  v <- log(S-tau_S)      # Compute v based on observations of LN2 to use for LBE equation
  var_v <- (1/n)*sum((v-mean(v))^2)   # Population variance of v
  
  # Estimator of correlation coefficient, r* (r1),  using Stedinger (1981) eqn 2
  # Stedinger (1981) eqn 3: estimator of variance of log of observations and simulations (CHECK!!!)
  s2_yOyS <- 1/n*sum((u-mean(u))*(v-mean(v))) 
  s2_yO <- 1/n*sum((u-mean(u))^2)   # Stedinger (1981) eqn 3: estimator of variance of log of observations
  s2_yS <- 1/n*sum((v-mean(v))^2)    # Stedinger (1981) eqn 3: estimator of variance of log of simulations
  r1 <- (exp(s2_yOyS)-1)/sqrt((exp(s2_yO)-1)*(exp(s2_yS)-1)) # Stedinger (1981) eqn 2: estimator of correlation coefficient, r* (r1)
  
  # LN3 estimators
  Co_LBE <- (sqrt(exp(2*mean(u)+var_u)*(exp(var_u)-1)))/(tau_O+exp(mean(u)+var_u/2)) # LBE estimate of coefficient of variation of observations
  theta_LBE <- sqrt((exp(2*mean(v)+var_v)*(exp(var_v)-1))/(exp(2*mean(u)+var_u)*(exp(var_u)-1))) # LBE estimate of theta
  delta_LBE <- 1-(tau_S+exp(mean(v)+var_v/2))/(tau_O+exp(mean(u)+var_u/2)) # LBE estimate of delta
  
  LBE <- 2*theta_LBE*r1-theta_LBE^2-delta_LBE^2/Co_LBE^2 # Barber et al. 2019 eq. 16
  LBEprime <- 1-sqrt((-1*delta_LBE)^2+(theta_LBE-1)^2+(r1-1)^2)
  
  effResults <- cbind(LBE, NSE, LogNSE, LBEprime, KGE, RNP)
  
  return(effResults)
}


################ Efficiency Parameters (LN3) ##################
# INPUT
# data should be a matrix with 2 columns, should be RAW data (not logged)
# column 1:observations, column 2: simulations

# OUTPUT
# returns matrix with efficiency parameters for a single site (LN3 estimators of parameters)
# mu_O, r1, Co_BLE, delta_BLE, theta_BLE

parms <- function(data){
  S <- data[,2]   # Specify similuations
  O <- data[,1]   # Specify observations
  # Get sample size of station data (how many observed and PRMS simulated streamflow data)
  n <- nrow(data)

  mu_O <- mean(O)  # Calculate mean of observations,mu_O (mean of column 1)

  ############## LBE ##############
  # Eq. 16
  # Calculate tau for observations and simulations
  if ((min(O)+max(O)-2*median(O))>0){
    tau_O <- (min(O)*max(O)-median(O)^2)/(min(O)+max(O)-2*median(O))     # Compute tau_O for each MC trial
  } else {
    tau_O <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }
  
  if ((min(S)+max(S)-2*median(S))>0){
    tau_S <- (min(S)*max(S)-median(S)^2)/(min(S)+max(S)-2*median(S))     # Compute tau_S for each MC trial
  } else {
    tau_S <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }
  
  #### ADDED 11.23.2019 to be consistent with Vogel's suggestion (11/1/2019) if tau values are negative set tau to zero. This means fitting LN2 at these sites.
  # If can't reproduce past results comment this out
  if (tau_O<0 | tau_S<0){
    tau_O <- 0
    tau_S <- 0
  }
  
  u <- log(O-tau_O)      # Compute u based on observations of LN2 to use for LBE equation
  mu_u <- mean(u)
  var_u <- (1/(n-1))*sum((u-mean(u))^2)   # Population variance of 
  sd_u <- sqrt(var_u)
  v <- log(S-tau_S)      # Compute v based on observations of LN2 to use for LBE equation
  mu_v <- mean(v)
  var_v <- (1/(n-1))*sum((v-mean(v))^2)   # Population variance of v
  sd_v <- sqrt(var_v)
  
  rho_log <- cor(u,v)
  
  # Estimator of correlation coefficient, r* (r1),  using Stedinger (1981) eqn 2
  # Stedinger (1981) eqn 3: estimator of variance of log of observations and simulations (CHECK!!!)
  s2_yOyS <- 1/n*sum((u-mean(u))*(v-mean(v))) 
  s2_yO <- 1/n*sum((u-mean(u))^2)   # Stedinger (1981) eqn 3: estimator of variance of log of observations
  s2_yS <- 1/n*sum((v-mean(v))^2)    # Stedinger (1981) eqn 3: estimator of variance of log of simulations
  r1 <- (exp(s2_yOyS)-1)/sqrt((exp(s2_yO)-1)*(exp(s2_yS)-1)) # Stedinger (1981) eqn 2: estimator of correlation coefficient, r* (r1)
  
  # LN3 estimators
  Co_LBE <- (sqrt(exp(2*mean(u)+var_u)*(exp(var_u)-1)))/(tau_O+exp(mean(u)+var_u/2)) # LBE estimate of coefficient of variation of observations
  theta_LBE <- sqrt((exp(2*mean(v)+var_v)*(exp(var_v)-1))/(exp(2*mean(u)+var_u)*(exp(var_u)-1))) # LBE estimate of theta
  delta_LBE <- 1-(tau_S+exp(mean(v)+var_v/2))/(tau_O+exp(mean(u)+var_u/2)) # LBE estimate of delta
  
  LBE <- 2*theta_LBE*r1-theta_LBE^2-delta_LBE^2/Co_LBE^2 # Barber et al. 2019 eq. 16
  LBEprime <- 1-sqrt((-1*delta_LBE)^2+(theta_LBE-1)^2+(r1-1)^2)
  
  LN3parms <- cbind(mu_O, r1, Co_LBE, delta_LBE, theta_LBE, tau_O, tau_S, mu_u, mu_v, sd_u, sd_v, rho_log)
  
  return(LN3parms)
}




################ Efficiency Parameters (LN3) EXPANDED ##################
# Not needed anymore
# INPUT
# data should be a matrix with 2 columns, should be RAW data (not logged)
# column 1:observations, column 2: simulations

# OUTPUT
# returns matrix with efficiency parameters for a single site (LN3 estimators of parameters)
# mu_O, r1, Co_BLE, delta_BLE, theta_BLE, tau_O, tau_S, mu_u, mu_v, sd_u, sd_v, rho_log

parms2 <- function(data){
  S <- data[,2]   # Specify similuations
  O <- data[,1]   # Specify observations
  # Get sample size of station data (how many observed and PRMS simulated streamflow data)
  n <- nrow(data)
  
  mu_O <- mean(O)  # Calculate mean of observations,mu_O (mean of column 1)
  mu_S <- mean(S)
  sd_O <- sd(O)
  sd_S <- sd(S)
  
  ############## LBE ##############
  # Eq. 16
  # Calculate tau for observations and simulations
  if ((min(O)+max(O)-2*median(O))>0){
    tau_O <- (min(O)*max(O)-median(O)^2)/(min(O)+max(O)-2*median(O))     # Compute tau_O for each MC trial
  } else {
    tau_O <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }
  
  if ((min(S)+max(S)-2*median(S))>0){
    tau_S <- (min(S)*max(S)-median(S)^2)/(min(S)+max(S)-2*median(S))     # Compute tau_S for each MC trial
  } else {
    tau_S <- 0   # Fit LN2 model instead of LN3 because Stedinger (1980) lower bound estimator cannot be computed
  }
  
  u <- log(O-tau_O)      # Compute u based on observations of LN2 to use for LBE equation
  mu_u <- log((mu_O-tau_O)/sqrt(1+(sd_O/(mu_O-tau_O))^2))
  var_u <- (1/n)*sum((u-mu_u)^2)   # Population variance of 
  sd_u <- sqrt(log(1+(sd_O/(mu_O-tau_O))^2))
  v <- log(S-tau_S)      # Compute v based on observations of LN2 to use for LBE equation
  mu_v <- log((mu_S-tau_S)/sqrt(1+(sd_S/(mu_S-tau_S))^2))
  var_v <- (1/n)*sum((v-mu_v)^2)   # Population variance of v
  sd_v <-  sqrt(log(1+(sd_S/(mu_S-tau_S))^2))
  
  rho_log <- cor(u,v)
  
  # Estimator of correlation coefficient, r* (r1),  using Stedinger (1981) eqn 2
  # Stedinger (1981) eqn 3: estimator of variance of log of observations and simulations (CHECK!!!)
  s2_yOyS <- 1/n*sum((u-mean(u))*(v-mean(v))) 
  s2_yO <- 1/n*sum((u-mean(u))^2)   # Stedinger (1981) eqn 3: estimator of variance of log of observations
  s2_yS <- 1/n*sum((v-mean(v))^2)    # Stedinger (1981) eqn 3: estimator of variance of log of simulations
  r1 <- (exp(s2_yOyS)-1)/sqrt((exp(s2_yO)-1)*(exp(s2_yS)-1)) # Stedinger (1981) eqn 2: estimator of correlation coefficient, r* (r1)
  
  # LN3 estimators
  Co_LBE <- (sqrt(exp(2*mean(u)+var_u)*(exp(var_u)-1)))/(tau_O+exp(mean(u)+var_u/2)) # LBE estimate of coefficient of variation of observations
  theta_LBE <- sqrt((exp(2*mean(v)+var_v)*(exp(var_v)-1))/(exp(2*mean(u)+var_u)*(exp(var_u)-1))) # LBE estimate of theta
  delta_LBE <- 1-(tau_S+exp(mean(v)+var_v/2))/(tau_O+exp(mean(u)+var_u/2)) # LBE estimate of delta
  
  LBE <- 2*theta_LBE*r1-theta_LBE^2-delta_LBE^2/Co_LBE^2 # Barber et al. 2019 eq. 16
  LBEprime <- 1-sqrt((-1*delta_LBE)^2+(theta_LBE-1)^2+(r1-1)^2)
  
  LN3parms <- cbind(mu_O, r1, Co_LBE, delta_LBE, theta_LBE, tau_O, tau_S, mu_u, mu_v, sd_u, sd_v, rho_log)
  
  return(LN3parms)
}


############### Custom Boxplots (1st and 99th percentiles) ########################
# not needed anymore
# After using this function, we want to plot the data using boxplot(newdata=boxplotQuants(data), range=0)
# NEED to specify range=0, or else it will plot quantiles as outliers
boxplotQuants <- function(data){
  # Calculate quantiles of the data
  quant <- quantile(data,probs=c(0.01,0.99))
  # Calculate the data used to make a boxplot, but don't plot it
  pnts <- boxplot(data, plot=FALSE)
  # Change the upper and lower whiskers to the quantile values
  pnts$stats[1] <- quant[1]
  pnts$stats[5] <- quant[2]
  # Calculate the data used to make the updated version of the boxplot
  boxquant <- boxplot(pnts$stats, range=0, plot=FALSE)
  whskdata <- as.vector(boxquant$stats)
  return(whskdata)
}

