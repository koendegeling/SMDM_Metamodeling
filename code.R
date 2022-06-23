# This script is part of the "Reducing the Computational Burden of Health Economic Models and Analyses 
# Using Metamodeling" workshop presented by Hendrik Koffijberg (University of Twente, the Netherlands) 
# and Koen Degeling (Lumen Value & Access, US/Netherlands) for Health Technology Assessment 
# International (HTAi). It illustrates the use of metamodels to address run time issues (i.e., the 
# computational burden) of a health economic simulation models to enable optimization.
#
# The hypothetical case study used in this script is based on the paper "Using Metamodeling to Identify  
# the Optimal Strategy for Colorectal Cancer Screening" in Value in Health (2021). The health economic
# model (HEM) generates the net monetary benefit (NMB) of a screening strategy compared to no screening
# at a willingness to pay of EUR 20,000 per quality-adjusted life year (QALY) gained. The screening 
# strategies are defined by four parameters:
#   - StartAge      integer defining the age (years) at which screening is started
#   - Interval      integer defining the interval (years) at which screening occurs
#   - NrScr         integer defining the number of screening rounds
#   - FITcutoff     integer defining the test threshold used
#
# The HEM and the function runHEM() to run the HEM are saved in the separate objects.RDA file.
#
# The script contains the following sections:
#   1. INITIALIZATION           preparing the session to run the rest of the script
#   2. RUN TIME EXPERIMENTS     experimenting with the run time of the HEM
#   3. SIMULATING DATASETS      simulating datasets to train and test the metamodels
#   4. FITTING METAMODELS       fitting the metamodels to the training dataset
#   5. ASSESSING PERFORMANCE    assessing the performance of the metamodels on the testing dataset
#   6. USING THE METAMODEL      using the metamodel to optimize the NMB
#
# The script was written en tested for R version 4.0.3 and requires the following packages:
#   - lhs     version 1.1.3
#   - gam     version 1.20
#   - GPfit   version 1.0-8
#
# Note that all code needs to be run from top to bottom, otherwise errors may occur.


### 1. INITIALIZATION ----

# Clear the Global Environment
rm(list = ls()); gc();

# Install packages if they are not already installed
if(!require(lhs))   install.packages('lhs')
if(!require(gam))   install.packages('gam')
if(!require(GPfit)) install.packages('GPfit')

# Load the required packages
library(lhs)        # latin hypercube sampling
library(gam)        # generalized additive models
library(GPfit)      # Gaussian process
library(parallel)   # multi-threading

# Load the HEM object and the function runHEM() to run the HEM from the objects.RDA file
# - ensure that the objects.RDA is located in the current working directory or update the file path (you
#   can use the getwd() function to see what the current working directory is)
load(file = 'objects.RDA')

# The runHEM() function simulates the NMB for a cancer screening strategy compared to no screening, at a
# willingness to pay of EUR 20,000 per QALY gained, using a hypothetical health economic model. The 
# screening strategy is defined using a set of parameter values defined through the 'pars' argument. This
# 'pars' argument should be a named vector of integers (defined in R, for example, as 1L), with following
# parameters/names:
# - StartAge      the start age of the screening strategy in years          range: 30L - 90L
# - Interval      the screening interval in years                           range: 1L - 60L
# - NrScr         the number of screening rounds                            range: 1L - 31L
# - FITcutoff     the threshold used for classifying the test results       50L, 75L, 100L, 150L
#
# For example, to simulate the outcomes for the current colorectal cancer screening program in the 
# Netherlands, where screening starts at age 55 and consists of 11 biennial screening rounds, using a 75
# ng/mL FIT cutoff, the corresponding 'pars' argument should be:
#   pars = c(StartAge = 55L, Interval = 2L, NrScr = 11L, FITcutoff = 75L)
#
# The function has a 'min_runtime' argument that can be used to specify the minimum run time of the model
# in seconds. This argument can be used to experiment with the impact of increasing run time of the model
# on the required time to perform, for example, a probabilistic analysis.
#
# The run time will be printed to the console if the 'print_runtime' argument is TRUE.

# Parameters according to the Dutch colorectal cancer screening program
pars_DutchScreening <- c(StartAge = 55L, Interval = 2L, NrScr = 11L, FITcutoff = 75L)

# Test the HEM
runHEM(pars = pars_DutchScreening, print_runtime = TRUE)

# Test the HEM with increased runtime
runHEM(pars = pars_DutchScreening, min_runtime = 3, print_runtime = TRUE)




### 2. RUN TIME EXPERIMENTS ----

# The code in this section can be used to experiment with the impact the run time and number of runs on
# the total time required to run an analysis. The number of runs can be defined using the n_runs 
# parameter and the run time (in seconds) can be defined through the min_runtime parameter.
n_runs <- 10
min_runtime <- 1

# The code below performs the runs one-by-one
system.time(expr = {
  runs_out <- sapply(1:n_runs, function(i) runHEM(pars = pars_DutchScreening, min_runtime = min_runtime))
})

# The code below performs the runs in parallel, using all the CPU cores available. The detectCores()
# function can be used to find out how many CPU cores your laptop/desktop has. The following steps are
# performed to run the code in parallel:
# - define the computing cluster 'cl' using the makeCluster() function
# - export the objects that need to be available on each node using the clusterExport() function
# - run the code in parallel using the parallel version of the sapply() function: parSapply()
# - stop the computing cluster using the stopCluster() function
detectCores()

system.time(expr = {
  cl <- makeCluster(detectCores())
  clusterExport(cl = cl, varlist = c('HEM', 'runHEM', 'pars_DutchScreening', 'min_runtime'))
  runs_out <- parSapply(cl, 1:n_runs, function(i) runHEM(pars = pars_DutchScreening, min_runtime = min_runtime))
  stopCluster(cl)
})
  



### 3. SIMULATING DATASETS ----

# This section defines the sets of experiments or datasets that will be used to train the metamodels and
# test their performance. Two designs of experiments will be considered: a random design and a Latin
# Hypercube design.

# Since some metamodeling methods require the input parameters to be defined on a 0-1 scale, functions
# are defined to normalize the parameters from their original scale to the normalized 0-1 scale.
fun_normalize <- function(org_pars) {
  
  # The minimum and maximum values of the parameters
  min_values <- c(30, 1, 1, 1)
  max_values <- c(90, 60, 31, 4)
  
  # Ensure the parameters are in the right order
  org_pars  <- org_pars[c('StartAge', 'Interval', 'NrScr', 'FITcutoff')]
  
  # Transform the FITcutoff parameter to a categorical variable
  org_pars['FITcutoff'] <- switch(as.character(org_pars['FITcutoff']), '50' = 1, '75' = 2, '100' = 3, '150' = 4)
  
  # Normalize the parameter values
  norm_pars <- (org_pars - min_values)/(max_values - min_values)

  # Check whether all values are on the 0-1 scale
  if(any(norm_pars < 0) | any(norm_pars > 1)) stop('Parameters outside the allowed range provided')
  
  return(norm_pars)
  
}
fun_denormalize <- function(norm_pars) {
  
  # Check whether all values are on the 0-1 scale
  if(any(norm_pars < 0) | any(norm_pars > 1)) stop('Parameters outside the allowed range provided')
  
  # The minimum and maximum values of the parameters
  min_values <- c(30, 1, 1, 1)
  max_values <- c(90, 60, 31, 4)
  
  # Ensure the parameters are in the right order
  norm_pars <- norm_pars[c('StartAge', 'Interval', 'NrScr', 'FITcutoff')]
  
  # Return the parameters to their original scale
  org_pars  <- min_values + (norm_pars * (max_values - min_values))
  
  # Ensure the parameters are integers and set the names
  org_pars  <- setNames(object = as.integer(org_pars), nm = names(org_pars))
  
  # Return the FITcutoff parameter to its continuous scale
  org_pars['FITcutoff'] <- switch(org_pars['FITcutoff'], 50L, 75L, 100L, 150L)
  
  return(org_pars)
  
}

# Illustration of the normalization and denormalization
pars_DutchScreening

(pars_DutchScreening_norm <- fun_normalize(pars_DutchScreening))

fun_denormalize(pars_DutchScreening_norm)

# Defining the number of experiments that are to be included in the training dataset used to fit the 
# metamodels (n_training) and testing dataset used to assess the performance of the models (n_testing)
n_training <- 100
n_testing  <- 50


## 3.1 Random design ----

# The random design can be simply generated by sampling from a Uniform distribution, which will generate
# normalized values that can then be transformed to the original scales of the parameters
# - the random number seed is set for reproducibility
set.seed(123)
pars_training_rand_norm <- matrix(
  data     = runif(n = n_training * 4), 
  ncol     = 4, 
  dimnames = list(NULL, c('StartAge', 'Interval', 'NrScr', 'FITcutoff'))
)

set.seed(456)
pars_testing_rand_norm <- matrix(
  data     = runif(n = n_testing * 4), 
  ncol     = 4, 
  dimnames = list(NULL, c('StartAge', 'Interval', 'NrScr', 'FITcutoff'))
)

# Transforming the values on 0-1 scale to the original parameter scale
pars_training_rand_org <- t(apply(pars_training_rand_norm, 1, fun_denormalize))
pars_testing_rand_org  <- t(apply(pars_testing_rand_norm, 1, fun_denormalize))


## 3.2 Latin hypercube design ----

# The Latin Hypercube design is generated using the maximinLHS() function of the lhs package for which 
# only the number of experiments/samples (n) and number of parameters (k) needs to be specified. The
# maximin Latin Hypercube design is an optimized version of the original Latin Hypercube for which the
# minimum distance between the different experiments is maximized. Similar to the random design, values
# are generated on a 0-1 scale.
# - again, the random number seed is set for reproducibility
set.seed(123)
pars_training_lhs_norm <- maximinLHS(n = n_training, k = 4)

set.seed(456)
pars_testing_lhs_norm  <- maximinLHS(n = n_testing, k = 4)

# Providing names to the generated matrices
colnames(pars_training_lhs_norm) <- colnames(pars_testing_lhs_norm) <- c('StartAge', 'Interval', 'NrScr', 'FITcutoff')
  
# Transforming the values on 0-1 scale to the original parameter scale
pars_training_lhs_org <- t(apply(pars_training_lhs_norm, 1, fun_denormalize))
pars_testing_lhs_org  <- t(apply(pars_testing_lhs_norm, 1, fun_denormalize))


## 3.3 Comparison ----

# Inspect the efficiency/coverage of the generated samples

par(mfrow = c(1, 2)) # two plots side-by-side

plot(x = pars_training_rand_org[ , 'StartAge'], y = pars_training_rand_org[ , 'Interval'],
     pch = 16, las = 1, xlab = 'StartAge', ylab = 'Interval', main = 'Random Design')

plot(x = pars_training_lhs_org[ , 'StartAge'], y = pars_training_lhs_org[ , 'Interval'],
     pch = 16, las = 1, xlab = 'StartAge', ylab = 'Interval', main = 'Latin Hypercube Design')




### 4. FITTING METAMODELS ----

# Three types of metamodels are fitted in this section: a linear regression model (lm), a generalized
# additive model (gam), and a Gaussian process (gp)
# - refer to the Help tab for more information on the functions used to fit the metamodels

# Before fitting the metamodels, the true values for the experiments need to be obtained from the HEM
out_training <- t(apply(pars_training_lhs_org, 1, runHEM))


## 4.1 Linear regression model ----

mm_lm <- lm(formula = NMB ~ StartAge + Interval + NrScr + FITcutoff,
            data = as.data.frame(out_training))

mm_lm


## 4.2 Generalized additive model ----

mm_gam <- gam(formula = NMB ~ s(StartAge) + s(Interval) + s(NrScr) + FITcutoff +  
                s(StartAge * Interval) + s(StartAge * NrScr) + s(StartAge * FITcutoff) + 
                s(Interval * NrScr) + s(Interval * FITcutoff) + s(NrScr * FITcutoff),
              data = as.data.frame(out_training))

mm_gam


## 4.3 Gaussian process ----

mm_gp <- GP_fit(X = pars_training_lhs_norm, Y = out_training[, 'NMB'])

mm_gp




### 5. ASSESSING PERFORMANCE ----

# In this section, the performance of the different metamodels is assessed for a range of error measures
# and in a calibration or quantile-quantile (Q-Q) plot. The considered error measures are the absolute
# error, relative error, relative absolute error, and the root mean squared error.

# Before the performance can be assessed, the outcomes need to be obtained for the testing datasets
out_testing <- t(apply(pars_testing_lhs_org, 1, runHEM))[, 'NMB']
out_lm      <- predict.lm(object = mm_lm, newdata = as.data.frame(pars_testing_lhs_org))
out_gam     <- predict.Gam(object = mm_gam, newdata = as.data.frame(pars_testing_lhs_org))
out_gp      <- predict.GP(object = mm_gp, xnew = pars_testing_lhs_norm)$Y_hat

## 5.1 Error (E) ----

getE <- function(predicted, observed) {
  
  diff <- predicted - observed
  
  out <- c(
    mean = mean(diff),
    min  = min(diff),
    max  = max(diff)
  )
  
  return(out)
  
}

getE(out_lm,  out_testing)
getE(out_gam, out_testing)
getE(out_gp,  out_testing)


## 5.2 Relative error (RE) ----

getRE <- function(predicted, observed) {
  
  rel_diff <- (predicted - observed) / observed
  
  out <- c(
    mean = mean(rel_diff),
    min  = min(rel_diff),
    max  = max(rel_diff)
  )
  
  return(out)
  
}

getRE(out_lm,  out_testing)
getRE(out_gam, out_testing)
getRE(out_gp,  out_testing)


## 5.3 Absolute error (AE) ----

getAE <- function(predicted, observed) {
  
  abs_diff <- abs(predicted - observed)
  
  out <- c(
    mean = mean(abs_diff),
    min  = min(abs_diff),
    max  = max(abs_diff)
  )
  
  return(out)
    
}

getAE(out_lm,  out_testing)
getAE(out_gam, out_testing)
getAE(out_gp,  out_testing)


## 5.4 Relative absolute error (RAE) ----

getRAE <- function(predicted, observed) {
  
  relabs_diff <- abs(predicted - observed) / observed
  
  out <- c(
    mean = mean(relabs_diff),
    min  = min(relabs_diff),
    max  = max(relabs_diff)
  )
  
  return(out)
  
}

getRAE(out_lm,  out_testing)
getRAE(out_gam, out_testing)
getRAE(out_gp,  out_testing)


## 5.5 Root mean squared error (RMSE) ----

getRMSE <- function(predicted, observed) {
  
  out <- sqrt(mean((predicted - observed)^2))
  
  return(out)
  
}

getRMSE(out_lm,  out_testing)
getRMSE(out_gam, out_testing)
getRMSE(out_gp,  out_testing)


## 5.6 Calibration plots ----

axis_lim <- c(-2000, 3000)  # limits for the axes
par(mfrow = c(1, 3))        # plot three plots side-by-side

plot(x = out_lm, y = out_testing, xlim = axis_lim, ylim = axis_lim, 
     las = 1, xlab = 'Predicted', ylab = 'Observed', main = 'Linear regression')
lines(x = axis_lim, y = axis_lim, col = 'red')

plot(x = out_gam, y = out_testing, xlim = axis_lim, ylim = axis_lim, 
     las = 1, xlab = 'Predicted', ylab = 'Observed', main = 'Generalized additive model')
lines(x = axis_lim, y = axis_lim, col = 'red')

plot(x = out_gp, y = out_testing, xlim = axis_lim, ylim = axis_lim, 
     las = 1, xlab = 'Predicted', ylab = 'Observed', main = 'Gaussian process')
lines(x = axis_lim, y = axis_lim, col = 'red')




### 6. USING THE METAMODEL ----

## 6.1 Optimization of StartAge ----

# The Gausian process model will be used to optimize the screening strategy with respect to the starting
# age given an interval of 2 years and for a FIT cut-off of 75. Based on the StartAge and Interval, the
# number of screening rounds (NrScr) is calculated to ensure no screening is performed after the age 90.

# Define the parameter values, except for the NrScr, which will be calculated subsequently
m_experiments_org <- cbind(
  StartAge  = seq(from = 30, to = 90, by = 1),
  Interval  = 2,
  NrScr     = NA,
  FITcutoff = 75
)

# Calculate the corresponding NrScr value
m_experiments_org[ , 'NrScr'] <- 1 + floor((90 - m_experiments_org[ , 'StartAge']) / 2)

# Normalize the experiments to 0-1 for the Gausian process
m_experiments_norm <- t(apply(m_experiments_org, 1, fun_normalize))

# Obtain the outcomes
out_experiments <- predict.GP(object = mm_gp, xnew = m_experiments_norm)$Y_hat

# Plot the NMB as the function of the StartAge
par(mfrow = c(1, 1))
plot(x = m_experiments_org[ , 'StartAge'], y = out_experiments, type = 'l',
     las = 1, xlab = 'StartAge', ylab = 'NMB (EUR 20k per QALY)',
     main = 'NMB as function of StartAge\nInterval = 2, FITcutoff = 75')


## 6.2 Optimization using "optim" ----

# To optimize the screening strategy with regard to all four parameters, optim, which is the standard 
# optimization function in R, will be used in combination with the Gaussian process metamodel.

# The optim function required a single function which it can call to evaluate the outcome for a set of
# parameter values specified through a (named) vector. Therefore, the fn_optim is defined to exactly do
# that by transforming the vector of parameters to a matrix and calling the predict.GP function and
# returning only the outcome estimate.
fn_optim <- function(x) predict.GP(object = mm_gp, xnew = matrix(data = x, ncol = 4))$Y_hat

# The optim function also requires start values for the parameters, for which the current Dutch 
# screening program is used. Bounds on the parameter are also defined, which is conveniently done by
# simplifying these to be 0 and 1, given that the GP metamodel requires normalized inputs. The argument
# "control = list(fnscale = -1)" is used to tell the algorithm that we want to maximize the outcome
# rather than minimize it, which is the default.
optim_out <- optim(
  par     = pars_DutchScreening_norm, 
  fn      = fn_optim, 
  lower   = c(0, 0, 0, 0), 
  upper   = c(1, 1, 1, 1), 
  method  = 'L-BFGS-B',
  control = list(fnscale = -1)
)

# Observing the results
optim_out

# Optimal screening strategy on the original parameter scale
fun_denormalize(optim_out$par)



