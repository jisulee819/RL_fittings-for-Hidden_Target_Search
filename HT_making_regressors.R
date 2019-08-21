########################################################################################
# For making the first regressors: probability of reward, lable(MAP or MIG)
#190508 Written by JISU
#190520 Modified by JISU
########################################################################################

library(MASS)
library(readr)
library(mvtnorm)
library(corpcor)
library(Hmisc)
library(reshape2)
library(ggplot2)
library(writexl)

### For Expected Probability 
set.radius <- 384/16
set.discount.rate <- 0.78

nx<-64
ny<-48
n=10

bayesian_update <- function(prior, inference, sigma, reward){
  
  likelihood <- matrix(data = NA, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
  posterior <- matrix(data = NA, nrow = nx, ncol = ny, byrow = FALSE, dimnames = NULL)
  
  for(i in 1:nx){
    for(j in 1:ny){
      mu <- matrix(c(i,j), nrow=2)
      likelihood[i,j] <- dmvnorm(t(inference), mu, sigma)}}
  maxl <- max(likelihood)
  if (reward == 0){
    likelihood <- maxl-likelihood}
  for(i in 1:nx){
    for(j in 1:ny){
      posterior[i,j] <- likelihood[i,j] * prior[i,j]}}
  return(posterior)
  
}

inference_bs <- function(prior){ ##
  inference <- which(prior == max(prior), arr.ind=TRUE) #
  inference <- inference[sample(1:nrow(inference), 1),]
  return(inference)
}

reward_bs <- function(inference, radius, treasure_location){ ##
  if (sqrt((treasure_location[1]-inference[1])^2 + (treasure_location[2]-inference[2])^2) < radius) return(1)
  else return(0)
}

modelfitting_BA <- function(s, lambda, roundn, subjID){
  rho <- 0 
  ratio <- 1024/nx
  
  subdata <- pilotdata[pilotdata$subjectID == subjID,] 
  subdata <- subdata[((roundn-1)*10+1):(roundn*10),]
  treasure_location <- c((subdata$treasureX[10]/ratio)+nx/2,(subdata$treasureY[10]/ratio)+ny/2)
  
  sigma <- matrix(c(s^2, 0, 0, s^2),nrow=2) #
  pr <- read.csv("/Volumes/clmnlab/HT/behavior_data/raw_data/priorreduced.csv")
  pr <- pr[,-1]
  prior <- matrix(NA, nrow=nx, ncol=ny)
  for (a in 1:nx){
    for (b in 1:ny){
      prior[a,b] <- pr[a,b]
    }
  }
  #prior.sigma <- matrix(c(priorsigma^2, 0, 0, priorsigma^2),nrow=2) #
  
  #for (a in c(1:nx)){
  #  for (b in c(1:ny)){
  #    prior[a,b] <-  dmvnorm(c(a,b), c(paraX[1],paraY[1]), prior.sigma)
  #  }}
  
  prior_history <- prior
  infer_history <- c(NA, NA) 
  infer_model_history <- c(NA, NA) 
  rew_history <- 0
  distfromtarget <- NA
  prob_reward_history <- NA
  entropy_prior <- NA
  fitting.error <- NA
  distfromMAP <- NA
  
  for (i in 1:10){
    
    infer <- c((subdata$inferX[i]/ratio)+nx/2, (subdata$inferY[i]/ratio)+ny/2)
    if (infer[1] < 1){
      infer[1] = 1
      } else if (infer[1] > 64){
      infer[1] = 64
    }
    
    if (infer[2] < 1){
      infer[2] = 1
    } else if (infer[2] > 48){
      infer[2] = 48
    }
    
    infer[1] <- round(infer[1])
    infer[2] <- round(infer[2])
    
    print(infer)
    infer_model <- inference_bs(prior)
    
    
    distfromMAP[i] <- sqrt((infer[1] - infer_model[1])^2 + (infer[2] - infer_model[2])^2)
    prob_reward <- 0
    for(k in 1:nx){for(l in 1:ny){if (sqrt((k-infer[1])^2 + (l-infer[2])^2) < sqrt(sigma[1,1])){
      prob_reward <- prob_reward + prior[k,l]
    }}}
    prob_reward_history[i] <- prob_reward
    
    entropy_prior[i] <- sum(prior*log2(1/prior), na.rm=TRUE)
    
    # when prior is uniform, entropy = 11.58496
    # when prior is 0s with only one 1, entropy = 0.
    
    infer_history <- rbind(infer_history, infer)
    infer_model_history <- rbind(infer_model_history, infer_model)
    
    if (i == 1){rad <- subdata$Radius[1]/16
    }else{rad <- subdata$Radius[i-1]/16}
    
    rew <- subdata$Smiley[i]
    rew_history[i] <- rew
    
    distfromtarget[i] <- subdata$Distance[i]/16
    posterior <- bayesian_update(prior, infer, sigma, rew)
    posterior <- posterior / sum(posterior)
    fitting.error[i] <- prior[infer[1],infer[2]]
    prior <- posterior
    prior_history <- rbind(prior_history, data.matrix(prior))
    
    s <- s*lambda^(rew)
    sigma <- matrix(c(s^2, 0, 0, s^2),nrow=2) 
    
  }  
  data = list('treasureX' = treasure_location[1],'treasureY' = treasure_location[2], "prior history" = prior_history, 'inference' = infer_history[-1,], 
              'inference of model' = infer_model_history[-1,], 'reward history' = rew_history, 'distance' = distfromtarget, 'distance from MAP' = distfromMAP,
              'prob reward history' = prob_reward_history, "sigma" = s, "lambda" = lambda, 'fitting error' = fitting.error, 'entropy_prior' = entropy_prior)
  
  return(data)
}




modelfitting.with.error <- function(s,lambda,subjID){
  model_err <- NA
  results <- NA
  roundn <- max(pilotdata$Round[pilotdata$subjectID == subjID])
  prior_history <- rep(NA,ny)
  round_error <- NA
  infer_tp <- NA
  infer_type <- NA
  for (i in 1:roundn){ 
    print(paste0("Round number ",i))
    data <- modelfitting_BA(s,lambda,i,subjID)
    inferences <- cbind(data$inference, data$`inference of model`)
    
    ################# criteria needs to be modified ####################
    criteria <- 5
    ####################################################################
    infer_tp <- ifelse(data$`distance from MAP`< criteria, "exploit", "explore")
    pilotdata$infer_type <- infer_tp
    infer_type[(10*(i-1)+1):(i*10)] <- infer_tp
    
    fitting_error <- (data$`fitting error`)#[data$`distance from MAP`< criteria])
    round_error <- -sum(log(fitting_error+0.00000001)) 
    model_err[i] <- round_error
    
    resultsA <- (data.frame("ID" = subjID,"round" = i,"trial" = c(rep(1:10)), inferences, "treasureX" = data$'treasureX', "treasureY" = data$'treasureY', "reward" = data$`reward history`,"distance" = data$distance, "prob_reward" = data$`prob reward history`, "entropy" = data$`entropy_prior`, "sigma" = s, "lambda" = lambda))
    results <- rbind(results, resultsA)
    prior_history <- rbind(prior_history, data$`prior history`)}
  rownames(results) <- NULL
  return(list("results" = cbind(results[-1,],infer_type), "model err" = sum(model_err), "prior history" = prior_history[-1,]))
  print(model_err)}






#setwd("/Users/jisu/Documents/Hidden Target/Behavioral Data/raw_data")


subjectpool <- c(23)

for (subjectID in subjectpool){
  
  if (subjectID < 10){
    filename <- paste0("HT0",subjectID)
  } else {filename <- paste0("HT",subjectID)
  }
  
print(filename)
#fileID <- paste0("0",subjectID)
setwd(paste0("/Volumes/clmnlab/HT/behavior_data/raw_data/",filename))
fileID <- paste0(subjectID)
pilotdata <- read.csv(paste0(filename ,".csv"), header=FALSE)
names(pilotdata)<- c("subjectID", "Trial", "Smiley", "treasureX", "treasureY",	"inferX",	"inferY",	"Radius", "Distance",	"RT", "Certainty")
pilotdata$Round <- c(rep(1:36, each = 10))
pilotdata$subjectID <- subjectID

data_opt_h <- modelfitting.with.error(set.radius, set.discount.rate, subjectID)
data <- data_opt_h$results

names(data) <- c("ID", "round", "trial", "inferX", "inferY", "infermodelX", "infermodelY", "treasureX", "treasureY", "reward", "distance","prob_reward", "entropy", "sigma","lambda","infer_type")
pilot_regressor <- data.frame("trial" = c(1:360), "RP" = data_opt_h$results$prob_reward, "nothing" = data$infer_type, 
                              "RP entropy" = (data_opt_h$results$prob_reward * log(1/data_opt_h$results$prob_reward) + (1-data_opt_h$results$prob_reward) * log(1/(1-data_opt_h$results$prob_reward)))
                              )


################################## to make regressors #####################################
write_xlsx(pilot_regressor, paste0("/Users/jisu/Documents/Hidden Target/Behavioral Data/regressor/", filename,"_regressor.xlsx"))
write_xlsx(pilot_regressor, paste0("/Volumes/clmnlab/HT/regressors/MAP_regressor/",filename,"_regressor.xlsx"))
###########################################################################################

}



#dim(posteriors)
#posteriors <- data_opt_h$`prior history`[(64*11*27+1):(64*11*28),]
#for (j in 1:11){
#post <- posteriors[(64*j-63):(64*j), ]
#ent <- sum((post+0.0000001) * log(1/(post+0.0000001)))
#print(ent)
#}
