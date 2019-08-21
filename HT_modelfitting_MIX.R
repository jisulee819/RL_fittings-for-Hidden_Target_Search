
# For fitting MAP/MIG and saving posterior heatmap
# Modified 190507 by JISU
# Modified 190817 by JISU (ver.13)


library(MASS)
library(readr)
library(mvtnorm)
library(corpcor)
library(Hmisc)
library(reshape2)
library(ggplot2)


#### Model Fitting : saving 1.inferences from the models 2.likelihood ####

modelfitting_MIX <- function(s,lambda, roundn, subjID){
  
  #### Prior Distribution ####
  # prior of the reward probability is analytically calculated in advance 
  pr <- read.csv("/Volumes/clmnlab/HT/behavior_data/raw_data/priorreduced2.csv")
  prior <<- data.matrix(pr[,-1])
  
  #### Settings ####
  
  subdata <- behavdata[((roundn-1)*10+1):(roundn*10),]
  treasure_location <- c((subdata$treasureX[10]/ratio)+nx/2,(subdata$treasureY[10]/ratio)+ny/2)
  
  sigma <- matrix(c(s^2, 0, 0, s^2),nrow=2)
  rho <- 0
  ratio <- 1024/nx
  
  prior_history <- prior
  entropy_diff_history <- rep(NA, ny)
  infer_history <- c(NA, NA)
  rew_history <- 0
  prob_reward_history <- NA
  likelihood <- NA
  lkh.ratio <- NA
  mig.history <- matrix(NA,nrow=10,ncol=2)
  migs.history <- matrix(0, nrow = 10*10, ncol=2)
  map.history <- matrix(NA,nrow=10,ncol=2)
  map.likelihood <- NA
  mig.likelihood <- NA
  z.map <- NA
  z.mig <- NA
  l.map <- NA
  l.mig <- NA
  f.map <- NA
  f.mig <- NA
  
  n = 10
  
  #### Iteration ####
  for (i in 1:10){
    print(paste0("Trial ",i))
    
    #### Radius and Reward ####
    if (i == 1){rad <- subdata$Radius[1]/16
    } else {rad <- subdata$Radius[i-1]/16}
    
    rew <- subdata$Smiley[i]
    rew_history[i] <- rew
    
    #### Truncate Inferences ####
    infer <- c((subdata$inferX[i]/ratio)+nx/2, (subdata$inferY[i]/ratio)+ny/2)
    if (infer[1] < 1){
      infer[1] = 1
    } else if (infer[1] > 32){
      infer[1] = 32
    }
    
    if (infer[2] < 1){
      infer[2] = 1
    } else if (infer[2] > 24) {
      infer[2] = 24
    }
    
    infer[1] <- round(infer[1])
    infer[2] <- round(infer[2])
    
    
    #### MIG fitting and saving ####
    mig_dataset <- inference_as(prior, sigma, rad)
    infer_model.migs <- mig_dataset$`inference mig`
    infer_model.mig <- infer_model.migs[1,]
    mig.history[i,] <- infer_model.mig
    migs.history[((i-1)*10+1):(i*10),] <- infer_model.migs
    
    entropy_diff_normalized <- mig_dataset$`entropy_difference`/sum(mig_dataset$`entropy_difference`)
    entropy_diff_history <- rbind(entropy_diff_history, entropy_diff_normalized)
    
    
    #### Probability of Reward ####
    prob_reward <- 0
    for(k in 1:nx){for(l in 1:ny){if (sqrt((k-infer[1])^2 + (l-infer[2])^2) < sqrt(sigma[1,1])){
      prob_reward <- prob_reward + prior[k,l]
    }}}
    prob_reward_history[i] <- prob_reward
    
    
    #### MAP fitting and saving ####
    infer_model.map <- inference_bs(prior) 
    map.history[i,] <- infer_model.map
    
    #migdist <- dist(rbind(infer, infer_model.mig))
    #mapdist <- dist(rbind(infer, infer_model.map))

    
    #### Bayesian Update of Prior ####
    posterior <- bayesian_update(prior, infer, sigma, rew)
    posterior <- posterior / sum(posterior)
    
    mat.prior <- melt(prior)
    mat.ent <- melt(entropy_diff_normalized)
    
    
    #### Calculate the Index of Model fitting ####
    # z.index: prior & entropy difference in the normalized space
    # l.index: prior & entropy difference linearly scaled btw [-1,1]
    # f.index: prior & entropy difference Fisher transformed
    
    z.prior <- (prior - mean(mat.prior[,3])) / sd(mat.prior[,3])
    z.ent <- (entropy_diff_normalized - mean(mat.ent[,3])) / sd(mat.ent[,3])
    
    z.map[i] <- z.prior[infer[1], infer[2]]
    z.mig[i] <- z.ent[infer[1], infer[2]]

    mat.z.prior <- melt(z.prior)
    mat.z.ent <- melt(z.ent)
    
    l.map[i] <- (2 / (max(mat.z.prior[,3]) - min(mat.z.prior[,3])))  * (z.map[i] - min(mat.z.prior[,3])) - 1
    l.mig[i] <- (2 / (max(mat.z.ent[,3]) - min(mat.z.ent[,3])))  * (z.mig[i] - min(mat.z.ent[,3])) - 1
    
    f.map[i] <- (1/2) * log((1+l.map[i])/(1-l.map[i])) 
    f.mig[i] <- (1/2) * log((1+l.mig[i])/(1-l.mig[i])) 
   
    map.lkh <- log(prior[infer[1],infer[2]]+0.0000001)
    mig.lkh <- log(entropy_diff_normalized[infer[1],infer[2]]+0.0000001)
    isEmpty <- function(x) {
    return(length(x)==0)
    }
    
    # index.lkh : prior & entropy difference in the original space
    
    map.likelihood[i] <- map.lkh
    mig.likelihood[i] <- mig.lkh
    
    if (!isEmpty(map.lkh - mig.lkh)){
    lkh.ratio[i] <- (map.lkh - mig.lkh) 
    } else {lkh.ratio[i] <- NA}
    
    prior <- posterior
    prior_history <- rbind(prior_history, data.matrix(prior))
    infer_history <- rbind(infer_history, infer)
    s <- s * lambda^(rew)
    sigma <- matrix(c(s^2, 0, 0, s^2),nrow=2) 
    
  }

  #### Save ####
  data = list('treasureX' = treasure_location[1],'treasureY' = treasure_location[2], "prior history" = prior_history, 'inference' = infer_history[-1,],
              'reward history' = rew_history, 'prob reward history' = prob_reward_history, "sigma" = s, "lambda" = lambda, 'entropy diff' = entropy_diff_history[-1,],
              'mig history' = mig.history, 'migs.history' = migs.history, 'map history' = map.history, 'lkh ratio' = lkh.ratio, 'map.likelihood' = map.likelihood, 'mig.likelihood' = mig.likelihood,
              'prior_history' = prior_history, 'entropy_diff_history' = entropy_diff_history, 'z.map' = z.map, 'z.mig' = z.mig, 'f.map' = f.map, 'f.mig' = f.mig)
  
  return(data)
}


