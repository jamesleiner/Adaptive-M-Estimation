args<-commandArgs(TRUE)
iter <- as.numeric(args[1])

source("funs.R")
library(sandwich)
library(nnet)
library(MASS)
library(mvtnorm)
library(reshape)
library(randomForest)
library(dplyr)
library(expm)
load("~/work/adaptive_inference/data/dataset_right.Rdata")
load("~/work/adaptive_inference/data/dataset_left.Rdata")


#Generate synthetic features
features <- dataset_right[,3:14]
features[,1:4] <- lapply(features[,1:4],as.factor)

ds <- features

ds$reward <- dataset_right[,"V06WOMADLR"] - dataset_right[,"V00WOMADLR"] 
model = randomForest(reward~.,data= ds,na.action = na.roughfix)

ds$var = (predict(model,ds) - ds$reward)**2
model_var = randomForest(var~.,data= ds,na.action = na.roughfix)




#Simulation parameters
sim_params <- expand.grid(tau=c(1.0),prob_sample=c(1.0),scenario = c(1,2,3,4,5,6),treat_method = c("active","epsilon_greedy","UCB","Thompson","Thompson-Clip"))
sim_params$UPDATE_NUM = 500
sim_params$BATCH_SIZE = 1 
sim_params$NUM_BATCHES = 2000/sim_params$prob_sample/sim_params$BATCH_SIZE

sim_params <- 
  do.call("rbind", replicate(n = 100, sim_params, simplify = FALSE)) 

sim_params$savename = sapply(rownames(sim_params),function(x) paste("/ocean/projects/mth230012p/jleiner/adaptive_inference/sims8_",x,".Rdata",sep=""))

num_per_job = 3
start = (iter-1)*num_per_job + 1

for(num in start:(start+num_per_job-1)){
  print(sim_params[num,])
  tau= sim_params[num,1]
  prob_sample = sim_params[num,2]
  scenario = sim_params[num,3]
  treat_method = sim_params[num,4]
  UPDATE_NUM  = sim_params[num,5]
  BATCH_SIZE = sim_params[num,6]
  NUM_BATCHES  = sim_params[num,7]
  savename = sim_params[num,8]
  INIT_BATCH = UPDATE_NUM
  NUM_EXTERNAL = 1000
  
  synth <- generate_synth(features,NUM_SYNTHETIC = NUM_BATCHES*BATCH_SIZE+UPDATE_NUM)
  synth$reward <- predict(model,synth)
  synth$var  <- predict(model_var,synth)
  feat_external <- generate_synth(features,NUM_SYNTHETIC = NUM_EXTERNAL)
  
  if(scenario == 1){
    treat_effect = c(0,0,1,2,2,3,4,4)
    treat_interact = c(1,-1,1,0,1,1,1,-3)/sd(synth$reward, na.rm=TRUE)
    num_treats = length(treat_effect)
    SD = 1.0
    
    erew = sapply(1:num_treats, function(x) treat_effect[x] + treat_interact[x]*synth$reward)
    vrew = sapply(1:num_treats, function(x) rep(SD, nrow(synth)))
    
  }
  
  
  
  if(scenario == 2){
    treat_effect = c(1,-1,1,0,1,1,1,-3)
    treat_interact = c(1,-1,1,0,1,1,1,-2)/sd(synth$reward)
    treat_var = c(1,1,1,1,1,1,1,1)*0.2
    num_treats = length(treat_effect)
    
    erew = sapply(1:num_treats, function(x) treat_effect[x] + treat_interact[x]*synth$reward)
    vrew = sapply(1:num_treats, function(x) treat_var[x]*synth$var)
  }
  
  
  
  if(scenario == 3){
    treat_effect = c(0,0,1,2,2,3,4,5)
    treat_interact = c(1,-1,1,0,1,1,1,-2)/sd(synth$reward)
    treat_var = c(1,2,3,4,5,5,5,5)*0.2
    num_treats = length(treat_effect)
    
    erew = sapply(1:num_treats, function(x) treat_effect[x] + treat_interact[x]*synth$reward)
    vrew = sapply(1:num_treats, function(x) treat_var[x]*synth$var)
  }
  
  if(scenario == 4){
    treat_effect = c(0,0,1,2,2,3,4,5)
    treat_var = c(1,1,1,1,1,1,1,1)*0.2
    num_treats = length(treat_effect)
    SD = 1.0
    
    erew = sapply(1:num_treats, function(x) treat_effect[x]+ 0*synth$reward)
    vrew = sapply(1:num_treats, function(x) rep(SD, nrow(synth)))
  }
  
  if(scenario == 5){
    treat_effect = c(0,0,1,2,2,3,5,5)
    treat_var = c(1,1,1,1,1,1,1,1)*0.2
    num_treats = length(treat_effect)
    SD = 1.0
    
    erew = sapply(1:num_treats, function(x) treat_effect[x]+ 0*synth$reward)
    vrew = sapply(1:num_treats, function(x) rep(SD, nrow(synth)))
  }
  if(scenario == 6){
    treat_effect = c(0,1,2,3,4,5,6,7)
    treat_var = c(1,1,1,1,1,1,1,1)*0.2
    num_treats = length(treat_effect)
    SD = 1.0
    
    erew = sapply(1:num_treats, function(x) treat_effect[x]+ 0*synth$reward)
    vrew = sapply(1:num_treats, function(x) rep(SD, nrow(synth)))
  }
  
  
  
  treat_labels = sapply(1:num_treats,function(x) paste("T",x,sep=""))
  
  noise = sqrt(vrew)*matrix(rnorm(ncol(erew)*nrow(erew),0,1),nrow(erew),ncol(erew))
  rewards = erew + noise
  colnames(rewards) = treat_labels
  colnames(erew) = treat_labels
  
  
  
  inf <- do_inference(synth[,1:12],feat_external,treat_labels,rewards,erew,
                      NUM_BATCHES=NUM_BATCHES,BATCH_SIZE=BATCH_SIZE,UPDATE_NUM=UPDATE_NUM,INIT_BATCH=INIT_BATCH,
                      prob_sample=prob_sample,tau=tau,savename=savename,
                      treat_method=treat_method,sample_method="active")
  
  
}
