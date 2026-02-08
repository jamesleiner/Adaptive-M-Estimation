
# Read job array index from command line (used to distribute work across cluster jobs)
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

# Define the real-data reward as a change from baseline in KL grade (follow-up minus baseline)
ds$reward <- dataset_right[,"V06WOMADLR"] - dataset_right[,"V00WOMADLR"] 
model = randomForest(reward~.,data= ds,na.action = na.roughfix)

# Calculate square residuals, we will predict this for outcome variance modelling.; 
ds$var = (predict(model,ds) - ds$reward)**2
model_var = randomForest(var~.,data= ds,na.action = na.roughfix)



###################################################
#Simulation parameters
# Grid over:
# - scenario: corresponds to different configurations of treatment effect, heteroskedasticity, and misspecification. Described in detail below
# - treat_method: adaptive treatment policy used inside do_inference:
#    * active corresponds to treating in proportion to model-based uncertainty
#    * epsilon-greedy, UCB and Thompson correspond to the standard Bandit sampling algorithms
#    * Thompson-Clip refers to Thompson sampling but with probabiltiies constrained to [0.05,0.95]
# - tau mixes the sampling rule with uniform sampling over treatments for the active method specifically. As we only choose tau = 1 here, active is just doing uniform sampling
# - UPDATE_NUM refers to the number of iteration before the predictive and uncertainty models are re-updated
# - BATCH_SIZE refers to the number of individuals treated with the same sampling rule. We only consider single batch experiments
# - NUM_BATCHES refers to the experiment size. We scale this based on the probability of sampling a data point as opposed to holding it out for variance calculations. We choose
#  the total number of samples based on the batch size and probability of including the data point in our sample 
###################################################
sim_params <- expand.grid(tau=c(1.0),prob_sample=c(1.0),scenario = c(1,2,3,4,5,6),treat_method = c("active","epsilon_greedy","UCB","Thompson","Thompson-Clip"))
sim_params$UPDATE_NUM = 500
sim_params$BATCH_SIZE = 1 
sim_params$NUM_BATCHES = 2000/sim_params$prob_sample/sim_params$BATCH_SIZE

sim_params <- 
  do.call("rbind", replicate(n = 100, sim_params, simplify = FALSE)) 

sim_params$savename = sapply(rownames(sim_params),function(x) paste("/ocean/projects/mth230012p/jleiner/adaptive_inference/sims8_",x,".Rdata",sep=""))

#We process three scenarios per job 
num_per_job = 3
start = (iter-1)*num_per_job + 1

for(num in start:(start+num_per_job-1)){
  #Load parameters
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
  
  #Generate synthetic data 
  synth <- generate_synth(features,NUM_SYNTHETIC = NUM_BATCHES*BATCH_SIZE+UPDATE_NUM)
  synth$reward <- predict(model,synth)
  synth$var  <- predict(model_var,synth)
  feat_external <- generate_synth(features,NUM_SYNTHETIC = NUM_EXTERNAL)
  
  
  # Scenario 1: Model Misspecification, Interaction Effects,
  if(scenario == 1){
    treat_effect = c(0,0,1,2,2,3,4,4)
    treat_interact = c(1,-1,1,0,1,1,1,-3)/sd(synth$reward, na.rm=TRUE)
    num_treats = length(treat_effect)
    SD = 1.0
    
    erew = sapply(1:num_treats, function(x) treat_effect[x] + treat_interact[x]*synth$reward)
    vrew = sapply(1:num_treats, function(x) rep(SD, nrow(synth)))
    
  }
  

  # Scenario 2: Model Misspecification, Interaction Effects
  if(scenario == 2){
    treat_effect = c(1,-1,1,0,1,1,1,-3)
    treat_interact = c(1,-1,1,0,1,1,1,-2)/sd(synth$reward)
    treat_var = c(1,1,1,1,1,1,1,1)*0.2
    num_treats = length(treat_effect)
    
    erew = sapply(1:num_treats, function(x) treat_effect[x] + treat_interact[x]*synth$reward)
    vrew = sapply(1:num_treats, function(x) treat_var[x]*synth$var)
  }
  
  
  # Scenario 3: Model Misspecification, Interaction Effects, Heteroskedasticity 
  if(scenario == 3){
    treat_effect = c(0,0,1,2,2,3,4,5)
    treat_interact = c(1,-1,1,0,1,1,1,-2)/sd(synth$reward)
    treat_var = c(1,2,3,4,5,5,5,5)*0.2
    num_treats = length(treat_effect)
    
    erew = sapply(1:num_treats, function(x) treat_effect[x] + treat_interact[x]*synth$reward)
    vrew = sapply(1:num_treats, function(x) treat_var[x]*synth$var)
  }
  
  # Scenario 4: Model Misspecification, No Interaction Effects, Homoskedasticity
  if(scenario == 4){
    treat_effect = c(0,0,1,2,2,3,4,5)
    treat_var = c(1,1,1,1,1,1,1,1)*0.2
    num_treats = length(treat_effect)
    SD = 1.0
    
    erew = sapply(1:num_treats, function(x) treat_effect[x]+ 0*synth$reward)
    vrew = sapply(1:num_treats, function(x) rep(SD, nrow(synth)))
  }
  
  # Scenario 5: Model Misspecification, No Interaction Effects, Homoskedasticity
  if(scenario == 5){
    treat_effect = c(0,0,1,2,2,3,5,5)
    treat_var = c(1,1,1,1,1,1,1,1)*0.2
    num_treats = length(treat_effect)
    SD = 1.0
    
    erew = sapply(1:num_treats, function(x) treat_effect[x]+ 0*synth$reward)
    vrew = sapply(1:num_treats, function(x) rep(SD, nrow(synth)))
  }
  # Scenario 5: Correct Specification, No Interaction Effects, Homoskedasticity
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
