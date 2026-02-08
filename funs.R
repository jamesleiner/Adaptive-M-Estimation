###################################################################################################
# Helper functions to run experiments
###################################################################################################

# -----------------------------------------
# pr_max: Compute per-arm probabilities of being the maximum under independent Normals
#
# Inputs:
#   - x: numeric vector of means for each arm
#   - y: numeric vector of variances for each arm (same length/order as x)
#
# Outputs:
#   - numeric vector p where p[i] = P(arm i attains the maximum)
# -----------------------------------------
pr_max <- function(x,y) sapply(1:length(x), function(i) prob_max(x,y,i))


# -----------------------------------------
# prob_max: Probability that a given arm is the argmax among independent Normal arms
#
# Inputs:
#   - means: numeric vector of means for X_1,...,X_n
#   - var: numeric vector of variances for X_1,...,X_n (assumes independence)
#   - ind: integer index (1..n) of the arm being tested as the maximum
#
# Outputs:
#   - scalar probability P(X_ind >= X_j for all j != ind)
# -----------------------------------------
prob_max <- function(means,var,ind){
  n<- length(means)
  diff<- -1*diag(n-1)


    # Build a (n-1) x n "difference" matrix that maps X -> (X_ind - X_j) for all j != ind.
  # Start from -I_{n-1} (for the "-X_j" pieces) and insert a +1 column for X_ind.
  if(ind!= n & ind!=1){
    diff = cbind(diff[,1:(ind-1)], rep(1,n-1),diff[,(ind+1):n-1])
  } else if( ind==1){
    diff = cbind( rep(1,n-1),diff)
  } else{
    diff = cbind(diff,rep(1,n-1))
  }
  # Differences D = diff %*% X are jointly Normal with:
  #   E[D] = diff %*% means
  #   Cov[D] = diff %*% diag(var) %*% t(diff)  (independence assumed)
  new_mu <- as.vector(diff %*% means)
  new_Sigma <- diff %*%diag(var) %*% t(diff)
  
  pr <- pmvnorm(lower=rep(0,n-1), mean = new_mu, sigma = new_Sigma)[[1]]
  return(pr)
}


# -----------------------------------------
# calc_st_lin: Sandwich-type variance estimator when a linear model for the rewards is assumed
#
# Inputs:
#   - theta: numeric vector of coefficients for linear model (typically length 2: intercept, slope)
#   - ref_feat: data.frame/matrix of covariates at which to evaluate predictions (n rows)
#   - model: fitted outcome model for E[Y | X, treat] supporting predict()
#   - model_uncertainty: fitted model for conditional variance/uncertainty proxy supporting predict()
#   - evaluation_policy: numeric vector of evaluation policy probabilities over treatments (length num_treats)
#   - prob_label: numeric vector of sampling/labeling probabilities per row 
#   - prob_treat: numeric matrix (n x num_treats) of treatment assignment probabilities
#
# Outputs:
#   - covariance matrix for theta (dimension matches length(theta))
# -----------------------------------------

calc_st_lin <- function(theta, ref_feat, model, model_uncertainty,
                        evaluation_policy, prob_label, prob_treat) {
  
  # Predict conditional mean for each (row, treatment)
  fmom_ext <- sapply(
    1:num_treats,
    function(x) predict(
      model,
      data.frame(ref_feat, treat = factor(treat_labels[x], levels = treat_labels))
    )
  )
  # Predict conditional variance/uncertainty proxy for each (row, treatment)
  var_ext <- sapply(
    1:num_treats,
    function(x) predict(
      model_uncertainty,
      data.frame(ref_feat, treat = factor(treat_labels[x], levels = treat_labels))
    )
  )
  
  colnames(fmom_ext) <- treat_labels
  colnames(var_ext)  <- treat_labels
  
  n <- nrow(fmom_ext)
  k <- ncol(fmom_ext)  # should equal num_treats
  
  ## ---------- Replace pivot_longer(fmom_ext) ----------
  # We want the same ordering as pivot_longer(everything()):
  # row 1 across all treatments, then row 2 across all treatments, etc.
  
  row_idx  <- rep(seq_len(n), each = k)       # 1,1,...,2,2,...,n,n
  col_idx  <- rep(seq_len(k), times = n)      # 1,2,...,k,1,2,...,k,...
  Treat    <- colnames(fmom_ext)[col_idx]
  FMOM     <- fmom_ext[cbind(row_idx, col_idx)]
  
  a <- data.frame(
    Treat = Treat,
    FMOM  = FMOM
  )
  
  # Extract numeric treatment index (e.g. "T1" -> 1)
  a$Treat <- as.numeric(sub("^[A-Za-z]+", "", a$Treat))
  
  X <- cbind(1, a$Treat)
  
  # term2: variance contribution scaled by evaluation_policy^2 / propensity
  # (variance/uncertainty divided by prob of observing that treatment)
  resid <- evaluation_policy * as.vector(a$FMOM - X %*% theta)
  
  ## ---------- term2 and its "long" version (replacing pivot_longer(term2)) ----------
  term2 <- t(apply(var_ext, 1, function(x) x * evaluation_policy^2 / prob_treat))
  colnames(term2) <- treat_labels
  
  Var <- term2[cbind(row_idx, col_idx)]
  b <- data.frame(
    Treat = Treat,
    Var   = Var
  )
  
  # Sandwich-type variance:
  # - cov(resid * X) captures empirical score variability
  # - second term adds model-based conditional variance contribution
  var_st <- cov(resid * X) + t(X) %*% (b$Var * X) / n
  
  return(var_st)
}



calc_st <- function(theta,ref_feat,model,model_uncertainty,evaluation_policy,prob_label,prob_treat) {
  
  fmom_ext = sapply(1:num_treats,function(x) predict(model,data.frame(ref_feat,treat=factor(treat_labels[x],levels=treat_labels))))
  var_ext = sapply(1:num_treats,function(x) predict(model_uncertainty,data.frame(ref_feat,treat=factor(treat_labels[x],levels=treat_labels))))
  
  term1 <- t(apply(fmom_ext, 1, function(x) evaluation_policy*(x - theta)))
  term2 <- t(apply(var_ext, 1, function(x) x*evaluation_policy**2/prob_treat)) 
  term3 <- t(apply(fmom_ext, 1, function(x) evaluation_policy**2*(x - theta)))
  
  var_st = cov(term1) + diag(colMeans(term2-term3))
  return(var_st)
}



# -----------------------------------------
# generate_synth: Generate synthetic covariates by MVN sampling + multinomial models for multi-class factors
#
# Inputs:
#   - covariates: data.frame of observed covariates (numeric + factor columns)
#   - NUM_SYNTHETIC: integer number of synthetic rows to generate
#
# Outputs:
#   - data.frame of synthetic covariates with the same column set/order as covariates
# -----------------------------------------

generate_synth <- function(covariates,NUM_SYNTHETIC=10000){
  # Identify factor columns
  cat <- names(which(sapply(covariates,is.factor)))
  cat_multi <- names(which(lapply(sapply(covariates,levels),length) >2))
  
  # Remove multi-class factors for the MVN step (avoid awkward dummy expansions here)
  covariates_uni <- covariates %>% 
    dplyr::select(-any_of(cat_multi))
  
  # Identify multi-class factor columns (levels > 2) handled via multinomial models
  covariates_numeric <- model.matrix(~.,covariates_uni)[,-1]

  # Fit multivariate Normal parameters on expanded covariates
  mu <- colMeans(covariates_numeric)
  Sigma <- cov(covariates_numeric)
  
  # Sample synthetic rows from MVN
  synth = mvrnorm(n=NUM_SYNTHETIC,mu,Sigma)
  colnames(synth) <- colnames(covariates_uni)
  synth[,which(colnames(synth) %in% cat)] <- (synth[,which(colnames(synth) %in% cat)] > 0.5)
  synth <- data.frame(synth)
  
   # Convert back to factors with original levels for binary factors
  for(var in colnames(synth)){
    if(var %in% cat){
      synth[,var] <- factor(synth[,var],labels=levels(covariates[,var]))
    }
  }
  

  # For multi-class factors: fit multinomial on real data (conditioning on covariates_uni) and predict for synth
  for(var in cat_multi){
    ds <- cbind(covariates_uni,covariates[,var])
    colnames(ds) <- c(colnames(covariates_uni),var)
    glm.fit=multinom(reformulate(colnames(covariates_uni),response=var), data=ds)
    synth[,var] = predict(glm.fit,synth)
  }
  synth <- synth[names(covariates)]
  return(synth)
}


# -----------------------------------------
# CI: Compute Wald confidence intervals using sandwich covariance and evaluate average width/coverage
#
# Inputs:
#   - reg: fitted regression object (e.g., lm) with $coefficients
#   - true_reg: fitted regression object providing "true" coefficients for coverage evaluation
#   - sig_level: numeric confidence level (e.g., 0.9 for 90% CI)
#
# Outputs:
#   - list with:
#       * vcov: sandwich covariance matrix
#       * stderr: standard errors (sqrt(diag(vcov)))
#       * interval: merged data.frame of intervals and true coefficients
#       * width: mean interval width
#       * cover: empirical coverage indicator averaged across coefficients
# -----------------------------------------
CI <- function(reg,true_reg,sig_level=0.9){
  vcov = sandwich(reg)
  theta= reg$coefficients
  intervals = cbind(theta + qnorm((1-sig_level)/2)*sqrt(diag(vcov)),theta-qnorm((1-sig_level)/2)*sqrt(diag(vcov)))
  width = mean(abs(intervals[,2] - intervals[,1]))

  # Merge with true coefficients to compute coverage
  intervals = merge(intervals,data.frame(true_reg$coefficients),by=0)
  cover = mean((intervals[,4] < intervals[,3])&(intervals[,4] > intervals[,2]))
  return(list(vcov=vcov,stderr=sqrt(diag(vcov)),vcov=vcov,interval=intervals,width=width,cover=cover))
}



# -----------------------------------------
# CI_cust: Wald confidence intervals from user-supplied (theta, vcov) and coverage vs true coefficients
#
# Inputs:
#   - theta: numeric vector of coefficients to intervalize
#   - vcov: covariance matrix for theta
#   - true_reg: fitted regression object providing "true" coefficients for coverage evaluation
#   - sig_level: numeric confidence level (e.g., 0.9 for 90% CI)
#   - projected: logical flag (interface; not used here)
#
# Outputs:
#   - list with:
#       * vcov: supplied covariance matrix
#       * stderr: standard errors (sqrt(diag(vcov)))
#       * interval: merged data.frame of intervals and true coefficients
#       * width: mean interval width
#       * cover: empirical coverage indicator averaged across coefficients
# -----------------------------------------
CI_cust <- function(theta,vcov,true_reg,sig_level=0.9,projected=FALSE){
  d <- length(theta)
  intervals = cbind(theta + qnorm((1-sig_level)/2)*sqrt(diag(vcov)),theta-qnorm((1-sig_level)/2)*sqrt(diag(vcov)))
  width = mean(abs(intervals[,2] - intervals[,1]))
  intervals = merge(intervals,data.frame(true_reg$coefficients),by=0)
  cover = mean((intervals[,4] < intervals[,3])&(intervals[,4] > intervals[,2]))
  return(list(vcov=vcov,stderr=sqrt(diag(vcov)),vcov=vcov,interval=intervals,width=width,cover=cover))
}



# -----------------------------------------
# CI_proj: Construct pointwise (Wald) or projected (chi-square) intervals and evaluate width/coverage
#
# Inputs:
#   - theta: numeric vector of coefficients
#   - vcov: covariance matrix for theta
#   - true_reg: fitted regression object providing "true" coefficients for coverage evaluation
#   - sig_level: numeric confidence level (e.g., 0.90 for 90%)
#   - projected: logical; if TRUE uses sqrt(chi^2_df) cutoff, else uses z cutoff
#
# Outputs:
#   - list with:
#       * vcov: supplied covariance matrix
#       * stderr: standard errors (sqrt(diag(vcov)))
#       * interval: merged data.frame with lower/upper and true coefficients
#       * width: mean interval width
#       * cover: empirical coverage indicator averaged across coefficients
# -----------------------------------------
CI_proj <- function(theta, vcov, true_reg, sig_level = 0.90, projected = TRUE) {
  
  theta <- as.numeric(theta)
  d <- length(theta)
  
  # Select cutoff depending on projected vs naive
  if (projected) {
    # Projected CI: sqrt(chi^2_d)
    cutoff <- sqrt(qchisq(sig_level, df = d))
  } else {
    # Naive Wald CI: z-score
    cutoff <- qnorm(1 - (1 - sig_level)/2)
  }
  
  # Radii for each coefficient
  radii <- cutoff * sqrt(diag(vcov))
  
  # Intervals: theta_j ± radius
  intervals <- cbind(theta - radii, theta + radii)
  colnames(intervals) <- c("lower", "upper")
  
  # Mean interval width
  width <- mean(intervals[, "upper"] - intervals[, "lower"])
  
  # Merge with true regression coefficients for coverage
  intervals_df <- merge(
    as.data.frame(intervals),
    data.frame(true_reg$coefficients),
    by = 0
  )
  
  # Coverage indicator
  cover <- mean(
    intervals_df[, "true_reg.coefficients"] < intervals_df[, "upper"] &
      intervals_df[, "true_reg.coefficients"] > intervals_df[, "lower"]
  )
  
  # Return EXACT same structure and names as your original
  return(list(
    vcov = vcov,
    stderr = sqrt(diag(vcov)),
    interval = intervals_df,
    width = width,
    cover = cover
  ))
}



# -----------------------------------------
# do_inference: Run batched adaptive treatment experiment
#
# Inputs:
#   - features: data.frame/matrix of covariates for all units (rows = time/units)
#   - feat_external: data.frame/matrix of external reference covariates for variance estimation
#   - treat_labels: character vector of treatment labels (e.g., c("T1","T2",...))
#   - rewards: matrix/data.frame of realized rewards, columns correspond to treatments
#   - erew: matrix/data.frame of expected rewards, columns correspond to treatments (used as source of truth)
#   - NUM_BATCHES: integer number of batches for each round of treamtnet
#   - BATCH_SIZE: number of individuals treated in each round
#   - UPDATE_NUM: integer frequency (in batches) of refitting outcome/uncertainty models
#   - INIT_BATCH: integer number of initial samples to fit first models
#   - prob_sample: baseline sampling/labeling probability
#   - sig_level: confidence level for CI evaluation
#   - tau: mixing/regularization parameter for probabilities (toward uniform / baseline)
#   - sig_level_UCB: confidence level used inside UCB treatment policy
#   - epsilon: exploration rate for epsilon-greedy / UCB policies
#   - savename: character path for saving intermediate results (.Rdata)
#   - treat_method: character specifying treatment assignment rule ("active","epsilon_greedy","UCB","Thompson","Thompson-Clip")
#   - sample_method: character specifying labeling rule ("active" uses uncertainty; otherwise reward-based)
#
# Outputs:
#   - list `res` containing width/coverage trajectories and run metadata; also saved to disk during the run
# -----------------------------------------


do_inference <- function(features,feat_external,treat_labels,rewards,erew,
                          NUM_BATCHES=10,BATCH_SIZE=100,UPDATE_NUM=10,INIT_BATCH=1000,
                          prob_sample=0.5,sig_level = 0.9,tau=0.5,sig_level_UCB = 0.9,epsilon=0.1,
                         savename="sim_results.Rdata",treat_method="epsilon_greedy",sample_method="active"){
  num_treats = length(treat_labels)
  treat_nums <- as.numeric(sub("^[A-Za-z]+", "", treat_labels))
  evaluation_policy = rep(1/num_treats, num_treats)
  prob_treat = matrix(rep(1/num_treats,num_treats*BATCH_SIZE),BATCH_SIZE,num_treats)
  prob_label = rep(prob_sample,BATCH_SIZE)
  BIN_FLAG =   length(unique(rewards[,1]))==2
  true_df = data.frame(melt(erew))[,-1]
  colnames(true_df) <- c("treat","reward")
  true_df$treat_num <- as.numeric(sub("^[A-Za-z]+", "", true_df$treat))
  
  #true models to get oracle parameters to compare against
  true_reg <- lm(reward~0+treat,data=true_df)
  true_reg2 <- lm(reward~treat_num,data=true_df)
  prob_treat = matrix(rep(1/num_treats,num_treats*INIT_BATCH),INIT_BATCH,num_treats)
  cov_batch = features[1:INIT_BATCH,]
  rew_batch = rewards[1:INIT_BATCH,]
  treat_batch = apply(t(apply(prob_treat,1, function(x) rmultinom(1,1,x))),1,which.max)
  rew_batch_actual <- as.vector(sapply(1:nrow(rew_batch),function(x) rew_batch[x,treat_batch[x]]))
  
  df_init <- data.frame(cov_batch, treat=treat_batch,reward=rew_batch_actual)
  df_init$treat <- as.factor(sapply(df_init$treat, function(x) treat_labels[x]))
  df_init$treat <- relevel(df_init$treat,treat_labels[1])
  model = randomForest(reward~.,data=df_init)
  df_init$resid = abs(df_init$reward - model$predicted)
  model_uncertainty = randomForest(resid~.,data=  df_init[ , which(names(df_init) %in% c(colnames(features),"treat","resid"))])
  
  
  
  for(i in 1:NUM_BATCHES){
    print(i)
    
    #get history
    cov_batch = features[(INIT_BATCH+BATCH_SIZE*(i-1)+1):(BATCH_SIZE*i+INIT_BATCH),]
    rew_batch = rewards[(INIT_BATCH+BATCH_SIZE*(i-1)+1):(BATCH_SIZE*i+INIT_BATCH),]
    erew_batch = erew[(INIT_BATCH+BATCH_SIZE*(i-1)+1):(BATCH_SIZE*i+INIT_BATCH),]
    vrew_batch = vrew[(INIT_BATCH+BATCH_SIZE*(i-1)+1):(BATCH_SIZE*i+INIT_BATCH),]
    
    
    pred_batch = sapply(1:num_treats,function(x) predict(model,data.frame(cov_batch,treat=factor(treat_labels[x],levels=treat_labels))))
    uncertainty_batch = sapply(1:num_treats,function(x) predict(model_uncertainty,data.frame(cov_batch,treat=factor(treat_labels[x],levels=treat_labels))))
    if(is.vector(uncertainty_batch)){
      pred_batch <- t(as.matrix(pred_batch))
      uncertainty_batch <- t(as.matrix(uncertainty_batch))
    }
    
    
    if(treat_method == "active"){
      prob_treat = uncertainty_batch/rowSums(uncertainty_batch)
      prob_treat = prob_treat*tau + (1-tau)/num_treats
      
    }
    if(treat_method == "epsilon_greedy"){
      prob_treat = rep(epsilon/(num_treats-1),num_treats)
      prob_treat <- matrix(rep(prob_treat,BATCH_SIZE),nrow=BATCH_SIZE,ncol=num_treats)
      maxes <- apply(pred_batch,1,which.max)
      for(k in 1:BATCH_SIZE){
        prob_treat[k,maxes[k]] = 1- epsilon
      }
    }
    if(treat_method == "UCB"){
      UCB = pred_batch - qnorm((1-sig_level_UCB)/2)*uncertainty_batch
      maxes <- apply(UCB,1,which.max)
      prob_treat = rep(epsilon/(num_treats-1),num_treats)
      prob_treat <- matrix(rep(prob_treat,BATCH_SIZE),nrow=BATCH_SIZE,ncol=num_treats)
      for(k in 1:BATCH_SIZE){
        prob_treat[k,maxes[k]] = 1- epsilon
      }
    }
    if(treat_method == "Thompson"){
      prob_treat <- t(sapply(1:BATCH_SIZE, function(x) pr_max(as.vector(pred_batch[x,]),as.vector(uncertainty_batch[x,]))))
    }
    if(treat_method == "Thompson-Clip"){
      prob_treat <- t(sapply(1:BATCH_SIZE, function(x) pr_max(as.vector(pred_batch[x,]),as.vector(uncertainty_batch[x,]))))
      prob_treat = (prob_treat < 0.01) * 0.01 + prob_treat
      prob_treat = prob_treat/sum(prob_treat)
    }
    
    
    avg_uncertainty = rowSums(uncertainty_batch*prob_treat)
    avg_reward = rowSums(pred_batch*prob_treat)
    sd_uncertainty <- sd(as.vector(uncertainty_batch))
    
    if((i ==1) |  (i %% UPDATE_NUM == 0)){
      treat_vec = prob_treat
    } else{
      treat_vec = rbind(treat_vec,prob_treat)
    }
    
    
    if(i==1){
      uncertainty_vec = avg_uncertainty
      reward_vec = avg_reward
      sd_vec = sd_uncertainty
    } else{
      uncertainty_vec = c(uncertainty_vec,avg_uncertainty)
      reward_vec = c(reward_vec,avg_reward)
      sd_vec = c(sd_vec,sd_uncertainty)
    }
    
    
    if(sample_method=="active"){
      prob_label = prob_sample*avg_uncertainty/mean(uncertainty_vec)
      prob_label = sapply(prob_label, function(x) min(x,1))
      
      #regualrize with uniform sampling rule 
      prob_label = prob_label*tau + (1-tau)*prob_sample
      
    } else{
      prob_label = prob_sample*exp(avg_reward)/mean(exp(reward_vec))
      prob_label = sapply(prob_label, function(x) min(x,1))
      
    }
    
    
    label_batch = rbinom(BATCH_SIZE,1,prob_label)
    
    if(i==1){
      prob_vec = prob_label
    } else{
      prob_vec = c(prob_vec, prob_label)
    }
    
    
    
    treat_batch = apply(t(apply(prob_treat,1, function(x) rmultinom(1,1,x))),1,which.max)
    label_batch = rbinom(BATCH_SIZE,1,prob_label)
    
    if(BATCH_SIZE > 1){
      rew_batch_actual <- as.vector(sapply(1:nrow(rew_batch),function(x) rew_batch[x,treat_batch[x]]))
    } else{
      rew_batch_actual <- rew_batch[treat_batch]
    }
    
    
    
    for(j in 1:num_treats){
      if(BATCH_SIZE >1){
        df <-data.frame(cov_batch, treat=treat_labels[j],reward = rew_batch[,j],erew=erew_batch[,j],vrew=vrew_batch[,j],in_samp=label_batch,in_treat= (j==treat_batch),samp_prob = prob_label,eval_prob=prob_treat[,j],treat_prob=prob_treat[,j],opt_rew = apply(erew_batch,1,max))
      }
      else{
        df <-data.frame(cov_batch, treat=treat_labels[j],reward = rew_batch[j],erew=erew_batch[j],vrew=vrew_batch[j],in_samp=label_batch,in_treat= (j==treat_batch),samp_prob = prob_label,eval_prob=evaluation_policy[j],treat_prob=prob_treat[j],opt_rew= min(which(erew_batch == max(erew_batch))))
      }
      if(j==1){
        df_all_init =df
      }
      else{
        df_all_init = rbind(df_all_init,df)
      }
    }
    df_all_init$treat_num <- as.numeric(sub("^[A-Za-z]+", "", df_all_init$treat))
    df_all_init$treat <- as.factor(df_all_init$treat)
    df_all_init$treat <- relevel(df_all_init$treat,treat_labels[1])
    df_all_init$pred <- predict(model,df_all_init)
    df_all_init$reward_dr <- df_all_init$pred + (df_all_init$reward - df_all_init$pred)*df_all_init$in_samp*df_all_init$in_treat/df_all_init$treat_prob
    df_all_init$resid <- (df_all_init$pred - df_all_init$reward)**2
    
    
    if(i == 1){
      df_all = df_all_init
      
    }
    if(i > 1){
      df_all = rbind(df_all,df_all_init)
    }
    
    
    
    df_actual <- df_all[df_all$in_samp&df_all$in_treat,]
    df_samples <- df_all[df_all$in_samp ==1,]
    
    metrics_naive <- list(width=NA,cover=NA)
    metrics_naive_proj <- list(width=NA,cover=NA)
    metrics_naive_cont <- list(width=NA,cover=NA)
    metrics_naive_contproj <- list(width=NA,cover=NA)
    metrics_ipw <- list(width=NA,cover=NA)
    metrics_ipw_proj <-  list(width=NA,cover=NA)
    metrics_ipw_cont <-  list(width=NA,cover=NA)
    metrics_ipw_contproj <-  list(width=NA,cover=NA)
    metrics_sqipw <- list(width=NA,cover=NA)
    metrics_sqipw_proj <-  list(width=NA,cover=NA)
    metrics_sqipw_cont <-  list(width=NA,cover=NA)
    metrics_sqipw_contproj <-  list(width=NA,cover=NA)
    
    
    if(length(unique(df_samples$treat)) > 1){
      reg_oracle <- lm(reward~0+treat,data=df_samples)
      metrics_oracle <- CI(reg_oracle,true_reg)
      
      reg_dr <- lm(reward_dr~0+treat,data=df_samples)
      regc_oracle <- lm(reward~treat_num,data=df_samples)
      
      regc_dr <- lm(reward_dr~treat_num,data=df_samples)
    }
    
    if(length(unique(df_actual$treat)) > 1){
      reg_naive <- lm(reward~0+treat,data=df_actual)
      regc_naive <- lm(reward~treat_num,data=df_actual)
      vcov_naive = vcovHC(reg_naive, type = "HC0")
      vcovc_naive = vcovHC(regc_naive, type = "HC0")
      metrics_naive <- CI_proj(reg_naive$coefficients,vcov_naive,true_reg,projected=FALSE,sig_level = sig_level)
      metrics_naive_proj <- CI_proj(reg_naive$coefficients,vcov_naive,true_reg,projected=TRUE,sig_level = sig_level)
      metrics_naive_cont <- CI_proj(regc_naive$coefficients,vcovc_naive,true_reg2,projected=FALSE,sig_level = sig_level)
      metrics_naive_contproj <- CI_proj(regc_naive$coefficients,vcovc_naive,true_reg2,projected=TRUE,sig_level = sig_level)
      
      reg_ipw <- lm(reward~0+treat,data=df_actual,weights = 1/df_actual$treat_prob)
      regc_ipw <- lm(reward~treat_num,data=df_actual,weights = 1/df_actual$treat_prob)
      vcov_ipw = vcovHC(reg_ipw, type = "HC0")
      vcovc_ipw = vcovHC(regc_ipw, type = "HC0")
      metrics_ipw <- CI_proj(reg_ipw$coefficients,vcov_ipw,true_reg,projected=FALSE,sig_level = sig_level)
      metrics_ipw_proj <- CI_proj(reg_ipw$coefficients,vcov_ipw,true_reg,projected=TRUE,sig_level = sig_level)
      metrics_ipw_cont <- CI_proj(regc_naive$coefficients,vcovc_ipw,true_reg2,projected=FALSE,sig_level = sig_level)
      metrics_ipw_contproj <- CI_proj(regc_naive$coefficients,vcovc_ipw,true_reg2,projected=TRUE,sig_level = sig_level)
      
      reg_sqipw <- lm(reward~0+treat,data=df_actual,weights = sqrt(1/df_actual$treat_prob))
      regc_sqipw <- lm(reward~treat_num,data=df_actual,weights = sqrt(1/df_actual$treat_prob))
      vcov_sqipw = vcovHC(reg_sqipw, type = "HC0")
      vcovc_sqipw = vcovHC(regc_sqipw, type = "HC0")
      metrics_sqipw <- CI_proj(reg_sqipw$coefficients,vcov_sqipw,true_reg,projected=FALSE,sig_level = sig_level)
      metrics_sqipw_proj <- CI_proj(reg_sqipw$coefficients,vcov_sqipw,true_reg,projected=TRUE,sig_level = sig_level)
      metrics_sqipw_cont <- CI_proj(regc_sqipw$coefficients,vcovc_sqipw,true_reg2,projected=FALSE,sig_level = sig_level)
      metrics_sqipw_contproj <- CI_proj(regc_sqipw$coefficients,vcovc_sqipw,true_reg2,projected=TRUE,sig_level = sig_level)
    }
    
    if(i==1){
      vcov1_i = calc_st(reg_dr$coefficients, feat_external, model,model_uncertainty,evaluation_policy,prob_label,prob_treat)
      vcov2_i = calc_st(reg_dr$coefficients, features[1:(BATCH_SIZE*i+INIT_BATCH),], model,model_uncertainty,evaluation_policy,prob_label,prob_treat)
      vcov3_i = calc_st(reg_dr$coefficients, rbind(df_init[,1:12],df_all[(!df_all$in_samp)&df_all$in_treat,1:12]), model,model_uncertainty,evaluation_policy,prob_label,prob_treat)
      vt_i = sqrtm(solve(vcov1_i))
      vt2_i = sqrtm(solve(vcov2_i))
      vt3_i = sqrtm(solve(vcov3_i))
      
      
      vcovc1_i = calc_st_lin(regc_dr$coefficients, feat_external, model,model_uncertainty,evaluation_policy,prob_label,prob_treat)
      vcovc2_i = calc_st_lin(regc_dr$coefficients, features[1:(BATCH_SIZE*i+INIT_BATCH),], model,model_uncertainty,evaluation_policy,prob_label,prob_treat)
      vcovc3_i = calc_st_lin(regc_dr$coefficients, rbind(df_init[,1:12],df_all[(!df_all$in_samp)&df_all$in_treat,1:12]), model,model_uncertainty,evaluation_policy,prob_label,prob_treat)
      vtc_i = sqrtm(solve(vcovc1_i))
      vtc2_i = sqrtm(solve(vcovc2_i))
      vtc3_i = sqrtm(solve(vcovc3_i))
      
      score = 0
      score2 = 0
      score3 = 0
      
      scorec = 0
      scorec2 = 0
      scorec3 = 0
      vt_seq =list(vt_i)
      vt2_seq =list(vt2_i)
      vt3_seq =list(vt3_i)
    }
    
    if(label_batch[1] == 1){
      s_t = diag(as.vector(evaluation_policy))
      score = score + vt_i %*% s_t
      score2 = score2 + vt2_i%*% s_t
      score3 = score3 + vt3_i%*% s_t
      
      s_t = matrix(c(1,sum(evaluation_policy*treat_nums),sum(evaluation_policy*treat_nums),sum(evaluation_policy*treat_nums**2)),2,2)
      scorec = scorec + vtc_i %*% s_t
      scorec2 = scorec2 + vtc2_i%*% s_t
      scorec3 = scorec3 + vtc3_i%*% s_t
      
      vt_seq =c(vt_seq,list(vt_i))
      vt2_seq =c(vt_seq,list(vt2_i))
      vt3_seq = c(vt_seq,list(vt3_i))
      
      
      a = solve(score/sqrt(nrow(df_actual)))
      vcov_ext = t(a) %*% a
      rownames(vcov_ext) <- names(reg_dr$coefficients)
      colnames(vcov_ext) <- names(reg_dr$coefficients)
      
      a = solve(score2/sqrt(nrow(df_actual)))
      vcov_split = t(a) %*% a
      rownames(vcov_split) <- names(reg_dr$coefficients)
      colnames(vcov_split) <- names(reg_dr$coefficients)
      
      a = solve(score3/sqrt(nrow(df_actual)))
      vcov_reuse = t(a) %*% a
      rownames(vcov_reuse) <- names(reg_dr$coefficients)
      colnames(vcov_reuse) <- names(reg_dr$coefficients)
      
      a = solve(scorec/sqrt(nrow(df_actual)))
      vcovc_ext = t(a) %*% a
      rownames(vcovc_ext) <- names(regc_dr$coefficients)
      colnames(vcovc_ext) <- names(regc_dr$coefficients)
      
      a = solve(scorec2/sqrt(nrow(df_actual)))
      vcovc_split = t(a) %*% a
      rownames(vcovc_split) <- names(regc_dr$coefficients)
      colnames(vcovc_split) <- names(regc_dr$coefficients)
      
      a = solve(scorec3/sqrt(nrow(df_actual)))
      vcovc_reuse = t(a) %*% a
      rownames(vcovc_reuse) <- names(regc_dr$coefficients)
      colnames(vcovc_reuse) <- names(regc_dr$coefficients)
      
      
      metrics_maipw <-   CI_proj(reg_dr$coefficients,vcov_ext ,true_reg,projected=FALSE,sig_level = sig_level)
      metrics_maipw_proj <-   CI_proj(reg_dr$coefficients,vcov_ext ,true_reg,projected=TRUE,sig_level = sig_level)
      
      metrics_maipw2 <-   CI_proj(reg_dr$coefficients,vcov_split ,true_reg,projected=FALSE,sig_level = sig_level)
      metrics_maipw2_proj <-   CI_proj(reg_dr$coefficients,vcov_split ,true_reg,projected=TRUE,sig_level = sig_level)
      
      metrics_maipw3 <-   CI_proj(reg_dr$coefficients,vcov_reuse,true_reg,projected=FALSE,sig_level = sig_level)
      metrics_maipw3_proj <-   CI_proj(reg_dr$coefficients,vcov_reuse,true_reg,projected=TRUE,sig_level = sig_level)
      
      
      metrics_cmaipw <-   CI_proj(regc_dr$coefficients,vcovc_ext ,true_reg2,projected=FALSE,sig_level = sig_level)
      metrics_cmaipw_proj <-   CI_proj(regc_dr$coefficients,vcovc_ext ,true_reg2,projected=TRUE,sig_level = sig_level)
      
      metrics_cmaipw2 <-   CI_proj(regc_dr$coefficients,vcovc_split ,true_reg2,projected=FALSE,sig_level = sig_level)
      metrics_cmaipw2_proj <-   CI_proj(regc_dr$coefficients,vcovc_split ,true_reg2,projected=TRUE,sig_level = sig_level)
      
      metrics_cmaipw3 <-   CI_proj(regc_dr$coefficients,vcovc_reuse,true_reg2,projected=FALSE,sig_level = sig_level)
      metrics_cmaipw3_proj <-   CI_proj(regc_dr$coefficients,vcovc_reuse,true_reg2,projected=TRUE,sig_level = sig_level)
      
    }
    
    width_i = c(metrics_oracle$width,metrics_naive$width,
                metrics_ipw$width,metrics_ipw_proj$width,
                metrics_sqipw$width,metrics_sqipw_proj$width,
                metrics_maipw$width,metrics_maipw_proj$width,
                metrics_maipw2$width,metrics_maipw2_proj$width,
                metrics_maipw3$width,metrics_maipw3_proj$width)
    
    cover_i = c(metrics_oracle$cover,metrics_naive$cover,
                metrics_ipw$cover,metrics_ipw_proj$cover,
                metrics_sqipw$cover,metrics_sqipw_proj$cover,
                metrics_maipw$cover,metrics_maipw_proj$cover,
                metrics_maipw2$cover,metrics_maipw2_proj$cover,
                metrics_maipw3$cover,metrics_maipw3_proj$cover)
    
    widthc_i = c(metrics_oracle$width,metrics_naive_cont$width,
                 metrics_ipw_cont$width,metrics_ipw_contproj$width,
                 metrics_sqipw_cont$width,metrics_sqipw_contproj$width,
                 metrics_cmaipw$width,metrics_cmaipw$width,
                 metrics_cmaipw2$width,metrics_cmaipw2$width,
                 metrics_cmaipw3$width,metrics_cmaipw3$width)
    
    coverc_i = c(metrics_oracle$cover,metrics_naive_cont$cover,
                 metrics_ipw_cont$cover,metrics_ipw_contproj$cover,
                 metrics_sqipw_cont$cover,metrics_sqipw_contproj$cover,
                 metrics_cmaipw$cover,metrics_cmaipw$cover,
                 metrics_cmaipw2$cover,metrics_cmaipw2$cover,
                 metrics_cmaipw3$cover,metrics_cmaipw3$cover)
    
    num_i = c(nrow(df_all),sum(df_all$in_treat *(!df_all$in_samp)),sum(df_all$in_treat *(df_all$in_samp)))
    
    if(i==1){
      num = num_i
      width = width_i
      cover = cover_i
      widthc = widthc_i
      coverc = coverc_i
    } else{
      num = rbind(num,num_i)
      width = rbind(width,width_i)
      cover = rbind(cover,cover_i)
      widthc = rbind(widthc,widthc_i)
      coverc = rbind(coverc,coverc_i)
    }
    
    
    if(i > 2){
      colnames(cover) <- c("Oracle","Naive","IPW","IPW (Projected)","SQ-IPW","SQ-IPW (Projected)",
                           "MAIPW-External","MAIPW-External (Projected)","MAIPW-Split","MAIPW-Split (Projected)","MAIPW-Reuse","MAIPW-Reuse (Projected)")
      colnames(width) <- c("Oracle","Naive","IPW","IPW (Projected)","SQ-IPW","SQ-IPW (Projected)",
                           "MAIPW-External","MAIPW-External (Projected)","MAIPW-Split","MAIPW-Split (Projected)","MAIPW-Reuse","MAIPW-Reuse (Projected)")
      
      rownames(width) <- NULL
      rownames(cover) <- NULL
      rownames(widthc) <- NULL
      rownames(coverc) <- NULL
      
      res = list(width=width,cover=cover,witdthc =widthc,coverc=coverc,num= num,tau=tau,prob_sample=prob_sample,
                 BATCH_SIZE=BATCH_SIZE,NUM_BATCHES=NUM_BATCHES,treat_method=treat_method,sample_method=sample_method)
      save(res,file=savename)
    }
    
    
    
    if((i %% UPDATE_NUM == 0)&(i < NUM_BATCHES)){
      start = Sys.time()
      df_samples = rbind(df_init,df_actual[ , which(names(df_actual) %in% c(colnames(features),"treat","reward","resid"))])
      print(nrow(df_samples))
      model = randomForest(reward~.,data=  df_samples[ , which(names(df_samples) %in% c(colnames(features),"treat","reward"))])
      model_uncertainty = randomForest(resid~.,data=  df_samples[ , which(names(df_samples) %in% c(colnames(features),"treat","resid"))])
      
      
      vcov1_i = calc_st(reg_dr$coefficients, feat_external, model,model_uncertainty,evaluation_policy,prob_label,prob_treat)
      vcov2_i = calc_st(reg_dr$coefficients, features[1:(BATCH_SIZE*i+INIT_BATCH),], model,model_uncertainty,evaluation_policy,prob_label,prob_treat)
      vcov3_i = calc_st(reg_dr$coefficients, rbind(df_init[,1:12],df_all[(!df_all$in_samp)&df_all$in_treat,1:12]), model,model_uncertainty,evaluation_policy,prob_label,prob_treat)
      
      vcovc1_i = calc_st_lin(regc_dr$coefficients, feat_external, model,model_uncertainty,evaluation_policy,prob_label,prob_treat)
      vcovc2_i = calc_st_lin(regc_dr$coefficients, features[1:(BATCH_SIZE*i+INIT_BATCH),], model,model_uncertainty,evaluation_policy,prob_label,prob_treat)
      vcovc3_i = calc_st_lin(regc_dr$coefficients, rbind(df_init[,1:12],df_all[(!df_all$in_samp)&df_all$in_treat,1:12]), model,model_uncertainty,evaluation_policy,prob_label,prob_treat)
      
      vt_i = sqrtm(solve(vcov1_i))
      vt2_i = sqrtm(solve(vcov2_i))
      vt3_i = sqrtm(solve(vcov3_i))
      vtc_i = sqrtm(solve(vcovc1_i))
      vtc2_i = sqrtm(solve(vcovc2_i))
      vtc3_i = sqrtm(solve(vcovc3_i))
      print(Sys.time()- start)
    }
    
  }
  
  
  return(res)
}

