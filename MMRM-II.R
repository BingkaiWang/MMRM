rm(list=ls())
library(tidyverse)
#' Title
#'
#' @param outcomes the outcome matrix where the last column is the final outcome, missing outcomes are NA
#' @param treatment a factor indicating treatment groups
#' @param baseline_variables a matrix with only binary or continuous outcomes. Factor must be converted to dummy variables (leaving out one level) using the model.matrix command.
#' @param stratification_variable_names the column names in baseline variables indicating stratificaiton factors. Each corresponding column must be a binary vector.
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
MMRM2 <- function(outcomes, treatment, baseline_variables, stratification_variable_names = NULL){
  n <- nrow(outcomes)
  treatment_levels <- levels(treatment)
  visits <- colnames(outcomes)
  K <- length(visits) 
  cov_x <- cov(baseline_variables)
  p <- 1 + ncol(baseline_variables)
  nonmissing_indi <- !is.na(outcomes)
  Y_K <- outcomes[,ncol(outcomes)]
  
  # point estimates of treatment effect
  d_K <- data.frame(Y = unlist(Y_K), treatment, baseline_variables)
  lm_fit <- lm(Y~., data = d_K)
  est <- lm_fit$coefficients[2:length(treatment_levels)]
  
  # variance estimation without adjusting for stratified randomization
  Beta <- cov(outcomes, baseline_variables, use = "pairwise.complete.obs") %*% solve(cov_x)
  residual <-scale(outcomes, scale = F) - scale(baseline_variables, scale = F) %*% t(Beta)
  simga <- cov(residual, use = "pairwise.complete.obs") * (length(residual)-1) / (length(residual) - p)
  sigma_j <- map(treatment_levels, function(j){
    residual <-scale(outcomes[treatment == j, ], scale = F) - scale(baseline_variables[treatment == j, ], scale = F) %*% t(Beta)
    cov(residual, use = "pairwise.complete.obs") * (length(residual)-1) / (length(residual) - p)
  })
  V <- map(1:n, function(i){
    m_i <- matrix(0, nrow = K, ncol = K)
    m_i[nonmissing_indi[i,],nonmissing_indi[i,]] <- solve(simga[nonmissing_indi[i,],nonmissing_indi[i,]])
    m_i
  })
  E_V <- Reduce("+", V)/n
  E_V_sigmaj_V <- map(sigma_j, function(mat){
    vv <- map(V,~(. %*% mat %*% .))
    Reduce("+", vv)/n
  })
  Var_tilde_mat <- matrix(0, nrow = length(treatment_levels), ncol = length(treatment_levels))
  for(j in 1: length(treatment_levels)){
    pi_j <- mean(treatment == treatment_levels[j])
    vvv <- solve(E_V) %*% E_V_sigmaj_V[[j]] %*% solve(E_V)
    Var_tilde_mat[j,j] <- vvv[K,K]/pi_j
  }
  
  # variance estimation with adjustment for stratified randomization
  if(is.null(stratification_variable_names)){
    Var_mat <- matrix(NA, nrow = length(treatment_levels), ncol = length(treatment_levels))
  }else{
    Var_mat <- matrix(0, nrow = length(treatment_levels), ncol = length(treatment_levels))
    # compute Var(E[X|S])=E[E[X|S]E[X|S]^T] - E[X]E[X]^T
    ps <- map_dbl(stratification_variable_names, function(j){
      mean(baseline_variables[,j])
    })
    ps <- c(1-sum(ps),ps)
    EXS <- matrix(NA, nrow = length(ps), ncol = ncol(baseline_variables))
    EXS[1,] <- colMeans(baseline_variables[rowSums(baseline_variables[,stratification_variable_names])==0,])
    for(j in 2:length(ps)){
      str_v_name <- stratification_variable_names[j-1]
      EXS[j,] <- colMeans(baseline_variables[baseline_variables[,str_v_name]==1, ])
    }
    Var_EXS <- -colMeans(baseline_variables) %*% t(colMeans(baseline_variables))
    for(j in 1: length(ps)){
      Var_EXS <- Var_EXS + ps[j] * EXS[j,] %*% t(EXS[j,])
    }
    # compute the variance difference term
    Y_K <- outcomes[,ncol(outcomes)]
    b_K <- cov(Y_K, baseline_variables, use = "pairwise.complete.obs") %*% solve(cov_x)
    bKj_bK <- map_dfr(treatment_levels,function(j){
      as.data.frame(cov(Y_K[treatment== j,], baseline_variables[treatment==j,], use = "pairwise.complete.obs") %*% solve(cov_x)-b_K)
    }) %>% as.matrix
    var_diff <- - bKj_bK %*% Var_EXS %*% t(bKj_bK)
    for(j in 1:nrow(var_diff)){
      pi_j <- mean(treatment == treatment_levels[j])
      var_diff[j,j] <- var_diff[j,j]+ bKj_bK[j,] %*% Var_EXS %*% bKj_bK[j,]/pi_j
    }
    Var_mat <- Var_tilde_mat - var_diff
  }
  summary_result <- matrix(NA, nrow = length(treatment_levels)-1, ncol = 3)
  for(j in 1:nrow(summary_result)){
    summary_result[j,1] <- est[j]
    summary_result[j,2] <- sqrt(c(-1, 1) %*% Var_tilde_mat[c(1,j+1), c(1,j+1)] %*% c(-1, 1)/(n-p))
    summary_result[j,3] <- sqrt(c(-1, 1) %*% Var_mat[c(1,j+1), c(1,j+1)] %*% c(-1, 1)/(n-p))
  }
  return(summary_result)
}


# generating data for testing
n <- 400
n_arm <- 4
n_stratum <- 4
n_addi_covariates <- 1
n_sim <- 1000
sim_result <- map(1:n_sim, function(kk){
  d <- data.frame(id = as.character(1:n),
                  treatment = as.factor(rep(1:n_arm, each = n/n_arm)),
                  stratum = as.factor(rep(rep(1:4, each = n/n_arm/n_stratum), n_arm)),
                  baseline = rnorm(n, mean = 0, sd = 1),
                  Y1 = rnorm(n, mean = 0, sd = 1),
                  Y2 = rnorm(n, mean = 0, sd = 1),
                  Y3 = rnorm(n, mean = 0, sd = 1))
  missing_pattern <- sample(1:4, size = n, replace = T)
  for(i in 1: n){
    if(missing_pattern[i] == 1){
      d$Y2[i] <- NA
      d$Y3[i] <- NA
    } else if(missing_pattern[i] == 2){
      d$Y3[i] <- NA
    } else if(missing_pattern[i] == 3){
      d$Y2[i] <- NA
    }
  }
  d <- d %>% 
    pivot_longer(Y1:Y3, names_to = "visit", values_to = "outcome") %>%
    .[complete.cases(.),]
  
  
  # function
  d_wide <- pivot_wider(d, names_from = visit, values_from = outcome)
  outcomes <- d_wide[,c("Y1", "Y2", "Y3")]
  treatment <- d_wide$treatment
  baseline_variables <- model.matrix(~ stratum + baseline, data = d_wide)[,-1]
  MMRM2(outcomes, treatment, baseline_variables, stratification_variable_names = c("stratum2", "stratum3", "stratum4"))
})

est1 <- map_dbl(sim_result, ~.[1,1])
sd_tilde <- map_dbl(sim_result, ~.[1,2])
sd <- map_dbl(sim_result, ~.[1,3])
sd(est1)
mean(sd_tilde)
mean(sd)

est2 <- map_dbl(sim_result, ~.[2,1])
sd_tilde <- map_dbl(sim_result, ~.[2,2])
sd <- map_dbl(sim_result, ~.[2,3])
sd(est2)
mean(sd_tilde)
mean(sd)
