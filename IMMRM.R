rm(list=ls())
library(tidyverse)
#' Title
#'
#' @param outcomes the outcome matrix where the last column is the final outcome, missing outcomes are NA
#' @param treatment a factor indicating treatment groups
#' @param baseline_variables a matrix with only binary or continuous outocmes. 
#' Factor must be coverted to dummy variables (leaving out one level) using the model.matrix command.
#'
#' @return
#' @export
#'
#' @examples
#' 
#' 
#' 
IMMRM <- function(outcomes, treatment, baseline_variables){
  n <- nrow(outcomes)
  treatment_levels <- levels(treatment)
  visits <- colnames(outcomes)
  K <- length(visits) 
  cov_x <- cov(baseline_variables)
  p <- 1 + ncol(baseline_variables)
  nonmissing_indi <- !is.na(outcomes)
  est_j <- map(treatment_levels, function(j){
    Beta <- cov(outcomes[treatment == j, ], baseline_variables[treatment == j, ], use = "pairwise.complete.obs") %*% solve(cov_x)
    colMeans(outcomes[treatment == j, ], na.rm = T) #- colMeans(baseline_variables[treatment == j, ]) %*% t(Beta)
  })
  sigma_j <- map(treatment_levels, function(j){
    Beta <- cov(outcomes[treatment == j, ], baseline_variables[treatment == j, ], use = "pairwise.complete.obs") %*% solve(cov_x)
    residual <-scale(outcomes[treatment == j, ], scale = F) - scale(baseline_variables[treatment == j, ], scale = F) %*% t(Beta)
    cov(residual, use = "pairwise.complete.obs") * (length(residual)-1) / (length(residual) - p)
  })
  r <- map_dfr(treatment_levels, function(j){
    data.frame(cov(outcomes[treatment == j, K], baseline_variables[treatment == j, ], use = "pairwise.complete.obs") %*% solve(cov_x))
  }) %>% as.matrix
  E_V_j <- map(sigma_j, function(mat){
    V <- map(1:n, function(i){
      m_i <- matrix(0, nrow = K, ncol = K)
      m_i[nonmissing_indi[i,],nonmissing_indi[i,]] <- solve(mat[nonmissing_indi[i,],nonmissing_indi[i,]])
      m_i
    })
    Reduce("+", V)/length(V)
  })
  Var_mat <- r %*% cov_x %*% t(r)
  for(j in 1:length(treatment_levels)){
    pi_j <- mean(treatment == treatment_levels[j])
    Var_mat[j,j] <- Var_mat[j,j] + solve(E_V_j[[j]])[K,K]/pi_j
  }
  est <- map_dbl(est_j, ~.[K])
  summary_result <- NULL
  for(j in 1:(length(treatment_levels)-1)){
    for(jj in (j+1):length(treatment_levels)){
      diff <- est[jj] - est[j]
      sd <- sqrt(c(-1, 1) %*% Var_mat[c(j,jj), c(j,jj)] %*% c(-1, 1)/n)
      row_name <- paste("treatment", jj, "vs treatment", j)
      summary_result <- rbind(summary_result, matrix(c(diff, sd), nrow = 1, ncol = 2, dimnames = list(row_name, c("est", "sd"))))
    }
  }
  return(summary_result)
}

# generating data for testing
n <- 300
n_arm <- 3
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
  IMMRM(outcomes, treatment, baseline_variables)
})


est1 <- map_dbl(sim_result, ~.[2,1])
sd <- map_dbl(sim_result, ~.[2,2])
sd(est1)
mean(sd)
