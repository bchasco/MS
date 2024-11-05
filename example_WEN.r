rm(list=ls())
library(dplyr)
library(tidyr)
data(WEN)
# WEN <- WEN %>%
#   mutate(stageID = ifelse(stageID=="sub",1,ifelse(stageID=="yr",2,3))) %>%
#   mutate(stageID = ifelse(is.na(stageID),3,stageID))

WEN <- WEN %>%
  mutate(stageID = ifelse(stageID=="sub" | stageID == "yr",1,stageID)) %>%
  mutate(stageID = ifelse(is.na(stageID),2,stageID))# %>%
  # filter(loc %in% c("Trib","First_Trap","MCJ"))

data <- list(
  n_fish = max(WEN$id),
  id = WEN$id,
  n_locs = length(unique(WEN$loc)),
  n_states = length(unique(WEN$stageID)),
  obs = as.integer(WEN$stageID),
  n = WEN$n
)

init_loc <- rep(1,data$n_fish)
init_state <- rep(1,data$n_fish)

for(i in 1:data$n_fish){
  tmp_obs <- data$obs[data$id == i]
  tmp_obs[is.na(tmp_obs)] <- 3
  init_loc[i] <- min(which(tmp_obs < data$n_states, arr.ind = TRUE))
  init_state[i] <- tmp_obs[init_loc[i]]
}


data$n_trans <-   (data$n_states * (data$n_states - 1)/2) - data$n_states + 1
data$init_state <- init_state
data$init_loc <- init_loc


# Create parameter list if not provided
parameters <- list(
  log_survival = base::matrix(-2, data$n_states-1, data$n_locs-2),
  log_transitions = base::matrix(0, data$n_trans, data$n_locs-2),
  logit_detection = base::matrix(-4, data$n_states-1, data$n_locs-2)
)

library(RTMB)

# Define the model function using R syntax
f <- function(params){

  RTMB::getAll(data,
               params)


  # Extract parameters
  log_survival <- params$log_survival         # (n_states-1) x (n_locs-1)
  log_transitions <- params$log_transitions   # (n_states*(n_states-1)/2 - n_states) x (n_locs-1)
  logit_detection <- params$logit_detection   # (n_states-1) x n_locs

  # Convert detection probabilities from logit scale
  detection <- RTMB::plogis(logit_detection)

  # Initialize negative log-likelihood
  nll <- rep(0,n_fish)

  # Create transition matrices for each location transition
  #There are n -1 transition locations because you don't need a transition for the
  #last location

  T <- list()
  for(l in 1:(n_locs-2)) {
    T[[l]] <- matrix(0, n_states, n_states)

    # For each state (except dead state)
    for(i in 1:(n_states-1)) {
      T_l <- rep(0,n_states)

      # Set survival probability
      T_l[i] <- log_survival[i,l]

      # Set transition probabilities
      trans_idx <- 1
      if((i+1)<n_states){
        for(j in (i+1):(n_states-1)) { #off-diagonal
          T_l[j] <- (log_transitions[trans_idx,l])
          trans_idx <- trans_idx + 1
        }
      }

      # Add baseline category (dead state) as 0 on log scale
      T_l[n_states] <- 1

      # Fill transition probabilities
      T[[l]][i,] <- exp(T_l[i:n_states])/sum(exp(T_l[i:n_states]))
    }

    # Dead state stays dead
    T[[l]][n_states,n_states] <- 1
    # print(l)
    # print(T[[l]])
  }
  T[[n_locs-1]] <- T[[n_locs-2]]

  # Create detection matrices for each location
  D <- list()
  for(l in 1:(n_locs-2)) {
    D[[l]] <- matrix(0, n_states, n_states)
  #
  #   # Fill detection probabilities
    for(i in 1:(n_states-1)) {
      D[[l]][i,i] <- detection[i,l]                    # Detected in correct state
      D[[l]][i,n_states] <- 1 - detection[i,l]         # Not detected
    }

    # Dead state is always "not detected"
    D[[l]][n_states,n_states] <- 1
  }

  D[[n_locs-1]] <- D[[n_locs-2]]
  # Forward algorithm with scaling for numerical stability
  for(i in 1:n_fish) {
    obs_i <- obs[id == i]
    n_i <- length(obs_i)

    delta <- rep(0,n_states)
    delta[init_state[i]] <- 1  # Initialize according to initial state vector

    # Forward algorithm through remaining locations
    prods <- diag(1,n_states)
    for(t in init_loc[i]:(n_i-1)) {
      # State transition & detection
      prods <- prods %*% T[[t]] %*% diag(D[[t]][,obs_i[t+1]])
    }

    # if(init_loc[i]==1){
      nll[i] <- -1. * log(sum(t(delta) %*% prods)) * n[i]
    # }
      print(obs_i)
      print(nll[i])
  }

  REPORT(T)
  REPORT(D)
  REPORT(nll)
  return(sum(nll))
}

obj <- RTMB::MakeADFun(f,
                       parameters,
                       silent = FALSE)

opt <- nlminb(obj$par, obj$fn, obj$gr)
rep <- obj$report()

print(rep$T)

