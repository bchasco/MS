library(RTMB)

# Define the model function using R syntax
f <- function(parms){

  RTMB::getAll(data,
               parms)

  # Declare data inputs
  n_fish <- max(data$id)      # Number of fish
  n_locs <- length(unique(data$n_locs))      # Number of locations
  n_states <- length(unique(data$stageID))  # Number of states
  obs <- data$stageID            # Observation matrix
  U <- rep(1,max(data$id))                # Initial state vector

  # Extract parameters
  log_survival <- params$log_survival         # (n_states-1) x (n_locs-1)
  log_transitions <- params$log_transitions   # (n_states*(n_states-1)/2 - n_states) x (n_locs-1)
  logit_detection <- params$logit_detection   # (n_states-1) x n_locs

  # # Convert detection probabilities from logit scale
  # detection <- plogis(logit_detection)
  #
  # # Initialize negative log-likelihood
  # nll <- 0
  #
  # # Create transition matrices for each location transition
  # T <- vector("list", n_locs-1)
  # for(l in 1:(n_locs-1)) {
  #   T[[l]] <- matrix(0, n_states, n_states)
  #
  #   # For each state (except dead state)
  #   for(i in 1:(n_states-1)) {
  #     log_probs <- numeric(n_states)
  #
  #     # Set survival probability
  #     log_probs[i] <- log_survival[i,l]
  #
  #     # Set transition probabilities
  #     trans_idx <- 1
  #     for(j in (i+1):(n_states-1)) {
  #       log_probs[j] <- log_transitions[trans_idx,l]
  #       trans_idx <- trans_idx + 1
  #     }
  #
  #     # Add baseline category (dead state) as 0 on log scale
  #     log_probs[n_states] <- 0
  #
  #     # Calculate probabilities using log-sum-exp trick
  #     max_log_prob <- max(log_probs)
  #     log_sum <- max_log_prob + log(sum(exp(log_probs - max_log_prob)))
  #
  #     # Fill transition probabilities
  #     T[[l]][i,] <- exp(log_probs - log_sum)
  #   }
  #
  #   # Dead state stays dead
  #   T[[l]][n_states,n_states] <- 1
  # }
  #
  # # Create detection matrices for each location
  # D <- vector("list", n_locs)
  # for(l in 1:n_locs) {
  #   D[[l]] <- matrix(0, n_states, n_states)
  #
  #   # Fill detection probabilities
  #   for(i in 1:(n_states-1)) {
  #     D[[l]][i,i] <- detection[i,l]                    # Detected in correct state
  #     D[[l]][i,n_states] <- 1 - detection[i,l]         # Not detected
  #   }
  #
  #   # Dead state is always "not detected"
  #   D[[l]][n_states,n_states] <- 1
  # }
  #
  # # Forward algorithm with scaling for numerical stability
  # for(i in 1:n_fish) {
  #   alpha <- numeric(n_states)
  #   alpha[U[i] + 1] <- 1  # Initialize according to initial state vector
  #
  #   # Initial observation
  #   alpha <- alpha * D[[1]][,obs[i,1]]
  #   log_scale_total <- log(sum(alpha))
  #   alpha <- alpha / sum(alpha)
  #
  #   # Forward algorithm through remaining locations
  #   for(t in 1:(n_locs-1)) {
  #     # State transition
  #     alpha <- as.vector(T[[t]] %*% alpha)
  #
  #     # Observation update
  #     alpha <- alpha * D[[t+1]][,obs[i,t+1]]
  #
  #     # Scale to prevent underflow
  #     scale <- sum(alpha)
  #     log_scale_total <- log_scale_total + log(scale)
  #     alpha <- alpha / scale
  #   }
  #
  #   # Add the total scaling factor to the negative log-likelihood
  #   nll <- nll - log_scale_total
  # }
  #
  # return(nll)
}


# #Obj
# obj <- RTMB::MakeADFun(f,
#                        parameters)
#
# # Optimize
# opt <- nlminb(obj$par, obj$fn, obj$gr)

