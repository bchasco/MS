#' Parse a Formula to Extract Covariate Terms
#'
#' This function extracts the terms (variables) from a formula object. It is used to
#' identify which covariates are involved in creating a state array for a hidden Markov model.
#'
#' @param formula A formula object representing the covariates.
#'
#' @return A character vector of the terms (covariates) in the formula.
#'
#' @examples
#' # Example formula
#' frm <- formula('~ LH + loc')
#'
#' # Extract terms
#' ftrms <- parseFormula(frm)
#'
#' #print terms of formula
#' print(ftrms)
#'
#' @export
parseFormula <- function(formula) {
  terms <- all.vars(formula)  # Extract the variables from the formula
  return(terms)
}

#' Extract Unique Values for Covariates from a Data Frame
#'
#' This function takes a vector of covariate names and a data frame, and returns
#' the unique values for each covariate. This is useful for determining the dimensions
#' of a state array for a hidden Markov model.
#'
#' @param terms A character vector of covariate names.
#' @param covars A data frame containing covariates with names matching the terms.
#'
#' @return A named list where each element contains the unique values of the corresponding covariate.
#'
#' @examples
#' # Example data
#' covars <- data.frame(LH = rep(c("Sub", "Yr", "Dead"),2),
#'                      loc = rep(c("Wenatchee", "McNary", "John Day"),2))
#'
#' # Extract unique values for LH and loc
#' unique_values <- extractUniqueValues(c("LH", "loc"), covars)
#'
#' #Print the unique values
#' print(unique_values)
#' @export
extractUniqueValues <- function(terms, covars) {
  unique_values <- lapply(terms, function(var) unique(covars[[var]]))
  names(unique_values) <- terms
  return(unique_values)
}


#' Create a State Array from a Formula and Covariate Data
#'
#' This function takes a formula and a data frame of covariates to generate
#' an array of possible states for a hidden Markov model. The dimensions of
#' the array correspond to the unique values of the covariates in the formula.
#'
#' @param formula A formula object representing the covariates to be used in the state array.
#' @param covars A data frame containing the covariates. The column names of the data frame
#' should match the terms in the formula.
#'
#' @return An array where the dimensions correspond to the unique values of the covariates
#' specified in the formula. If there is one covariate term, the function returns a 2D array (NxN),
#' and if there are multiple covariates, it returns a multi-dimensional array.
#'
#' @examples
#' # Example data
#' state_vars <- data.frame(LH = rep(c("sub", "yr", "unk"),2),
#'                      loc = rep(c("A", "B", "C"),2))
#'
#' #You have to explicitly define the levels for the states
#' state_vars$LH <- factor(state_vars$LH, levels = c("sub", "yr", "unk"))
#' state_vars$loc <- factor(state_vars$loc, levels = c("C", "A", "B"))
#'
#' # Example formula
#' frm <- formula('~ LH + loc')
#'
#' # Create the state array
#' state_array <- createStateArray(frm, state_vars)
#'
#' @export
createStateArray <- function(formula, covars) {
  # Check that formula and covars are valid
  if (!inherits(formula, "formula")) stop("`formula` must be a formula object.")
  if (!is.data.frame(covars)) stop("`covars` must be a data frame.")

  # Extract terms from the formula and ensure all are in covars
  terms <- all.vars(formula)
  missing_terms <- setdiff(terms, names(covars))
  if (length(missing_terms) > 0) stop("Missing terms in covars: ", paste(missing_terms, collapse = ", "))

  # Extract and sort unique values for each variable in the specified order
  unique_values <- lapply(terms, function(term) sort(unique(covars[[term]])))
  names(unique_values) <- terms


  # Define array dimensions based on unique values of each covariate
  dims <- sapply(unique_values, length)
  array_dim <- if (length(dims) == 1) c(dims, dims) else c(dims[1], dims[1], dims[-1])

  # Initialize the array and set the dimension names
  state_array <- array(NA, dim = array_dim)

  # Assign dimension names based on sorted unique values, preserving order of terms
  if (length(dims) == 1) {
    dimnames(state_array) <- list(unique_values[[1]], unique_values[[1]])
  } else {
    # Set dimension names for each level based on unique_values
    dimnames_list <- c(list(unique_values[[1]], unique_values[[1]]), unique_values[-1])
    dimnames(state_array) <- dimnames_list
  }

  # print(state_array)
  return(state_array)
}

test_f <- function(parms){

  RTMB::getAll(tmb.data,
               parms)

  nll <- 0
  nll2 <- 0

  ni <- length(unique(data$id))
  # pred_ni <- length(unique(data$id))

  gamma <- array(0,dim = c(ns,ns, length(levels(data$loc)), ni), #),#,ni,
                 dimnames = list(next_state = levels(state),
                                 current_state = levels(state),
                                 j = levels(data$loc),
                                 i = 1:ni))
  omega <- array(0,dim = c(ns,ns,length(levels(data$loc)), ni), #),#,ni
                 dimnames = list(next_state = levels(state),
                                 current_state = levels(state),
                                 j = levels(data$loc),
                                 i = 1:ni))



  #Fill the survival list
  phi_s <- list()
  icnt <- 0
  if(length(phi_dim)>1){
    for(i in seq_along(phi_dim)){
      phi_s[[i]] <- phi[(icnt+1):(icnt+phi_dim[[i]])] #phi is a vector
      icnt <- icnt + sum(phi_dim[[i]])
    }
  }else{
    for(i in 1:(ns-1)){
      phi_s[[i]] <- phi #phi is a vector
    }
  }

  #Fill the survival list
  phi_re_s <- list()
  icnt <- 0
  if(length(phi_re_dim)>1){
    for(i in seq_along(phi_re_dim)){
      if(phi_re_dim[[i]]>0){
        phi_re_s[[i]] <- phi_re[(icnt+1):(icnt+phi_re_dim[[i]])] #phi is a vector
      }else{
        phi_re_s[[i]] <- NULL
      }
      icnt <- icnt + sum(phi_re_dim[[i]])
    }
  }else{
    for(i in 1:(ns-1)){
      phi_re_s[[i]] <- phi_re #phi is a vector
    }
  }

  #Fill the detection list
  p_s <- list()
  icnt <- 0
  if(length(p_dim)>1){
    for(i in seq_along(p_dim)){
      p_s[[i]] <- p[(icnt+1):(icnt + p_dim[[i]])] #phi is a vector
      icnt <- icnt + sum(p_dim[[i]])
    }
  }else{
    for(i in 1:(ns-1)){
      p_s[[i]] <- p #p is a vector
    }
  }

  #Fill the survival list
  p_re_s <- list()
  icnt <- 0
  if(length(p_re_dim)>1){
    for(i in seq_along(p_re_dim)){
      if(p_re_dim[[i]]>0){
        p_re_s[[i]] <- p_re[(icnt+1):(icnt+p_re_dim[[i]])] #phi is a vector
      }else{
        p_re_s[[i]] <- NULL
      }
      icnt <- icnt + sum(p_re_dim[[i]])
    }
  }else{
    for(i in 1:(ns-1)){
      p_re_s[[i]] <- p_re #phi is a vector
    }
  }

  #Fill the last survival location list
  lam_s <- list()
  icnt <- 0
  if(length(lam_dim)>1){
    for(i in seq_along(lam_dim)){
      lam_s[[i]] <- lam[(icnt+1):(icnt + lam_dim[[i]])] #phi is a vector
      icnt <- icnt + sum(lam_dim[[i]])
    }
  }else{
    for(i in 1:(ns-1)){
      lam_s[[i]] <- lam #lam is a vector
    }
  }

  #Fill the survival list
  lam_re_s <- list()
  icnt <- 0
  if(length(lam_re_dim)>1){
    for(i in seq_along(lam_re_dim)){
      if(lam_re_dim[[i]]>0){
        lam_re_s[[i]] <- lam_re[(icnt+1):(icnt+lam_re_dim[[i]])] #phi is a vector
      }else{
        lam_re_s[[i]] <- NULL
      }
      icnt <- icnt + sum(lam_re_dim[[i]])
    }
  }else{
    for(i in 1:(ns-1)){
      lam_re_s[[i]] <- lam_re #phi is a vector
    }
  }

  phi_tmp <- c(RTMB::plogis(phi),RTMB::plogis(lam))
  p_tmp <- c(RTMB::plogis(p),RTMB::plogis(lam))

  id_cnt <- 1 #fish group increment
  ii <- 1 #design matrix increment
  eta_ii <- 1 #transition matrix increment

  for(i in unique(data$id)){
    #fill the gamma matrix,
    #Observation for the ith "fish"
    i_obs <- as.integer(data$state[data$id==i])
    t_obs <- data$time[data$id==i] * 1.

    #sample size
    n_i <- max(data$n[data$id == i])

    #Initial location of the ith fish
    init_loc <- min(which(i_obs < ns,arr.ind=TRUE))

    #initial state
    init_state <- i_obs[init_loc]

    last_loc <- max(which(i_obs < ns,arr.ind=TRUE))


    if(MR_settings$mod=="CJS"){
      chi_i <- 1
      if(last_loc<length(i_obs)){
        for(j in (length(i_obs)-1):last_loc){
          chi_i <- (1 - phi_tmp[j]) + phi_tmp[j]*(1 - p_tmp[j]) * chi_i
        }
      }
      if(last_loc>1){
        for(j in (init_loc+1):last_loc){
          y_i <- ifelse(i_obs[j]<ns, 1., 0.)
          nll2 <- nll2 - sum(RTMB::dbinom(1., 1., phi_tmp[j-1], log=TRUE)) * n_i
          nll2 <- nll2 - sum(RTMB::dbinom(y_i, 1, p_tmp[j-1], log=TRUE)) * n_i
        }
      }
      nll2 <- nll2 - sum(RTMB::dbinom(1, 1, chi_i, log=TRUE)) * n_i
    }

    if(MR_settings$mod=="MS" | MR_settings$mod=="MSt"){

      delta <- rep(0,ns)
      delta[init_state] <- 1

      prods <- list(RTMB::diag(1, ns))

      for(j in (init_loc+1):length(i_obs)){
        # if(MR_settings$mod=="MS"){
        #   if(j < length(i_obs)){
        #     for(s in 1:(ns-1)){
        #       if(phi_dim[[s]]>0){
        #         gamma[s,s,j,id_cnt] <- t(dm$phi[[s]]$X[ii,]) %*% phi_s[[s]]
        #         if(phi_re_dim[[s]]>0){
        #           gamma[s,s,j,id_cnt] <- gamma[s,s,j,id_cnt] + t(dm$phi[[s]]$reTrms$Zt[,ii]) %*% phi_re_s[[s]]
        #         }
        #       }else{
        #         gamma[s,s,j,id_cnt] <- -Inf
        #       }
        #
        #       omega[s,s,j,id_cnt] <- t(dm$p[[s]]$X[ii,]) %*% p_s[[s]]
        #       if(p_re_dim[[s]]>0){
        #         omega[s,s,j,id_cnt] <- omega[s,s,j,id_cnt] + t(dm$p[[s]]$reTrms$Zt[,ii]) %*% p_re_s[[s]]
        #       }
        #
        #       if(ns>2 & s==1){
        #         gamma[1,2,j,id_cnt] <- dm$eta[[1]]$X[eta_ii,] %*% eta
        #       }
        #     }
        #
        #     ii <- ii + 1 #location incrementor. The design matrix is referenced by location ii
        #     eta_ii <- eta_ii + 1
        #
        #   }
        #
        #   if(j == length(i_obs)){
        #     for(s in 1:(ns-1)){
        #       gamma[s,s,j,id_cnt] <- t(dm$lam[[s]]$X[id_cnt,]) %*% lam_s[[s]]
        #       omega[s,s,j,id_cnt] <- t(dm$lam[[s]]$X[id_cnt,]) %*% lam_s[[s]]
        #       if(ns>2 & s==1){
        #         gamma[1,2,j,id_cnt] <- dm$eta[[1]]$X[eta_ii,] %*% eta
        #       }
        #     }
        #     eta_ii <- eta_ii + 1
        #   }
        #
        #   gamma[,ns,j,id_cnt] <- 1
        #   # gamma_tmp[,,j,id_cnt] <- gamma[,,j,id_cnt]
        #   omega[,ns,j,id_cnt] <- 1
        #
        #   for(s in 1:(ns-1)){
        #     gamma[s,s:ns,j,id_cnt] <- exp((gamma[s,s:ns,j,id_cnt]))/sum(exp(gamma[s,s:ns,j,id_cnt]))
        #     omega[s,c(s,ns),j,id_cnt] <- exp(omega[s,c(s,ns),j,id_cnt])/sum(exp(omega[s,c(s,ns),j,id_cnt]))
        #   }
        #   prods[[j - init_loc + 1]] <-  prods[[j - init_loc]] %*% gamma[,,j,id_cnt] %*% RTMB::diag(omega[,i_obs[j],j,id_cnt])
        # }


        # if(MR_settings$mod == "MSt"){
          if(j < length(i_obs)){
            for(s in 1:(ns-1)){
              if(phi_dim[[s]]>0){
                gamma[s,s,j,id_cnt] <- t(dm$phi[[s]]$X[ii,]) %*% phi_s[[s]]
                if(phi_re_dim[[s]]>0){
                  gamma[s,s,j,id_cnt] <- -1. * exp(gamma[s,s,j,id_cnt]) * exp(t(dm$phi[[s]]$reTrms$Zt[,ii]) %*% phi_re_s[[s]])
                }else{
                  gamma[s,s,j,id_cnt] <- -1. * exp(gamma[s,s,j,id_cnt])
                }

              }else{
                gamma[s,s,j,id_cnt] <- 0
              }

              omega[s,s,j,id_cnt] <- t(dm$p[[s]]$X[ii,]) %*% p_s[[s]]
              if(p_re_dim[[s]]>0){
                omega[s,s,j,id_cnt] <- -1. * exp(omega[s,s,j,id_cnt]) * exp(t(dm$p[[s]]$reTrms$Zt[,ii]) %*% p_re_s[[s]])
              }
              omega[s,s,j,id_cnt] <- -1. * exp(omega[s,s,j,id_cnt])

              if(ns>2 & s==1){
                gamma[1,2,j,id_cnt] <- - 1. * (dm$eta[[1]]$X[eta_ii,] %*% eta)
              }
            }

            if(!is.na(t_obs[j])){
              nll2 <- nll2 - RTMB::dnorm(t_obs[j], exp(tau[j-1]), exp(fsig_tau[j-1]), log = TRUE) * n_i
            }

            ii <- ii + 1 #location incrementor. The design matrix is referenced by location ii
            eta_ii <- eta_ii + 1
          }

          if(j == length(i_obs)){
            for(s in 1:(ns-1)){
              lam <- t(dm$lam[[s]]$X[id_cnt,]) %*% lam_s[[s]]
              if(lam_re_dim[[s]]>0){
                lam <- -1. * exp(lam) * exp(t(dm$lam[[s]]$reTrms$Zt[,id_cnt]) %*% lam_re_s[[s]])
              }else{
                lam <- -1. * exp(lam)
              }
              gamma[s,s,j,id_cnt] <- lam
              omega[s,s,j,id_cnt] <- lam

              if(ns>2 & s==1){
                gamma[1,2,j,id_cnt] <- dm$eta[[1]]$X[eta_ii,] %*% eta
              }
            }
            eta_ii <- eta_ii + 1
          }

          for(s in 1:(ns-1)){
            gamma[s,ns,j,id_cnt] <- -1 * sum(gamma[s,s:(ns-1),j,id_cnt])
            omega[s,ns,j,id_cnt] <- -1 * sum(omega[s,s:(ns-1),j,id_cnt])
          }

          if(MR_settings$mod == "MSt"){
            if(is.na(t_obs[j])){
              prods[[j - init_loc + 1]] <-  prods[[j - init_loc]] %*% Matrix::expm((gamma[,,j,id_cnt]) * exp(tau[j-1])) %*% RTMB::diag(Matrix::expm(omega[,,j,id_cnt])[,i_obs[j]])
            }else{
              prods[[j - init_loc + 1]] <-  prods[[j - init_loc]] %*% Matrix::expm((gamma[,,j,id_cnt]) * t_obs[j]) %*% RTMB::diag(Matrix::expm(omega[,,j,id_cnt])[,i_obs[j]])
            }
          }else{
            prods[[j - init_loc + 1]] <-  prods[[j - init_loc]] %*% Matrix::expm((gamma[,,j,id_cnt])) %*% RTMB::diag(Matrix::expm(omega[,,j,id_cnt])[,i_obs[j]])
          }
        # }

      }
      pp <- (sum(t(delta) %*% prods[[length(prods)]]))
      nll2 <- nll2 - RTMB::dbinom(1,1,pp, log = TRUE) * n_i
      id_cnt <- id_cnt + 1
    }
  }
  nll2 <- nll2 - sum(RTMB::dnorm(phi_re_s[[1]],
                                 0,
                                 exp(phi_re_sig[dm$phi[[1]]$reTrms$Lind]),
                                 log = TRUE))

  nll2 <- nll2 - sum(RTMB::dnorm(p_re_s[[1]],
                                 0,
                                 exp(p_re_sig[dm$p[[1]]$reTrms$Lind]),
                                 log = TRUE))

  nll2 <- nll2 - sum(RTMB::dnorm(lam_re_s[[1]],
                                 0,
                                 exp(lam_re_sig[dm$lam[[1]]$reTrms$Lind]),
                                 log = TRUE))

  pd <- tmb.data$pred.grid
  pdm <- tmb.data$pred$dm
  id_cnt <- 1 #fish group increment
  ii <- 1 #design matrix increment
  eta_ii <- 1 #transition matrix increment
  pni <- length(unique(pd$id))
  gamma_pred <- array(0,dim = c(ns,ns, length(levels(data$loc)), pni), #),#,ni,
                      dimnames = list(next_state = levels(state),
                                      current_state = levels(state),
                                      j = levels(data$loc),
                                      i = 1:pni))
  omega_pred <- array(0,dim = c(ns,ns, length(levels(data$loc)), pni), #),#,ni,
                      dimnames = list(next_state = levels(state),
                                      current_state = levels(state),
                                      j = levels(data$loc),
                                      i = 1:pni))
  for(i in unique(pd$id)){
    #fill the gamma matrix,
    #Observation for the ith "fish"
    i_obs <- as.integer(pd$state[pd$id==i])
    #Initial location of the ith fish
    init_loc <- 1#min(which(i_obs < ns,arr.ind=TRUE))

    #initial state
    init_state <- 2#i_obs[init_loc]

    last_loc <- 3#max(which(i_obs < ns,arr.ind=TRUE))

    for(j in (init_loc+1):length(i_obs)){
      if(j < length(i_obs)){
        for(s in 1:(ns-1)){
          if(phi_dim[[s]]>0){
            gamma_pred[s,s,j,id_cnt] <- exp(t(pdm$phi[[s]]$X[ii,]) %*% phi_s[[s]])
            if(phi_re_dim[[s]]>0){
              gamma_pred[s,s,j,id_cnt] <- gamma_pred[s,s,j,id_cnt] * exp(t(pdm$phi[[s]]$reTrms$Zt[,ii]) %*% phi_re_s[[s]])
            }
            gamma_pred[s,s,j,id_cnt] <- -1. * gamma_pred[s,s,j,id_cnt]
          }else{
            gamma_pred[s,s,j,id_cnt] <- 0
          }

          omega_pred[s,s,j,id_cnt] <- exp(t(pdm$p[[s]]$X[ii,]) %*% p_s[[s]])
          if(p_re_dim[[s]]>0){
            omega_pred[s,s,j,id_cnt] <- omega_pred[s,s,j,id_cnt] * exp(t(pdm$p[[s]]$reTrms$Zt[,ii]) %*% p_re_s[[s]])
          }
          omega_pred[s,s,j,id_cnt] <- -1. * exp(omega_pred[s,s,j,id_cnt])

          if(ns>2 & s==1){
            gamma_pred[1,2,j,id_cnt] <- - 1. * (pdm$eta[[1]]$X[eta_ii,] %*% eta)
          }
        }
        ii <- ii + 1 #location incrementor. The design matrix is referenced by location ii, but only from init + 1 : last - 1
        eta_ii <- eta_ii + 1
      }

      if(j == length(i_obs)){
        for(s in 1:(ns-1)){
          gamma_pred[s,s,j,id_cnt] <- -1. * exp(t(pdm$lam[[s]]$X[id_cnt,]) %*% lam_s[[s]])
          omega_pred[s,s,j,id_cnt] <- -1. * exp(t(pdm$lam[[s]]$X[id_cnt,]) %*% lam_s[[s]])
          if(ns>2 & s==1){
            gamma_pred[1,2,j,id_cnt] <- pdm$eta[[1]]$X[eta_ii,] %*% eta
          }
          gamma_pred[s,s,j,id_cnt] <- -1. * exp(gamma_pred[s,s,j,id_cnt])
        }
        eta_ii <- eta_ii + 1
      }

      for(s in 1:(ns-1)){
        gamma_pred[s,ns,j,id_cnt] <- -1 * sum(gamma_pred[s,s:(ns-1),j,id_cnt])
        omega_pred[s,ns,j,id_cnt] <- -1 * sum(omega_pred[s,s:(ns-1),j,id_cnt])
      }
      if(is.na(pd$time[ii-1])){
        tmp_time <- 1
      }else{
        tmp_time <- pd$time[ii-1]
      }
      # print(tmp_time)

    }
    id_cnt <- id_cnt + 1
  }

  RTMB::REPORT(ii)
  RTMB::REPORT(id_cnt)
  RTMB::REPORT(ii)

  RTMB::REPORT(phi)
  RTMB::REPORT(phi_s)
  RTMB::REPORT(p)
  RTMB::REPORT(p_s)
  RTMB::REPORT(phi_re_s)
  RTMB::REPORT(p_re_s)
  RTMB::REPORT(eta)
  RTMB::REPORT(lam)
  RTMB::REPORT(lam_s)
  # RTMB::REPORT(lam_re_s)
  RTMB::REPORT(gamma)
  RTMB::REPORT(omega)
  RTMB::REPORT(nll2)
  RTMB::REPORT(prods)
  if(MR_settings$mod=="CJS"){
    RTMB::REPORT(phi_tmp)
    RTMB::REPORT(p_tmp)
  }
  RTMB::REPORT(gamma_pred)
  RTMB::REPORT(omega_pred)
  # RTMB::ADREPORT(gamma_pred)
  # RTMB::ADREPORT(omega_pred)
  return(nll2)
}


create_design_matrices <- function(settings, data) {
  dm <- list()

  #What is that state column in the original data
  state <- settings$state #model state from MR_settings

  #Get rid of the last state
  states <- unique(data$state)
  states <- states[states != "unk"]

  # Helper function for generating design matrices
  generate_design_matrix <- function(formula_list, data_subset, model_type) {
    design_matrices <- list()
    pred_matrices <- list()
    reTrms <- list()

    for (i in seq_along(states)) { #Removed the last state ("unk")
      if(length(formula_list)==1){
        pi <- 1
      }else{
        pi <- i
      }
      current_formula <- formula_list[[pi]]
      tmp_data <- data_subset()
      tryCatch({
        if(min(attr(gregexpr("\\|",current_formula)[[1]],"match.length")>0)){
          lF <- lme4::lFormula(formula(paste(state,current_formula)), tmp_data)
          design_matrices[[i]] <- lF
        }else{
          design_matrices[[i]] <- list()
          design_matrices[[i]]$X <- model.matrix(formula(current_formula), tmp_data)
          design_matrices[[i]]$reTrms <- list()
          design_matrices[[i]]$fr <- tmp_data
        }

      }, error = function(e) {
        design_matrices[[i]] <- NULL
        # pred_matrices[[i]] <- NULL
      })
    }
    list(design = design_matrices)

  }

  # Generate design matrices for phi
  phi_data_fn <- function() {
    data %>%
      mutate(tag_site = init_site, last_site = last_site) %>%
      filter(loc != tag_site & loc != last_site) %>%
      droplevels()
  }
  phi_matrices <- generate_design_matrix(settings$frm$phi, phi_data_fn, "phi")
  dm$phi <- phi_matrices$design

  # Generate design matrices for p
  p_data_fn <- function() {
    data %>%
      mutate(tag_site = init_site, last_site = last_site) %>%
      filter(loc != tag_site & loc != last_site) %>%
      droplevels()
  }
  p_matrices <- generate_design_matrix(settings$frm$p, p_data_fn, "p")
  dm$p <- p_matrices$design

  # Generate design matrices for lam
  lam_data_fn <- function() {
    data %>%
      mutate(tag_site = init_site, last_site = last_site) %>%
      filter(loc == last_site) %>%
      droplevels()
  }
  lam_matrices <- generate_design_matrix(settings$frm$lam, lam_data_fn, "lam")
  dm$lam <- lam_matrices$design

  # Generate design matrices for eta
  eta_data_fn <- function() {
    data %>%
      mutate(
        loc = factor(loc, levels = locs),
        tag_site = init_site,
        last_site = locs[length(locs)]
      ) %>%
      filter(loc != last_site) %>%
      droplevels()
  }
  eta_matrices <- generate_design_matrix(settings$frm$eta, eta_data_fn, "eta")
  dm$eta <- eta_matrices$design  # Only one formula for eta


  # Generate design matrices for eta
  time_data_fn <- function() {
    data %>%
      mutate(
        loc = factor(loc, levels = locs),
        tag_site = init_site,
        last_site = locs[length(locs)]
      ) %>%
      filter(loc != tag_site) %>%
      droplevels()
  }
  time_matrices <- generate_design_matrix(settings$frm$time, time_data_fn, "time")
  dm$time <- time_matrices$design  # Only one formula for eta

  return(list(dm = dm, data = data))

}

create_pred_matrices <- function(object, basis, new_grid) {
  dm <- list()

  #What is that state column in the original data
  state <- object@MR_settings$state #model state from MR_settings

  #Get rid of the last state
  states <- unique(object@data$state)
  states <- states[states != "unk"]

  # Helper function for generating design matrices
  generate_design_matrix <- function(formula_list, data_subset, model_type, basis) {
    design_matrices <- list()
    reTrms <- list()

    for (i in seq_along(states)) { #Removed the last state ("unk")
      if(length(formula_list)==1){
        pi <- 1
      }else{
        pi <- i
      }
      current_formula <- formula_list[[pi]]
      tmp_data <- data_subset()
      # print(model_type)
      # print(current_formula)
      # print("basis")
      # print(basis$phi)
      # print('gregexpr')
      # print(attr(gregexpr("\\|",current_formula)[[1]],"match.length"))
      # print('useBytes')
      # print(attr(gregexpr("\\|",current_formula)[[1]],"useBytes"))

      tryCatch({
        if(max(attr(gregexpr("\\|",current_formula)[[1]],"match.length"))>0){
          print(paste(model_type,"lF"))
            lF <- lme4::lFormula(formula(paste(state,current_formula)), tmp_data)
            design_matrices[[i]] <- lF

          if(attr(gregexpr("poly", current_formula)[[1]],"match.length")>(-1)){
            poly_vars <- names(basis[[model_type]][[i]])
            for(p in poly_vars){
              design_matrices[[i]]$X[grep(p, colnames(design_matrices[[i]]$X), value = TRUE)] <- basis$phi[[p]]
            }
          }
        }else{
          design_matrices[[i]] <- list()
          #you need to get the original poly basis transformation for consistency
          if(attr(gregexpr("poly", current_formula)[[1]],"match.length")>(-1)){
            design_matrices[[i]]$X <- model.matrix(formula(current_formula), tmp_data)
            poly_vars <- names(basis[[model_type]][[i]])
            for(p in poly_vars){
              design_matrices[[i]]$X[grep(p, colnames(design_matrices[[i]]$X), value = TRUE)] <- basis$phi[[p]]
            }
            #replace poly terms
            design_matrices[[i]]$reTrms <- list()
          }else{
            # print("test")

            design_matrices[[i]]$X <- model.matrix(formula(current_formula), tmp_data)
            design_matrices[[i]]$reTrms <- list()
          }
        }

      }, error = function(e) {
        design_matrices[[i]] <- NULL
      })
    }
    list(design = design_matrices)

  }

  # Generate design matrices for phi
  phi_data_fn <- function() {
    new_grid %>%
      mutate(tag_site = init_site, last_site = last_site) %>%
      filter(loc != tag_site & loc != last_site) %>%
      droplevels()
  }

  phi_matrices <- generate_design_matrix(object@MR_settings$frms$phi, phi_data_fn, "phi", basis)
  dm$phi <- phi_matrices$design

  # Generate design matrices for p
  p_data_fn <- function() {
    new_grid %>%
      mutate(tag_site = init_site, last_site = last_site) %>%
      filter(loc != tag_site & loc != last_site) %>%
      droplevels()
  }
  p_matrices <- generate_design_matrix(object@MR_settings$frms$p, p_data_fn, "p", basis)
  dm$p <- p_matrices$design
  # print("test create p pred matrices")

  # Generate design matrices for lam
  lam_data_fn <- function() {
    new_grid %>%
      mutate(tag_site = init_site, last_site = last_site) %>%
      filter(loc == last_site) %>%
      droplevels()
  }
  lam_matrices <- generate_design_matrix(object@MR_settings$frms$lam, lam_data_fn, "lam", basis)
  dm$lam <- lam_matrices$design

  # Generate design matrices for eta
  eta_data_fn <- function() {
    new_grid %>%
      mutate(
        loc = factor(loc, levels = locs),
        tag_site = init_site,
        last_site = locs[length(locs)]
      ) %>%
      filter(loc != last_site) %>%
      droplevels()
  }
  eta_matrices <- generate_design_matrix(object@MR_settings$frm$eta, eta_data_fn, "eta", basis)
  dm$eta <- eta_matrices$design  # Only one formula for eta

  # Generate design matrices for eta
  time_data_fn <- function() {
    new_grid %>%
      mutate(
        loc = factor(loc, levels = locs),
        tag_site = init_site,
        last_site = locs[length(locs)]
      ) %>%
      filter(loc != tag_site) %>%
      droplevels()
  }
  time_matrices <- generate_design_matrix(object@MR_settings$frm$time, time_data_fn, "time", basis)
  dm$time <- time_matrices$design  # Only one formula for eta

  return(list(dm = dm))
}
