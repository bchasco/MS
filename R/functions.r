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

  # Print unique values for debugging
  print(unique_values)

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

  print(state_array)
  return(state_array)
}

test_f <- function(parms){

  RTMB::getAll(tmb.data,
               parms)

  nll <- 0
  nll2 <- 0

  gamma <- array(0,dim = c(ns,ns, length(levels(data$loc))), #),#,ni,
                 dimnames = list(next_state = levels(state),
                                 current_state = levels(state),
                                 j = levels(data$loc)))
  omega <- array(0,dim = c(ns,ns,length(levels(data$loc))), #),#,ni
                 dimnames = list(next_state = levels(state),
                                 current_state = levels(state),
                                 j = levels(data$loc)))

  #Fill the survival list
  phi_s <- list()
  icnt <- 0
  for(i in seq_along(phi_dim)){
    if(length(phi_dim[[i]])>1){
      phi_s[[i]] <- phi[(icnt+1):phi_dim[[i]]] #phi is a vector
      icnt <- sum(phi_dim[[i]])
    }else{
      phi_s[[i]] <- phi #phi is a vector
    }
  }

  #Fill the survival list
  p_s <- list()
  icnt <- 0
  for(i in seq_along(p_dim)){
    if(length(p_dim[[i]])>1){
      p_s[[i]] <- p[(icnt+1):p_dim[[i]]] #phi is a vector
      icnt <- sum(p_dim[[i]])
    }else{
      p_s[[i]] <- p #phi is a vector
    }
  }

  #Fill the survival list
  lam_s <- list()
  icnt <- 0
  for(i in seq_along(lam_dim)){
    if(length(lam_dim[[i]])>1){
      lam_s[[i]] <- lam[(icnt+1):lam_dim[[i]]] #phi is a vector
      icnt <- sum(lam_dim[[i]])
    }else{
      lam_s[[i]] <- lam #phi is a vector
    }
  }

  ii <- 1

  phi_tmp <- c(RTMB::plogis(phi),RTMB::plogis(lam))
  p_tmp <- c(RTMB::plogis(p),RTMB::plogis(lam))

  id_cnt <- 1
  for(i in unique(data$id)){

    #fill the gamma matrix,
    #Observation for the ith "fish"
    i_obs <- as.integer(data$state[data$id==i])
    # i_obs[is.na(i_obs)] <- 2

    #sample size
    n_i <- max(data$n[data$id == i])

    #Initial location of the ith fish
    init_loc <- min(which(i_obs < ns,arr.ind=TRUE))

    #initial state
    init_state <- i_obs[1]
    last_loc <- max(which(i_obs == 1,arr.ind=TRUE))

    if(MR_settings$mod=="CJS"){
      chi_i <- 1
      if(last_loc<length(i_obs)){
        for(j in (length(i_obs)-1):last_loc){
          chi_i <- (1 - phi_tmp[j]) + phi_tmp[j]*(1 - p_tmp[j]) * chi_i
        }
      }
      if(last_loc>1){
        for(j in (init_loc+1):last_loc){
          y_i <- ifelse(i_obs[j]==1, 1., 0.)
          nll2 <- nll2 - sum(RTMB::dbinom(1., 1., phi_tmp[j-1], log=TRUE)) * n_i
          nll2 <- nll2 - sum(RTMB::dbinom(y_i, 1, p_tmp[j-1], log=TRUE)) * n_i
        }
      }
      nll2 <- nll2 - sum(RTMB::dbinom(1, 1, chi_i, log=TRUE)) * n_i
    }

    if(MR_settings$mod=="MS"){

      delta <- rep(0,ns)
      delta[1] <- 1


      prods <- list(RTMB::diag(1, ns))
      for(j in (init_loc+1):length(i_obs)){
        if(j < length(i_obs)){
          for(s in 1:(ns-1)){
            gamma[s,s,j] <- t(dm$phi[[s]][ii,]) %*% phi_s[[s]] #phi[s,j-1]
            omega[s,s,j] <- t(dm$p[[s]][ii,]) %*% p_s[[s]]
            if(ns>2){
              gamma[1,2] <- dm$eta[[1]][ii,] %*% eta
            }
          }
          ii <- ii + 1
        }
        # print(dm$lam[[s]][id_cnt,])
        if(j == length(i_obs)){
          for(s in 1:(ns-1)){
            gamma[s,s,j] <- t(dm$lam[[s]][id_cnt,]) %*% lam_s[[s]]
            omega[s,s,j] <- t(dm$lam[[s]][id_cnt,]) %*% lam_s[[s]]
            if(ns>2){
              gamma[1,2,j] <- data$eta[[1]][ii,] %*% eta
            }
          }
        }
        gamma[1:ns,ns,j] <- 1
        omega[1:ns,ns,j] <- 1
        for(s in 1:(ns-1)){
          gamma[s,s:ns,j] <- exp((gamma[s,s:ns,j]))/sum(exp((gamma[s,s:ns,j])))
          omega[s,s:ns,j] <- exp(omega[s,s:ns,j])/sum(exp(omega[s,s:ns,j]))
        }
        prods[[j - init_loc + 1]] <-  prods[[j - init_loc]] %*% gamma[,,j] %*% RTMB::diag(omega[,i_obs[j],j])
        # print(gamma[,])

      }
      pp <- (sum(t(delta) %*% prods[[length(prods)]]))
      nll2 <- nll2 - RTMB::dbinom(1,1,pp, log = TRUE) * n_i
      id_cnt <- id_cnt + 1
    }
  }

    rel_site <- unique(data$ReleaseSite)
    #Predicted survival
    gamma_pred <- array(0,dim = c(ns,ns, length(levels(data$loc)), length(rel_site)),
                        dimnames = list(next_state = levels(state),
                                        current_state = levels(state),
                                        Locations = levels(data$loc),
                                        rel_site = rel_site
                        ))
    cum_gamma_pred <- array(1,dim = c(ns,ns, length(levels(data$loc)), length(rel_site)),
                        dimnames = list(next_state = levels(state),
                                        current_state = levels(state),
                                        Locations = levels(data$loc),
                                        rel_site = rel_site
                        ))

    icnt <- 1
    pred_x <- rep(0,nrow(MR_settings$dm$pred_m))
    cum_pred_x <- rep(1,nrow(MR_settings$dm$pred_m))
    for(i in 1:length(pred_x)){
      cur_s <- as.integer(MR_settings$pred_data$current_state[i])
      next_s <- as.integer(MR_settings$pred_data$next_state[i])
      if(next_s<cur_s){
        pred_x[i] <- 0
      }
      if(cur_s<=next_s){
        pred_x[i] <- MR_settings$dm$pred_m[i,] %*% phi_s[[s]]
      }
      if(next_s==ns){
        pred_x[i] <- 1
      }
    }

    for(site in levels(MR_settings$pred_data$ReleaseSite)){
      lcnt <- 1
      for(loc in levels(MR_settings$pred_data$loc)){
        for(cur_s in levels(MR_settings$pred_data$current_state)){

          next_states <- pred_x[MR_settings$pred_data$ReleaseSite == site &
                                  MR_settings$pred_data$loc == loc &
                                  MR_settings$pred_data$current_state == cur_s &
                                  as.integer(MR_settings$pred_data$current_state) <= as.integer(MR_settings$pred_data$next_state)]


          pred_x[MR_settings$pred_data$ReleaseSite == site &
                   MR_settings$pred_data$loc == loc &
                   MR_settings$pred_data$current_state == cur_s &
                   as.integer(MR_settings$pred_data$current_state) <= as.integer(MR_settings$pred_data$next_state)] <- exp(next_states)/sum(exp(next_states))

          # tmp_cum <- pred_x[MR_settings$pred_data$ReleaseSite == site &
          #                                 MR_settings$pred_data$loc == loc &
          #                                 MR_settings$pred_data$current_state == cur_s &
          #                                 MR_settings$pred_data$next_state == levels(MR_settings$pred_data$next_state)[1]]

          if(loc ==levels(MR_settings$pred_data$loc)[1]){
            cum_pred_x[MR_settings$pred_data$ReleaseSite == site &
                       MR_settings$pred_data$loc == loc &
                       MR_settings$pred_data$current_state == cur_s &
                       MR_settings$pred_data$next_state == cur_s] <- pred_x[MR_settings$pred_data$ReleaseSite == site &
                                                                          MR_settings$pred_data$loc == loc &
                                                                          MR_settings$pred_data$current_state == cur_s &
                                                                          MR_settings$pred_data$next_state == cur_s]
          }
        }

        if(loc !=levels(MR_settings$pred_data$loc)[1]){
          for(s in levels(MR_settings$pred_data$next_state)[1]){
            cum_pred_x[MR_settings$pred_data$ReleaseSite == site &
                           MR_settings$pred_data$loc == loc &
                           MR_settings$pred_data$current_state == s &
                           MR_settings$pred_data$next_state == s] <- pred_x[MR_settings$pred_data$ReleaseSite == site &
                                                                                          MR_settings$pred_data$loc == loc &
                                                                                          MR_settings$pred_data$current_state == s &
                                                                                          MR_settings$pred_data$next_state == s] *
              cum_pred_x[MR_settings$pred_data$ReleaseSite == site &
                       MR_settings$pred_data$loc == levels(MR_settings$pred_data$loc)[lcnt - 1] &
                       MR_settings$pred_data$current_state == s &
                       MR_settings$pred_data$next_state == s]



          }
        }
        lcnt <- lcnt + 1
      }
    }

    cum_pred_x <- RTMB::qlogis(cum_pred_x)

  #   for(i in 1:length(pred_x)){
  #     cur_s <- MR_settings$pred_data$current_state[i]
  #     next_s <- MR_settings$pred_data$next_state[i]
  #     loc <- MR_settings$pred_data$loc[i]
  #     init_loc <- MR_settings$pred_data$tag_site[i]
  #   }
  #
  #
    nj <- length(i_obs)
    icnt <- 1
    for(i in 1:nrow(MR_settings$dm$pred_m)){
      site <- as.integer(MR_settings$pred_data$ReleaseSite[i])
      loc <- as.integer(MR_settings$pred_data$loc[i])
      cur_s <- as.integer(MR_settings$pred_data$current_state[i])
      next_s <- as.integer(MR_settings$pred_data$next_state[i])
      # for(j in 2:(nj-1)){
      #   for(s in 1:(ns-1)){
      if(next_s<cur_s){
        gamma_pred[cur_s,next_s,loc,site] <- MR_settings$dm$pred_m[i,] %*% phi_s[[s]]
        # icnt <- icnt + 1
      }
      if(next_s==ns){
        gamma_pred[cur_s,next_s,loc,site] <- 1
        # icnt <- icnt + 1
      }
        # for(s in 1:(ns-1)){
        #   gamma_pred[next_s,cur_s:ns,j,site] <- exp((gamma_pred[s,s:ns,j,site]))/sum(exp((gamma_pred[s,s:ns,j,site])))
        #   cum_gamma_pred[next_s,cur_s:ns,j,site] <- cum_gamma_pred[s,s:ns,j-1,site] * gamma_pred[s,s:ns,j,site]
        # }
      }
      #Last location
      # gamma_pred[s,s,nj,site] <- lam_s[[s]][site]
      # gamma_pred[1:ns,ns,nj,site] <- 1
      # for(s in 1:(ns-1)){
      #   gamma_pred[s,s:ns,nj,site] <- exp((gamma_pred[s,s:ns,nj,site]))/sum(exp((gamma_pred[s,s:ns,nj,site])))
      # }
    # }
  # }

  # pred_est <-  MR_settings$dm$pred_m %*% phi_s[[s]]
  # for(i in unique(MR_settings$pred_data$id)){
  #   for(j in unique(MR_settings$pred_data$current_state)){
  #     est <- rep(0,ns)
  #     cum_est <- rep(0,ns)
  #     jcnt <- 1
  #     for(j2 in unique(MR_settings$pred_data$next_state)){
  #       if(j2 == "unk"){
  #         est[jcnt] <- 1
  #       }else{
  #         est[jcnt] <- pred_est[MR_settings$pred_data$id == i & MR_settings$pred_data$current_state == j & MR_settings$pred_data$next_state == j2, ]
  #       }
  #       jcnt <- jcnt + 1
  #     }
  #     pred_est[MR_settings$pred_data$id == i & MR_settings$pred_data$current_state == j] <- exp(est)/sum(exp(est))
  #   }
  # }
  # for(i in unique(MR_settings$pred_data$id)){
  #       for(j in unique(MR_settings$pred_data$current_state)){
  #         pred_est[MR_settings$pred_data$id == i &
  #                    MR_settings$pred_data$loc == unique(MR_settings$pred_data$loc)[1] &
  #                    MR_settings$pred_data$current_state == unique(MR_settings$pred_data$current_state)[j] &
  #                    MR_settings$pred_data$current_state == unique(MR_settings$pred_data$current_state)[j]] <-
  #           for(l in 2:length(unique(MR_settings$pred_data$loc))){
  #             # pred_est[MR_settings$pred_data$id == i &
  #       #            MR_settings$pred_data$loc == l &
  #       #            MR_settings$pred_data$current_state == unique(MR_settings$pred_data$current_state)[j] &
  #       #            MR_settings$pred_data$current_state == unique(MR_settings$pred_data$current_state)[j]] <-
  #       #   pred_est[MR_settings$pred_data$id == i &
  #       #              MR_settings$pred_data$loc == l &
  #       #              MR_settings$pred_data$current_state == unique(MR_settings$pred_data$current_state)[j-1]] *
  #       #   pred_est[MR_settings$pred_data$id == i &
  #       #              MR_settings$pred_data$loc == l &
  #       #              MR_settings$pred_data$current_state == unique(MR_settings$pred_data$current_state)[j-1]] *
  #       }
  #   }
  # }
  gamma_pred <- RTMB::qlogis(gamma_pred)
  cum_gamma_pred <- RTMB::qlogis(cum_gamma_pred)

  RTMB::REPORT(phi)
  RTMB::REPORT(phi_s)
  RTMB::REPORT(p)
  RTMB::REPORT(eta)
  RTMB::REPORT(lam)
  RTMB::REPORT(gamma)
  RTMB::REPORT(pred_x)
  RTMB::REPORT(gamma_pred)
  RTMB::REPORT(cum_gamma_pred)
  # RTMB::REPORT(pred_est)
  RTMB::REPORT(cum_pred_x)
  RTMB::REPORT(omega)
  RTMB::REPORT(nll2)
  # RTMB::ADREPORT(pred_est)
  RTMB::ADREPORT(cum_pred_x)
  return(nll2)
}


create_design_matrices <- function(settings, data){
  dm <- list()
  phi_m <- list()
  pred_m <- list()
  eta_m <- list()
  lam_m <- list()
  p_m <- list()

  for(i in c("yr")){
    tmp_WEN <- data %>%
      mutate(state = ifelse(is.na(stage),"unk","yr")) %>%
      mutate(state = factor(state, levels = c("yr","unk"))) %>%
      mutate(loc = factor(loc, levels = c("Trib","First_Trap","RIS_RIA","WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
      droplevels() %>%
      # mutate(state = !!i) %>%
      mutate(tag_site = "Trib",
             last_site = last_site) %>%
      filter(loc != tag_site,
             loc != last_site) %>%
      droplevels()


    #For each state
    phi_m[[i]] <-  unlist(tryCatch(Matrix::sparse.model.matrix(formula(MR_settings$frm$phi), tmp_WEN),
                                   error = function(e) e,
                                   warning = function(w) w))

    names_phi <- colnames(phi_m[[i]])

    pred_data <- MR_settings$pred_data %>%
      mutate(tag_site = "Trib", last_site = locs[length(locs)]) %>%
      filter(loc != tag_site & loc != last_site)

    pred_m[[i]] <-  unlist(tryCatch(Matrix::sparse.model.matrix(formula(MR_settings$frm$phi), pred_data),
                                    error = function(e) e,
                                    warning = function(w) w))
    names_pred <- colnames(pred_m[[i]])

    pred_m <- pred_m[[i]][,names_pred %in% names_phi]

  }
  MR_settings$pred_data <- pred_data
  dm$phi <- phi_m
  dm$pred_m <- pred_m

  for(i in c("yr")){
    tmp_WEN <- data %>%
      mutate(state = ifelse(is.na(stage),"unk","yr")) %>%
      mutate(state = factor(state, levels = c("yr","unk"))) %>%
      mutate(loc = factor(loc, levels = c("Trib","First_Trap","RIS_RIA","WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
      # mutate(loc_phi = factor(loc, levels = c("Trib","First_Trap","WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
      droplevels() %>%
      # mutate(state = !!i) %>%
      mutate(tag_site = "Trib",
             last_site = last_site) %>%
      filter(loc != tag_site,
             loc != last_site) %>%
      droplevels()


    #For each state
    p_m[[i]] <-  unlist(tryCatch(Matrix::sparse.model.matrix(formula(MR_settings$frm$p), tmp_WEN),
                                 error = function(e) e,
                                 warning = function(w) w))


  }
  dm$p <- p_m

  for(i in c("yr")){
    tmp_WEN <- data %>%
      mutate(state = ifelse(is.na(stage),"unk","yr")) %>%
      mutate(state = factor(state, levels = c("yr","unk"))) %>%
      mutate(loc = factor(loc, levels = c("Trib","First_Trap","RIS_RIA","WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
      # mutate(loc_phi = factor(loc, levels = c("Trib","First_Trap","WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
      droplevels() %>%
      # mutate(state = !!i) %>%
      mutate(tag_site = "Trib",
             last_site = last_site) %>%
      filter(loc == last_site) %>%
      droplevels()


    #For each state
    lam_m[[i]] <-  unlist(tryCatch(Matrix::sparse.model.matrix(formula(MR_settings$frm$lam), tmp_WEN),
                                   error = function(e) e,
                                   warning = function(w) w))


  }
  dm$lam <- lam_m

  for(i in c("sub")){
    tmp_WEN <- data %>%
      mutate(state = ifelse(is.na(stage), "unk", stage)) %>%
      mutate(state = factor(state, levels = c("yr","unk"))) %>%
      mutate(loc = factor(loc, levels = c("Trib","First_Trap","RIS_RIA", "WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
      mutate(loc_eta = case_when(
        i == "sub" & loc %in% c("Wen","MCJ","JDJ", "BON") ~ "Wen",
        i == "sub" ~ as.character(loc)
      )) %>%
      mutate(loc_eta = factor(loc_eta, levels = c("Trib","First_Trap","WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
      droplevels() %>%
      mutate(state = !!i) %>%
      mutate(tag_site = "Trib",
             last_site = locs[length(locs)])

    tmp_WEN <- tmp_WEN[tmp_WEN$loc!=tmp_WEN$tag_site &
                         tmp_WEN$loc!=tmp_WEN$last_site,] %>%
      droplevels()

    eta_m[[i]] <-  tryCatch(Matrix::sparse.model.matrix(formula('state ~ loc_eta'), tmp_WEN),
                            error = function(e) e,
                            warning = function(w) w)
  }
 dm$eta <- eta_m

 return(dm)
}
