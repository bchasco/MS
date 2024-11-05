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

  gamma <- array(0,dim = c(ns,ns,ni, length(levels(data$loc))),
                 dimnames = list(next_state = levels(state),
                                 current_state = levels(state),
                                 id = unique(data$id),
                                 j = levels(data$loc)))
  omega <- array(0,dim = c(ns,ns,ni,length(levels(data$loc))),
                 dimnames = list(next_state = levels(state),
                                 current_state = levels(state),
                                 id = unique(data$id),
                                 j = levels(data$loc)))

  #Fill the survival list
  phi_s <- list()
  icnt <- 0
  for(i in seq_along(phi_dim)){
    # print((icnt+1):sum(phi_dim[1:i]))
    # print(phi)
    # print(phi[(icnt+1):sum(phi_dim[1:i])])

    phi_s[[i]] <- phi[(icnt+1):sum(phi_dim[1:i])] #phi is a vector
    icnt <- sum(phi_dim[1:i])
  }

  ii <- 1
  for(i in unique(data$id)){

    #fill the gamma matrix,
    #Observation for the ith "fish"
    print(i)
    i_obs <- as.integer(data$state[data$id==i])
    i_obs[is.na(i_obs)] <- 2

    #Initial location of the ith fish
    init_loc <- min(which(i_obs < ns,arr.ind=TRUE))
    #initial state
    init_state <- i_obs[init_loc]

    #sample size
    n_i <- max(data$n[data$id == i])

    delta <- rep(0,ns)
    delta[init_state] <- 1


    prods <- list(RTMB::diag(1, ns))
    print(i_obs)
    print((init_loc+1):length(i_obs))
    for(j in (init_loc+1):length(i_obs)){
      if(j < length(i_obs)){
        for(s in 1:(ns-1)){
          gamma[s,s,ii,j] <- dm$phi[[s]][ii,] %*% phi_s[[s]] #phi[s,j-1]
          omega[s,s,ii,j] <- p[s,j-1]
          if(ns>2){
            gamma[1,2,ii,j] <- dm$eta[[1]][ii,] %*% eta
          }
        }
      }
      if(j == length(i_obs)){
        for(s in 1:(ns-1)){
          gamma[s,s,ii,j] <- lam
          omega[s,s,ii,j] <- lam
          if(ns>2){
            gamma[1,2,ii,j] <- data$eta[[1]][ii,] %*% eta
          }
        }
      }
      gamma[1:ns,ns,ii,j] <- 1
      omega[1:ns,ns,ii,j] <- 1

      for(s in 1:(ns-1)){
        gamma[s,s:ns,ii,j] <- exp((gamma[s,s:ns,ii,j]))/sum(exp((gamma[s,s:ns,ii,j])))
        omega[s,s:ns,ii,j] <- exp(omega[s,s:ns,ii,j])/sum(exp(omega[s,s:ns,ii,j]))
      }
      prods[[j - init_loc + 1]] <-  prods[[j - init_loc]] %*% gamma[,,ii,j] %*% RTMB::diag(omega[,i_obs[j],ii,j])
    }
    pp <- (sum(t(delta) %*% prods[[length(prods)]]))
    nll2 <- nll2 - RTMB::dbinom(1,1,pp, log = TRUE) * n_i
    ii <- ii + 1
  }

  # RTMB::REPORT(gamma)
  RTMB::REPORT(phi)
  RTMB::REPORT(phi_s)
  RTMB::REPORT(p)
  RTMB::REPORT(eta)
  RTMB::REPORT(lam)
  RTMB::REPORT(gamma)
  RTMB::REPORT(nll2)
  return(nll2)
}


