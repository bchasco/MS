# Generic build method
setGeneric("MStmb", function(object) standardGeneric("MStmb"))
# Define the generic function first


#' Estimate
#'
#' Build the tmb object and estimate that parameters of the model
#'
#' @param object A tmb_list object.
#' @return The estimate tmb object
#' @export
#' @examples
#' # Define MR_settings (as empty list here; customize as needed)
#' MR_settings_example <- list()
#'
#' # Create an instance of the tmb_list class using the WEN data set
#' input <- new("tmb_list",
#'       MR_settings = MR_settings_example,
#'       data = WEN)
#'
#' # Run the model
#' fit <- MStmb(input)
#'
# Method building and runnning the tmb model
setMethod("MStmb", "tmb_list", function(object) {

  design_matrices <- create_design_matrices(object@MR_settings, object@data)
  poly_basis <- create_poly_basis(object)


  tmb.data <<- list()
  tmb.data[['data']] <<- object@data
  object@MR_settings$dm <- design_matrices$dm
  tmb.data[['dm']] <<- object@MR_settings$dm

  # pred_grid <- create_expanded_grid(object)
  # pred_matrices <- create_pred_matrices(object, poly_basis, pred_grid)
  # tmb.data[['pred']] <<- pred_matrices
  # tmb.data[['pred.grid']] <<- pred_grid

  #You are going to have multiple recaptures
  #You have to create a state transition matrix for each recapture
  #Except you don't want to estimate a parameter for the first capture event
  #and the parameters for detection probability and survival for the last capture event are the same

  tmb.data[['state']] <<- unique(tmb.data$data$state)

  tmb.data[['ns']] <<- ns <- length(levels(tmb.data$state))
  tmb.data[['ni']] <<- ni <- length(unique(tmb.data$data$id))

  phi_dim <- list()
  phi_re_dim <- list()
  phi_states <- length(tmb.data[['dm']][['phi']])
  for(i in 1:phi_states){
    if(length(object@MR_settings$frms$phi)>1){
      phi_dim[[i]] <- ncol(tmb.data[['dm']][['phi']][[i]]$X)
    }else{
      phi_dim[[1]] <- ncol(tmb.data[['dm']][['phi']][[1]]$X)
    }
    if(length(tmb.data[['dm']][['phi']][[i]][['reTrms']])>0){
      phi_re_dim[[i]] <- sum(tmb.data[['dm']][['phi']][[i]][['reTrms']]$nl)
    }else{
      phi_re_dim[[i]] <- 0
    }
  }
  tmb.data[['phi_dim']] <<- phi_dim
  tmb.data[['phi_re_dim']] <<- phi_re_dim

  p_dim <- list()
  p_re_dim <- list()
  p_states <- length(tmb.data[['dm']][['p']])
  for(i in 1:p_states){
    if(length(object@MR_settings$frms$p)>1){
      p_dim[[i]] <- ncol(tmb.data[['dm']][['p']][[i]]$X)
    }else{
      p_dim[[1]] <- ncol(tmb.data[['dm']][['p']][[1]]$X)
    }
    if(length(tmb.data[['dm']][['p']][[i]][['reTrms']])>0){
      p_re_dim[[i]] <- sum(tmb.data[['dm']][['p']][[i]][['reTrms']]$nl)
    }else{
      p_re_dim[[i]] <- 0
    }
  }
  tmb.data[['p_dim']] <<- p_dim
  tmb.data[['p_re_dim']] <<- p_re_dim

  lam_dim <- list()
  lam_states <- length(tmb.data[['dm']][['lam']])
  for(i in 1:lam_states){
    if(length(object@MR_settings$frms$lam)>1){
      lam_dim[[i]] <- ncol(tmb.data[['dm']][['lam']][[i]]$X)
    }else{
      lam_dim[[1]] <- ncol(tmb.data[['dm']][['lam']][[1]]$X)
    }
  }
  tmb.data[['lam_dim']] <<- lam_dim

  eta_dim <- sapply(tmb.data[['dm']][['eta']],function(x){ncol(x$X)})
  eta_dim <- eta_dim[1]
  tmb.data[['eta_dim']] <<- eta_dim
  if(ns < 3){
    eta_dim <- numeric(0)
    tmb.data[['eta_dim']] <<- eta_dim
  }

  tmb.data$MR_settings <<- object@MR_settings

  parameters <<- list(phi = rep(-3 , sum(unlist(phi_dim))),
                      phi_re = rep(0 , sum(unlist(phi_re_dim))),
                     p = rep(-2, sum(unlist(p_dim))),
                     p_re = rep(0 , sum(unlist(p_re_dim))),
                     lam = rep(-2, sum(unlist(lam_dim))),
                     eta = rep(0 , sum(unlist(eta_dim))),
                     sig = 0)

  if(object@MR_settings$mod == "MSt"){
    parameters$tau <<- rep(0, length(unique(tmb.data$data$loc)) - 1)
    parameters$fsig_tau <<- rep(0, length(unique(tmb.data$data$loc)) - 1)
  }

  environment(test_f) <- .GlobalEnv

  obj <- RTMB::MakeADFun(test_f,
                         parameters,
                         silent = FALSE,
                         # map = list(sig = as.factor(NA)),
                         random = c("phi_re", "p_re"))
  #
  object@TMB$obj <- obj
  object@TMB$obj$tmb.data <- tmb.data

  object@TMB$opt <- nlminb(obj$par,
                           obj$fn,
                           obj$gr
                           )
  object@TMB$rep <- obj$report()

  object@TMB$parameters <- parameters
  invisible(return(object))  # Return the updated object with results

})

