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

  tmb.data <<- list()


  design_matrices <- create_design_matrices(object@MR_settings, object@data)
  poly_basis <- create_poly_basis(object)

  tmb.data[['data']] <<- object@data
  object@MR_settings$dm <- design_matrices$dm
  tmb.data[['dm']] <<- object@MR_settings$dm

  # print(names(object@MR_settings))
  unique_cols <- extract_unique_vars(object@MR_settings$frms)
  unique_cols <- unique_cols[!(unique_cols%in%c("loc"))] #Need this because you add a time formula that may include location

  pred_grid <<- expanded_result <- expand_grid(
    data = object@data,
    cols = unique_cols,
    include_time = object@MR_settings$mod=="MSt"
  )
  # print(pred_grid)
  #
  pred_matrices <- create_pred_matrices(object, poly_basis, pred_grid)
  tmb.data[['pred']] <<- pred_matrices
  tmb.data[['pred.grid']] <<- pred_grid

  #You are going to have multiple recaptures
  #You have to create a state transition matrix for each recapture
  #Except you don't want to estimate a parameter for the first capture event
  #and the parameters for detection probability and survival for the last capture event are the same

  tmb.data[['state']] <<- unique(tmb.data$data$state)

  tmb.data[['ns']] <<- ns <- length(levels(tmb.data$state))
  tmb.data[['ni']] <<- ni <- length(unique(tmb.data$data$id))

  phi_dim <- list()
  phi_re_dim <- list()
  phi_re_sig <- list()

  phi_states <- length(tmb.data[['dm']][['phi']])
  for(i in 1:(ns-1)){
    if(length(object@MR_settings$frms$phi)>1){
      phi_dim[[i]] <- ncol(tmb.data[['dm']][['phi']][[i]]$X)

      if(length(tmb.data[['dm']][['phi']][[i]][['reTrms']])>0){
          phi_re_dim[[i]] <- sum(tmb.data[['dm']][['phi']][[i]][['reTrms']]$nl)
          phi_re_sig[[i]] <- length(unique(tmb.data[['dm']][['phi']][[i]][['reTrms']]$Lind))
      }else{
        phi_re_dim[[i]] <- 0
        phi_re_sig[[i]] <- NULL
      }
    }else{
      phi_dim[[1]] <- ncol(tmb.data[['dm']][['phi']][[1]]$X)
      if(length(tmb.data[['dm']][['phi']][[1]][['reTrms']])>0 & i == 1){
          phi_re_dim[[i]] <- sum(tmb.data[['dm']][['phi']][[1]][['reTrms']]$nl)
          phi_re_sig[[i]] <- length(unique(tmb.data[['dm']][['phi']][[i]][['reTrms']]$Lind))
      }else{
        phi_re_dim[[i]] <- 0
        phi_re_sig[[i]] <- NULL
      }
    }

  }
  tmb.data[['phi_dim']] <<- phi_dim
  tmb.data[['phi_re_dim']] <<- phi_re_dim

  p_dim <- list()
  p_re_dim <- list()
  p_re_sig <- list()
  p_states <- length(tmb.data[['dm']][['p']])
  for(i in 1:(ns-1)){
    if(length(object@MR_settings$frms$p)>1){
      p_dim[[i]] <- ncol(tmb.data[['dm']][['p']][[i]]$X)

      if(length(tmb.data[['dm']][['p']][[i]][['reTrms']])>0){
        p_re_dim[[i]] <- sum(tmb.data[['dm']][['p']][[i]][['reTrms']]$nl)
        p_re_sig[[i]] <- length(unique(tmb.data[['dm']][['p']][[i]][['reTrms']]$Lind))
      }else{
        p_re_dim[[i]] <- 0
        p_re_sig[[i]] <- NULL
      }
    }else{
      p_dim[[1]] <- ncol(tmb.data[['dm']][['p']][[1]]$X)
      if(length(tmb.data[['dm']][['p']][[1]][['reTrms']])>0 & i == 1){
        p_re_dim[[i]] <- sum(tmb.data[['dm']][['p']][[1]][['reTrms']]$nl)
        p_re_sig[[i]] <- length(unique(tmb.data[['dm']][['p']][[i]][['reTrms']]$Lind))
      }else{
        p_re_dim[[i]] <- 0
        p_re_sig[[i]] <- NULL
      }
    }

  }
  tmb.data[['p_dim']] <<- p_dim
  tmb.data[['p_re_dim']] <<- p_re_dim

  lam_dim <- list()
  lam_re_dim <- list()
  lam_re_sig <- list()
  lam_states <- length(tmb.data[['dm']][['lam']])
  for(i in 1:lam_states){
    if(length(object@MR_settings$frms$lam)>1){
      lam_dim[[i]] <- ncol(tmb.data[['dm']][['lam']][[i]]$X)
    }else{
      lam_dim[[1]] <- ncol(tmb.data[['dm']][['lam']][[1]]$X)
    }
    if(length(tmb.data[['dm']][['lam']][[i]][['reTrms']])>0){
      lam_re_dim[[i]] <- sum(tmb.data[['dm']][['lam']][[i]][['reTrms']]$nl)
      lam_re_sig[[i]] <- length(unique(tmb.data[['dm']][['lam']][[i]][['reTrms']]$Lind))
    }else{
      lam_re_dim[[i]] <- 0
    }
  }
  tmb.data[['lam_dim']] <<- lam_dim
  tmb.data[['lam_re_dim']] <<- lam_re_dim


  eta_dim <- sapply(tmb.data[['dm']][['eta']],function(x){ncol(x$X)})[1]
  eta_re_dim <- list()
  eta_re_sig <- list()
  eta_states <- length(tmb.data[['dm']][['eta']])
  tmb.data[['eta_dim']] <<- eta_dim
  if(ns < 3){
    eta_dim <- numeric(0)
    tmb.data[['eta_dim']] <<- eta_dim
  }
  if(length(tmb.data[['dm']][['eta']][[1]][['reTrms']])>0){
    eta_re_dim[[1]] <- sum(tmb.data[['dm']][['eta']][[1]][['reTrms']]$nl)
    eta_re_sig[[1]] <- length(unique(tmb.data[['dm']][['eta']][[i]][['reTrms']]$Lind))
  }else{
    eta_re_dim[[1]] <- 0
  }
  tmb.data[['eta_re_dim']] <<- eta_re_dim

  tau_dim <- list()
  tau_re_dim <- list()
  tau_re_sig <- list()
  tau_states <- length(tmb.data[['dm']][['tau']])
  for(i in 1:(ns-1)){
    if(length(object@MR_settings$frms$tau)>1){
      tau_dim[[i]] <- ncol(tmb.data[['dm']][['tau']][[i]]$X)

      if(length(tmb.data[['dm']][['tau']][[i]][['reTrms']])>0){
        tau_re_dim[[i]] <- sum(tmb.data[['dm']][['tau']][[i]][['reTrms']]$nl)
        tau_re_sig[[i]] <- length(unique(tmb.data[['dm']][['tau']][[i]][['reTrms']]$Lind))
      }else{
        tau_re_dim[[i]] <- 0
        tau_re_sig[[i]] <- NULL
      }
    }else{
      tau_dim[[1]] <- ncol(tmb.data[['dm']][['tau']][[1]]$X)
      if(length(tmb.data[['dm']][['tau']][[1]][['reTrms']])>0 & i == 1){
        tau_re_dim[[i]] <- sum(tmb.data[['dm']][['tau']][[1]][['reTrms']]$nl)
        tau_re_sig[[i]] <- length(unique(tmb.data[['dm']][['tau']][[i]][['reTrms']]$Lind))
      }else{
        tau_re_dim[[i]] <- 0
        tau_re_sig[[i]] <- NULL
      }
    }

  }
  tmb.data[['tau_dim']] <<- tau_dim
  if(is.null(tmb.data$data$time)){
    tmb.data[['tau_dim']] <<- 0
  }
  tmb.data[['tau_re_dim']] <<- tau_re_dim

  tmb.data$MR_settings <<- object@MR_settings

  parameters <<- list(phi = rep(-0.04 , sum(unlist(phi_dim))),
                      phi_re = rep(0 , sum(unlist(phi_re_dim))),
                      phi_re_sig = rep(0 , sum(unlist(phi_re_sig))),
                      p = rep(-0.4, sum(unlist(p_dim))),
                      p_re = rep(0 , sum(unlist(p_re_dim))),
                      p_re_sig = rep(0 , sum(unlist(p_re_sig))),
                      lam = rep(-3, sum(unlist(lam_dim))),
                      lam_re = rep(0 , sum(unlist(lam_re_dim))),
                      lam_re_sig = rep(1 , sum(unlist(lam_re_sig))),
                      eta = rep(0 , sum(unlist(eta_dim))),
                      eta_re = rep(0 , sum(unlist(eta_re_dim))),
                      eta_re_sig = rep(0 , sum(unlist(eta_re_sig))),
                      tau = rep(2, sum(unlist(tau_dim))),
                      tau_sig = rep(1,length(tau_dim)),
                      tau_re = rep(0 , sum(unlist(tau_re_dim))),
                      tau_re_sig = rep(2 , sum(unlist(tau_re_sig)))
                      )


  environment(test_f) <- .GlobalEnv

  print("test")
  map <- list()
  if(MR_settings$mod=="MS"){
    map = append(map,
                 list(tau_sig = as.factor(rep(NA,length(tau_dim)))
                      ,tau = as.factor(rep(NA, sum(unlist(tau_dim))))
                      ,tau_re = as.factor(rep(NA , sum(unlist(tau_re_dim))))
                      ,tau_re_sig = as.factor(rep(NA , sum(unlist(tau_re_sig))))
                 )
                 )
  }
  print(map)

  obj <- RTMB::MakeADFun(test_f,
                         parameters,
                         silent = FALSE,
                         map = map,
                         random = c("phi_re", "p_re","lam_re","tau_re","eta_re"))
  #
  # if(MR)
  object@TMB$obj <- obj
  object@TMB$obj$tmb.data <- tmb.data

  print("test")
  object@TMB$opt <- nlminb(obj$par,
                           obj$fn,
                           obj$gr
                           )
  object@TMB$rep <- obj$report()
  opt <- object@TMB$opt
  object@TMB$AICc <- opt$objective +
    2 * length(opt$par)
  object@TMB$parameters <- parameters
  invisible(return(object))  # Return the updated object with results

})

