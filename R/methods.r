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

  object@MR_settings$dm <- create_design_matrices(object@MR_settings, object@data)

  tmb.data <<- list()

  tmb.data[['data']] <<- object@data %>%
    mutate(state = ifelse(is.na(stage), "unk", "yr")) %>%
    mutate(state = factor(state, levels = c("yr","unk"))) %>%
    mutate(loc = factor(loc, levels = c("Trib","First_Trap","RIS_RIA","WEN","MCJ","JDJ","BON","TWX_EST"))) %>%
    droplevels()

  #
  # tmp_WEN <- tmb.data[['data']] %>%
  #   filter(loc!=tag_site &
  #                  loc!=last_site) %>%
  #   droplevels()
  #
  #
  # m <-  tryCatch(model.matrix(formula('state ~ loc'), tmp_WEN),
  #                error = function(e) e,
  #                warning = function(w) w)
  #
  #
  # tmb.data[['m']] <<- m

  tmb.data[['dm']] <<- object@MR_settings$dm

  #You are going to have multiple recaptures
  #You have to create a state transition matrix for each recapture
  #Except you don't want to estimate a parameter for the first capture event
  #and the parameters for detection probability and survival for the last capture event are the same

  tmb.data[['state']] <<- unique(tmb.data$data$state)

  tmb.data[['ns']] <<- ns <- length(levels(tmb.data$state))
  tmb.data[['ni']] <<- ni <- length(unique(tmb.data$data$id))

  phi_dim <- list()
  phi_states <- length(tmb.data[['dm']][['phi']])
  for(i in 1:phi_states){
    phi_dim[[i]] <- ncol(tmb.data[['dm']][['phi']][[i]])
  }
  tmb.data[['phi_dim']] <<- phi_dim

  p_dim <- list()
  p_states <- length(tmb.data[['dm']][['p']])
  for(i in 1:p_states){
    p_dim[[i]] <- ncol(tmb.data[['dm']][['p']][[i]])
  }
  tmb.data[['p_dim']] <<- p_dim

  lam_dim <- list()
  lam_states <- length(tmb.data[['dm']][['lam']])
  for(i in 1:lam_states){
    lam_dim[[i]] <- ncol(tmb.data[['dm']][['lam']][[i]])
  }
  tmb.data[['lam_dim']] <<- lam_dim

  eta_dim <- sapply(tmb.data[['dm']][['eta']],function(x){ncol(x)})
  tmb.data[['eta_dim']] <<- eta_dim
  if(ns < 3){
    eta_dim <- numeric(0)
    tmb.data[['eta_dim']] <<- eta_dim
  }

  tmb.data$MR_settings <<- object@MR_settings

  parameters <<- list(phi = rep(0 , sum(unlist(phi_dim))),
                     p = rep(0 , sum(unlist(p_dim))),
                     lam = rep(0 , sum(unlist(lam_dim))),
                     eta = rep(0 , sum(eta_dim)))

  environment(test_f) <- .GlobalEnv

  obj <- RTMB::MakeADFun(test_f,
                         parameters)

  object@TMB$obj <- obj
  object@TMB$obj$tmb.data <- tmb.data

  object@TMB$opt <- nlminb(obj$par,
                           obj$fn,
                           obj$gr
                           )
  # object@TMB$frm <- frm
  # object@TMB$m <- m
  object@TMB$rep <- obj$report()

  object@TMB$parameters <- parameters
  invisible(return(object))  # Return the updated object with results

})

