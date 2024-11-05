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

  tmb.data[['data']] <<- object@data %>%
    mutate(state = ifelse(is.na(stage), "unk", stage)) %>%
    mutate(state = factor(state, levels = c("yr","unk"))) %>%
    mutate(loc = factor(loc, levels = c("Trib","First_Trap","Wen","MCJ","JDJ","BON","TWX_EST"))) %>%
    mutate(tag_site = "Trib",
           last_site = "TWX_EST") %>%
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

  phi_dim <- sapply(tmb.data[['dm']][['phi']],function(x){ncol(x)})
  tmb.data[['phi_dim']] <<- phi_dim

  eta_dim <- sapply(tmb.data[['dm']][['eta']],function(x){ncol(x)})
  tmb.data[['eta_dim']] <<- eta_dim
  if(ns < 3){
    eta_dim <- numeric(0)
    tmb.data[['eta_dim']] <<- eta_dim
  }

  parameters <<- list(phi = rep(0,sum(phi_dim)),
                     p = matrix(0, ns-1, length(levels(tmb.data$data$loc))-2),
                     lam = rep(0),
                     eta = rep(0,sum(eta_dim)))

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

