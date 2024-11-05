# Define the general Model class
setClass(
  "Model",
  slots = list(
    data = "ANY",
    model_type = "character"
  )
)
#' tmb_list Class
#'
#' This class holds the data needed for TMB-based mark-recapture analysis.
#'
#' @slot MR_settings A list of mark-recapture settings, expected to contain ...
#' @slot data A list of TMB data. Must contain `id`, `loc`, and `n` columns in its data frames.
#' @export
setClass("tmb_list",
         slots = list(
           MR_settings = "list",
           data = "data.frame",
           TMB = 'list'
         ),
         validity = function(object) {
           # Check MR_settings is a list with required elements
           if (!is.list(object@MR_settings)) {
             return("MR_settings must be a list.")
           }

           # Check data contains required columns
           if (!is.data.frame(object@data) || !all(c("id", "loc", "n") %in% names(object@data))) {
             return("data must be a data.frame containing 'id' and 'state' elements.")
           }
           if (!is.numeric(object@data$n)) {
             return("\n***data$n must be numeric***")
           }
           if(!('time') %in% names(object@data)){
             print("Without a column for time, the temporal aspect of the model, \nwith be set to one for all observations.")
           }


           TRUE
         }
)

#'TMBmodel class
#' @slot TMB The TMB obj
#' @export
setClass(
  "tmb_model",
  contains = "Model",
  slots = list(
    TMB = "list"
  )
)
