#' Wenatchee example data
#'
#' A dataset containing mark recapture information in the Wenatchee Basin
#'
#' @format A data frame with X rows and Y columns:
#' \describe{
#'   \item{loc}{Detection location}
#'   \item{cum_time}{cumulative time between detections}
#'   \item{init_year}{Year when the fish was first tagged}
#'   \item{event_day}{Day of year for the detection}
#'   \item{stageID}{The stage of the fish based on the its migration day}
#'   \item{n}{The number of fish observed to have the same unique combination over all locations}
#'   \item{id}{An index of the unique combination}
#'
#'   ...
#' }
#' @source Look at the build_wenatchee_data.r file in the data folder. V2 does not work.
"WEN"
