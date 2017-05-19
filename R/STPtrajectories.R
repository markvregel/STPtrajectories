#' STPtrajectories.
#'
#' @name STPtrajectories
#' @docType package
#' @description  Package for handling Space-Time Prism(STP) trajectories.
#'It contains functions to calculate Potential Path Areas(PPAs), create random
#'trajectories and to test for possible encounters by applying the alibi query.
#'It also provides functions to visulize the STPs treajectories in 3D.
#' @section Background:
#' A trajectory consists of successive points in space and time.
#' The location of the individual between two successive points is unknown,
#' but based on a maximum speed a space-time prism can be calculated.
#' A space-time prism(STP) is defined as the collection of space-time locations the individual could reach,
#' given a speed limitation. The \link{STP_Track} class can be used to handle space-time prism trajectories.
#'
#' The \link{alibi_query} uses this concept to test if two individuals could have met each other.
#' The package provides related functions to help users in their analysis.
#' of their trajectories and take into account the uncertainty about the location of an individual.
#' These include methods to visualise space-time prisms(\link{STP_plot}),
#' calculate the potential path area(\link{calculate_PPA}) and a random trajectory generator(\link{RTG}).
#' @section help:
#'  \strong{need help:}
#' \itemize{
#' \item manual: \url{https://github.com/markvregel/STPtrajectories/blob/master/STPtrajectories.pdf}
#' \item vignette: \url{http://htmlpreview.github.io/?https://raw.githubusercontent.com/markvregel/STPtrajectories/master/vignettes/STP_Tracks.html}
#' \item check functions help and examples
#' }
#' @seealso
#' \itemize{
#' \item github: \url{https://github.com/markvregel/STPtrajectories}
#' \item trajecotries package \url{https://cran.rstudio.com/web/packages/trajectories/index.html}
#' }
#' @author Mark ten Vregelaar
#'
NULL
