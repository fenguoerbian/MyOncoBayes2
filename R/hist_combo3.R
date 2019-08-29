#' Dataset: historical and concurrent data on a three-way combination
#'
#' This dataset involves a hypothetical dose-escalation study of combination
#' therapy with three treatment components. From two previous studies
#' \code{HistAgent1} and \code{HistAgent2}, historical data is available on each
#' of the treatments as single-agents, as well as two of the two-way
#' combinations. However, due to a difference in treatment schedule between the
#' \code{Combo} study and the historical studies, a stratification (through \code{stratum})
#' is made between the groups to allow differential discounting of the
#' alternate-schedule data.
#'
#' @format A data frame with 18 rows and 7 variables:
#' \describe{
#'   \item{group_id}{study}
#'   \item{DosesAdm1}{dose of Drug A}
#'   \item{DosesAdm2}{dose of Drug B}
#'   \item{DosesAdm3}{dose of Drug C}
#'   \item{Npat}{number of patients}
#'   \item{Ntox}{number of DLTs}
#'   \item{stratum}{stratum for \code{group_id}'s used for differential discounting}
#' }
#'
#' @template example-start
#' @template example-combo3
#' @template example-stop
#'
"hist_combo3"
