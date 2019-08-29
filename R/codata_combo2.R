#' Dataset: historical and concurrent data on a two-way combination
#'
#' One of two datasets from the application described in Neuenschwander et al
#' (2016). In the study \code{trial_AB}, the risk of DLT was studied as a function of
#' dose for two drugs, drug A and drug B. Historical information on the toxicity
#' profiles of these two drugs is available from single agent trials \code{trial_A}
#' and \code{trial_B}. Another study \code{IIT} was run concurrently to \code{trial_AB}, and
#' studies the same combination. A second dataset \code{hist_combo2} is
#' available from this example, which includes only the data from the single
#' agent studies, prior to the initiation of \code{trial_AB} and \code{IIT}.
#'
#' @format A data frame with 20 rows and 5 variables:
#' \describe{
#'   \item{group_id}{study}
#'   \item{DosesAdm1}{dose of Drug A}
#'   \item{DosesAdm2}{dose of Drug B}
#'   \item{Npat}{number of patients}
#'   \item{Ntox}{number of DLTs}
#' }
#' @template ref-mac
#'
#' @template example-start
#' @template example-combo2
#' @template example-stop
#'
"codata_combo2"
