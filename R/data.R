#' Example nodelist file of a bipartite network.
#'
#' A dataset consisting of anonymized patients (n=798) and symptoms (d=8), where each patient has one or more symptoms.
#'
#' @format A data frame with 806 rows and 5 columns:
#' \describe{
#'   \item{Label}{Node ID of patients and symptoms.}
#'   \item{Cluster}{Cluster membership found by BipartiteModularityMaximization. Can be changed to any customized cluster membership.}
#'   \item{X}{X coordinates found by Fruchterman & Reingold layout. Can be changed to any customized coordinates. }
#'   \item{Y}{Y coordinates found by Fruchterman & Reingold layout. Can be changed to any customized coordinates. }
#'   \item{Entity}{Indicating whether a node is a patient (assigned 1) or a symptom (assigned 2).}
#' }
"example_nodelist"


#' Sample dataset of a bipartite network..
#'
#' The sample dataset consists of a simulated data consisting of patients (n=798) and symptoms (d=8), where each patient has one or more symptoms referred to here as an incidence matrix.
#'
#' @format A data frame with 798 rows and 8 binary variables:
#' \describe{
#'   \item{Symptom_1}{}
#'   \item{Symptom_2}{}
#'   \item{Symptom_3}{}
#'   \item{Symptom_4}{}
#'   \item{Symptom_5}{}
#'   \item{Symptom_6}{}
#'   \item{Symptom_7}{}
#'   \item{Symptom_8}{}
#' }
"example_incidmat"
