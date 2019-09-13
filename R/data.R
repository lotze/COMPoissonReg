#' Couple dataset
#' 
#' A dataset investigating the impact of education level and level of anxious
#' attachment on unwanted pursuit behaviors in the context of couple separation.
#' 
#' @usage
#' data(couple)
#' 
#' @format
#' \describe{
#' \item{UPB}{number of unwanted pursuit behavior perpetrations.}
#' \item{EDUCATION}{1 if at least bachelor's degree; 0 otherwise.}
#' \item{ANXIETY}{continuous measure of anxious attachment.}
#' }
#' 
#' @references 
#' Loeys, T., Moerkerke, B., DeSmet, O., Buysse, A., 2012. The analysis of
#' zero-inflated count data: Beyond zero-inflated Poisson regression. British
#' J. Math. Statist. Psych. 65 (1), 163-180.
#' @name couple
#' @docType data
"couple"

#' Freight dataset
#' 
#' A set of data on airfreight breakage (breakage of ampules filled with some
#' biological substance are shipped in cartons).
#' 
#' @usage
#' data(freight)
#' 
#' @format
#' \describe{
#' \item{broken}{number of ampules found broken upon arrival.}
#' \item{transfers}{number of times carton was transferred from one aircraft to
#' another.}
#' }
#' 
#' @references
#' Kutner MH, Nachtsheim CJ, Neter J (2003). Applied Linear Regression Models,
#' Fourth Edition. McGraw-Hill.
#' 
#' @name freight
#' @docType data
"freight"
