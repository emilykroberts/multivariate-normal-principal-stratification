% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/final_results.R
\name{final_results}
\alias{final_results}
\title{results}
\usage{
final_results(params_matrix, write, holdR, holdS, res, array_id)
}
\arguments{
\item{params_matrix}{results from mcmc}

\item{write}{boolean of writing results to csv}

\item{holdR}{holdR correlations}

\item{holdS}{holdS variances}

\item{res}{results}

\item{array_id}{identifier related to random seed for simulation replication}
}
\value{
final results
}
\description{
creates output of final results
}
\examples{
example(final_results(params_matrix = params, holdR = holdR, holdS = holdS,
res = results, write = F))

}
