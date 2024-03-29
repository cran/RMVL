% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RMVL.R
\name{mvl_compute_repeats}
\alias{mvl_compute_repeats}
\title{Find stretches of repeated rows among vectors}
\usage{
mvl_compute_repeats(L)
}
\arguments{
\item{L}{list of vector like MVL_OBJECTs}
}
\value{
partition describing repeated rows
}
\description{
This function is passed a list of vector like MVL_OBJECTs which are considered as columns in a table.
It returns a vector V starting with 1 and ending with number of rows plus 1, so that stretches of repeated rows can be found as V[i]:V[i+1]
}
