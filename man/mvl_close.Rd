% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RMVL.R
\name{mvl_close}
\alias{mvl_close}
\title{Close MVL file}
\usage{
mvl_close(MVLHANDLE)
}
\arguments{
\item{MVLHANDLE}{handle to opened MVL file as generated by mvl_open()}
}
\value{
None
}
\description{
Closes MVL file releasing all resources.
For read-only files the memory is unmapped, reducing the virtual memory footprint.
For files opened for writing the directory is written out, so it is important to call \code{mvl_close} or the newly written file will be corrupt.
After \code{mvl_close()} all previously obtained MVL_OBJECT's with this handle become invalid.
}
\seealso{
\code{\link{mvl_open}}, \code{\link{mvl_remap}}
}
