% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RMVL.R
\name{mvl_write_object}
\alias{mvl_write_object}
\title{Write R object into MVL file}
\usage{
mvl_write_object(MVLHANDLE, x, name = NULL, drop.rownames = FALSE)
}
\arguments{
\item{MVLHANDLE}{a handle to MVL file produced by mvl_open()}

\item{x}{a suitable R object (vector, array, list, data.frame) or a vector-like MVL_OBJECT}

\item{name}{if specified add a named entry to MVL file directory}

\item{drop.rownames}{set to TRUE to prevent rownames from being written}
}
\value{
an object of class MVL_OFFSET that describes an offset into this MVL file. MVL offsets are vectors and can be concatenated. They can be written to MVL file directly, or as part of another object such as list.
}
\description{
Write R object into MVL file
}
\examples{
\dontrun{
Mtmp<-mvl_open("tmp_a.mvl", append=TRUE, create=TRUE)
mvl_write_object(Mtmp, runif(100), "vec1")
L<-list()
L[["x"]]<-mvl_write_object(Mtmp, 1:5)
L[["y"]]<-mvl_write_object(Mtmp, c("a", "b"))
L[["df"]]<-mvl_write_object(Mtmp, data.frame(x=1:100, z=runif(100)))
mvl_write_object(Mtmp, L, "L")
Mtmp<-mvl_remap(Mtmp)
print(Mtmp$L)
}
}
\seealso{
\code{\link{mvl_indexed_copy}}, \code{\link{mvl_merge}}
}
