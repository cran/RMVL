require("RMVL")

M<-mvl_open("test_checksums1.mvl", append=TRUE, create=TRUE)

mvl_write_object(M, rnorm(10000), "x")
mvl_write_object(M, 1:100, "y")
mvl_write_object(M, as.character(1:100), "s")
mvl_write_object(M, as.list(1:100), "L")
mvl_write_object(M, data.frame(i=1:100, x=rnorm(100), s=as.character(1:100)), "df")

M<-mvl_remap(M)

mvl_write_extent_index(M, list(M$df[,"i", ref=TRUE]), "extent_index")

mvl_close(M, checksums=TRUE)


M<-mvl_open("test_checksums1.mvl")
mvl_verify(M)
mvl_verify(M["x", ref=TRUE])
mvl_verify(M["y", ref=TRUE])
mvl_verify(M["s", ref=TRUE])
mvl_verify(M["L", ref=TRUE])
mvl_verify(M["df", ref=TRUE])
mvl_verify(M["extent_index", ref=TRUE])
mvl_close(M)

unlink("test_checksums1.mvl")
