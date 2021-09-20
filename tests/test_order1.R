require("RMVL")
 
M3<-mvl_open("test2.mvl", append=TRUE, create=TRUE)
L<-list()
df<-data.frame(x=1:1e5, y=rnorm(1e5), s=rep(c("a", "b"), 5e4))
L[["x"]]<-mvl_write_object(M3, df)
L[["description"]]<-"Example of large data frame"
mvl_write_object(M3, L, "test_object")
mvl_close(M3)

M3<-mvl_open("test2.mvl")
L2<-M3["test_object"]

N<-dim(df)[1]

idx0<-50:100

idx1a<-idx0[order(df[idx0, "x"])]
idx1b<-mvl_order_vectors(list(M3$test_object$x[,"x"]), idx0)
if(!isTRUE(all.equal(idx1a, idx1b)))cat("test1 failed\n")

idx1a<-idx0[order(df[idx0, "y"])]
idx1b<-mvl_order_vectors(list(M3$test_object$x[,"y"]), idx0)
if(!isTRUE(all.equal(idx1a, idx1b)))cat("test2 failed\n")

idx1a<-idx0[order(df[idx0, "s"])]
idx1b<-mvl_order_vectors(list(M3$test_object$x[,"s"]), idx0)
if(!isTRUE(all.equal(idx1a, idx1b)))cat("test3 failed\n")

idx1a<-idx0[order(df[idx0, "s"], df[idx0, "y"])]
idx1b<-mvl_order_vectors(list(M3$test_object$x[,"s"], M3$test_object$x[,"y"]), idx0)
if(!isTRUE(all.equal(idx1a, idx1b)))cat("test4 failed\n")

idx1a<-order(df[, "y"])
idx1b<-mvl_order_vectors(list(M3$test_object$x[,"y"]))
if(!isTRUE(all.equal(idx1a, idx1b)))cat("test5 failed\n")

idx1a<-order(df[, "s"], df[, "y"])
idx1b<-mvl_order_vectors(list(M3$test_object$x[,"s"], M3$test_object$x[,"y"]))
if(!isTRUE(all.equal(idx1a, idx1b)))cat("test6 failed\n")

idx1a<-order(df[, "s"], df[, "y"], decreasing=TRUE)
idx1b<-mvl_order_vectors(list(M3$test_object$x[,"s"], M3$test_object$x[,"y"]), decreasing=TRUE)
if(!isTRUE(all.equal(idx1a, idx1b)))cat("test7 failed\n")

idx1a<-idx0[order(df[idx0, "s"], df[idx0, "y"], decreasing=TRUE)]
idx1b<-mvl_order_vectors(list(M3$test_object$x[,"s"], M3$test_object$x[,"y"]), idx0, decreasing=TRUE)
if(!isTRUE(all.equal(idx1a, idx1b)))cat("test8 failed\n")


mvl_close(M3)

unlink("test2.mvl")

