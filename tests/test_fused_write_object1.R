require("RMVL")

M1<-mvl_open("test1.mvl", append=TRUE, create=TRUE)

L1<-list(1:10, 12:41)
mvl_fused_write_objects(M1, L1, "L1")

L2<-list(runif(10), rnorm(12))
mvl_fused_write_objects(M1, L2, "L2")

L3<-list(paste("test", 1:10, sep=""), as.character(12:41))
mvl_fused_write_objects(M1, L3, "L3")

L4<-list(matrix(1:6, 2, 3), matrix(5:12, 2, 4))
mvl_fused_write_objects(M1, L4, "L4")

L5<-list(data.frame(x=1:10, b=paste("test", 10:1)), data.frame(x=5:2, b=paste("test", 2:5)))
mvl_fused_write_objects(M1, L5, "L5")

mvl_close(M1)

M2<-mvl_open("test1.mvl")

if(any(do.call(c, L1)!=M2$L1))cat("L1 error 1\n")
if(any(do.call(c, L2)!=M2$L2))cat("L2 error 1\n")
if(any(do.call(c, L3)!=M2$L3))cat("L3 error 1\n")
if(any(do.call(cbind, L4)!=M2$L4))cat("L4 error 1\n")
if(any(do.call(rbind, L5)!=M2$L5[,]))cat("L5 error 1\n")


M1<-mvl_open("test1b.mvl", append=TRUE, create=TRUE)

flatten<-function(L) {
	return(lapply(L, function(x) { 
		if(!inherits(x, "MVL_OBJECT"))return(x)
		if(mvl_inherits(x, "data.frame"))return(x[,])
		return(x[])
		}))
	}

L1<-list(10:5, M2["L1", ref=TRUE], 7:3)
mvl_fused_write_objects(M1, L1, "L1")

L2<-list(runif(10), M2["L2", ref=TRUE], rnorm(12))
mvl_fused_write_objects(M1, L2, "L2")

L3<-list(paste("test", 10:5, sep=""), M2["L3", ref=TRUE], as.character(7:3))
mvl_fused_write_objects(M1, L3, "L3")

L4<-list(matrix(2:7, 2, 3), M2["L4", ref=TRUE], matrix(9:16, 2, 4))
mvl_fused_write_objects(M1, L4, "L4")

L5<-list(data.frame(x=1:9, b=paste("test", 9:1)), M2["L5", ref=TRUE], data.frame(x=4:2, b=paste("test", 2:4)))
mvl_fused_write_objects(M1, L5, "L5")

L6<-list(M2["L5", ref=TRUE], data.frame(x=1:9, b=paste("test", 9:1)), M2["L5", ref=TRUE], data.frame(x=4:2, b=paste("test", 2:4)))
mvl_fused_write_objects(M1, L6, "L6")

mvl_write_object(M1, M2["L4", ref=TRUE], "L4copy")

mvl_close(M1)

M1<-mvl_open("test1b.mvl")

if(any(do.call(c, flatten(L1))!=M1$L1))cat("L1 error 2\n")
if(any(do.call(c, flatten(L2))!=M1$L2))cat("L2 error 2\n")
if(any(do.call(c, flatten(L3))!=M1$L3))cat("L3 error 2\n")
if(any(do.call(cbind, flatten(L4))!=M1$L4))cat("L4 error 2\n")
if(any(do.call(rbind, flatten(L5))!=M1$L5[,]))cat("L5 error 2\n")
if(any(do.call(rbind, flatten(L6))!=M1$L6[,]))cat("L6 error 2\n")

if(any(M2$L4[,]!=M1$L4copy[,]))cat("L4copy error 2\n")

mvl_close(M1)
mvl_close(M2)

unlink("test1.mvl")
unlink("test1b.mvl")




