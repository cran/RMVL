require("RMVL")
 
M3<-mvl_open("test2.mvl", append=TRUE, create=TRUE)
L<-list()
df<-data.frame(x=(1:1e5) %% 11, y=round(cos(1:1e5)*1000), s=rep(c("a", "b"), 5e4))
mvl_write_object(M3, df, "test_object")

M3<-mvl_remap(M3)

mvl_write_extent_index(M3, list(M3$test_object[,"x", ref=TRUE]), "test_object_extent_index_x")
mvl_write_extent_index(M3, list(M3$test_object[,"s", ref=TRUE]), "test_object_extent_index_s")

mvl_close(M3)

M3<-mvl_open("test2.mvl")
L2<-M3["test_object", ref=TRUE]

N<-dim(df)[1]

idx0<-50:100

S1<-sum(df[df[,"x"] %in% idx0, "y"])
S2<-sum(unlist(mvl_extent_index_lapply(M3["test_object_extent_index_x", ref=TRUE], list(as.integer(idx0)), function(i, idx) {
	return(sum(M3$test_object[idx, "y"]))
	})))
if(abs(S1-S2)>1e-5)cat("test1 failed: S1=", S1, "S2=", S2, "\n")

S2<-sum(unlist(mvl_extent_index_lapply(M3["test_object_extent_index_x", ref=TRUE], list(as.numeric(idx0)), function(i, idx) {
	return(sum(M3$test_object[idx, "y"]))
	})))
if(abs(S1-S2)>1e-5)cat("test2 failed: S1=", S1, "S2=", S2, "\n")

S1<-sum(df[df[,"s"]=="a", "y"])
S2<-sum(unlist(mvl_extent_index_lapply(M3["test_object_extent_index_s", ref=TRUE], list("a"), function(i, idx) {
	return(sum(M3$test_object[idx, "y"]))
	})))
if(abs(S1-S2)>1e-5)cat("test2 failed: S1=", S1, "S2=", S2, "\n")

mvl_close(M3)

unlink("test2.mvl")

