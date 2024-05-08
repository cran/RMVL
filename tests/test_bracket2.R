require("RMVL")

M3<-mvl_open("test_bracket2a.mvl", append=TRUE, create=TRUE)

L<-list()

df<-data.frame(x=1:1e5, y=rnorm(1e5), s=rep(c("a", "b"), 5e4), b=rnorm(1e5)<0.5)
L[["x"]]<-mvl_write_object(M3, df)

aa<-array(rnorm(10000), c(10, 50, 20))
L[["y"]]<-aa

mm<-matrix(rnorm(10000), 10, 1000)
L[["z"]]<-mm

LL2<-as.list(rnorm(10000))
names(LL2)<-paste("x", 1:10000, sep="")
L[["LL2"]]<-LL2

L[["description"]]<-"Example of large data frame"
mvl_write_object(M3, L, "test_object")

LM1<-lm(rnorm(100)~runif(100))
mvl_write_serialized_object(M3, LM1, "LM1")

mvl_close(M3)


M3<-mvl_open("test_bracket2a.mvl")
print(names(M3))

L2<-M3["test_object", ref=TRUE]

N<-dim(df)[1]

compare_df<-function(x, y) {
	if(length(dim(x))!=length(dim(y)))return(FALSE)
	if(any(dim(x)!=dim(y)))return(FALSE)
	if(any(names(x)!=names(y)))return(FALSE)
	if(dim(x)[2]>0) {
		for(i in 1:(dim(x)[2])) {
			if(any(x[,i]!=y[,i]))return(FALSE)
			}
		}
	return(TRUE)
	}

if(!compare_df(df, mvl2R(L2[["x"]]))) {
	cat("test1a failed\n")
	print(attributes(df))
	print(attributes(mvl2R(L2[["x"]])))
	cat("-----------\n")
	}
	
if(!isTRUE(all.equal(aa, mvl2R(L2[["y"]])))) {
	cat("test1b failed\n")
	print(all.equal(aa, mvl2R(L2[["y"]])))
	print(attributes(aa))
	print(attributes(mvl2R(L2[["y"]])))
	cat("-----------\n")
	}
	
if(!compare_df(mm, mvl2R(L2[["z"]]))) {
	cat("test1c failed\n")
	print(all.equal(mm, mvl2R(L2[["z"]])))
	print(attributes(mm))
	print(attributes(mvl2R(L2[["z"]])))
	cat("-----------\n")
	}
	
if(!isTRUE(all.equal(LL2, mvl2R(L2[["LL2"]])))) {
	cat("test1d failed\n")
	print(all.equal(LL2, mvl2R(L2[["LL2"]])))
	print(attributes(LL2))
	print(attributes(mvl2R(L2[["LL2"]])))
	cat("-----------\n")
	}
	
if(!isTRUE(all.equal("Example of large data frame", L2[["description"]]))) {
	cat("test1e failed\n")
	print(all.equal("Example of large data frame", L2[["description"]]))
	print(attributes("Example of large data frame"))
	print(attributes(L2[["description"]]))
	cat("-----------\n")
	}

# # R behaviour is mixed in this situation
# # For lists R returns empty list, but (1:5)[[NA]] throws an exception
# # It would not be unreasonable to think that vec[[NA]] should be NA
# # On the other hand, subscripting with NA is inefficient, and throwing an exception
# # forces to filter out NAs first
# # For now, we throw an exception and bypass the test
# if(!isTRUE(all.equal(L[[NA]], mvl2R(L2[[NA]])))) {
# 	cat("test1e failed\n")
# 	print(all.equal(L[[NA]], mvl2R(L2[[NA]])))
# 	print(attributes(L[[NA]]))
# 	print(attributes(mvl2R(L2[[NA]])))
# 	cat("-----------\n")
# 	}

if(!compare_df(df, mvl2R(L2[[1]]))) {
	cat("test1f failed\n")
	print(attributes(df))
	print(attributes(mvl2R(L2[[1]])))
	cat("-----------\n")
	}
	
if(!isTRUE(all.equal(aa, mvl2R(L2[[2]])))) {
	cat("test1g failed\n")
	print(all.equal(aa, mvl2R(L2[[2]])))
	print(attributes(aa))
	print(attributes(mvl2R(L2[[2]])))
	cat("-----------\n")
	}
	
if(!compare_df(mm, mvl2R(L2[[3]]))) {
	cat("test1h failed\n")
	print(all.equal(mm, mvl2R(L2[[3]])))
	print(attributes(mm))
	print(attributes(mvl2R(L2[[3]])))
	cat("-----------\n")
	}
	
if(!isTRUE(all.equal(LL2, mvl2R(L2[[4]])))) {
	cat("test1i failed\n")
	print(all.equal(LL2, mvl2R(L2[[4]])))
	print(attributes(LL2))
	print(attributes(mvl2R(L2[[4]])))
	cat("-----------\n")
	}
	
if(!isTRUE(all.equal("Example of large data frame", L2[[5]]))) {
	cat("test1j failed\n")
	print(all.equal("Example of large data frame", L2[[5]]))
	print(attributes("Example of large data frame"))
	print(attributes(L2[[5]]))
	cat("-----------\n")
	}

# We need check.attributes=FALSE because on R 4.1.x there is a mismatch in array attributes - one has class and the other does not
if(!isTRUE(all.equal(L[c(2, 3)], L2[c(2,3), recurse=TRUE], check.attributes=FALSE))) {
	cat("test2a failed\n")
	print(all.equal(L[c(2, 3)], L2[c(2, 3), recurse=TRUE], check.attributes=FALSE))
	print(attributes(L[c(2, 3)]))
	print(attributes(L2[c(2, 3), recurse=TRUE]))
	cat("-----------\n")
	}

# Some of the names are NA and all.equal() does not handle this properly
if(!isTRUE(all.equal.list(L[c(2, NA, 3)], L2[c(2, NA, 3), recurse=TRUE], use.names=FALSE, check.attributes=FALSE))) {
	cat("test2b failed\n")
	print(all.equal.list(L[c(2, NA, 3)], L2[c(2, NA, 3), recurse=TRUE], use.names=FALSE, check.attributes=FALSE))
	print(attributes(L[c(2, NA, 3)]))
	print(attributes(L2[c(2, NA, 3), recurse=TRUE]))
	cat("-----------\n")
	}
	
if(!isTRUE(all.equal(L[c("y", "z")], L2[c("y", "z"), recurse=TRUE], check.attributes=FALSE))) {
	cat("test2c failed\n")
	print(all.equal(L[c("y", "z")], L2[c("y", "z"), recurse=TRUE], check.attributes=FALSE))
	print(attributes(L[c("y", "z")]))
	print(attributes(L2[c("y", "z"), recurse=TRUE]))
	cat("-----------\n")
	}

if(!isTRUE(all.equal(L[c("W", "y", "z")], L2[c("W", "y", "z"), recurse=TRUE], check.attributes=FALSE))) {
	cat("test2d failed\n")
	print(all.equal(L[c("W", "y", "z")], L2[c("W", "y", "z"), recurse=TRUE], check.attributes=FALSE))
	print(attributes(L[c("W", "y", "z")]))
	print(attributes(L2[c("W", "y", "z"), recurse=TRUE]))
	cat("-----------\n")
	}

if(!isTRUE(all.equal(L[c("W", "y", NA, "z")], L2[c("W", "y", NA, "z"), recurse=TRUE], check.attributes=FALSE))) {
	cat("test2e failed\n")
	print(all.equal(L[c("W", "y", NA, "z")], L2[c("W", "y", NA, "z"), recurse=TRUE], check.attributes=FALSE))
	print(attributes(L[c("W", "y", NA, "z")]))
	print(attributes(L2[c("W", "y", NA, "z"), recurse=TRUE]))
	cat("-----------\n")
	}
	
