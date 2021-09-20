require("RMVL")
 
M3<-mvl_open("test2.mvl", append=TRUE, create=TRUE)
L<-list()
df<-data.frame(x=1:1e5, y=rnorm(1e5), s=rep(c("a", "b"), 5e4))
L[["x"]]<-mvl_write_object(M3, df)
L[["description"]]<-"Example of large data frame"
mvl_write_object(M3, L, "test_object")
mvl_close(M3)

M3<-mvl_open("test2.mvl")
print(names(M3))
L2<-M3["test_object"]

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

if(!compare_df(df[1:20,], L2[["x"]][1:20,]))cat("test1 failed\n")

if(!compare_df(df[(1:N) %% 5==0,], L2[["x"]][(1:N) %% 5==0,]))cat("test2 failed\n")

if(!compare_df(df[(1:N) %% 5==0,], L2[["x"]][(1:N)[(1:N) %% 5==0],]))cat("test3 failed\n")

if(!compare_df(df[(1:N) %% 5==0, c("y", "s")], L2[["x"]][(1:N)[(1:N) %% 5==0], c("y", "s")]))cat("test4 failed\n")

if(!identical(df[(1:N) %% 5==0, c("s")], L2[["x"]][(1:N)[(1:N) %% 5==0], c("s")]))cat("test5 failed\n")

if(!compare_df(df[(1:N) %% 5==0, 2:3], L2[["x"]][(1:N)[(1:N) %% 5==0], 2:3]))cat("test6 failed\n")

print(mvl_object_stats(L2[["x"]])[c("length", "type")])

mvl_close(M3)

unlink("test2.mvl")

