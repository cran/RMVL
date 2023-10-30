require("RMVL")

M1<-mvl_open("test1.mvl", append=TRUE, create=TRUE)

df1<-data.frame(x=runif(5000), i=1:5000, n=as.character(1:5000))
mvl_write_object(M1, df1, "df1")

df2<-data.frame(x=runif(100), i=1:100, n=as.character(1:100))
mvl_write_object(M1, df2, "df2")

x<-runif(5000)
mvl_write_object(M1, x, "x")

i<-1:5000
mvl_write_object(M1, i, "i")

n<-as.character(1:5000)
mvl_write_object(M1, n, "n")

P1<-rev(1:5000)
mvl_write_object(M1, P1, "P1")

P2<-rev(1:100)
mvl_write_object(M1, P2, "P2")

M1<-mvl_remap(M1)

mvl_indexed_copy(M1, M1["df1", ref=TRUE], P1, name="df1_p1")
mvl_indexed_copy(M1, M1["x", ref=TRUE], P1, name="x_p1")
mvl_indexed_copy(M1, M1["i", ref=TRUE], P1, name="i_p1")
mvl_indexed_copy(M1, M1["n", ref=TRUE], P1, name="n_p1")

mvl_indexed_copy(M1, M1["df1", ref=TRUE], M1["P1", ref=TRUE], name="df1_mp1")
mvl_indexed_copy(M1, M1["x", ref=TRUE], M1["P1", ref=TRUE], name="x_mp1")
mvl_indexed_copy(M1, M1["i", ref=TRUE], M1["P1", ref=TRUE], name="i_mp1")
mvl_indexed_copy(M1, M1["n", ref=TRUE], M1["P1", ref=TRUE], name="n_mp1")

mvl_indexed_copy(M1, M1["df2", ref=TRUE], P2, name="df2_p2")
mvl_indexed_copy(M1, M1["df2", ref=TRUE], M1["P2", ref=TRUE], name="df2_mp2")

mvl_indexed_copy(M1, df2, P2, name="df2r_p2")
mvl_indexed_copy(M1, df2, M1["P2", ref=TRUE], name="df2r_mp2")

mvl_indexed_copy(M1, df1, P1, name="df1r_p1")
mvl_indexed_copy(M1, df1, M1["P1", ref=TRUE], name="df1r_mp1")

M1<-mvl_remap(M1)

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

if(!compare_df(df1[P1,,drop=FALSE], mvl2R(M1$df1_p1))) {
	cat("test1a failed\n")
	print(attributes(df1))
	print(attributes(mvl2R(M1$df1_p1)))
	cat("----------------\n")
	}

if(!compare_df(df1[P1,,drop=FALSE], mvl2R(M1$df1_mp1))) {
	cat("test1b failed\n")
	print(attributes(df1))
	print(attributes(mvl2R(M1$df1_mp1)))
	cat("----------------\n")
	}
	
if(!isTRUE(all.equal(x[P1], mvl2R(M1$x_p1)))) {
	cat("test2a failed\n")
	print(all.equal(x[P1], mvl2R(M1$x_p1)))
	print(attributes(x[P1]))
	print(attributes(mvl2R(M1$x_p1)))
	cat("----------------\n")
	}

if(!isTRUE(all.equal(x[P1], mvl2R(M1$x_mp1)))) {
	cat("test2b failed\n")
	print(all.equal(x[P1], mvl2R(M1$x_mp1)))
	print(attributes(x[P1]))
	print(attributes(mvl2R(M1$x_mp1)))
	cat("----------------\n")
	}
	
if(!isTRUE(all.equal(i[P1], mvl2R(M1$i_p1)))) {
	cat("test3a failed\n")
	print(all.equal(i[P1], mvl2R(M1$i_p1)))
	print(attributes(i[P1]))
	print(attributes(mvl2R(M1$i_p1)))
	cat("----------------\n")
	}

if(!isTRUE(all.equal(i[P1], mvl2R(M1$i_mp1)))) {
	cat("test3b failed\n")
	print(all.equal(i[P1], mvl2R(M1$i_mp1)))
	print(attributes(i[P1]))
	print(attributes(mvl2R(M1$i_mp1)))
	cat("----------------\n")
	}
	
if(!isTRUE(all.equal(n[P1], mvl2R(M1$n_p1)))) {
	cat("test4a failed\n")
	print(all.equal(n[P1], mvl2R(M1$n_p1)))
	print(attributes(n[P1]))
	print(attributes(mvl2R(M1$n_p1)))
	cat("----------------\n")
	}

if(!isTRUE(all.equal(n[P1], mvl2R(M1$n_mp1)))) {
	cat("test4b failed\n")
	print(all.equal(n[P1], mvl2R(M1$n_mp1)))
	print(attributes(n[P1]))
	print(attributes(mvl2R(M1$n_mp1)))
	cat("----------------\n")
	}
	
if(!compare_df(df2[P2,,drop=FALSE], mvl2R(M1$df2_p2))) {
	cat("test5a failed\n")
	print(attributes(df2))
	print(attributes(mvl2R(M1$df2_p2)))
	cat("----------------\n")
	}

if(!compare_df(df2[P2,,drop=FALSE], mvl2R(M1$df2_mp2))) {
	cat("test5b failed\n")
	print(attributes(df2))
	print(attributes(mvl2R(M1$df2_mp2)))
	cat("----------------\n")
	}

if(!compare_df(df2[P2,,drop=FALSE], mvl2R(M1$df2r_p2))) {
	cat("test6a failed\n")
	print(attributes(df2))
	print(attributes(mvl2R(M1$df2r_p2)))
	cat("----------------\n")
	}

if(!compare_df(df2[P2,,drop=FALSE], mvl2R(M1$df2r_mp2))) {
	cat("test6b failed\n")
	print(attributes(df2))
	print(attributes(mvl2R(M1$df2r_mp2)))
	cat("----------------\n")
	}

if(!compare_df(df1[P1,,drop=FALSE], mvl2R(M1$df1r_p1))) {
	cat("test7a failed\n")
	print(attributes(df1))
	print(attributes(mvl2R(M1$df1r_p1)))
	cat("----------------\n")
	}

if(!compare_df(df1[P1,,drop=FALSE], mvl2R(M1$df1r_mp1))) {
	cat("test7b failed\n")
	print(attributes(df1))
	print(attributes(mvl2R(M1$df1r_mp1)))
	cat("----------------\n")
	}
	
unlink("test1.mvl")

