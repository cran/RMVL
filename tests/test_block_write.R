require("RMVL")

mvl_block_write<-function(MVLHANDLE, N, f, name=NULL, block=1000, type=NULL, ...) {
	obj<-NULL
	offset<-NULL
	for(i in seq(1.0, N, by=block)) {
		idx<-i-1+(1:block)
		idx<-idx[idx<=N]
		
		v<-f(idx, ...)
		
		if(is.null(obj)) {
			if(!is.null(type))attr(v, "MVL_TYPE")<-type
			offset<-mvl_start_write_vector(MVLHANDLE, v, expected.length=N, name=name)
			MVLHANDLE<-mvl_remap(MVLHANDLE)
			obj<-MVLHANDLE[offset,ref=TRUE]
			} else {
			mvl_rewrite_vector(obj, i, v)
			}
		}
	return(invisible(offset))
	}

test_write<-function(MVLHANDLE, vec, name, block=1000) {
	f<-function(idx) { return(vec[idx]) }
	
	mvl_block_write(M, length(vec), f, name=name, block=block)
	
	MVLHANDLE<-mvl_remap(MVLHANDLE)
	
	if(any(vec!=mvl2R(MVLHANDLE[name]))) {
		cat("Block write test", name, "failed\n")
		print(all.equal(vec, mvl2R(MVLHANDLE[name])))
		}
	}
	
M<-mvl_open("test1.mvl", append=TRUE, create=TRUE)

test_write(M, rnorm(1e6), "double1")
test_write(M, as.integer(1:1e6), "int1")
test_write(M, as.integer(1:1e6) %% 2==0, "logical1")
test_write(M, as.Date("2001-01-01")+as.integer(1:1e6), "date1")

mvl_close(M)
unlink("test1.mvl")



