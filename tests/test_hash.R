require("RMVL")

N<-10000

df<-data.frame(x=(1:N) %% 5, y=((1:N) %% 6)+0.5, ab=paste("ab", ((1:N) %% 10)))

Mtmp<-mvl_open("tmp.mvl", append=TRUE, create=TRUE)
mvl_write_object(Mtmp, df, "df")
mvl_close(Mtmp)

Mtmp<-mvl_open("tmp.mvl")

h<-mvl_hash_vectors(list(Mtmp$df[,"x",ref=TRUE]))
LL<-split(h, Mtmp$df[,"x"][])
LL0<-unlist(lapply(LL, function(x){if(length(x)<1)return(x);return(max(abs(diff(x))))}))
if(any(LL0!=0))cat("test1 failed\n")

h<-mvl_hash_vectors(list(Mtmp$df[,"y",ref=TRUE]))
LL<-split(h, Mtmp$df[,"y"][])
LL0<-unlist(lapply(LL, function(x){if(length(x)<1)return(x);return(max(abs(diff(x))))}))
if(any(LL0!=0))cat("test2 failed\n")

h<-mvl_hash_vectors(list(Mtmp$df[,"ab",ref=TRUE]))
LL<-split(h, Mtmp$df[,"ab"][])
LL0<-unlist(lapply(LL, function(x){if(length(x)<1)return(x);return(max(abs(diff(x))))}))
if(any(LL0!=0))cat("test3 failed\n")

h<-mvl_hash_vectors(list(Mtmp$df[,"x",ref=TRUE], Mtmp$df[,"y",ref=TRUE], Mtmp$df[,"ab",ref=TRUE]))
LL<-split(h, Mtmp$df[,])
LL0<-unlist(lapply(LL, function(x){if(length(x)<1)return(x);return(max(abs(diff(x))))}))
if(any(LL0!=0))cat("test4 failed\n")

mvl_close(Mtmp)
unlink("tmp.mvl")
