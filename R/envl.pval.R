envl.pval <-
function(n, alpha=0.05, k=10000){
   result=.C("alpha_envelope",as.integer(n),as.double(alpha),as.integer(k),pval=double(1),PACKAGE='envelope');
   return(result$pval)
}

