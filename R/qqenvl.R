qqenvl <-
function(x,datax=FALSE,conf=95,k=10000){
   x=x[complete.cases(x)];
   mx=mean(x);
   sdx=sd(x);
   n=length(x);
   alpha=1-conf/100;
   if(conf!=95){
      cat("Approximating value of p for ",conf,"% confidence...\n\n",sep="");
      p=envl.pval(n,1-conf/100,k);
   } else {
      p=ifelse(n<=75,.8025685*n^(-.7546436),.3675769*n^(.573086));
   }
   p_i=qnorm(ppoints(n));
   p_l=mx+sdx*qnorm(qbeta(p/2,(1:n),n-(1:n)+1));
   p_h=mx+sdx*qnorm(qbeta(1-p/2,(1:n),n-(1:n)+1));
   if(datax==FALSE){
      lines(p_i,p_l,'l',lty=2);
      lines(p_i,p_h,'l',lty=2);
   } else {
      lines(p_l,p_i,'l',lty=2);
      lines(p_h,p_i,'l',lty=2);
   }
}

