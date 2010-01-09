envl.plot <-
function(x,conf=95,plot.it=TRUE,mc=FALSE,k=10000){
   x=x[complete.cases(x)];
   dataname=deparse(substitute(x));
   is.between=function(p,l,h){
      if(l<p & p<h)
        return(TRUE);
      return(FALSE);
   }
   isnormal=TRUE;
   n=length(x);
   alpha=1-conf/100;
   if(conf!=95 | mc==TRUE){
      cat("Approximating value of p for ",conf,"% confidence...\n\n",sep="");
      p=envl.pval(n,alpha,k);
   } else {
      p=ifelse(n<=75,.8025685*n^(-.7546436),.3675769*n^(-.573086));
   }
   
   x=sort((x-mean(x))/sd(x));
   p_i=qnorm(ppoints(n))
   
   p_l=qnorm(qbeta(p/2,(1:n),n-(1:n)+1));
   p_h=qnorm(qbeta(1-p/2,(1:n),n-(1:n)+1));
   for(i in 1:n){
      isnormal=isnormal&is.between(x[i],p_l[i],p_h[i]);
   }
   if(plot.it==TRUE){
      plot.new();
      plot.window(xlim=range(p_i),ylim=c(min(p_l,x),max(p_h,x)));
      axis(1); axis(2); box();
      points(p_i,x,pch=20,col='red');
      lines(p_i,p_l,'l',col='blue');
      lines(p_i,p_h,'l',col='blue');
      mt=paste("Normal probability plot with ",conf,"% confidence band",sep="");
      if(isnormal==TRUE){
         title(main=mt,sub=list('Not enough evidence to reject normality',col='blue',font=3),xlab='Theoretical Quantiles',ylab='Sample Quantiles');
      } else {
         title(main=mt,sub=list('Data was not likely sampled from a normal population',col='red',font=3),xlab='Theoretical Quantiles',ylab='Sample Quantiles');
      }
   }
   
   TRES=list(quantiles=p_i,
    obs=x,
    low.band=p_l,
    high.band=p_h,
    method="Envelope test",
    data.name=dataname,
    alternative="Data was not sampled from a Normal Distribution",
    confidence.level=conf,
    n=n,
    criterion="Reject if at least one observed values lies outside of the confidence band",
    outcome=!isnormal);
   class(TRES)="htest";
   invisible(TRES);
}

