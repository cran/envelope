#include <R.h>
#include <Rmath.h>
#define is_between(x,l,h) (((l<x)&(x<h))?true:false)

extern "C" {
double mean(double *&sample,int n){
   double x=0;
   for(int i=0; i<n; i++)
      x=x+sample[i];
   return x/n;
}
}
extern "C"{
double sd(double *&sample, int n){
   double x=0;
   for(int i=0; i<n; i++)
      x=x+pow(sample[i],2);
   return sqrt((x-n*pow(mean(sample,n),2))/(n-1));
}
}
extern "C"{
void normalize(double *&sample, int n){
   double smean=mean(sample,n);
   double ssd=sd(sample,n);
   for(int i=0; i<n; i++)
      sample[i]=(sample[i]-smean)/ssd;
}
}
extern "C"{
   void alpha_envelope(int *n, double *alpha, int *k, double *pval){
   double error=1;
   double pl=0;
   double ph=1;
   int tot=0;
   double psi=0;
   bool interval_check=true;
   double *low_limit=new double[*n];
   double *high_limit=new double[*n];
   double *sample=new double[*n];
   
   while(error>0.0001){
      *pval=(ph+pl)/2;
      for(int i=0; i<*n; i++){
	 low_limit[i]=qnorm(qbeta(*pval/2,i+1,*n-i,1,0),0,1,1,0);
	 high_limit[i]=qnorm(qbeta(1-*pval/2,i+1,*n-i,1,0),0,1,1,0);
      }
      for(int i=0; i<*k;i++){
	 //get random sample of size n, sort it and normalize it
	 GetRNGstate();
	 for(int i=0; i<*n; i++)
	    sample[i]=rnorm(0,1);
	 PutRNGstate();
	 //Rprintf("%f\n",sample[1]);
	 //sort(sample,0,*n-1);
	 R_rsort(sample,*n);
	 normalize(sample,*n);
	 for(int i=0; i<*n; i++){
	    interval_check=interval_check && is_between(sample[i],low_limit[i],high_limit[i]);
	    if(interval_check==false)
	       i=*n;
	 }
	 if(interval_check==true)
	    tot++;
	 interval_check=true;
      }
      psi=(double) tot/ *k;
      if(psi>1-*alpha)
	 pl=*pval;
      else
	 ph=*pval;
      tot=0;
      error=ph-pl;
      //Rprintf("%f\n",*pval);
   }
   
   delete [] low_limit;
   delete [] high_limit;
   delete [] sample;
}
}
