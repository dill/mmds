var.Np<-function(pars,mix.terms,width,z,zdim,pt,n,H.inv,N,pa.vec,pa){
   # calculate the variance of N and the average p

   # most of this is taken from the summary.ds method in mrds
   # by Jeff Laake.

   ret<-list()

   # function to calculate N, to be passed to DeltaMethod()
   Nfct<-function(pars,mix.terms,width,z,zdim,pt,n){
      mu<-mu.calc(pars,mix.terms,width,z,zdim,pt)

      if(length(mu)!=n){
         mu<-rep(mu,n)
      }
      return(width*sum(1/mu))
   }


   ### calculate the variance, se and CV of N in the covered region
   Nhatvar.list<-DeltaMethod(pars,Nfct,H.inv,0.001,mix.terms=mix.terms,
                             width=width,z=z,zdim=zdim,pt=pt,n=n)
   Nhatvar<-Nhatvar.list$variance + sum((1-pa.vec)/pa.vec^2)
   cvN<-sqrt(Nhatvar)/N
   # put them in the results list
   ret$N.var<-Nhatvar
   ret$N.se<-sqrt(Nhatvar)
   ret$N.CV<-cvN



   ### calculate variance, se and CV for average detectability

   # function for calculating pa, for DeltaMethod()
   pafct<-function(pars,mix.terms,width,z,zdim,pt,n){
      mu<-mu.calc(pars,mix.terms,width,z,zdim,pt)

      if(length(mu)!=n){
         mu<-rep(mu,n)
      }
      return((n/width)*1/(sum(1/mu)))
   }
   vc1.list<-DeltaMethod(pars,pafct,H.inv,0.001,mix.terms=mix.terms,
                      width=width,z=z,zdim=zdim,pt=pt,n=n)
   vc1<-vc1.list$variance
   vc2<-sum(pa.vec^2*(1 - pa.vec)/pa.vec^2)
   covar<-sum(pa.vec*(1 - pa.vec)/pa.vec^2)
   var.pbar.list<-list(var=vc1+vc2,partial=vc1.list$partial,covar=covar)

   #var.pbar.list=prob.se(model,avgp,H.inv)

   covar<-t(Nhatvar.list$partial)%*%H.inv%*%var.pbar.list$partial+
                        var.pbar.list$covar
   var.pbar<-pa^2*(cvN^2 + var.pbar.list$var/n^2-2*covar/(n*N))
   # put the results in the return list...
   ret$pa.var<-var.pbar
   ret$pa.se<-sqrt(var.pbar)
   ret$pa.CV<-sqrt(var.pbar)/pa

   return(ret)
}
