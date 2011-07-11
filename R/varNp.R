var.Np<-function(pars,mix.terms,width,z,zdim,pt,n,H.inv,N,pa){
   # calculate the variance of N and the average p

   # return list
   ret<-list()

   if(all(zdim==0)){
   # Calculate the standard error of P_a using sandwich estimator
      # derivatives of mu/nu
      muderivs<-muderiv(pars,mix.terms,width,z,zdim,pt)

      pa.se<-(t(muderivs)%*%H.inv%*%muderivs)

      # now split based on point-transect vs. not pt
      if(pt){
         pa.se<-sqrt((1/width^4)*pa.se)
      }else{
         pa.se<-sqrt((1/width^2)*pa.se)
      }

      ret$pa<-pa.se/sqrt(n)

   }else{
   # calculate variance of p_a in the covariate case 

      # using numerical derivatives....
      parlen<-length(pars)
      eps<-rep((.Machine$double.eps^(1/3)),parlen)

      mdx<-rep(0,parlen)
      pdx<-rep(0,parlen)

      # identity matrix for picking out elements
      picker<-diag(parlen)

      for(i in 1:parlen){
         mdx[i]<-1/sum(1/mu.calc(pars-picker[i,]*eps,
                                 mix.terms,width,z,zdim,pt))
         pdx[i]<-1/sum(1/mu.calc(pars+picker[i,]*eps,
                                 mix.terms,width,z,zdim,pt))
      }

      paderivs<-(pdx-mdx)/(2*eps)
################# HACK
      if(!pt){
         mu<-mu.calc(pars,mix.terms,width,z,zdim,pt)
         varpa<-sum((mu/width)^2 - mu/width)

      }else{
         h.at.zero<-2*pi*detfct(0,pars,mix.terms,zdim,z)/
                        mu.calc(pars,mix.terms,width,z,zdim,pt)

         varpa<-(4/width^4)*sum((1/h.at.zero)^2)-pa

      }

      for(j in 1:parlen){
         for(m in 1:parlen){
            varpa<-varpa+paderivs[j]*paderivs[m]*H.inv[j,m]
         }
      }
################# </HACK>
      ret$pa<-varpa
   }

      ret$pa<-0


   ### calculate the standard error of N
   # always do this numerically
   parlen<-length(pars)
   eps<-rep((.Machine$double.eps^(1/3)),parlen)

   mdx<-rep(0,parlen)
   pdx<-rep(0,parlen)

   # identity matrix for picking out elements
   picker<-diag(parlen)


   # covars
   if(all(zdim!=0)){
      mult<-width
      if(pt) mult<-width^2*pi

      for(i in 1:parlen){
         mdx[i]<-sum(mult/mu.calc(pars-picker[i,]*eps,mix.terms,
                               width,z,zdim,pt))
         pdx[i]<-sum(mult/mu.calc(pars+picker[i,]*eps,mix.terms,
                               width,z,zdim,pt))
      }
      Nderivs<-(pdx-mdx)/(2*eps)
      # using the formula on p. 42 of Buckland etal 2004
      if(!pt){
         mu<-mu.calc(pars,mix.terms,width,z,zdim,pt)
         varN<-sum((width/mu)^2 - width/mu)

      }else{
         # NB. _no_ r! see p.44 of ADS
         h.at.zero<-2*pi*detfct(0,pars,mix.terms,zdim,z)/
                        mu.calc(pars,mix.terms,width,z,zdim,pt)

         varN<-(width^4/4)*sum(h.at.zero^2)-N

      }

      for(j in 1:parlen){
         for(m in 1:parlen){
            varN<-varN+Nderivs[j]*Nderivs[m]*H.inv[j,m]
         }
      }
   # no covariates
   }else{
      # just use sandwich estimator
      for(i in 1:parlen){
         mdx[i]<-1/mu.calc(pars-picker[i,]*eps,mix.terms,
                               width,z,zdim,pt)
         pdx[i]<-1/mu.calc(pars+picker[i,]*eps,mix.terms,
                               width,z,zdim,pt)
      }
      Nderivs<-(pdx-mdx)/(2*eps)
      varN<-(n*width)^2*t(Nderivs)%*%H.inv%*%Nderivs

      if(pt){
         varN<-varN*pi^2*width^2
      }
   }

   # put the variance of N in the results list
   ret$N<-varN

   return(ret)

}
