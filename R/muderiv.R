# derivatives of mu
muderiv<-function(par,mix.terms,width,z=NULL,zdim=0,pt=FALSE){
   got.pars<-getpars(par,mix.terms,zdim,z)
   key.scale<-got.pars$key.scale
   key.shape<-got.pars$key.shape
   mix.prop<-got.pars$mix.prop

   # storage
   dbeta<-c()

   if(pt){
      intfcn<-integrate.hn.pt
   }else{
      intfcn<-integrate.hn
   }

   ### derivatives with respect to beta first...
   # intercept terms
   for(j in 1:mix.terms){

      if(is.list(z)|all(zdim>0)){
         keysc<-key.scale[j,]
      }else{
         keysc<-key.scale[j]
      }

      # calculate muj
      muj<-intfcn(keysc,width)

      if(pt){
         mid<-mix.prop[j]*(2*pi*(2*muj+2*pi*width^2*exp(-width^2/(2*keysc^2))))
      }else{
         mid<-mix.prop[j]*(muj-width*exp(-(width^2)/(2*keysc^2)))
      }
      
      dbeta<-cbind(dbeta,mid)
   }

   # if we have covariates - taken from flt.gr
   K<-length(par)-mix.terms-(mix.terms-1)
   if(K>0){
      # loop over covariates
      for(k in 1:K){

         n<-nrow(z[[1]])
         sb<-rep(0,n)
         for(j in 1:mix.terms){
            if(is.list(z)|all(zdim>0)){
               if(length(z)==1){
                  zz<-z[[1]]
               }else{
                  zz<-z[[j]]
               }
               keysc<-key.scale[j,]
            }else{
               keysc<-key.scale[j]
               zz<-matrix(1,length(x$distance),1)
            }

            intj<-intfcn(keysc,width)

            # second bit -- different for lt and pt
            if(pt){
               sb<-sb+mix.prop[j]*2*(zz[,k+1]*intj-
                        zz[,k+1]*pi*width^2*keyfct.hn(width,keysc))
            }else{
               sb<-sb+mix.prop[j]*zz[,k+1]*
                      (intj-width*keyfct.hn(width,keysc))
            }
         }
         dbeta<-cbind(dbeta,sb)
      }
   }


   ### then wrt alpha...
   if(mix.terms>1){
      dpi<-rep(0,(mix.terms-1))
      alphas<-par[(length(par)-mix.terms+2):length(par)]

      for(j in 1:(mix.terms-1)){
         midterms<-0
         for(js in 1:(mix.terms-1)){
            if(is.list(z)|all(zdim>0)){
               keysc<-key.scale[js,]
            }else{
               keysc<-key.scale[js]
            }

            # alpha premultiplier
            alphaterms<-0
            if(j<= js){
               esum1<-dgamma(sum(exp(alphas[1:js])),3,2)
               if((js-1)<1 | (js-1)<j){
                  esum2<-0
               }else{
                  esum2<-dgamma(sum(exp(alphas[1:(js-1)])),3,2)
               }
               alphaterms<-esum1-esum2
            }

            midterms<-midterms+alphaterms*intfcn(keysc,width)
         }

         if(is.list(z)|all(zdim>0)){
            keyscJ<-key.scale[mix.terms,]
         }else{
            keyscJ<-key.scale[mix.terms]
         }

         # ESW for the Jth component
         intJ<-intfcn(keyscJ,width)
         # alpha premultiplier
         alphaterms<- -dgamma(sum(exp(alphas[1:(mix.terms-1)])),3,2)
         intJ<-intJ*alphaterms

         #put it all together....   
         dpi[j]<-exp(alphas[j])*sum(midterms+intJ)
      }

      if(is.list(z)|all(zdim>0)){
         dpar<-cbind(dbeta,rep(dpi,n))
      }else{
         dpar<-c(dbeta,dpi)
      }

   }else{
      dpar<-dbeta
   }
   return(dpar)
}
