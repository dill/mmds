inv.reparam.pi<-function(mix.prop,lastpar=FALSE){
   # re-parameterise the pis so that when there is more
   # than a 2 point mixture
   # INVERSE VERSION!
   # pis->thetas

   # functions to use for the transforms
   F<-function(x) pgamma(x,3,2)
   Finv<-function(x) qgamma(x,3,2)
   #F<-function(x){exp(x)/(1+exp(x))}
   #Finv<-function(x){log(x/(1-x))}

   thetas<-c()

   if(lastpar){
      looper<-length(mix.prop)
   }else{
      looper<-(length(mix.prop)-1)
   }

   for(i in 1:looper){
      if(i==1){
         theta<-log(Finv(mix.prop[i]))
      }else{
         theta<-log(
                  Finv(mix.prop[i]+F(sum(exp(thetas))))-
                  sum(exp(thetas)))
      }
      thetas<-c(thetas,theta)
   }

   # return the thetas
   return(thetas)
}
