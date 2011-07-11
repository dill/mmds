reparam.pi<-function(thetas){
   # re-parameterise the pis 

   # function to use for the transform
   F<-function(x) pgamma(x,3,2) # gamma cdf
   #F<-function(x){exp(x)/(1+exp(x))} # logistic

   # so we don't have to constrain to be >0
   ethetas<-exp(thetas)

   mix.props<-rep(0,length(thetas))

   for(i in 1:length(ethetas)){

      lastbit<-F(sum(ethetas[1:(i-1)]))
      if(i==1) lastbit<-0

      mix.props[i]<-F(sum(ethetas[1:i]))-lastbit
   }

   # return the pis
   return(c(mix.props,1-sum(mix.props)))
}
