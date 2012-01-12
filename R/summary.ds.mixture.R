summary.ds.mixture<-function(object,...){
# summary function for ds.mixture model object
#
# Arguments:
#  model  a ds.mixture model object
#
# Value: summary object
#
# See print.summary.ds for how this gets printed

   model<-object

   ans <- list()
   # Number of observations
   ans$n <- length(model$distance)

   # number of mixture components
   ans$mix.terms<-model$mix.terms

   ans$coeff<-list()

   # Parameter estimates and their standard errors
   ans$coeff$pars<-data.frame(estimate=model$pars,se=model$pars.se)

   # separate mixture proportions
   gp<-getpars(model$pars,model$mix.terms,model$zdim,model$z)
   ans$coeff$mix.prop<-gp$mix.prop
   names(ans$coeff$mix.prop)<-paste("pi_",1:model$mix.terms,sep="")

   # AIC
   ans$aic <- model$aic

   # Truncation distance
   ans$width <- model$width
  
   # point or line transects?
   if(model$pt){
      ans$ttype<-"point"
   }else{
      ans$ttype<-"line"
   }

   # average p
   ans$average.p<-model$pa
   ans$average.p.se<-model$pa.se
   ans$average.p.cv<-model$pa.se/model$pa
   # abundance
   ans$Nhat <- model$N
   ans$Nhat.se <- model$N.se
   ans$Nhat.cv <- model$N.se/model$N

   # set the class and return
   class(ans) <- "summary.ds.mixture"
   return(ans)
}
