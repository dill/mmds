#' Simulate data from a mixture model detection function
#'
#' Simulate data from a (line or point transect) mixture model detection 
#' function with or without covariates using rejection sampling.
#'
#' This routine uses rejection sampling, so may be rather slow of large sample 
#' sizes. Direct sampling will be available soon.
#'
#' @param pars Parameters of the model to fit. See \code{\link{mmds.pars}} 
#'    for details.
#' @param mix.terms Number of mixture components.
#' @param n Number of data to generate.
#' @param width Truncation distance.
#' @param zdim Number of columns of \code{z}. Defaults to 0.
#' @param z Covariate data. Defaults to NULL. See details for more information.
#' @param pt Should point transect data be generated? Defaults to FALSE.
#' @param showit Print the acceptance rate. Defaults to FALSE.
#' @return a \code{data.frame} with the following columns:
#'    \tabular{ll}{
#'    observed \tab Whether the object was observed, always \code{n} 1s. Kept for \code{mmds} compatability.\cr
#'    object \tab Object identifier, numbered 1 to \code{n}. Kept for \code{mmds} compatability.\cr
#'    distance \tab Observed distances.\cr
#'   \tab Then follows as many columns as there are columns as \code{z}, named as in \code{z}.}
#'
#' @export
#'
#' @author David L. Miller, Len Thomas
#' @examples
#' library(mmds)
#' set.seed(0)
#' ## simulate some line transect data from a 2 point mixture
#' sim.dat<-sim.mix(c(-0.223,-1.897,inv.reparam.pi(0.3)),2,100,1)
#' hist(sim.dat$distance)
sim.mix<-function(pars,mix.terms,n,width,zdim=0,z=NULL,pt=FALSE,showit=FALSE,key="hn",hr.shape=5){

  # need to integrate these into the rest of the code?
  detfct.tmp<-function(x,pars,mix.terms,zdim=0,z=NULL,key="hn",hr.shape=NULL){
    # evaluate the total detection function for either hn or hr keys - no error
    # checking on inputs, assumes that's been done earlier
  
    gp<-getpars(pars,mix.terms,zdim,z)
    sigmas<-gp$key.scale
    pis<-gp$mix.prop
  
    res<-rep(0,length(x))
  
    if(key=="hr"){
      if(length(hr.shape)==1){
        hr.shape <- rep(hr.shape,mix.terms)
      }
    }

    for(j in 1:mix.terms){
      if(is.null(z)){
         keysc<-sigmas[j]
      }else{
         keysc<-sigmas[j,]
      }
      res<- res+pis[j] * switch(key,
        "hn"=keyfct.hn(x,keysc),
        "hr"=keyfct.hr(x,keysc,hr.shape[j]),
        0)
    }
    return(res)
  }
  
  keyfct.hr<-function(x,keysc,hr.shape){
    return(1-exp(-(x/keysc)^-hr.shape))
  }

# hist(sim.dat<-sim.mix(c(log(0.1),log(0.8),inv.reparam.pi(0.5)),2,100000,1,pt=F,key="hr",hr.shape=5)$distance)


   #check input on "key" parameter
   if(key=="hr"){
     #Check no covarriates and no pt
     if(pt) stop("point transects not supported for hazard rate mixtures")
     if(zdim>0) stop("covariates not supported for hazard rate mixtures")
   
   } else {
     if (key!="hn") stop ("argument to 'key' parameter not recognized")
   }

   out<-rep(0,n)
   counter<-0
   total.sim<-0

   mult<-1

   while(counter<n){

      # for covariate models, set z for each observation
      if(!is.null(z)){
         z.obj<-list(matrix(z[[1]][counter+1,],nrow=1))
      }else{
         z.obj<-NULL
      }

      total.sim<-total.sim+1
      proposal<-runif(1,0,width)
      U<-runif(1)

      if(pt){
         # for point transects we don't know what the maximum of f(x) is
         # so, calculate the max(f(x)) using optimize()

         # dummy function with x first
         eval.pdf2<-function(x,fpar,width,mix.terms,showit=0,ftype="hn",
                              z=NULL,zdim=0,pt=TRUE){
            sum(eval.pdf(fpar,data.frame(distance=x),width,mix.terms,
                           showit,ftype,z,zdim,pt))
         }
               
         M<-optimize(eval.pdf2,interval=c(0,width),maximum=TRUE,fpar=pars,
                     mix.terms=mix.terms,pt=pt,z=z.obj,zdim=zdim,
                     width=width)$objective

         nu<-mu.calc(pars,mix.terms,width,z.obj,zdim,pt)

         mult<-(2*pi*proposal/nu)/M
      }

      # accept/reject
      if(U<=(mult*width)*detfct.tmp(proposal,pars,mix.terms,
                     zdim=zdim,z=z.obj,key=key,hr.shape=hr.shape)){
         counter<-counter+1
         out[counter]<-proposal
      }

   }

   out<-data.frame(observed=rep(1,n),object=1:n,distance=out)

   if(!is.null(z)){
      out<-cbind(out,z[[1]])
   }

   if(showit){
      cat("Acceptance rate=",counter/total.sim,"\n")
   }

   return(out)
}


