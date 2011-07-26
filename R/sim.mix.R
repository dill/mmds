# simulate a mixture model detection function
sim.mix<-function(pars,mix.terms,n,width,zdim=0,z=NULL,pt=FALSE,showit=FALSE){

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
      if(U<=(mult*width)*detfct(proposal,pars,mix.terms,
                     zdim=zdim,z=z.obj)){
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
