# fitmix - Fit a mixture detection function to data.
"fitmix"<-function(data,width,mix.terms=1,pt=FALSE,model.formula="~1",initialvalues=NULL,showit=0,ctrl.options=c(maxit=10000),opt.method="BFGS+SANN",usegrad=TRUE,ftype="hn"){
   # Arguments
	#  data            - data of format as in mrds
	#  ftype           - function type (at the moment only half-normal)
   #  model.formula   - c(~1) by default, covariate formula 
	#  mix.terms       - number of mixture terms in the model
	#  initialvalues   - starting values
	#  width           - trunction width for the data
	#  showit          - debugging level
	#  ctrl.options    - options to pass to optim
	#  opt.method      - optimisation method -- BFGS+SANN, BFGS or EM
   #  usegrad         - use analytic gradients?
   #  pt              - is the data from point transects?

   # truncate the distances to width
   data<-data[data$distance<=width,]
   n<-length(data$distance)
   
   # use analytic gradient?
   if(usegrad==TRUE){
      gradfunc<-flt.gr
   }else{
      gradfunc<-NULL
   }

   # first check if we just have one formula for all parts of the mixture
   if(!is.list(model.formula)){
      # is the formula just ~1? -- mixture or non-mixture but with no covars
      if(model.formula=="~1"){
         zdim<-0
         z<-NULL
      # one formula - keep covars across mixtures, vary intercept
      }else{
         # this just picks out those columns from data that are
         # specified in the formula
         z<-list(model.matrix(as.formula(model.formula),data))
         zdim<-dim(z[[1]])[2]
      } 

   ### WE DON'T THINK THIS WORKS ANY MORE!
   }else{
      # different covars for each part of the mixture
      z<-list()
      zdim<-c()
      i<-1
   
      # create the covariate matrices and record their sizes
      for(modform in model.formula){
         z[[i]]<-model.matrix(as.formula(modform),data)
         zdim<-c(zdim,dim(z[[i]])[2])
         i<-i+1
      }
   }

   # set the initial values for optim
   if(is.null(initialvalues)){
      initialvalues<-setinitialvalues(data,mix.terms,model.formula,z,zdim,
                                      width,showit)
   # check initialvalues
   }else{
      initialvalues<-checkinitialvalues(initialvalues,zdim,mix.terms)
   }

   # Debug
   if(showit>=1){
      cat("Initial values =",initialvalues,"\n")
   }

   # if we have a 1-point mixture, use BFGS+SANN
   if(opt.method=="EM" & mix.terms==1){
      message("EM algorithm selected with 1 mixture: switching to BFGS.\n")
      opt.method<-"BFGS"
   }

   ########################
   # Optimise

   if(opt.method=="BFGS+SANN"){
      # store the models
      run.models<-list()
      aics<-c()
      iter<-1
      nit<-1
 
      # run "many" times (or at least until 1 works or there are 
      # too many its) then find best AIC
      while(length(aics)<5 & nit<20){
  
         # first run SANN for a bit...
         lt1 <- try(optim(initialvalues, flt, method="SANN", 
                   control=c(lmm=100,fnscale=-1,ctrl.options,maxit=500),
                   mix.terms=mix.terms,width=width,showit=showit,x=data,
                   ftype=ftype,hessian=TRUE,z=z,zdim=zdim,pt=pt))
      
      
         # if SANN died or if we can't evaluate the likelihood at 
         # this point then don't do the BFGS step
         if(class(lt1)=="try-error"){
            if(showit>1){
               cat("Try error in SANN step... skipping...\n")
            }
            iter<-(iter-1)
         }else{

            # debug output
            if(showit>1){
               cat("SANN step finished, starting values for BFGS=",lt1$par,"\n")
            }

            # then do BFGS
            lt <- try(optim(lt1$par,flt,gr=gradfunc,method="BFGS",control=c(lmm=1000,
                  fnscale=-1,ctrl.options),
                  mix.terms=mix.terms, width=width, showit=showit, x=data,ftype=ftype, 
                  hessian=TRUE,z=z,zdim=zdim,pt=pt))

            # if nothing bad happened...
            if(class(lt)!="try-error" & all(!is.nan(lt$hessian))){
               if((lt$convergence==0) & (lt$value<1e15)){
                  if(all(diag(lt$hessian)<0)){
                  # store this model
                  run.models[[iter]]<-lt
                  # and it's AIC
                  aics<-c(aics,-(2*lt$value)+2*length(lt$par))
                  iter<-iter+1
                  if(showit>1){
                     cat("BFGS step finished, parameters are=",lt$par,"\n")
                  }
                  }
               }
            }
         }
         nit<-nit+1
      } # end of while
   

      if(length(aics)==0){
         cat(paste(mix.terms,"-point model could not be fitted! Try changing initial values\n"))
         lt<-list(par=NA)
         class(lt)<-"try-error"
      }else{

         # pick the best model
         lt<-run.models[[which.min(aics)]]

         # run for 1 it, with pars in correct order
         # put things in the right order
         fixed.pars<-switchpars(lt$par,mix.terms,zdim,z)
         pars<-fixed.pars$fpar
         z<-fixed.pars$z
         zdim<-fixed.pars$zdim

         # run optim for 1 iteration to calulate the hessian etc needed later
         lt<-try(optim(pars,flt,gr=gradfunc,method="BFGS", 
                  control=c(lmm=1000,fnscale=-1,maxit=1),
                  mix.terms=mix.terms,
                  width=width, showit=showit, x=data,ftype=ftype,
                  hessian=TRUE,z=z,zdim=zdim,pt=pt))

         if(showit>=1){
            cat("BFGS+SANN model AICs=",aics,"\n")
            cat("Best model had AIC=",min(aics),"\n")
         }
      }

   # just use BFGS
   }else if(opt.method=="BFGS"){

      lt <- try(optim(initialvalues,flt,gr=gradfunc,method="BFGS",control=c(lmm=1000,
            fnscale=-1,ctrl.options),
            mix.terms=mix.terms, width=width, showit=showit, x=data,ftype=ftype, 
            hessian=TRUE,z=z,zdim=zdim,pt=pt))

   # use the E-M algorithm
   }else if(opt.method=="EM"){
      lt<-em(data,initialvalues,mix.terms,zdim=zdim,z=z,ftype="hn",width,
                 showit=showit,grad=!is.null(gradfunc),pt=pt)
   }else{
      stop("Not an optimisation method!\n")
   }

   ### At this point the model has been fit.
   # Now to format output and calculate other things

   # When everything goes well
   if(class(lt)!="try-error"){

      # name the parameters
      names(lt$par)<-namepars(lt$par,mix.terms,zdim,z)

      # create the return object
      ret.list<-list(distance=data$distance,likelihood=lt$value,pars=lt$par,
                     mix.terms=mix.terms,width=width,z=z,zdim=zdim,
                     hessian=lt$hessian,pt=pt,data=data,ftype=ftype,
                     ctrl.options=ctrl.options,showit=showit,
                     opt.method=opt.method,usegrad=usegrad,
                     model.formula=model.formula)

      # calculate mu
      ret.list$mu<-mu.calc(lt$par,mix.terms,width,z,zdim,pt)
      names(ret.list$mu)<-NULL

      # calculate P_a
      if(pt){
         ret.list$pa.vec<-ret.list$mu/(pi*width^2)
      }else{
         ret.list$pa.vec<-ret.list$mu/width
      }

      # estimate N and average P_a if necessary
      if(is.null(z)){
         ret.list$N<-n/ret.list$pa.vec
         ret.list$pa<-ret.list$pa.vec
      }else{
         ret.list$N<-sum(1/ret.list$pa.vec)
         ret.list$pa<-n/ret.list$N
      }


      ## invert the Hessian to get the var/covar matrix
      H.inv<-solvecov(-lt$hessian)$inv
      
      # store the standard errors of the parameters

      ### Calculate the standard error of P_a and N
      varnp<-var.Np(lt$par,mix.terms,width,z,zdim,pt,n,
                    H.inv,ret.list$N,ret.list$pa.vec,ret.list$pa,data)
      ret.list$N.se<-varnp$N.se
      ret.list$pa.se<-varnp$pa.se
      ret.list$pars.se<-sqrt(diag(varnp$vcov))

      # Calculate the AIC
      ret.list$aic<- -(2*lt$value)+(2*length(lt$par))

      # K-S test & C-vM test
      gofres<-goftests(x=data,width=width,pars=lt$par,
                       mix.terms=mix.terms,z=z,zdim=zdim,pt=pt)
      ret.list$cvm<-gofres$cvm
      ret.list$ks<-gofres$ks

      # set the class
      class(ret.list)<-c("ds.mixture")

      # return it
      return(ret.list)
   }
   
   # if we did't return then...
   # What to do if optim() fails
   stop("Optimisation has failed or ran out of iterations.\n")
   return(NA)
}
