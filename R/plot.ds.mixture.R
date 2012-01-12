#plot.ds.mixture<-function(x,...){
#plot.ds.mixture<-function(fit.object,style="",main="",breaks="Sturges",ylim=NULL,xlim=NULL,pdf=FALSE,plot.formula=NULL,hide.hist=FALSE,nomf=FALSE,x.axis=NULL){

plot.ds.mixture<-function(x,style="",main="",breaks="Sturges",ylim=NULL,xlim=NULL,pdf=FALSE,plot.formula=NULL,hide.hist=FALSE,nomf=FALSE,x.axis=NULL,...){
   fit.object<-x

   # todo:
   #  * level legend

   #  fit.object
   #  style
   #  main
   #  breaks="Sturges"
   #  ylim=NULL
   #  xlim=NULL
   #  pdf=FALSE
   #  plot.formula=NULL
   #  hide.hist=FALSE, hide the histogram bars
   #  nomf=FALSE - don't alter mfrow, useful in big plots
   #  x.axis=NULL specify x axis points  

   width<-fit.object$width
   data<-fit.object$data
   pars<-fit.object$pars
   mix.terms<-fit.object$mix.terms
   zdim<-fit.object$zdim
   z<-fit.object$z
   pt<-fit.object$pt

   # plot resolution
   plot.res<-1000

   # set up the parameters, shortcuts
   swpars<-switchpars(pars,mix.terms,z=z,zdim=zdim)
   gp<-getpars(swpars$fpar,mix.terms,swpars$zdim,swpars$z)
   sigma<-gp$key.scale
   pis<-gp$mix.prop

   # use the correct function for effective strip width
   intfcn<-integrate.hn
   if(pt)
      intfcn<-integrate.hn.pt

   # don't plot the histogram if asked
   if(hide.hist){
      hist.col<-"white"
   }else{
      hist.col<-"black"
   }

   if(is.null(z)){
   ##### no covariates

      # calculate mu
      mu<-c()
      for(j in 1:mix.terms)
         mu<-c(mu,intfcn(sigma[j],width))
      
      # create a sequence to evaluate the detection function, and evaluate
      x.seq<-seq(0,width,len=plot.res)
      plotvals<-matrix(0,mix.terms,plot.res)
      for(j in 1:mix.terms)
         plotvals[j,]<-keyfct.hn(x.seq,sigma[j])

      # what is the value of the detection function at 0 distance? 
      p.at.zero<-c()
      for(j in 1:mix.terms)
         p.at.zero<-c(p.at.zero,pis[j]*keyfct.hn(0,sigma[j]))

      p.at.zero<-sum(p.at.zero)

      # make the histogram object
      a<-hist(data$distance,plot=FALSE,breaks=breaks)
      ylabel<-"Probability of detection"

      # point transect data
      if(pt){
         # do pdf plotting
         if(pdf){
            plotvals<-matrix(0,mix.terms,plot.res)
            for(j in 1:mix.terms)
               plotvals[j,]<-2*pi*pis[j]*x.seq*keyfct.hn(x.seq,sigma[j])/
                                         intfcn(sigma[j],width)
   
            a$density<-a$density*(1/sum(a$density*diff(a$breaks)))
            # this isn't really p at zero, it's the max on the y axis
            p.at.zero<-max(plotvals,a$density)
      
            ylabel<-"PDF of detected distances"
         # not pdf plotting
         }else{
            # rescale hist by distance
            rsf<-a$mids
            musum<-sum(mu*pis)/(2*pi)
            a$density<-a$density*(musum/sum(a$density*diff(a$breaks)))
            a$density<-a$density/rsf
         }
      }else{
         musum<-sum(mu*pis)
         hist.area<-sum(a$density*diff(a$breaks))
         a$density<-a$density*(musum/hist.area)
      }      

      par(las=1)
      # actually do the plotting
      plot(a,freq=FALSE,ylab=ylabel,axes=FALSE,
           ylim=c(0,max(a$density,p.at.zero)),xlim=c(0,width),
           xlab="Distance",border=hist.col,
           main=main)
      if(!is.null(x.axis)){
         axis(1,x.axis)
      }else{
         axis(1)
      }
      axis(2,at=c(0,p.at.zero/2,p.at.zero),labels=c(0,0.5,1))
      if(mix.terms>1){
         lines(x.seq,pis%*%plotvals)
      }else{
         lines(x.seq,plotvals)
      }

      if(style == "comp" & mix.terms>1){
         for(i in 1:mix.terms){
            lines(x.seq,plotvals[i,],lty=2)
         }
      }

      box()

   }else{
   ##### covariate code

      ylabel<-"Probability of detection"

      # calculate mu per observation
      mu<-apply(sigma,1,intfcn,width)
      mus<-mu%*%matrix(pis,length(pis),1)

      ##### evaluate the detection function
      # create a sequence of xs to evaluate the detection function at
      x.seq<-seq(0,width,len=1000)
      # quick detection function function
      eval.detfct<-function(x,sigma,pis){
         matrix(pis,1,length(pis))%*%apply(as.matrix(x),1,keyfct.hn,key.scale=sigma)
      }
      # evaluate the detection function at x.seq for each sigma
      # 1000 x # samples matrix results
      plotvals.mat<-apply(sigma,2,eval.detfct,pis=pis,x=x.seq)

      if(pt & pdf){
         plotvals.mat<-2*pi*x.seq*plotvals.mat/as.vector(mus)
      }

      # add in a weighting
      pas<-mus/width
      ws<-(1/pas)/sum(1/pas)
      plotvals.mat.w<-t(t(plotvals.mat)*c(ws))
      # average according to the weighting
      plotvals<-rowSums(plotvals.mat.w)

      ##### what is the value of the detection function at 0 distance? 
      # evaluate the detection function at 0 distance
      paz.mat<-matrix(pis,1,length(pis))%*%apply(sigma,2,keyfct.hn,distance=0)

      if(pt & pdf){
         paz.mat<-paz.mat/as.vector(mus)
      }

      # put in the weights
      paz.mat<-t(t(paz.mat)*c(ws))
      # average
      p.at.zero<-rowSums(paz.mat)

      ##### calculate the overall mu (ie. area under the detection function)
      mu<-sum(c(mus)*c(ws))

      ##### make the histogram object
      a<-hist(data$distance,plot=FALSE,breaks=breaks)

      if(pt){
         # do pdf plotting
         if(pdf){
            a$density<-a$density*(1/sum(a$density*diff(a$breaks)))
            # this isn't really p at zero, it's the max on the y axis
            p.at.zero<-max(plotvals,a$density)
      
            ylabel<-"PDF of detected distances"

            style<-""
         # not pdf plotting
         }else{
            # rescale hist by distance
            rsf<-a$mids
            musum<-sum(mu*pis)/(2*pi)
            a$density<-a$density*(musum/sum(a$density*diff(a$breaks)))
            a$density<-a$density/rsf
         }
      }else{
         # re-weight the heights to the area under histogram == area under curve
         a$density<-a$density*(mu/sum(a$density*diff(a$breaks)))
      }


      # now store plotvals into a matrix
      # the covariates are cbinded to this...
      pv<-matrix(plotvals,length(plotvals),1)
      colnames(pv)[ncol(pv)]<-"Average detection function"

      plot.names<-"Average detection function"

      # store the plot sequence
      # vector of number of "levels" for each covar so
      # the graphics line up...
      plot.seq<-c(1)
      
      if(style=="comp"){
         plotvals.comp<-c()
         for(j in 1:mix.terms){
            curr<-apply(matrix(sigma[j,],1,length(sigma[j,])),2,keyfct.hn,distance=x.seq)
            
            # add in a weighting
            pas<-mus/width
            ws<-(1/pas)/sum(1/pas)
            plotvals.mat.w<-t(t(curr)*c(ws))
            # average according to the weighting
            plotvals.comp<-rbind(plotvals.comp,rowSums(plotvals.mat.w))

         }
         plot.seq[1]<-plot.seq[1]+mix.terms
         pv<-cbind(pv,t(plotvals.comp))
      }

      ##### plot covariate levels

      # default to all
      if(is.null(plot.formula))
         plot.formula<-fit.object$model.formula
   
      # do nothing if the user asks
      if((as.character(plot.formula) != "none") & !(pt&pdf)){

         # make sure it's a formula
         plot.formula<-as.formula(plot.formula)


         # work out levels etc...
         for(curr.term in attr(terms.formula(plot.formula),"term.labels")){
   
            # if we have a factor covariate
            if(grepl("^as.factor",curr.term)){
               # indicator for column that we want
               ind<-colnames(data)==gsub("^as.factor\\((\\w+)\\)$","\\1",curr.term)
               # pull that column out
               this.col<-data[,ind]
               # find its levels
               this.col.levels<-levels(as.factor(this.col))

               for(tl in this.col.levels){
                  ind0<-this.col==tl

                  fac0<-sum(a$density*diff(a$breaks))/mean(mus[ind0])
fac0<-1

                  plotvals.mat0<-plotvals.mat[,ind0]
                  pas0<-mus[ind0]/width
                  ws0<-(1/pas0)/sum(1/pas0)
ws0<-1
                  plotvals.mat0<-t(t(plotvals.mat0)*c(ws0))
                  plotvals0<-rowSums(plotvals.mat0)*fac0
                  pv<-cbind(pv,plotvals0/plotvals0[1])
                  colnames(pv)[ncol(pv)]<-paste(curr.term,"==",tl,sep="")
         
               }
               plot.names<-c(plot.names,
                             paste("Levels of ",curr.term,sep=""))

            # continuous covariates
            # estimate the 0.25, 0.5 and 0.75 quantiles
            }else{
               ind<-colnames(data)==curr.term
               # pull that column out
               this.col<-data[,ind]
               # find its "levels" - the 0.25, 0.5, 0.75 quantiles
               this.col.levels<-quantile(this.col,c(0.25,0.5,0.75),names=FALSE)

               for(tl in this.col.levels){
                  ind0<-this.col<=tl

                  fac0<-sum(a$density*diff(a$breaks))/mean(mus[ind0])
fac0<-1
                  plotvals.mat0<-plotvals.mat[,ind0]
                  pas0<-mus[ind0]/width
                  ws0<-(1/pas0)/sum(1/pas0)
ws0<-1
                  plotvals.mat0<-t(t(plotvals.mat0)*c(ws0))
                  plotvals0<-rowSums(plotvals.mat0)*fac0
                  pv<-cbind(pv,plotvals0/plotvals0[1])
                  colnames(pv)[ncol(pv)]<-paste(curr.term," ",tl," quantile",sep="")
               }
               plot.names<-c(plot.names,
                             paste(curr.term," 0.25/0.5/0.75 quantiles",sep=""))
            }

            plot.seq<-c(plot.seq,length(this.col.levels))
         }
         # how big does mfrow have to be?
         if(length(plot.seq)>4){
            mfrows<-c(ceiling(length(plot.seq)/4),4)
         }else{
            mfrows<-c(1,length(plot.seq))
         }
      }else{
         mfrows<-c(1,1)
      }

      # the user can specify the names if they want...
      plot.names[1:length(main)]<-main

      # save before changing!
      if(!nomf){
         op<-par()
         par(mfrow=mfrows)
      }
      par(las=1)
      
      k<-1
      for(i in 1:length(plot.seq)){


         # plot the histogram
         plot(a,freq=FALSE,ylab=ylabel,axes=FALSE,
              ylim=c(0,max(a$density,p.at.zero)),xlab="Distance",
              main=plot.names[i],xlim=c(0,width),las=1,border=hist.col)

         if(!is.null(x.axis)){
            axis(1,x.axis)
         }else{
            axis(1)
         }
         
         # get the labels right
         if(pt){
            axis(2,at=c(0,p.at.zero/2,p.at.zero),labels=round(c(0,p.at.zero/2,p.at.zero),3))
         }else{
            axis(2,at=c(0,p.at.zero/2,p.at.zero),labels=c(0,0.5,1))
         }
         box()

         if(i==1){
            line.style<-c(1,rep(2,mix.terms))
            line.colours<-rep(1,plot.seq[i])
         }else{
            line.style<-rep(1,plot.seq[i])
            line.colours<-grey(seq(0,0.5,len=plot.seq[i]))
         }

         for(j in 1:plot.seq[i]){
            lines(x.seq,pv[,k],col=line.colours[j],lty=line.style[j])
            k<-k+1
         }

      }

      # restore pars
      if(!nomf){
         par(op)
      }

   }
}
