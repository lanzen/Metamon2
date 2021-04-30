 ###################################################################################
 # Quantile spline regression spline using polynomial models                       # 
 # Follows the approach as closely as possible of Anderson (2008) "Animal-sediment #
 # relationships re-visited: Characterising species' distributions along  an       #
 # environmental gradient using canonical analysis and quantile regression splines"#
 #                                                                                 # 
 # Anders Lanzen Oct 2019                                                          #
 ###################################################################################

require(quantreg)
require(splines)
require(irr)

splinePredict = function(goodTITANInds, t, otus.train, md.train, imageOutDir = NA) {
  
  groups = data.frame(row.names=names(otus.train),value=rep(0,dim(otus.train)[2]))
  
  for (ti in c(1:dim(goodTITANInds)[1])){
    
    sp = row.names(goodTITANInds)[ti]
    titanGrp = goodTITANInds$maxgrp[ti]
    if (t=="Redox" | t=="Redox5yAvg" | t=="Redox10yAvg"| t=="NSI"| t=="ISI") titanGrp = 3-titanGrp
    
    ## Model distribution using 4th degree polynomial splines based on 95% quantiles 
    ## (Andersson et al 2008 and Nigel, except not sure how to check AIC and compare to 3 and 5 df)
    mt = data.frame(ab=otus.train[,sp],t=md.train[,t])
    
    bestAIC = 1E6
    bestModel = NA
    bestDF = NA
    for (df in c(2:5)){
      try = tryCatch({
        bsp <- rq(ab ~ bs(t,degree=df),data=mt,tau=.9)
        #print(AIC(bsp))
        aic = AIC(bsp, k=df)[1]
        if (aic<bestAIC-2){ 
          bestAIC = aic
          bestModel=bsp
          bestDF = df
        }
      },
      error= function(what_condition_my_condition_is_in){
        print(paste("Warning: degree",df,"spline failed for",sp))
        try=NULL
      })
    }
    
    print(paste("Best prediction for",df,"d.f. AIC =",bestAIC))
    
    if (bestAIC < 1E6 ){
      st = seq(min(md.train[,t]),max(md.train[,t]),length.out=1000)
      pred<-predict(bestModel,data.frame(t=st))
      maxT = st[pred==max(pred)]
      ecoGroup = getBIGroupFromValue(maxT,bi=t)
      
      ## Plot distribution and check indicator status manually if needed 
      
      if (!is.na(imageOutDir) & !is.na(ecoGroup)){
        
        change=NA
        if (goodTITANInds$maxgrp[ti]==2) {
          change <- "+"
        } else {
          if (goodTITANInds$maxgrp[ti]==1) change <- "-"
        }
        
        png(paste(imageOutDir,"/",sp,"_",t,".png",sep=""),width=400,height=400)
        
        plot(otus.train[,sp]~md.train[,t],
             main=paste(sp,change,"@",round(goodTITANInds$ienv.cp[ti],1)),
             sub=paste("Ecogroup",ecoGroup,"max @",round(maxT,1)),
             xlab=t,ylab="Abundance")
        abline(v=goodTITANInds$ienv.cp[ti], col="blue")
        abline(v=maxT, col="red", lty=2)
        lines(st,pred,col="orange")
        dev.off()
      }
      
      if ((titanGrp==1 & ecoGroup>3) | (titanGrp==2 & ecoGroup<3)){
        print(paste("Warning: conflict for",sp,"TITAN =",titanGrp,"but predicted EC =",ecoGroup))
        print("Setting EC to zero")
        ecoGroup = 0
      }
    }
    else {
        print(paste("Warning: cannot predict spline for",sp))
        ecoGroup = 0
    } 
    
    groups[sp,"value"] = ecoGroup
    
  }
  return (groups)
}
