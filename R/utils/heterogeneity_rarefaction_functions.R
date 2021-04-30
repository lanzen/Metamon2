### Functions for repeated subsampling to evaluate effect of technical replication 
### at constant read depth

setClass(Class="Resampled",
         representation(
           total="array",
           nonmet="array",
           met="array"
         )
)

setClass(Class="Rarefied",
         representation(
           rarefied="array",
           stdev="array"
         )
)

repeatedSubsamplingPlate1 = function (otus, res, res_p=NULL, res_m=NULL, noSamples, noReps, iterations,reads,
                                classification=NA,PMSamples=4, PCRRepSamples=10,metadata=md.n.18SR) {
  
    for (i in 1:noReps){
    print (i)
    for (j in 1:iterations){
      #print(j)
      for (s in 1:noSamples){
        if (!(s==4 & i>PMSamples) &!(s==5 & i>PCRRepSamples)){ #We stop Vortex_PM (s==4) at 4 since no more reps are available
          #print(grp[s])
          red=rrarefy(otus[metadata$KitHom==grp[s] & (metadata$Group=="ExtractionRep" | metadata$Group=="PCRRepX58"),]
                      ,sample=reads/i)
          # get an i repeated random selection of all reps for sample s
          toSum = sample(dim(red)[1],i)
          if (i>1){
            res[s,i,j,] = colSums(red[toSum,])
            if(!is.null(res_p)){
              res_p[s,i,j,] = colSums(red[toSum,-grep("Metazoa",classification)])
            }
            if(!is.null(res_m)){
              res_m[s,i,j,] = colSums(red[toSum,grep("Metazoa",classification)])
            }
          }
          else{ #i==1
            
            res[s,i,j,] = red[toSum,]
            if(!is.null(res_p)){
              res_p[s,i,j,] = red[toSum,-grep("Metazoa",classification)]
            }
            if(!is.null(res_m)){
                res_m[s,i,j,] = red[toSum,grep("Metazoa",classification)]
            }
          }
        }
      }
    }
    }
  if (is.null(res_p)) res_p = array(dim=0)
  if (is.null(res_m)) res_m = array(dim=0)
  return(new("Resampled", total=res, nonmet=res_p, met=res_m))
}


repeatedSubsampling = function (otus, res, metadata, by="Group", res_p=NA, res_m=NA, noSamples, noReps, iterations,reads,
                                classification=NA) {
  
  grp = unique(metadata[[by]])
  
  for (i in 1:noReps){
    print (i)
    for (j in 1:iterations){
      for (s in 1:noSamples){
        #print(grp[s])
        red=rrarefy(otus[metadata[[by]]==grp[s],] ,
                    sample=reads/i)
        # get an i repeated random selection of all reps for sample s
        toSum = sample(dim(red)[1],i)
        if (i>1){
          res[s,i,j,] = colSums(red[toSum,])
          if(!is.na(res_p)){
            res_p[s,i,j,] = colSums(red[toSum,-grep("Metazoa",classification)])
          }
          if(!is.na(res_m)){
            res_m[s,i,j,] = colSums(red[toSum,grep("Metazoa",classification)])
          }
        }
        else{
          res[s,i,j,] = red[toSum,]
          if(!is.na(res_p)){
            res_p[s,i,j,] = red[toSum,-grep("Metazoa",classification)]
          }
          if(!is.na(res_m)){
            res_m[s,i,j,] = red[toSum,grep("Metazoa",classification)]
          }
        }
      }
    }
  }
  if (is.na(res_p)) res_p=array(dim=c(0,0,0,0))
  if (is.na(res_m)) res_m=array(dim=c(0,0,0,0))
  
  return(new("Resampled", total=res, nonmet=res_p, met=res_m))
}

ppPlot = function(svgFile, resampled, ylims, colour, grp=grp, 
                  noReps=10, noSamples=5, customOrder=c(1:noSamples)){
  # Richness vs reps same read depth calc + graph
  avgS = matrix(nrow=noSamples,ncol=noReps)
  sdS = matrix(nrow=noSamples,ncol=noReps)
  
  for (i in 1:noReps){
    for (s in 1:noSamples){
      avgS[s,i] = mean(specnumber(resampled[s,i,,]))
      sdS[s,i] = sd(specnumber(resampled[s,i,,]))
    }
  }
  
  svg(svgFile,height=6,width=6)
  
  plot(1:noReps,avgS[1,],xlab="Pooled replicates",ylab="OTU richness",
       main=paste(reads,"reads"),type="b",col=colour[1],ylim=ylims,log="x")
  arrows(1:noReps,avgS[1,]-sdS[1,], 1:noReps, avgS[1,]+sdS[1,], 
         length=0.05, angle=90, code=3,col=colour[1])
  
  for (s in c(2:noSamples)) {
    lines(1:noReps,avgS[s,],type="b",col=colour[s],pch=s)
    arrows(1:noReps,avgS[s,]-sdS[s,], 1:noReps, avgS[s,]+sdS[s,], 
           length=0.05, angle=90, code=3,col=colour[s])
  }
  
  legend("topleft",legend=grp[customOrder],col=colour[customOrder],pch=customOrder,box.lwd=0)
  dev.off()

}


ppHPlot = function(svgFile, resampled, ylims, colour, grp=grp, 
                   noReps=10, noSamples=5,customOrder=c(1:noSamples)) {
  # Shannon vs reps same read depth calc + graph
  avgH = matrix(nrow=noSamples,ncol=noReps)
  sdH = matrix(nrow=noSamples,ncol=noReps)
  
  for (i in 1:noReps){
    for (s in 1:noSamples){
      avgH[s,i] = mean(diversity(resampled[s,i,,]))
      sdH[s,i] = sd(diversity(resampled[s,i,,]))
    }
  }
  svg(svgFile,height=6,width=6)
  
  plot(1:noReps,avgH[1,],xlab="Pooled replicates",ylab="Shannon diversity",
       main=paste(reads,"reads"),type="b",col=colour[1],ylim=ylims, log="x")
  arrows(1:noReps,avgH[1,]-sdH[1,], 1:noReps, avgH[1,]+sdH[1,], 
         length=0.05, angle=90, code=3,col=colour[1])
  
  for (s in c(2:noSamples)) {
    lines(1:noReps,avgH[s,],type="b",col=colour[s],pch=s)
    arrows(1:noReps,avgH[s,]-sdH[s,], 1:noReps, avgH[s,]+sdH[s,], 
           length=0.05, angle=90, code=3,col=colour[s])
  }
  
  legend("bottomright",legend=grp[customOrder],col=colour[customOrder],pch=customOrder,box.lwd=0)
  
  dev.off()
}


plotRare = function(rarefied, stdev, file, height, width, ylims=ylims, xlims=xlims,
                    diff=diff, colour=colour){
  
  # Rarefaction - now including PCR rep for clarity
  svg(file=file, height=height, width=width)
  plot(rarefied[,1],rarefied[,2],type="l",ylim=ylims, xlim=xlims,
       xlab="Reads per replicate",ylab="OTU richness",col=colour[1],
       lwd=2)#,log="x")
  arrows(rarefied[,1]-diff, rarefied[,2]-stdev[,2], 
         rarefied[,1]-diff, rarefied[,2]+stdev[,2], 
         length=0.05, angle=90, code=0,col=colour[1]) 
  
  # Make right order of s iterations and take direct from rspooled and rsis!
  
  for (s in c(2:5)) {
    lines(rarefied[,1],rarefied[,s+1],col=colour[s],lwd=2)
    arrows(rarefied[,1]+diff*(s-2), rarefied[,s+1]-stdev[,s+1], 
           rarefied[,1]+diff*(s-2), rarefied[,s+1]+stdev[,s+1], 
           length=0.05, angle=90, code=0,col=colour[s])
  }
  rpCorr = c(1,3,4,2)
  rsCorr = c(1,3,4,2,5)
  for (s in c(1:4)) {
    lines(rarefied[,1],rarefied[,s+6],lty=2,col=colour[rpCorr[s]],lwd=2)
  }
  for (s in c(1:5)) {
    lines(rarefied[,1],rarefied[,s+10],lty=3,col=colour[rsCorr[s]],lwd=2)
  }
  
  legend("topleft",
         legend=c(as.character(grp),"Technical rep.","In vitro pooled","In silico pooled"),
         col=c(colour,rep("black",3)),lwd=2,box.lwd=0,cex=.95,lty=c(rep(1,6),2,3))
  dev.off()
}