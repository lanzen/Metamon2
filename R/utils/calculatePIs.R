###################################################################################
# Calculates abiotic pressure indices based on Aylagas et al (2017), but without  #
# including redox potential. Limits taken from Basque regional values (Menchaca   #
# 2014)                                                                           #
#                                                                                 # 
# Anders Lanzen 2019-11-12                                                        #
###################################################################################

#TODO get limit values from table allowing different ones such as EU or Basque

calculatePIs = function(metadata){
  
  metadata$pi = 0
  metadata$piMetals = 0
  metadata$piHC = 0
  metadata$THC_PI = NA
  metadata$PAH_PI = NA
  
  baLM = lm(log(md.18S.p$Ba)~md.18S.p$NSI)
  plot(log(md.18S.p$Ba)~md.18S.p$NSI)
  abline(baLM,col="grey")
  summary(baLM,x=20)
  exp(6.14)
  
  thcLM = lm(log(md.18S.p$THC)~md.18S.p$NSI)
  plot(log(md.18S.p$THC)~md.18S.p$NSI)
  abline(thcLM, col="grey")
  summary(thcLM,col="grey")
  predict(thcLM)
  exp(10.96-.334*20)
  # Limit for Ba is done using fitting to NSI and expected value at NSI==20
  # Others «Klassifisering-av-miljotilstand-i-vann-02-2018» p 208--211 AA EQS
  metalLimits = data.frame(row.names=c("Cu", "Ba"), #"Cr","Ni",Zn","Pb","Hg","Cd"
                           limits = c(84,464)) # was 20, 100

  # PAH in ppb, HC in ppm, l
  # imit for PAH is is done using fitting to NSI and expected value at NSI==20
  hcLimits = data.frame(row.names=c("TotalPAH","TotalHC"), #THC based on article from Zeina
                         limits = c(2000,72)) #was 300,538
  
  #omLimit = 2.0
  allLimits = as.data.frame(cbind(t(metalLimits), t(hcLimits)))#, data.frame(OM=omLimit))
  
  #geo = c(0,1,2,5,10,100,1000) <- old series was not geometric and was missing upper limit
  geo = c(0,.5,1,2,4,8,16,32)
  for (s in row.names(metadata)){
    metalData = metadata[s,row.names(metalLimits)[row.names(metalLimits) %in% names(metadata)]]
    hcData = metadata[s,row.names(hcLimits)[row.names(hcLimits) %in% names(metadata)]]
    allPIData = cbind(metalData, hcData)#, metadata[s,"OM"])
    
    mLim = metalLimits[row.names(metalLimits) %in% names(metadata),]
    hcLim = hcLimits[row.names(hcLimits) %in% names(metadata),]
    aLim = allLimits[,names(allLimits) %in% names(metadata)]
    
    hasPI = !is.na(allPIData)
    hasMData = !is.na(metalData)
    hasHCData = !is.na(hcData)
    
    # Get k-values (1/number of parameters)
    ka = 1/sum(hasPI)
    km = 1/sum(hasMData)
    khc = 1/sum(hasHCData)
    # 
    # print(hasHCData)
    # print(khc)
    # 
    # Assign NA if no data
    if (sum(hasPI)==0) metadata[s,"pi"] = NA
    if (sum(hasMData)==0) metadata[s,"piMetals"] = NA
    if (sum(hasHCData)==0) metadata[s,"piHC"] = NA
    
    # Calculate total PI
    for (i in c(1:length(allPIData))[hasPI]){
      v = allPIData[,i]
      l = aLim[,i]
      for (j in c(1:6)){
        #print(paste(v,"<",l*geo[j+1],"?"))
        if (v<l*geo[j+1]){
          metadata[s,"pi"] = metadata[s,"pi"] + ka*(v-l*geo[j])/(l*geo[j+1] - l*geo[j])
          break
        } else {
          metadata[s,"pi"] = metadata[s,"pi"] + ka
        }
      }
    }
      
    # Calculate metal PI
    for (i in c(1:length(metalData))[hasMData]){
      v = metalData[,i]
      l = mLim[i]
      for (j in c(1:6)){
        #print(paste(v,"<",l*geo[j+1],"?"))
        if (v<l*geo[j+1]){
          metadata[s,"piMetals"] = metadata[s,"piMetals"] + km*(v-l*geo[j])/(l*geo[j+1] - l*geo[j])
          break
        } else {
          metadata[s,"piMetals"] = metadata[s,"piMetals"] + km
        }
      }
    }
      
    
    # Calculate HC PI
    for (i in c(1:length(hcData))[hasHCData]){
      v = hcData[,i]
      l = hcLim[i]
      for (j in c(1:6)){
        #print(paste(v,"<",l*geo[j+1],"?"))
        if (v<l*geo[j+1]){
          metadata[s,"piHC"] = metadata[s,"piHC"] + khc*(v-l*geo[j])/(l*geo[j+1] - l*geo[j])
          break
        } else {
          metadata[s,"piHC"] = metadata[s,"piHC"] + khc
        }
      }
    }
    
    
    # Calculate Partial hydrocarbon PIs (TotalPAH and TotalHC)
    if (!is.na(metadata[s,"TotalPAH"])){
      pah = metadata[s,"TotalPAH"]
      j=1
      PAH_PI = 0
      l = hcLimits$limits[1]
      while(pah>l*geo[j] & j<8){
        PAH_PI = PAH_PI + min(1,(pah-l*geo[j])/(l*geo[j+1] - l*geo[j]))  
        j=j+1
      }
      metadata[s,"PAH_PI"] = PAH_PI
    }
    if (!is.na(metadata[s,"TotalHC"])){
      pah = metadata[s,"TotalHC"]
      j=1
      THC_PI = 0
      l = hcLimits$limits[2]
      while(pah>l*geo[j] & j<7){
        THC_PI = THC_PI + min(1,(pah-l*geo[j])/(l*geo[j+1] - l*geo[j]))  
        j=j+1
      }
      metadata[s,"THC_PI"] = THC_PI
    }
    # 
    # # Calculate OM PI
    # if (!is.na(metadata[s,"OM"])){
    #   om = metadata[s,"OM"]
    #   j=1
    #   piOM = 0
    #   while(om>omLimit*geo[j] & j<6){
    #     piOM = piOM + min(1,(om-omLimit*geo[j])/(omLimit*geo[j+1] - omLimit*geo[j]))  
    #     j=j+1
    #   }
    #   metadata[s,"piOM"] = piOM
    # }
  }
  return(metadata)
}
