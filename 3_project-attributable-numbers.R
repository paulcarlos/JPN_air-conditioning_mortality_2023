library(lubridate)
library(dlnm)
library(MASS)

# knots used for natural cubic B-spline
vk <- c(0.5,0.9)

###################### 7 Kansai Cities ###########################

# upload data list
dlist1 <- readRDS("data/cities/split_time/cities_52_japan_1976-2009.rds")
loc <- c("Osaka","Sakai","Kobe","Kyoto","Nara","Otsu","Wakayama")
dlist2 <- dlist1[loc]
dlist3 <- lapply(dlist2,function(z)do.call(rbind,z))
mort <- lapply(dlist3,function(z)z[grep(paste0(2000:2009,"-08-",collapse="|"),z$date),c("date","mort")])
rm(dlist1,dlist2,dlist3)

# second-stage outputs
load("cities52_2nd-stage_predicted-AC_coef_vcov.rda")

# temperatures
ctrl <- readRDS("data/exp_aist/daily_tmean_ctrl_calib_scaling_v2.rds")
noah <- readRDS("data/exp_aist/daily_tmean_noah_calib_scaling_v2.rds")

loc1 <- c("osaka","sakai","kobe","kyoto","nara","otsu","waka")
tm <- readRDS("data/exp_aist/daily_tmean_ctrl_noah_imp2004_coords.rds")
ctrl <- list()
for (i in seq(loc1)) {
  #sdat <- tm[tm$cat1=="ctrl",c("cat2","day",paste0(loc1[i],".tmean"))]; sdat$loc <- loc1[i]
  sdat <- tm[tm$cat1=="ctrl",c("cat2","day",loc1[i])]; sdat$loc <- loc1[i] #for coordinates data
  colnames(sdat) <- c("gw","date","t2m","loc")
  sdat1 <- sdat[,c("loc","gw","date","t2m")]
  ctrl[[loc[i]]] <- sdat1
}
rm(i,sdat,sdat1)
noah <- list()
for (i in seq(loc1)) {
  #sdat <- tm[tm$cat1=="noah",c("cat2","day",paste0(loc1[i],".tmean"))]; sdat$loc <- loc1[i]
  sdat <- tm[tm$cat1=="noah",c("cat2","day",loc1[i])]; sdat$loc <- loc1[i] #for coordinates data
  colnames(sdat) <- c("gw","date","t2m","loc")
  sdat1 <- sdat[,c("loc","gw","date","t2m")]
  noah[[loc[i]]] <- sdat1
}
rm(i,sdat,sdat1)


# create objects
#gw <- factor(c("none",paste0("+",seq(0.5,3,0.5),"k")),levels=c("none",paste0("+",seq(0.5,3,0.5),"k")))
gw <- factor(c("none","+0.5k","+1.0k","+1.5k","+2.0k","+2.5k","+3.0k"),
             levels=c("none","+0.5k","+1.0k","+1.5k","+2.0k","+2.5k","+3.0k"))
sim <- c(paste0("noah_",names(pred)[1]),paste0("ctrl_",names(pred)[2]),paste0("noah_",names(pred)[2]))
nsim <- 1000
yr <- 2000:2010
# everything is RCP 8.5, assumed no population growth

# create arrays and vectors
anloc <- afloc <- array(0,dim=c(length(loc),3,2,3,length(gw),length(sim)),
                        dimnames=list(loc,c("est","ci.l","ci.u"),c("abs","rel"),c("tot","cold","heat"),gw,sim))
antot <- aftot <- array(0,dim=c(3,2,3,length(gw),length(sim)),
                        dimnames=list(c("est","ci.l","ci.u"),c("abs","rel"),c("tot","cold","heat"),gw,sim))
dperiod <- rep(0,length(loc)); names(dperiod) <- loc

# loop
for (i in seq(sim)) {
  # print
  cat("\n\n",sim[i],"\n")
  # store uncertainty
  anlocsim <- array(0,dim=c(length(loc),2,3,length(gw),nsim+1),
                    dimnames=list(loc,c("abs","rel"),c("tot","cold","heat"),gw,c("est",paste0("sim",seq(nsim)))))
  # get projection list according to simulation
  tlist <- get(substr(sim[i],1,4)) #gets the first four letters
  
  # loop by locations or cities
  for (j in seq(loc)){
    # print
    cat(j,"")
    
    # tmean projection
    tproj <- tlist[[loc[j]]]
    
    # outcome projection
    ddoy <- tapply(mort[[loc[j]]]$mort,day(mort[[loc[j]]]$date),mean)
    #ddoy <- tapply(dlist[[loc[j]]]$mort,day(dlist[[loc[j]]]$date),mean)
    dproj <- rep(ddoy,length=nrow(tproj))
    speriod <- factor(rep(gw,each=31*length(yr))) #no period, just global warming categories
    dperiod[j] <- sum(ddoy)*length(yr)
    
    # model specs
    #lim <- grep(paste0(loc[j],"_",c(2000,2005),collapse="|"),names(temp[[loc[j]]]))
    #var1 <- temp[[loc[j]]][lim]
    var1 <- tmean[[loc[j]]]
    argvar <- list(fun="ns",knots=quantile(var1,vk),Bound=range(var1))
    coef1 <- pred[[substr(sim[i],6,nchar(sim[i]))]]$fit
    vcov1 <- pred[[substr(sim[i],6,nchar(sim[i]))]]$vcov
    cen1 <- cen[loc[j]]
    
    # multivariate normal distribution
    set.seed(1515)
    coefsim <- mvrnorm(nsim,coef1,vcov1)
    
    # centred basis
    bvar <- do.call(onebasis,c(list(x=tproj$t2m),argvar))
    cenvec <- do.call(onebasis,c(list(x=cen1),argvar))
    bvarcen <- scale(bvar,center=cenvec,scale=FALSE)
    
    # indicator for heat
    indheat <- tproj$t2m>cen1
    
    # daily attributable outcomes
    an <- (1-exp(-bvarcen%*%coef1))*dproj
    
    # sum by range and period
    anlocsim[j,1,1,,1] <- tapply(an,speriod,sum)
    anlocsim[j,1,2,unique(speriod[!indheat]),1] <- tapply(an[!indheat],factor(speriod[!indheat]),sum)
    anlocsim[j,1,3,unique(speriod[indheat]),1] <- tapply(an[indheat],factor(speriod[indheat]),sum)
    
    # loop by simulation
    for (s in seq(nsim)){
      # daily attributable outcomes
      an <- (1-exp(-bvarcen%*%coefsim[s,]))*dproj
      
      # simulations
      anlocsim[j,1,1,,s+1] <- tapply(an,speriod,sum)
      anlocsim[j,1,2,unique(speriod[!indheat]),s+1] <- tapply(an[!indheat],factor(speriod[!indheat]),sum)
      anlocsim[j,1,3,unique(speriod[indheat]),s+1] <- tapply(an[indheat],factor(speriod[indheat]),sum)
    }
  }
  # relative to no global warming
  anlocsim[,2,,,] <- anlocsim[,1,,,] - anlocsim[,1,,rep(1,length(gw)),] #no global warming or "none"
  
  # aggregate 
  antotsim <- apply(anlocsim[,,,,],2:length(dim(anlocsim)),sum)
  
  # attributable numbers 
  anloc[,1,,,,i] <- anlocsim[,,,,1] # mean
  anloc[,2,,,,i] <- apply(anlocsim[,,,,-1],1:4,quantile,0.025) # lower eCI
  anloc[,3,,,,i] <- apply(anlocsim[,,,,-1],1:4,quantile,0.975) # upper eCI
  antot[1,,,,i] <- apply(antotsim[,,,1],1:3,mean) # mean by location
  antot[2,,,,i] <- apply(antotsim[,,,-1],1:3,quantile,0.025) # lower eCI across locations
  antot[3,,,,i] <- apply(antotsim[,,,-1],1:3,quantile,0.975) # upper eCI across locations
  
  # attributable fractions (net effect)
  afloc[,,,,,i] <- anloc[,,,,,i]/dperiod*100
  aftot[,,,,i] <- antot[,,,,i]/sum(dperiod)*100
}

# save
save(cen,anloc,afloc,antot,aftot,file="proj_outputs_kansai_ac_all-cause-ages_calib.rda")