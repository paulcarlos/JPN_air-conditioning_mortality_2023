library(splines)
library(dlnm)
library(readxl)

# import data
dat <- readRDS("cities_52_japan_1976-2009.rds")

# list of cities and years
city <- read_excel("city_code.xlsx",sheet=1)
y1 <- list(1976:1979,1980:1984,1985:1989,1990:1994,1995:1999,2000:2004,2005:2009)
y2 <- c("1976_1979","1980_1984","1985_1989","1990_1994","1995_1999","2000_2004","2005_2009")
y3 <- c("1976-1979","1980-1984","1985-1989","1990-1994","1995-1999","2000-2004","2005-2009")

# import AC prevalence
ac <- read.csv("ac_prev_1972-2009.csv",stringsAsFactors = FALSE)
ac <- data.frame(t(ac)); colnames(ac) <- ac[1,]; ac <- ac[-1,]; ac$year <- 1972:2009; ac <- ac[ac$year %in% 1976:2009,]

# model specifications using Sera et al 2020 (doi:10.1097/EDE.0000000000001241)
tlag <- 2 # max lag days
vk <- c(0.5,0.9) # variable knots

# storage for coefficients and variance-covariance
coef1 <- matrix(data=NA,nrow=0,ncol=length(vk)+1)
vcov1 <- list()

# storage for meta-regressors
mreg <- data.frame(matrix(NA,nrow=0,ncol=6)) ; colnames(mreg) <- c("id","year","av_ac","av_temp","iqr_temp","loc")

# storage daily mean temperatures
tmean <- list()

# pairs from locations and time periods
g1 <- data.frame("loc"=rep(1:length(dat),each=length(y1)),"period"=rep(1:length(y1),times=length(dat)))

# first stage models with plots
for (i in 1:length(dat)){
  cat(i," ")
  # select data
  d1 <- dat[[i]]
  # store mean temperatures
  x1 <- unlist(lapply(d1,FUN=function(x)x$tmean))
  tmean[[i]] <- x1
  names(tmean)[i] <- names(dat)[i]
  # loop coefficients and vcov
  for (j in 1:length(d1)) {
    # sub data
    s1 <- d1[[j]]
    var <- s1$tmean
    # create crossbasis
    cb1 <- crossbasis(var,lag=c(0,tlag),argvar=list(fun="ns",knots=quantile(var,vk)),arglag=list(fun="integer"),group=s1$year)
    # centring in August
    int1 <- ((s1$doy-75)/120)*cb1
    # model
    mod1 <- glm(mort~cb1+ns(doy,4):factor(year)+factor(dow)+int1,family=quasipoisson(),data=s1)
    # reduce
    red1 <- crossreduce(cb1,mod1,cen=quantile(var,0.5))
    #plot(red1)
    # store coefficients and variance-covariance
    r <- as.numeric(row.names(g1[g1$loc==i & g1$period==j,]))
    coef1 <- rbind(coef1,coef(red1))
    dimnames(coef1)[[1]][r] <- paste0(names(dat)[i],"_",y2[j])
    vcov1[[r]] <- vcov(red1)
    names(vcov1)[r] <- paste0(names(dat)[i],"_",y2[j])
    #meta regressors
    mreg[r,] <- c(r,round(mean(y1[[j]])),mean(as.numeric(ac[ac$year%in%y1[[j]],city$Prefecture[i]])),mean(var),IQR(var),names(dat)[i])
  }
}
rm(i,d1,j,s1,var,cb1,int1,mod1,red1,r)

# convert data types
mreg$id <- as.integer(mreg$id)
mreg$year <- as.integer(mreg$year)
mreg$av_ac <- as.numeric(mreg$av_ac)
mreg$av_temp <- as.numeric(mreg$av_temp)
mreg$iqr_temp <- as.numeric(mreg$iqr_temp)
mreg$loc <- factor(mreg$loc,levels=city$City)
mreg$region <- city$Region[charmatch(mreg$loc,city$City)]
mreg$region <- factor(mreg$region,levels=c("Kanto","Hokkaido","Tohoku","Chubu","Kansai","Chugoku","Shikoku","Kyushu"))

# save RDA file
save(coef1,vcov1,mreg,tmean,file="cities52_1st-stage_loc-time_coef_vcov.rda")