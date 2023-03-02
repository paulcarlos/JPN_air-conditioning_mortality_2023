library(splines)
library(dlnm)
library(mixmeta)

# WALD TEST
fwald <- function(model,var) {
  ind <- grep(var,names(coef(model)),fixed=TRUE)
  coef <- coef(model)[ind]
  vcov <- vcov(model)[ind,ind]
  waldstat <- coef%*%solve(vcov)%*%coef
  df <- length(coef)
  return(1-pchisq(waldstat,df))
}

# percentile of temperature
ptile <- function(data,value) ecdf(data)(value)

# Kansai cities
loc <- c("Osaka","Sakai","Kobe","Kyoto","Nara","Otsu","Wakayama")

# import data
load("cities52_1st-stage_loc-time_coef_vcov.rda")

# meta-regression
fixed <- formula(coef1~av_ac+av_temp+iqr_temp+ns(year,knots=1992))
mod <- mixmeta(fixed,vcov1,random=~1|loc,control=list(showiter=TRUE),data=mreg,method="reml")
summary(mod)
blup1 <- blup(mod,vcov=TRUE)

# test effects
fwald(mod,"av_ac") # air conditioning

# predict
dat <- data.frame("av_temp"=NA,"iqr_temp"=NA,"year"=2007,"av_ac"=rep(c(0,100),times=length(loc)),loc=rep(loc,each=2))
for (i in loc) {
  dat$av_temp[dat$loc==i] <- mreg$av_temp[mreg$year==2007 & mreg$loc==i]
  dat$iqr_temp[dat$loc==i] <- mreg$iqr_temp[mreg$year==2007 & mreg$loc==i]
}
rm(i)
pred <- predict(mod,dat,vcov=TRUE,format="list")
names(pred) <- paste0(dat$loc,dat$av_ac)

# temp
x1 <- unname(rowMeans(do.call(cbind,tmean)))
bvar <- onebasis(x1,fun="ns",knots=quantile(x1,c(0.5,0.9)))
# center
pctl <- c(1:99)
predvar <- quantile(x1,pctl/100,na.rm=T)
argvar1 <- list(x=predvar,fun="ns",knots=quantile(x1,c(0.5,0.9)),Bound=range(x1,na.rm=T))
bvar1 <- do.call(onebasis,argvar1)
mperc <- (pctl)[which.min((bvar1%*%pred[[2]]$fit))]
cen1 <- round(quantile(x1,mperc/100,na.rm=T),1)
# difference
#cp1 <- crosspred(bvar,coef=pred[[1]]$fit,vcov=pred[[1]]$vcov,model.link="log",by=0.1,cen=cen1)
#cp2 <- crosspred(bvar,coef=pred[[2]]$fit,vcov=pred[[2]]$vcov,model.link="log",by=0.1,cen=cen1)
#plot(cp1$allRRfit-cp2$allRRfit,xaxt="n")
#axis(1,at=1:length(cp1$allRRfit),labels=names(cp1$allRRfit))
#abline(h=0)
# plot
plot(0,type="n",xlim=c(min(x1),max(x1)),xaxt="n",xlab="temperature",ylim=c(0.9,1.5),ylab="relative risk",
     main=paste0("Overall - 52 cities"," (I^2=",round(summary(mod)$i2stat[1]),"%)"))
axis(1,at=seq(floor(min(x1)),ceiling(max(x1)),1))
abline(h=1)
legend("topleft",legend=c(paste0(dat$year[1],"|AC=",round(dat$av_ac[1]),"%"),
                          paste0(dat$year[2],"|AC=",round(dat$av_ac[2]),"%")),
       lty=1,lwd=3,cex=1.2,col=1:3,bty="n")
for (i in 1:length(pred)) {
  cp1 <- crosspred(bvar,coef=pred[[i]]$fit,vcov=pred[[i]]$vcov,model.link="log",by=0.1,cen=cen1)
  lines(cp1,ci="area",lwd=3,col=i,ci.arg=list(density=10+(10*i),angle=60,col=i))
}
abline(v=cen1,lty=3)
text(cen1,1.4,paste0(cen1," - ",names(cen1)),cex=1.2)
rm(i,cp1,x1,bvar,pctl,predvar,argvar1,bvar1,mperc,cen1)

# select center or MMT
cen <- rep(0,length(loc)); names(cen) <- loc
pctl <- c(1:99)
par(mfrow=c(2,3),oma=c(4,4,2,2),mar=c(2,3,3,2))
for (i in 1:length(loc)){
  # temperatures
  lim <- grep(paste0(loc[i],"_",c(2000:2009),collapse="|"),names(tmean[[loc[i]]]))
  x1 <- tmean[[loc[i]]][lim] #temperatures 2000-2009
  #plot(x1)
  x2 <- split(x1,rep(2000:2009,each=122))
  x1 <- rowMeans(do.call(cbind,x2))
  # center
  predvar <- quantile(x1,pctl/100,na.rm=T)
  argvar1 <- list(x=predvar,fun="ns",knots=quantile(x1,c(0.5,0.9)),Bound=range(x1,na.rm=T))
  bvar <- do.call(onebasis,argvar1)
  mperc <- (pctl)[which.min((bvar%*%pred[[2]]$fit))]
  cen1 <- round(quantile(x1,mperc/100,na.rm=T),1)
  cen[loc[i]] <- cen1
  # crosspredict
  bvar1 <- onebasis(x1,fun="ns",knots=quantile(x1,c(0.5,0.9)))
  p1 <- paste0(loc[i],100); p2 <- paste0(loc[i],0)
  cp1 <- crosspred(bvar1,coef=pred[[p1]]$fit,vcov=pred[[p1]]$vcov,model.link="log",by=0.1,cen=cen1)
  cp2 <- crosspred(bvar1,coef=pred[[p2]]$fit,vcov=pred[[p2]]$vcov,model.link="log",by=0.1,cen=cen1)
  # sample plot
  plot(0,type="n",xlim=c(floor(min(x1)),ceiling(max(x1))),xaxt="n",xlab="",ylim=c(0.8,1.6),ylab="relative risk",main=loc[i])
  axis(1,at=seq(floor(min(x1)),ceiling(max(x1)),1))
  abline(h=1)
  lines(cp1,ci="area",lwd=3,col="blue",ci.arg=list(density=20,angle=60,col="blue"))
  lines(cp2,ci="area",lwd=3,col="red",ci.arg=list(density=20,angle=30,col="red"))
  abline(v=cen1,lty=3,col="blue")
  text(cen1,1.4,paste0(cen1," - ",names(cen1)),cex=1.2)
  legend("topleft",legend=c(paste0("AC=",round(dat$av_ac[1]),"%"),paste0("AC=",round(dat$av_ac[2]),"%")),lty=1,lwd=3,cex=1.2,col=c("red","blue"),bty="n")
}
#mtext("mean temperatures °C",side=1,line=1.5,outer=TRUE)
#mtext("relative risk",side=2,line=1.5,outer=TRUE)
rm(i,lim,x1,predvar,argvar1,bvar,mperc,cen1,bvar1,cp1,cp2)

# save
#temp <- tmean[loc]
#cen <- lapply(temp,quantile,mperc/100)
save(tmean,cen,pred,file="data/rda/cities52_2nd-stage_predicted-AC_coef_vcov.rda")

