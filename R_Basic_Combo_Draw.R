setwd('c:/users/user/desktop/R')
library(openxlsx)
uptemp <- read.xlsx('Rh.xlsx',colNames=T,cols=c(1,2,3,4))
downtemp <- read.xlsx('Rh.xlsx',colNames=T,cols=c(5,6,7,8))
colors <- c('Blue','Gray','green','Yellow','Orange','Red')
uck <- uptemp[uptemp$treatment1=='CK',]
uhn <- uptemp[uptemp$treatment1=='HN',]
uln <- uptemp[uptemp$treatment1=='LN',]
uw <- uptemp[uptemp$treatment1=='W',]
uwhn <- uptemp[uptemp$treatment1=='WHN',]
uwln <- uptemp[uptemp$treatment1=='WLN',]
dck <- downtemp[downtemp$treatment2=='CK',]
dhn <- downtemp[downtemp$treatment2=='HN',]
dln <- downtemp[downtemp$treatment2=='LN',]
dw <- downtemp[downtemp$treatment2=='W',]
dwhn <- downtemp[downtemp$treatment2=='WHN',]
dwln <- downtemp[downtemp$treatment2=='WLN',]
#png('test.png',width=674,height=359)
opar <- par(no.readonly=T)
par(fig=c(0,0.585,0,0.9),mar=c(4,4,2.1,0))
plot(uck$temperature1,uck$Rh1,type='p',pch=20,bty='l',col='Blue',
     ylab='Rh',xlab=' ',ylim=c(0,48),xlim=c(4,26),mgp=c(2,0.5,0))
lines(lowess(uck$temperature1,uck$Rh1,f=0.5),col='Blue')
points(uhn$temperature1,uhn$Rh1,pch=20,col='Gray')
lines(lowess(uhn$temperature1,uhn$Rh1,f=0.5),col='Gray')
points(uln$temperature1,uln$Rh1,pch=20,col='green')
lines(lowess(uln$temperature1,uln$Rh1,f=0.5),col='green')
points(uw$temperature1,uw$Rh1,pch=20,col='Yellow')
lines(lowess(uw$temperature1,uw$Rh1,f=0.5),col='Yellow')
points(uwhn$temperature1,uwhn$Rh1,pch=20,col='Orange')
lines(lowess(uwhn$temperature1,uwhn$Rh1,f=0.5),col='Orange')
points(uwln$temperature1,uwln$Rh1,pch=20,col='Red')
lines(lowess(uwln$temperature1,uwln$Rh1,f=0.5),col='Red')
arrows(x0=5,y0=45,x1=24.5,y1=45,col='black',length=0.15,angle=12)
text(15,47.5,'Rising',col='black')
par(fig=c(0.51,1,0,0.9),mar=c(4,0,2.1,0.5),new=T)
plot(dck$temperature2,dck$Rh2,type='p',pch=20,bty='l',col='Blue',
     axes=F,xlab=' ',ylim=c(0,48),xlim=c(26,4))
lines(lowess(dck$temperature2,dck$Rh2,f=0.5),col='Blue')
axis(1,at=c(25,20,15,10,5),labels=c(' ',20,15,10,5),mgp=c(2,0.5,0))
points(dhn$temperature2,dhn$Rh2,pch=20,col='Gray')
lines(lowess(dhn$temperature2,dhn$Rh2,f=0.5),col='Gray')
points(dln$temperature2,dln$Rh2,pch=20,col='green')
lines(lowess(dln$temperature2,dln$Rh2,f=0.5),col='green')
points(dw$temperature2,dw$Rh2,pch=20,col='Yellow')
lines(lowess(dw$temperature2,dw$Rh2,f=0.5),col='Yellow')
points(dwhn$temperature2,dwhn$Rh2,pch=20,col='Orange')
lines(lowess(dwhn$temperature2,dwhn$Rh2,f=0.5),col='Orange')
points(dwln$temperature2,dwln$Rh2,pch=20,col='Red')
lines(lowess(dwln$temperature2,dwln$Rh2,f=0.5),col='Red')
arrows(x0=24.5,y0=45,x1=5,y1=45,col='black',length=0.12,angle=12)
text(15,47.5,'Cooling',col='black')
legend(x=9,y=43,legend=c('CK','HN','LN','W','WHN','WLN'),col=colors,
      cex=0.6,pt.cex=1,pch=20)
par(fig=c(0,1,0,0.9),mar=c(4,4,2.1,0.5),new=T)
title(main='xiao wang',xlab='temperature',mgp=c(2,0.5,0))