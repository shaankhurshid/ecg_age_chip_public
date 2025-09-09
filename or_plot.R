library(Cairo)
library(data.table)

# Load overall results
overall <- fread(file='exercise_all_chip_smoke3_070825.csv')
overall$outcome[1] <- 'All CHIP'; overall$outcome[4] <- 'Other CHIP'

CairoPDF(file='or_plot_predicted_adjusted.pdf',height=4,width=6)
par(oma=c(2,2.5,1,1),mar=c(2,2.5,1,1))

x_var <- overall$both_age_or
y_var <- 12:1
col <- ifelse(overall$both_age_lower >= 1,'#a50026',
              ifelse(overall$both_age_lower >= 0.97,'#fe9929','darkgray')); col[1] <- 'black'
sym <- c(18,rep(16,12))

plot(x=0.1,y=0.1,col='white',log='x',
     bty='n',xaxt='n',yaxt='n',xlim=c(0.5,16),
     ylim=c(0.5,12.5),xlab='',ylab='')

segments(1,0,1,12.5,lwd=0.75,lty=5,col='darkgray')

segments(ifelse(overall$both_age_lower < 0.5,0.5,overall$both_age_lower),y_var,ifelse(overall$both_age_upper > 16,16,overall$both_age_upper),y_var,col=col)
arrows(overall$both_age_lower[6],7,16,7,col='darkgray',code=2,lwd=1,angle=45,length=0.1)
arrows(0.5,4,overall$both_age_upper[9],4,col='darkgray',code=1,lwd=1,angle=45,length=0.1)
arrows(0.5,3,overall$both_age_upper[10],3,col='darkgray',code=1,lwd=1,angle=45,length=0.1)
arrows(0.5,2,overall$both_age_upper[11],2,col='darkgray',code=1,lwd=1,angle=45,length=0.1)
arrows(0.5,1,overall$both_age_upper[12],1,col='darkgray',code=1,lwd=1,angle=45,length=0.1)

par(new=TRUE)
plot(x=x_var,y=y_var,col=col,log='x',
     bty='n',xaxt='n',yaxt='n',pch=sym,xlim=c(0.5,16),
     ylim=c(0.5,12.5),xlab='',ylab='',cex=c(1.8,rep(1,16)))

axis(1,at=c(0.5,1,2,4,8,16),cex.axis=0.8,pos=0.3,
     labels=c(paste0(c('0.5','1','2','4','8','16'))),padj=-1)
axis(2,at=12:1,cex.axis=0.8,pos=0.49,
     labels=paste0(overall$outcome),las=2)

mtext("Odds ratio per 10-year increase",1,line=1.5,cex=1)

dev.off()