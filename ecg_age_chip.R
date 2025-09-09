# Depends
library(data.table)
library(stringr)
library(plyr)
library(survival)
source('./functions/functions.R')

# Load age inferences
## Vector of true instance 2 ages
true_age <- fread(file='/phenotypes/2024_01/instances_all_202401.csv')[instance==0]
censor <- fread(file='/phenotypes/2024_01/censor_202401.csv')
censor_old <- fread(file='/phenotypes/2021_09/censor_202109.csv')
setkey(censor,sample_id); setkey(true_age,sample_id); setkey(censor_old,sample_id)
true_age[censor,birth_date := as.Date(i.birthdate,format='%Y-%m-%d')]
true_age[censor_old,sex := i.sex]
true_age[,actual_age := (as.numeric(value)-as.numeric(birth_date))/365.25]

## Remove withdrawals
withdrawals <- fread(file='/phenotypes/withdrawals/w7089_20241216.csv')
true_age <- true_age[!(sample_id %in% withdrawals$V1)] #502250 - 117 = 502133

## ECG
ecg_age0 <- fread(file='/age_chip/inference_age_0.tsv')
setnames(ecg_age0,c('21003_Age-when-attended-assessment-centre_0_0_prediction','21003_Age-when-attended-assessment-centre_0_0_actual'),
         c('predicted_age0','actual_age0'))

ecg_age1 <- fread(file='inference_age_1.tsv')
setnames(ecg_age0,c('21003_Age-when-attended-assessment-centre_0_0_prediction','21003_Age-when-attended-assessment-centre_0_0_actual'),
         c('predicted_age0','actual_age0'))

setkey(ecg_age,sample_id); setkey(true_age,sample_id)
ecg_age[true_age,':='(actual_age = i.actual_age, sex = i.sex)]

## CHIP calls
chip <- fread(file='/chip_calls/ukb_chip_450k_oct2022_small.txt')

## PCS
pcs <- fread(file='/central_qc/ukb_sqc_v2_7089.tsv')

## Race
race <- fread(file='/phenotypes/2020_06/race_0.csv')

## Remove no instance 0 ecg
ecg_bike_method = fread(file='ecg_bike_method.csv')
instance0_ecg = ecg_bike_method[instance==0]
ecg_age <- ecg_age[sample_id %in% instance0_ecg$sample_id]

## Remove missing actual age
ecg_age <- ecg_age[!is.na(actual_age)] # 502133 - 431073 = 71060
ecg_age[,delta_age := predicted_age - actual_age]
ecg_age[,delta_age_abs := abs(predicted_age - actual_age)]
ecg_age[,':='(actual_age10 = actual_age/10,
              predicted_age10 = predicted_age/10)]

## Merge
setkey(ecg_age,sample_id); setkey(chip,id)
chip_analysis <- ecg_age[chip,nomatch=0] # 71060 - 5764 = 65296

## Add PCS
setkey(pcs,eid); setkey(chip_analysis,sample_id)
chip_analysis[pcs,':='(PC1 = i.PC1,PC2 = i.PC2,PC3 = i.PC3,PC4 = i.PC4,PC5 = i.PC5)]
chip_analysis <- chip_analysis[!is.na(PC1) & !is.na(PC2) & !is.na(PC3) & !is.na(PC4) & !is.na(PC5)] #65296 - 64 = 65232

## Add race
setkey(race,sample_id); setkey(chip_analysis,sample_id)
chip_analysis[race,':='(race_binary = ifelse(i.value %in% c(1,1001,1002,1003),1,0))]

# Corr plot
pdf(file='ecg_age_ukb_trained0.pdf',height=4,width=4,pointsize=5)
par(mar=c(5,5,2,2),oma=c(1,1,1,1))
plot(x=chip_analysis$actual_age,y=chip_analysis$predicted_age,
     xaxt='n',yaxt='n',xlab='',ylab='',bty='n',xlim=c(40,70),ylim=c(40,70),
     pch=19,col='#2c7fb84D')
axis(1,at=seq(30,70,10),cex.axis=1.6)
axis(2,at=seq(30,70,10),cex.axis=1.6,las=2)
mtext('Actual age',1,line=3.2,cex=1.8)
mtext("Predicted age",2,line=3.5,cex=1.8)
segments(-0.1,-0.1,71,71,lty=5)
dev.off()

cor.test(chip_analysis$actual_age,y=chip_analysis$predicted_age)

# MAE
boot_mae <- ci_mae(y='actual_age',x1='predicted_age',
                   data=chip_analysis,runs=1000)
raw_mae <- mean(abs(chip_analysis$actual_age - chip_analysis$predicted_age))
ci <- c(raw_mae,raw_cimae - 1.96*sd(boot_mae),raw_mae + 1.96*sd(boot_mae))

# BA plot
chip_analysis[,mean_pair := apply(.SD,FUN=mean,MARGIN=1),.SDcols=c('actual_age','predicted_age')]
chip_analysis[,diff_pair := predicted_age - actual_age]
setkey(chip_analysis,mean_pair)

png('ecg_age_ukb_trained_ba0.png',pointsize=6,res=300,
    height=1200,width=1200)
par(oma=c(1,1,1,1))
par(mar=c(4,3,1,1))
col <- "#2c7fb84D"

# Plot
plot(x=chip_analysis$mean_pair,y=chip_analysis$diff_pair,
     pch=19,col=col,xlab='',ylab='',xaxt='n',yaxt='n',frame=F,
     xlim=c(40,90),ylim=c(-40,40),cex=0.8)

# CI Lines
mean_diff <- mean(chip_analysis$diff_pair)
upper <- mean(chip_analysis$diff_pair)+1.96*sd(chip_analysis$diff_pair)
lower <- mean(chip_analysis$diff_pair)-1.96*sd(chip_analysis$diff_pair)
segments(40,upper,91,upper,lty=5,lwd=1,col='black')
segments(40,lower,91,lower,lty=5,lwd=1,col='black')
segments(40,mean_diff,91,mean_diff,lty=1,lwd=1,col='#f03b20')

# Axes
axis(1,cex.axis=1.6,pos=-40,at=seq(0,90,10))
axis(2,cex.axis=1.6,las=2,pos=40,at=seq(-40,40,10))
mtext(side=1,"Mean of ages",cex=1.8,line=2.8)
mtext(side=2,"Predicted minus actual age",cex=1.8,line=2.5)

# Agreement
par(xpd=TRUE)
agree <- nrow(chip_analysis[c((diff_pair >= lower) & (diff_pair <= upper))])/nrow(chip_analysis)
text(x=70,y=-30,labels=paste0("Limits of Agreement: ",round(lower,2),'-',round(upper,2)),cex=1.4)
text(x=70,y=-32,labels=paste0("% Within Limits: ",round(agree*100,1)),cex=1.4)

dev.off()

#### Primary model
outcomes <- c('CHIP','DNMT3A','TET2','ASXL1','JAK2','PPM1D',
              'TP53','SF3B1','SRSF2','DTA','DDR','Splice','other')

out <- list(); n <- 1
for (j in outcomes){
  mod_chrono <- glm(chip_analysis[,get(j)] ~ chip_analysis$actual_age10 + chip_analysis$sex + chip_analysis$PC1 + chip_analysis$PC2 + chip_analysis$PC3 + chip_analysis$PC4 + chip_analysis$PC5,family='binomial')
  coeff_chrono <- exp(c(coefficients(mod_chrono)[2],confint(mod_chrono)[2,]))
  mod_ecg <- glm(chip_analysis[,get(j)] ~ chip_analysis$predicted_age10 + chip_analysis$sex + chip_analysis$PC1 + chip_analysis$PC2 + chip_analysis$PC3 + chip_analysis$PC4 + chip_analysis$PC5,family='binomial')
  coeff_ecg <- exp(c(coefficients(mod_ecg)[2],confint(mod_ecg)[2,]))
  mod_both <- glm(chip_analysis[,get(j)] ~ chip_analysis$actual_age10 + chip_analysis$predicted_age10 + chip_analysis$sex + chip_analysis$PC1 + chip_analysis$PC2 + chip_analysis$PC3 + chip_analysis$PC4 + chip_analysis$PC5,family='binomial')
  coeff_both <- exp(c(coefficients(mod_both)[3],confint(mod_both)[3,]))
  out[[n]] <- data.table(outcome=j,
                         actual_age_or=coeff_chrono[1],actual_age_lower=coeff_chrono[2],actual_age_upper=coeff_chrono[3],
                         pred_age_or=coeff_ecg[1],pred_age_lower=coeff_ecg[2],pred_age_upper=coeff_ecg[3],
                         both_age_or=coeff_both[1],both_age_lower=coeff_both[2],both_age_upper=coeff_both[3])
  n <- n + 1
}

all_chip <- do.call(rbind,out)
write.csv(all_chip,file='exercise_all_chip0_112624.csv',row.names=F)

#### Age diff
outcomes <- c('CHIP','DNMT3A','TET2','ASXL1','JAK2','PPM1D',
              'TP53','SF3B1','SRSF2','DTA','DDR','Splice','other')

chip_analysis[,age_diff := predicted_age - actual_age]

out <- list(); n <- 1
for (j in outcomes){
  mod_diff <- glm(chip_analysis[,get(j)] ~ chip_analysis$age_diff + chip_analysis$sex + chip_analysis$PC1 + chip_analysis$PC2 + chip_analysis$PC3 + chip_analysis$PC4 + chip_analysis$PC5,family='binomial')
  coeff_diff <- exp(c(coefficients(mod_diff)[2],confint(mod_diff)[2,]))
  mod_both <- glm(chip_analysis[,get(j)] ~ chip_analysis$actual_age + chip_analysis$age_diff + chip_analysis$sex + chip_analysis$PC1 + chip_analysis$PC2 + chip_analysis$PC3 + chip_analysis$PC4 + chip_analysis$PC5,family='binomial')
  coeff_both <- exp(c(coefficients(mod_both)[3],confint(mod_both)[3,]))
  out[[n]] <- data.table(outcome=j,
                         diff_or=coeff_diff[1],diff_lower=coeff_diff[2],diff_upper=coeff_diff[3],
                         diff_adj_or=coeff_both[1],diff_adj_lower=coeff_both[2],diff_adj_upper=coeff_both[3])
  n <- n + 1
}

diff_chip <- do.call(rbind,out)
write.csv(diff_chip,file='exercise_diff_chip0_112624.csv',row.names=F)

#### Expanded
outcomes <- c('expandedCHIP','expandedDNMT3A','expandedTET2','expandedASXL1','expandedJAK2','expandedPPM1D',
              'expandedTP53','expandedSF3B1','expandedSRSF2','expandedDTA','expandedDDR','expandedSplice')

out <- list(); n <- 1
for (j in outcomes){
  mod_chrono <- glm(chip_analysis[,get(j)] ~ chip_analysis$actual_age10 + chip_analysis$sex + chip_analysis$PC1 + chip_analysis$PC2 + chip_analysis$PC3 + chip_analysis$PC4 + chip_analysis$PC5,family='binomial')
  coeff_chrono <- exp(c(coefficients(mod_chrono)[2],confint(mod_chrono)[2,]))
  mod_ecg <- glm(chip_analysis[,get(j)] ~ chip_analysis$predicted_age10 + chip_analysis$sex + chip_analysis$PC1 + chip_analysis$PC2 + chip_analysis$PC3 + chip_analysis$PC4 + chip_analysis$PC5,family='binomial')
  coeff_ecg <- exp(c(coefficients(mod_ecg)[2],confint(mod_ecg)[2,]))
  mod_both <- glm(chip_analysis[,get(j)] ~ chip_analysis$actual_age10 + chip_analysis$predicted_age10 + chip_analysis$sex + chip_analysis$PC1 + chip_analysis$PC2 + chip_analysis$PC3 + chip_analysis$PC4 + chip_analysis$PC5,family='binomial')
  coeff_both <- exp(c(coefficients(mod_both)[3],confint(mod_both)[3,]))
  out[[n]] <- data.table(outcome=j,
                         actual_age_or=coeff_chrono[1],actual_age_lower=coeff_chrono[2],actual_age_upper=coeff_chrono[3],
                         pred_age_or=coeff_ecg[1],pred_age_lower=coeff_ecg[2],pred_age_upper=coeff_ecg[3],
                         both_age_or=coeff_both[1],both_age_lower=coeff_both[2],both_age_upper=coeff_both[3])
  n <- n + 1
}

expanded_chip <- do.call(rbind,out)
write.csv(expanded_chip,file='exercise_expanded_chip0_112624.csv',row.names=F)

#C-index
mod_basic <- glm(chip_analysis$CHIP ~ chip_analysis$actual_age10  + chip_analysis$sex + chip_analysis$PC1 + chip_analysis$PC2 + chip_analysis$PC3 + chip_analysis$PC4 + chip_analysis$PC5,family='binomial')
mod_plus <- glm(chip_analysis$CHIP ~ chip_analysis$predicted_age10 + chip_analysis$sex + chip_analysis$PC1 + chip_analysis$PC2 + chip_analysis$PC3 + chip_analysis$PC4 + chip_analysis$PC5,family='binomial')
concordance(mod_basic)
concordance(mod_plus)

mod_tet2 <- glm(chip_analysis$TET2 ~ chip_analysis$actual_age10 + chip_analysis$predicted_age10 + chip_analysis$sex + chip_analysis$PC1 + chip_analysis$PC2 + chip_analysis$PC3 + chip_analysis$PC4 + chip_analysis$PC5,family='binomial')
mod_dta <- glm(chip_analysis$DTA ~ chip_analysis$actual_age10 + chip_analysis$predicted_age10 + chip_analysis$sex + chip_analysis$PC1 + chip_analysis$PC2 + chip_analysis$PC3 + chip_analysis$PC4 + chip_analysis$PC5,family='binomial')


