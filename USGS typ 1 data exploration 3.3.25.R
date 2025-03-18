# The purpose of this document is to analyze trends form the 2025 data 
#update put out by Tristan Mohs 


setwd("C:/Users/brown/OneDrive/Desktop/USGS CMP/USGS")

data <- read.csv("Summary_Metrics_Deliverable.csv")
colnames(data)

d1 <- read.csv("Site_Totals_Density_Biomass_Richness2.csv")
colnames(d1)

df <- cbind(d1,data, na.rm=TRUE)

plot(df$Streambank_Lateral_Root_Exposure_Distance_Mean, df$Streambank_Lateral_Erosion_Rate_Mean)

rtlm <- lm(Streambank_Lateral_Erosion_Rate_Mean~Streambank_Lateral_Root_Exposure_Distance_Mean, data=df)
summary(rtlm)
abline(rtlm)

#generate residuals plot 
fit_vals <- rtlm$fitted
stud_resid <- rstudent(rtlm)

#residuals vs fitted plot 
plot(fit_vals,stud_resid)
abline(h=0, col="red"
       pch = 16,
       main = 'Studentized Residuals vs. Fitted Values',
       xlab = 'Fitted (a.k.a Predicted) (a.k.a Y-hat) Values',
       ylab = 'Studentized Residuals')
abline(h = 0, col = 'grey')
abline(h = c(-2.5,2.5),col ='blue', lty = 2)

#qqplot to check normality of residuals 
qqnorm(stud_resid)
abline(0,1,col='red')
#the residuals are not normally distributed 


#Biomass plots
plot(df$Streambank_Lateral_Erosion_Rate_Mean, df$Mean_Biomass)
#density plots
plot(df$Streambank_Lateral_Erosion_Rate_Mean, df$Mean_Density)
#Taxarichness plots
plot(df$Embed_Median, df$Mean_Taxa_Richness)

#biomass model 
Bmdl <- lm(Mean_Biomass~Streambank_Lateral_Erosion_Rate_Mean, data=df)
summary(Bmdl)

#density model
Dmdl <- lm(Mean_Density~Streambank_Lateral_Erosion_Rate_Mean, data=df)
summary(Dmdl)

#Richness model 
Rmdl <- lm(Mean_Taxa_Richness~Embed_Median, data=df)
summary(Rmdl)
#reg plot
plot(df$Embed_Median, df$Mean_Taxa_Richness, 
     xlab="Median Embeded Score",
     ylab="Mean Taxa Richness")
abline(Rmdl, col="red",
       pch = 16)
library (ggplot2)

library(ggplot2) 
ggplot(Rmdl$model, aes_string(x = names(Rmdl$model)[2], y = names(Rmdl$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  geom_label(aes(x = 65, y = 70), hjust = 0, 
             label = paste("Adj R2 = ",signif(summary(Rmdl)$adj.r.squared, 5),
                           "\nIntercept =",signif(Rmdl$coef[[1]],5 ),
                           " \nSlope =",signif(Rmdl$coef[[2]], 5),
                           " \nP =",signif(summary(Rmdl)$coef[2,4], 5)))




#generate residuals plot 
fit_vals <- Rmdl$fitted
stud_resid <- rstudent(Rmdl)

#residuals vs fitted plot 
plot(fit_vals,stud_resid)
abline(h=0, col="red"
       pch = 16,
       main = 'Studentized Residuals vs. Fitted Values',
       xlab = 'Fitted (a.k.a Predicted) (a.k.a Y-hat) Values',
       ylab = 'Studentized Residuals')
abline(h = 0, col = 'grey')
abline(h = c(-2.5,2.5),col ='blue', lty = 2)

#qqplot to check normality of residuals 
  qqnorm(stud_resid)
abline(0,1,col='red')
#the residuals are not normally distributed 

#full mdl for adjusted lateral erosion 
full_mdl <- lm(Adjusted_Streambank_Lateral_Erosion_Rate_Mean~., data=df)


#Adjusted Streambank erosion
plot(df$ws_pasture__prop, df$Streambank_Lateral_Erosion_Rate_Mean)


#BMP Score tests 
plot(df$BMP.Score, df$Streambank_Lateral_Erosion_Rate_Mean)
#interesting negative relationship here, lots of clustering 
# near zero and 1
bpElm <- lm(Streambank_Lateral_Erosion_Rate_Mean~BMP.Score, data=df)
summary(bpElm)


plot(df$BMP.Score, df$Streambank_Lateral_Root_Exposure_Distance_Mean)
#intersting negative relationship here as well, with clustering 
# near zero and 1 
bpRlm <- lm(Streambank_Lateral_Root_Exposure_Distance_Mean~BMP.Score, data=df)
summary(bpRlm)


plot(df$BMP.Score, df$Instream_Wood_Volume)
wblm<- lm(Instream_Wood_Volume ~ BMP.Score, data=df)
summary(wblm)
ggplot(wblm$model, aes_string(x = names(wblm$model)[2], y = names(wblm$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  geom_label(aes(x = 3, y = 2), hjust = 0, 
             label = paste("Adj R2 = ",signif(summary(wblm)$adj.r.squared, 5),
                           "\nIntercept =",signif(wblm$coef[[1]],5 ),
                           " \nSlope =",signif(wblm$coef[[2]], 5),
                           " \nP =",signif(summary(wblm)$coef[2,4], 5)))


#BMP score and bug metrics 

#richenss and bmp score
plot(df$BMP.Score, df$Mean_Taxa_Richness)
BMPR <- lm(Mean_Taxa_Richness ~ BMP.Score, data=df)
summary(BMPR)
#non significant

#biomass and BMP socre 
plot(df$BMP.Score, df$Mean_Biomass)
BMPB <- lm(Mean_Biomass ~ BMP.Score, data=df)
summary(BMPB)

#Densty and BMP score 
plot(df$BMP.Score, df$Mean_Density)
BMPD <- lm(Mean_Density ~ BMP.Score, data=df)
summary(BMPD)

#pasture prop and Bug metrics 
plot(df$ws_pasture__prop, df$Mean_Taxa_Richness)


