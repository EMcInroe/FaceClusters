###############################################################################################
## FACE FIT 
## Cleaner Version of PCA_explor.R
## Date Created: October 30, 2023
## Melissa McInroe
## Last Updated January 5, 2024
##
## This script uses penalized regression, boosting and best subset selection to find variables
## that will be used to crate clusters for the subjects, based on craniofacial measurements.
## The selected variables are then used in a PCA. 
## 99 Participants
## Updates include changing the box and whisker plots.
################################################################################################

library(data.table)
library(ggplot2)
library(corrplot)
library(stringr)
library(dplyr)
library(tidyr)
library(e1071)
library(factoextra)
library(cluster)
library(boot)
library(gbm)
library(leaps)
library(glmnet)
library(caret)
library(rstatix)
library(RColorBrewer)


###############################################################################
## Load Data
## Demographics
## Craniometric Measurements
## Mask Performance (FFE)
###############################################################################
setwd("")

demographics<-fread("R_scripts_input/demographics.csv")
names(demographics)
# "Study #"        "Date"           "Age"            "Sex"            "Race/Ethnicity" "Height (cm)"    "Weight (kg)"   
# "BMI" 

measurements<-fread("PythonOutput/total_table_scaled_8-14-2023.csv")
names(measurements)
# "subject"               "Condition"             "Bending Mean"          "Bending SD"            "Reading Mean"         
# "Reading SD"            "LR Mean"               "LR SD"                 "UD Mean"               "UD SD"                
# "Overall Mean"          "Overall SD"            "Bitragion Chin Arc"    "Bitragion Coronal Arc" "Bitragion Frontal Arc"
# "Head Circumference"    "Neck Circumference"    "Ear Breadth"           "Ear Length"            "Lip Length"           
# "Upper Facial Breadth"  "Menton-Sellion Length" "Nose Breadth"          "Nose Length"           "Bigonial Breadth"     
# "Bizygomatic Breadth"   "Head Breadth"          "Head Length"           "nose_gap_area"         "nose_gap_curve_length"
# "Bizyg_1_11"            "MenSel_6_89"           "EarLength_2_9"         "LipLength_59_65"       "sel-earb"             
# "nose-earb"             "chin-earb"             "lip-earb"              "sphenomaxillary_angle" "ratio_sel-nose_tip"   
# "ratio_sel-lip"         "ratio_sel-chin"  

demographics$subject<-str_sub(demographics$`Study #`, -3,-1)    ## Change subject number to match other datasets 
measurements$subject<-str_pad(measurements$subject,3,pad="0")
summary(measurements)
FFE<-fread("R_scripts_input/FFE.csv")
names(FFE)
# "Subject"       "Date"          "Sex"           "Age"           "BMI"           "N95"           "KN95"          "KN95_clip"    
# "Surgical"      "Surgical_clip" "KF94"          "KF94_clip"     "MKF94"         "MKF94_clip" 

## Calculate difference between mask with clip and mask without clip 
FFE$KN95_diff<-FFE$KN95_clip-FFE$KN95
FFE$Surgical_diff<-FFE$Surgical_clip-FFE$Surgical
FFE$KF94_diff<-FFE$KF94_clip-FFE$KF94
FFE$MKF94_diff<-FFE$MKF94_clip-FFE$MKF94

#merge all data
all_stats<-left_join(measurements, demographics, by="subject", multiple = "all")

#select only measurement to be included in PCA for exploratory analysis
#response variables excluded
select_stats<-data.table(all_stats[,c(1,13:29,41,45,46,48:50)], key="subject")
select_stats<-select_stats[,head(.SD,1), key(select_stats)]
names(select_stats)
# "subject"               "Bitragion Chin Arc"    "Bitragion Coronal Arc" "Bitragion Frontal Arc" "Head Circumference"   
# "Neck Circumference"    "Ear Breadth"           "Ear Length"            "Lip Length"            "Upper Facial Breadth" 
# "Menton-Sellion Length" "Nose Breadth"          "Nose Length"           "Bigonial Breadth"      "Bizygomatic Breadth"  
# "Head Breadth"          "Head Length"           "nose_gap_area"         "ratio_sel-lip"         "Age"                  
# "Sex"                   "Height (cm)"           "Weight (kg)"           "BMI"    

## Remove larger data sets that aren't needed

rm(demographics, measurements, all_stats)

## distributions
select_stats[,-c(1,21)] %>% gather() %>% head()
ggplot(gather(select_stats[,-c(1,21)]), aes(value)) +
  geom_histogram(bins = 10) +
  facet_wrap(~key, scales = "free_x")

# scale demographics
scale_stats<-data.frame(select_stats[,c(1,21)], apply(select_stats[,-c(1,21)], 2, scale))

scale_stats[,-c(1,2)] %>% gather() %>% head()
ggplot(gather(scale_stats[,-c(1,2)]), aes(value)) +
  geom_histogram(bins = 10) +
  facet_wrap(~key, scales = "free_x")

## Correlation
#remove row 39, subject 42 because they do not have nose_gap measurements (look into it)
cormatrix<-cor_mat(scale_stats[-39,-c(1,2,22,23)], method = "pearson")
cor_pvals<-cor_gather(cormatrix)
cor_pvals<-as.data.frame(cor_pvals)
corrplot(cor(scale_stats[-39,-c(1,2,22,23)], method = "pearson"),col = COL2('BrBG'),addCoef.col = "black",
         number.cex = 0.6,tl.cex = .7)
#write.csv(cor_pvals,"Craniometrics Paper/cor_matrix_data.csv")

scale_stats<-subset(scale_stats, !is.na(nose_gap_area))

## remove outliers
#outliers<-c("001", "003","009","011","025","032","081")
#scale_stats_subset<-subset(scale_stats, !(subject %in% outliers))
FFE$Subject<-str_sub(FFE$Subject, -3,-1)
stats_FFE<-left_join(scale_stats, FFE, by=c("subject"="Subject"))
stats_FFE<-stats_FFE[,-c(25:28)]

stats_ffe_cor<-na.omit(stats_FFE[,-c(3:6,9:13,15,17,18,20:24)])
corrplot(cor(stats_ffe_cor[,-c(1,2)], method = "pearson"), method = "color", type = "upper",
         diag = FALSE, addCoef.col = "black", tl.cex = .8, tl.col = "black",col = COL2('BrBG'),
         number.cex = 0.6)
################################################################################
## Boosting
## gbm package
################################################################################
##KN95
set.seed(5)
boost_mod<-gbm(KN95~., data = stats_FFE[-c(1,2,25,27:37)], distribution = "gaussian",
               n.trees = 5000, interaction.depth = 4)
summary(boost_mod)

stats_FFE_long<-pivot_longer(stats_FFE, cols = c(25,26,28,30), names_to ="Condition", values_to = "FFE")
set.seed(5)
boost_mod<-gbm(FFE~., data = stats_FFE_long[-c(1,2,21:34)], distribution = "gaussian",
               n.trees = 5000, interaction.depth = 4)
summary(boost_mod)

plot(boost_mod, i="Ear.Breadth")
plot(boost_mod, i="nose_gap_area")
plot(boost_mod, i="Nose.Breadth")

stats_FFE$Sex<-ifelse(stats_FFE$Sex =="F",1,0)
stats_FFE$Sex<-as.numeric(stats_FFE$Sex)
## KF94
boost_mod<-gbm(KF94~., data = stats_FFE[-c(1,2,25:29,31:37)], distribution = "gaussian",
               n.trees = 5000, interaction.depth = 4)
summary(boost_mod)

####  Demonstrate Boosting with Sex
# complete_stats$Sex<-as.factor(complete_stats$Sex)
# boost_mod<-gbm(KF94~., data = complete_stats[,-c(1, 10:16,20:24,26:32)], distribution = "gaussian",
#                n.trees = 5000, interaction.depth = 4)
# summary(boost_mod)

## Surgical
boost_mod<-gbm(Surgical~., data = stats_FFE[-c(1,2,25:27,29:37)], distribution = "gaussian",
               n.trees = 5000, interaction.depth = 4)
summary(boost_mod)

## MFK94
stat_mkf94<-na.omit(stats_FFE[-c(1,2,25:31,33:37)])
boost_mod<-gbm(MKF94~., data = stat_mkf94, distribution = "gaussian",
               n.trees = 5000, interaction.depth = 4)
summary(boost_mod)

#################################################################################
## Elastic Net Regression 
## glmnet package
#################################################################################

## Elastic
model<-train(KN95~., data=stats_FFE[-c(1,2,25,27:37)], method = "glmnet",
             trControl = trainControl("cv", number = 10),
             tuneLength = 10)
model$bestTune
# alpha   lambda
# 10   0.1 3.704994

model<-train(FFE~., data = stats_FFE_long[-c(1,2,21:34)], method = "glmnet",
             trControl = trainControl("cv", number = 10),
             tuneLength = 10)
model$bestTune

coef(model$finalModel, model$bestTune$lambda)
# s1
# (Intercept)           69.427926563
# Bitragion.Chin.Arc     0.868310870 
# Head.Circumference     0.250662861
# Neck.Circumference     1.100548082
# Ear.Breadth            1.847755691
# Lip.Length             1.694117054
# Upper.Facial.Breadth  -0.785949160
# Menton.Sellion.Length  1.015436809
# Nose.Breadth           0.562156594
# Nose.Length            0.673638269
# Bigonial.Breadth       0.430640694
# Bizygomatic.Breadth    0.934639906
# Head.Breadth           0.003943001 
# nose_gap_area         -3.011890165 
# Age.x                  0.404811635 
# BMI.x                 -0.598810227


#KF94
x<-model.matrix(KF94~., stats_FFE[-c(1,2,25:29,31:37)])[,-1]
model<-train(KF94~., data=stats_FFE[-c(1,2,25:29,31:37)], method = "glmnet",
             trControl = trainControl("cv", number = 10),
             tuneLength = 10)
model$bestTune
# alpha   lambda
# 0.1     3.727006
coef(model$finalModel, model$bestTune$lambda)

#Surgical
model<-train(Surgical~., data=stats_FFE[-c(1,2,25:27,29:37)], method = "glmnet",
             trControl = trainControl("cv", number = 10),
             tuneLength = 10)
model$bestTune
# alpha    lambda
# 67   0.7 0.3715513
coef(model$finalModel, model$bestTune$lambda)

#MKF94
stats_FFE_mkf94<-na.omit(stats_FFE[-c(1,2,25:31,33:37)])
model<-train(MKF94~., data=stats_FFE_mkf94, method = "glmnet",
             trControl = trainControl("cv", number = 10),
             tuneLength = 10)
model$bestTune
# alpha   lambda
# 39   0.4 1.654041
coef(model$finalModel, model$bestTune$lambda)


#################################################################################
## Best Subset Selection
## leaps package
#################################################################################

#KN95
regfit<-regsubsets(KN95~.,data = stats_FFE[-c(1,2,25,27:37)])
summary(regfit)
reg.sum<-summary(regfit)
reg.sum$rsq

par(mfrow=c(2,2))
plot(reg.sum$rss, xlab = "Number of Variables", ylab = "RSS", type = "l")
plot(reg.sum$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq", type = "l")
which.max(reg.sum$adjr2)
plot(reg.sum$cp, xlab = "Number of Variables", ylab = "Cp", type = "l")
plot(reg.sum$bic, xlab = "Number of Variables", ylab = "BIC", type = "l")

coef(regfit,4)
coef(regfit,6)

#KF94
regfit<-regsubsets(KF94~., stats_FFE[-c(1,2,25:29,31:37)])
summary(regfit)
reg.sum<-summary(regfit)
reg.sum$rsq

par(mfrow=c(2,2))
plot(reg.sum$rss, xlab = "Number of Variables", ylab = "RSS", type = "l")
plot(reg.sum$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq", type = "l")
plot(reg.sum$cp, xlab = "Number of Variables", ylab = "Cp", type = "l")
plot(reg.sum$bic, xlab = "Number of Variables", ylab = "BIC", type = "l")

coef(regfit,4)

# Surgical
regfit<-regsubsets(Surgical~., stats_FFE[-c(1,2,25:27,29:37)])
summary(regfit)
reg.sum<-summary(regfit)
reg.sum$rsq

par(mfrow=c(2,2))
plot(reg.sum$rss, xlab = "Number of Variables", ylab = "RSS", type = "l")
plot(reg.sum$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq", type = "l")
plot(reg.sum$cp, xlab = "Number of Variables", ylab = "Cp", type = "l")
plot(reg.sum$bic, xlab = "Number of Variables", ylab = "BIC", type = "l")

coef(regfit,5)
coef(regfit,6)

#MKF94
regfit<-regsubsets(MKF94~., stats_FFE_mkf94)
summary(regfit)
reg.sum<-summary(regfit)
reg.sum$rsq

par(mfrow=c(2,2))
plot(reg.sum$rss, xlab = "Number of Variables", ylab = "RSS", type = "l")
plot(reg.sum$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq", type = "l")
plot(reg.sum$cp, xlab = "Number of Variables", ylab = "Cp", type = "l")
plot(reg.sum$bic, xlab = "Number of Variables", ylab = "BIC", type = "l")

coef(regfit,5)
## All
regfit<-regsubsets(FFE~.,data = stats_FFE_long[-c(1,2,21:34)])
summary(regfit)
reg.sum<-summary(regfit)
reg.sum$rsq

par(mfrow=c(2,2))
plot(reg.sum$rss, xlab = "Number of Variables", ylab = "RSS", type = "l")
plot(reg.sum$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq", type = "l")
plot(reg.sum$cp, xlab = "Number of Variables", ylab = "Cp", type = "l")
plot(reg.sum$bic, xlab = "Number of Variables", ylab = "BIC", type = "l")

coef(regfit,4)
coef(regfit,5)

#######################################################################################
## Clustering with selected variables 
## Bizygotmatic Breath
## Ear Breadth
## nose length
## nose gap
## neck circumference
#######################################################################################
select_var<-c("Bizygomatic.Breadth","Nose.Length","nose_gap_area","Neck.Circumference","Ear.Breadth","Age","Height..cm.")
select_var_stats<-scale_stats[,c("subject", "Sex", select_var)]
colnames(select_var_stats)[9]<-"Height"

set.seed(91)
wss<-function(k) {
  kmeans(select_var_stats[,-c(1,2,8,9)], k, nstart=20)$tot.withinss
}
k.values<-1:10

wss_values<-purrr::map_dbl(k.values,wss)

plot(k.values,wss_values, type="b", pch = 19,
     frame = FALSE, xlab = "Number of Clusters",
     ylab = "Total within-cluster sum of squares")
## k-means
#k=4
set.seed(91)
km.out<-kmeans(select_var_stats[,-c(1,2,8,9)], 4, nstart=20)
km.out

col_pal<-c("#78160C","#C4830A","#164C45", "#508C9B")  
plot(select_var_stats[,-c(8,9)], col = col_pal[km.out$cluster], pch = 20, cex = 2)

#k=5
# set.seed(91)
# km.out<-kmeans(select_var_stats[,-c(1,2)], 5, nstart=20)
# km.out
# 
# plot(select_var_stats, col = (km.out$cluster +1), pch = 20, cex = 2)

FFE_subset<-subset(FFE, FFE$Subject %in% select_var_stats$subject)
FFE_subset$k_clust<-as.factor(km.out$cluster)
stats_FFE$k_cluster<-as.factor(km.out$cluster)

ggplot(FFE_subset, aes(x = KN95, y = KN95_clip, color = k_clust)) +
  geom_point(size = 3)

ggplot(FFE_subset, aes(x = KF94, y = KF94_clip, color = k_clust)) +
  geom_point(size = 3)

ggplot(FFE_subset, aes(x = Surgical, y = Surgical_clip, color = k_clust)) +
  geom_point(size = 3)

ggplot(FFE_subset, aes(x = MKF94, y = MKF94_clip, color = k_clust)) +
  geom_point(size = 3)

## make scatter plot of baseline and clip with clusters as colors.
############################### PCA ##################################
pr.out<-prcomp(select_var_stats[,-c(1,2,8,9)])
biplot(pr.out)

pv_dat<-as.data.frame(pr.out$x)

pr.var<-(pr.out$sdev)
eigenval<-pr.var^2
pve<-pr.var/sum(pr.var)

plot(pve, xlab="Principal Component", ylab="Proportion of Variance Explained", type = "b")
plot(cumsum(pve), xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", type="b")
cumsum(pve)

pv_dat$k_clust<-as.factor(km.out$cluster)

ggplot(pv_dat, aes(x = PC1, y=PC2, color = k_clust)) +
  geom_point(size = 3)

## Plot for poster
col_pal<-c("#78160C","#C4830A","#164C45", "#508C9B")        
ggplot(pv_dat, aes(x = PC1, y = PC2, color = k_clust)) +
  geom_point(size = 4) +
  scale_color_manual(values = col_pal) +
  labs(x = "PC1", y = "PC2", color = "Cluster") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title = element_text(face = "bold", size = 18),
        axis.text = element_text(face = "bold",size = 18),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 
ggsave("L:/PRIV/EPHD_CRB/FACEFIT/DATA/Craniometrics Paper/PCA_cluster_plot.png", width = 1400, height = 1000, units = "px")
#########################################################################
## scatterplots of classes
## FFE baseline vs FFE with clip
#########################################################################
KN95<-ggplot(FFE_subset, aes(x = KN95, y = KN95_clip, color = k_clust)) +
  geom_point(size = 4) +
  scale_color_manual(values = col_pal) +
  annotate(geom="text", x=80, y=30, label = "KN95", color = "black", size =12, fontface = "bold")+
  labs(title = "", x = "Baseline FFE", y = "FFE with Clip", color = "Cluster") +
  scale_x_continuous(limits = c(15,105), breaks =c(20,30,40,50,60,70,80,90,100), expand = c(0,0)) +
  scale_y_continuous(limits = c(15,105), breaks =c(20,30,40,50,60,70,80,90,100), expand = c(0,0)) +
  theme(axis.title.y = element_text(size = 18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "gray80"),
        legend.position = "top")

surgical<- ggplot(FFE_subset, aes(x = Surgical, y = Surgical_clip, color = k_clust)) +
  geom_point(size = 4) +
  scale_color_manual(values = col_pal) +
  annotate(geom="text", x=80, y=30, label = "Surgical", color = "black", size =12, fontface = "bold")+
  labs(title = "", x = "Baseline FFE", y = "FFE with Clip", color = "Cluster") +
  scale_x_continuous(limits = c(15,105), breaks =c(20,30,40,50,60,70,80,90,100), expand = c(0,0)) +
  scale_y_continuous(limits = c(15,105), breaks =c(20,30,40,50,60,70,80,90,100), expand = c(0,0)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "gray80"),
        legend.position = "top")

KF94<-ggplot(FFE_subset, aes(x = KF94, y = KF94_clip, color = k_clust)) +
  geom_point(size = 4) +
  scale_color_manual(values = col_pal) +
  annotate(geom="text", x=80, y=30, label = "KF94", color = "black", size =12, fontface = "bold")+
  labs(title = "", x = "Baseline FFE", y = "FFE with Clip", color = "Cluster") +
  scale_x_continuous(limits = c(15,105), breaks =c(20,30,40,50,60,70,80,90,100), expand = c(0,0)) +
  scale_y_continuous(limits = c(15,105), breaks =c(20,30,40,50,60,70,80,90,100), expand = c(0,0)) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        axis.text = element_text(size = 18, face = "bold"),
        panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "gray80"),
        legend.position = "top")



MKF94<-ggplot(FFE_subset, aes(x = MKF94, y = MKF94_clip, color = k_clust)) +
  geom_point(size = 4) +
  scale_color_manual(values = col_pal) +
  annotate(geom="text", x=80, y=30, label = "MKF94", color = "black", size =12, fontface = "bold")+
  labs(title = "", x = "Baseline FFE", y = "FFE with Clip", color = "Cluster") +
  scale_x_continuous(limits = c(15,105), breaks =c(20,30,40,50,60,70,80,90,100), expand = c(0,0)) +
  scale_y_continuous(limits = c(15,105), breaks =c(20,30,40,50,60,70,80,90,100), expand = c(0,0)) +
  theme(axis.title.x = element_text(size = 18, face = "bold"),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major = element_line(color = "gray80"),
        legend.position = "top")


ggpubr::ggarrange(KN95, surgical, KF94, MKF94, nrow =2, ncol = 2, widths = c(1,.85), common.legend = TRUE)

## sex breakdown
stat_per<-stats_FFE %>%
  group_by(k_cluster, Sex.x) %>%
  summarise(count = n())


############################################################################
## Box and Whisker plots
## separated by cluster and condition (unmodified, modified/with clip)
############################################################################
FFE_long<-pivot_longer(FFE_subset, col = 7:14, names_to="Condition", values_to="FFE")

#KN95 
FFE_long_KN95<-subset(FFE_long, Condition == "KN95" | Condition == "KN95_clip")

KN95<-ggplot(FFE_long_KN95, aes(x = k_clust, y = FFE, color = k_clust, fill = Condition)) +
  geom_point(aes(group = Condition), position = position_dodge(width = .5)) +
  geom_boxplot(color = "black", alpha = 0.7, width = 0.5) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c("gray10","gray70"), labels = c("Unmodified", "Modified")) +
  annotate(geom="text", x="2", y=30, label = "KN95", color = "black", size =5, fontface = "bold", hjust=0.05)+
  geom_hline(yintercept = 80, linetype = 2) +
  geom_hline(yintercept = 65, linetype = 2) +
  scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
     theme(panel.grid.major.x = element_blank(),
           panel.grid.minor.y = element_blank(),
           legend.position = "bottom",                                                                   #"bottom"
           #legend.text = element_text(size = 14),
           #legend.title = element_text(size = 14),
           axis.title.y = element_text(size = 18, face = "bold"),
           axis.text.y = element_text(size = 18, face = "bold"),
           axis.title.x = element_blank(),
           axis.text.x = element_blank(),
           axis.ticks.x = element_blank(),
           plot.margin = unit(c(0,0,0,0),"cm"),
           panel.background = element_blank(),
           panel.border = element_rect(color = "black", fill = NA)) +
  labs(y = "FFE (%)", color = "Cluster", fill = "Condition")

#Surgical
FFE_long_surgical<-subset(FFE_long, Condition == "Surgical" | Condition == "Surgical_clip")

surgical<-ggplot(FFE_long_surgical, aes(x = k_clust, y = FFE, color = k_clust, fill = Condition)) +
  geom_point(aes(group = Condition), position = position_dodge(width = .5)) +
  geom_boxplot(color = "black", alpha = 0.7, width = 0.5) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c("gray10","gray70"), labels = c("Unmodified", "Modified")) +
  annotate(geom="text", x="2", y=30, label = "Surgical", color = "black", size =5, fontface = "bold", hjust = 0.2)+
  geom_hline(yintercept = 80, linetype = 2) +
  geom_hline(yintercept = 65, linetype = 2) +
  scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom",                                                                   #"bottom"
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.border = element_rect(color = "black", fill = NA)) +
  labs(color = "Cluster", fill = "Condition")

#KF94
FFE_long_KF94<-subset(FFE_long, Condition == "KF94" | Condition == "KF94_clip")

KF94<-ggplot(FFE_long_KF94, aes(x = k_clust, y = FFE, color = k_clust, fill = Condition)) +
  geom_point(aes(group = Condition), position = position_dodge(width = .5)) +
  geom_boxplot(color = "black", alpha = 0.7, width = 0.5) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c("gray10","gray70"), labels = c("Unmodified", "Modified")) +
  annotate(geom="text", x="2", y=30, label = "KF94", color = "black", size =5, fontface = "bold", hjust=0.08)+
  geom_hline(yintercept = 80, linetype = 2) +
  geom_hline(yintercept = 65, linetype = 2) +
  scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size=18, face = "bold"),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "bottom",                                                                   #"bottom"
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = unit(c(0,0,0,0),"cm")) +
  labs(y = "FFE (%)", x = "", color = "Cluster", fill = "Condition")

#MKF94
FFE_long_MKF94<-subset(FFE_long, Condition == "MKF94" | Condition == "MKF94_clip")
FFE_long_MKF94<-subset(FFE_long_MKF94, !is.na(FFE))

MKF94<-ggplot(FFE_long_MKF94, aes(x = k_clust, y = FFE, color = k_clust, fill = Condition)) +
  geom_point(aes(group = Condition), position = position_dodge(width = .5)) +
  geom_boxplot(color = "black", alpha = 0.7, width = 0.5) +
  scale_color_manual(values = col_pal) +
  scale_fill_manual(values = c("gray10","gray70"), labels = c("Unmodified", "Modified")) +
  geom_hline(yintercept = 80, linetype = 2) +
  geom_hline(yintercept = 65, linetype = 2) +
  scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
  annotate(geom="text", x="2", y=30, label = "MKF94", color = "black", size =5, fontface = "bold", hjust = 0.15)+
  theme(legend.position = "bottom",                                                                   #"bottom"
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.border = element_rect(color = "black", fill = NA)) +
  labs(color = "Cluster", fill = "Condition")

ggpubr::ggarrange(KN95, surgical, KF94, MKF94, nrow =2, ncol = 2, widths = c(1,.85), common.legend = TRUE)

#################################################################################
## tukey of FFE by cluster
#################################################################################
FFE_long$Condition<-factor(FFE_long$Condition, levels = c("KN95","Surgical","KF94","MKF94","KN95_clip",
                                                          "Surgical_clip","KF94_clip","MKF94_clip"),
                           labels =c("KN95","Surgical","KF94","MKF94","Mod-KN95",
                                     "Mod-Surgical","Mod-KF94","Mod-MKF94"))

## cluster 1
model1<-lm(FFE~Condition , data = FFE_long, subset = c(k_clust == "1"))
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ Condition,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

clust_1<-ggplot(data = mod_means, aes(x = Condition, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, linewidth = 1.5, color = "#FFAE00") +
  geom_point(size = 4, color = "#FFAE00") +
  annotate(geom="text", x="Surgical", y=93, label = "Cluster 1", color = "black", size =8, fontface = "bold")+
  geom_text(aes(label = gsub(" ","", .group)),
            position = position_nudge(x = -0.3), color = "black", size = 8) +
  ylim(c(35,100)) +
  labs(y = "") +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 1.5),
        panel.grid.major.y = element_line(color = "gray80"),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.text = element_text(size=18, face = "bold"),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(colour="black"))

## cluster 2
model1<-lm(FFE~Condition , data = FFE_long, subset = c(k_clust == "2"))
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ Condition,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

clust_2<-ggplot(data = mod_means, aes(x = Condition, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, linewidth = 1.5, color = "#0C7D74") +
  geom_point(size = 4, color = "#0C7D74") +
  annotate(geom="text", x="Surgical", y=93, label = "Cluster 2", color = "black", size =8, fontface = "bold")+
  geom_text(aes(label = gsub(" ","", .group)),
            position = position_nudge(x = -0.3), color = "black", size = 8) +
  labs(y = "") +
  ylim(c(35,100)) +
  theme(axis.title = element_blank(),
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        # = element_line(color = "black", linewidth = 1.5),
        panel.grid.major.y = element_line(color = "gray80"),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.text.x = element_text(size=18, face = "bold", angle = 90),
        axis.text.y = element_text(size=18, face = "bold"),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(colour="black"))

## Cluster 3
model1<-lm(FFE~Condition , data = FFE_long, subset = c(k_clust == "3"))
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ Condition,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

clust_3<-ggplot(data = mod_means, aes(x = Condition, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, linewidth = 1.5, color = "#0F353D") +
  geom_point(size = 4, color = "#0F353D") +
  annotate(geom="text", x="Surgical", y=93, label = "Cluster 3", color = "black", size =8, fontface = "bold")+
  geom_text(aes(label = gsub(" ","", .group)),
            position = position_nudge(x = -0.4), color = "black", size = 8) +
  labs(y = "") +
  ylim(c(35,100)) +
  theme(axis.title = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_line(color = "black", linewidth = 1.5),
        panel.grid.major.y = element_line(color = "gray80"),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.text.x = element_text(size=18, face = "bold", angle = 90, hjust = 1),
        axis.title.y = element_blank(),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(colour="black"))

## cluster 4
model1<-lm(FFE~Condition , data = FFE_long, subset = c(k_clust == "4"))
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ Condition,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

clust_4<-ggplot(data = mod_means, aes(x = Condition, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, linewidth = 1.5, color = "#D55E00") +
  geom_point(size = 4, color = "#D55E00") +
  annotate(geom="text", x="Surgical", y=93, label = "Cluster 4", color = "black", size =8, fontface = "bold")+
  geom_text(aes(label = gsub(" ","", .group)),
            position = position_nudge(x = -0.2), color = "black", size = 8) +
  labs(y = "") +
  ylim(c(35,100)) +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text= element_blank(),
        axis.line = element_line(color = "black", linewidth = 1.5),
        panel.grid.major.y = element_line(color = "gray80"),
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank())

ggpubr::ggarrange(clust_1,clust_4,clust_2,clust_3, ncol =2, nrow = 2, heights = c(1.5,2))

###########################################################################
# Correlation with FFE
###########################################################################
all_stats_FFE<-left_join(select_stats, FFE, by = c("subject" = "Subject"))
all_stats_FFE<-all_stats_FFE[,-c(20:28)]

cormatrix<-cor_mat(all_stats_FFE[,-1], method = "pearson")
cor_pvals<-cor_gather(cormatrix)
cor_pvals<-as.data.frame(cor_pvals)
masks<-colnames(all_stats_FFE)[20:32]
cor_pvals<-subset(cor_pvals, var2 %in% masks)
cor_pvals<-subset(cor_pvals, !(var1 %in% masks))
#write.csv(cor_pvals,"Craniometrics Paper/cor_matrix_data.csv")

sel_vars<-c("Bizygomatic Breadth","Nose Length","nose_gap_area","Neck Circumference","Ear Breadth")

cor_pvals<-subset(cor_pvals, var1 %in% sel_vars & !(var2 %like% "_diff"))
cor_pvals$p<-round(cor_pvals$p,2)
cor_pvals$cor<-round(cor_pvals$cor,2)
cor_pvals$var2<-factor(cor_pvals$var2, levels =c("N95","KN95","KN95_clip",
                         "Surgical","Surgical_clip","KF94","KF94_clip", "MKF94",
                         "MKF94_clip"))

ggplot(cor_pvals, aes(x = var2, y = var1, fill = cor)) + 
  geom_tile() +
  scale_fill_gradientn(colors = c("chocolate4","lightyellow","turquoise4"))+
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 90, face = "bold", size = 15),
        axis.text.y = element_text(face = "bold", size = 15)) +
  geom_text(aes(label = paste("\np-value:\n",p))) +
  labs(fill = "Correlation") 
  




#############################################################################
## Summary Stats by Cluster
#############################################################################
summary(FFE_subset$k_clust)
summary(FFE_subset$k_clust[which(FFE_subset$Sex=="F")])
summary(FFE_subset$k_clust[which(FFE_subset$Sex=="M")])

##############################################################################
## Models
##############################################################################
# FFE_long$Clip<-as.factor(ifelse(str_sub(FFE_long$Condition,-4,-1)=="clip", 1,0))
# FFE_long$Condition<-as.factor(ifelse(str_sub(FFE_long$Condition,-4,-1)=="clip", gsub("_clip","", FFE_long$Condition), FFE_long$Condition))

FFE_long_clip<-FFE_long

FFE_long_clip$modified<-as.factor(ifelse(str_sub(FFE_long_clip$Condition,-4,-1)=="clip",1,0))
FFE_long_clip$Condition<-as.factor(ifelse(str_sub(FFE_long_clip$Condition,-4,-1)=="clip", 
                                          gsub("_clip","", FFE_long_clip$Condition), 
                                          FFE_long_clip$Condition))

model<-lm(FFE~Condition:k_clust, data = FFE_long_clip)
summary(model)
TukeyHSD(aov(model))
#--------------------------------------------------------------------------------------------------
## KN95
#baseline
model1<-lm(FFE~k_clust, data = FFE_long, subset=(Condition == "KN95"))
summary(model1)
TukeyHSD(aov(model1))

#clip
model2<-lm(FFE~k_clust, data = FFE_long, subset=(Condition == "KN95_clip"))
summary(model2)
TukeyHSD(aov(model2))
#------------------------------------------------------------------------------------------

## Surgical
#baseline
model3<-lm(FFE~k_clust, data = FFE_long, subset = (Condition == "Surgical"))
summary(model3)
TukeyHSD(aov(model3))

#clip
model4<-lm(FFE~k_clust, data = FFE_long, subset = (Condition == "Surgical_clip"))
summary(model4)
TukeyHSD(aov(model4))
#-------------------------------------------------------------------------------------------

## KF94
#baseline
model3<-lm(FFE~k_clust, data = FFE_long, subset = (Condition == "KF94"))
summary(model3)
TukeyHSD(aov(model3))

#clip
model3<-lm(FFE~k_clust, data = FFE_long, subset = (Condition == "KF94_clip"))
summary(model3)
TukeyHSD(aov(model3))


model4<-lm(FFE~k_clust, data = FFE_long, subset = (Condition == "MKF94"))
summary(model4)
TukeyHSD(aov(model4))

model4<-lm(FFE~k_clust, data = FFE_long, subset = (Condition == "MKF94_clip"))
summary(model4)
TukeyHSD(aov(model4))
#------------------------------------------------------------------------------------------

####################  Models By Cluster Comparing Conditions ######################
df_all<-data.frame(matrix(nrow = 28))
k_clusters<-levels(FFE_long$k_clust)

for(i in k_clusters) {
model5<-lm(FFE~Condition, data = FFE_long, subset = (k_clust == i))
summary(model5)
df<-as.data.frame(TukeyHSD(aov(model5))$Condition)
df<-setNames(df,paste0("K",i,"_",names(df)))
df_all<-cbind(df_all,df)
}

df_all<-df_all[,-1]

write.csv(df_all, "Craniometrics Paper/pairwise_tukey_mask_bycluster.csv")
##############################################################################
## Average Overall FFE By Cluster
##############################################################################
avg_FFE<-FFE_subset %>% group_by(k_clust) %>%
  summarise(N95 = mean(N95), 
            KN95 = mean(KN95),
            KN95_clip = mean(KN95_clip),
            Surgical = mean(Surgical),
            Surgical_clip = mean(Surgical_clip),
            KF94 = mean(KF94),
            KF94_clip  = mean(KF94_clip),
            MKF94 = mean(MKF94, na.rm = TRUE),
            MKF94_clip = mean(MKF94_clip, na.rm = TRUE))

summary(FFE_subset)
FFE_sd<-lapply(FFE_subset[which(FFE_subset$k_clust=="4"),6:14],sd, na.rm = TRUE)
FFE_sd

### Add tests to check significance of difference of means.
##############################################################################
## T-tests by clip and cluster
#############################################################################
## Overall
result<-by(FFE_long_clip, FFE_long_clip$Condition, function(x) 
  t.test(x$FFE[which(x$modified=="0")], x$FFE[which(x$modified=="1")], mu=0, alt="two.sided", paired = TRUE, var.equal == TRUE))

output<-type.convert(as.data.frame(do.call(rbind,result), as.is=TRUE))
print(output) 

## By Cluster
result<-by(FFE_long_clip[which(FFE_long_clip$k_clust=="1"),], FFE_long_clip$Condition[which(FFE_long_clip$k_clust=="1")], function(x) 
  t.test(x$FFE[which(x$modified=="0")], x$FFE[which(x$modified=="1")], mu=0, alt="two.sided", paired = TRUE, var.equal == TRUE))

output<-type.convert(as.data.frame(do.call(rbind,result), as.is=TRUE))
print(output) 

result<-by(FFE_long_clip[which(FFE_long_clip$k_clust=="2"),], FFE_long_clip$Condition[which(FFE_long_clip$k_clust=="2")], function(x) 
  t.test(x$FFE[which(x$modified=="0")], x$FFE[which(x$modified=="1")], mu=0, alt="two.sided", paired = TRUE, var.equal == TRUE))

output<-type.convert(as.data.frame(do.call(rbind,result), as.is=TRUE))
print(output) 

result<-by(FFE_long_clip[which(FFE_long_clip$k_clust=="3"),], FFE_long_clip$Condition[which(FFE_long_clip$k_clust=="3")], function(x) 
  t.test(x$FFE[which(x$modified=="0")], x$FFE[which(x$modified=="1")], mu=0, alt="two.sided", paired = TRUE, var.equal == TRUE))

output<-type.convert(as.data.frame(do.call(rbind,result), as.is=TRUE))
print(output) 

result<-by(FFE_long_clip[which(FFE_long_clip$k_clust=="4"),], FFE_long_clip$Condition[which(FFE_long_clip$k_clust=="4")], function(x) 
  t.test(x$FFE[which(x$modified=="0")], x$FFE[which(x$modified=="1")], mu=0, alt="two.sided", paired = TRUE, var.equal == TRUE))

output<-type.convert(as.data.frame(do.call(rbind,result), as.is=TRUE))
print(output) 

#############################################################################
## Average Raw Measurements by Cluster
#############################################################################
measurements<-fread("PythonOutput/total_table_8-14-2023.csv")
measurements<-as.data.frame(measurements)
#measurements<-na.omit(measurements)
summary(measurements)
lapply(measurements, sd)
measurements$subject<-str_pad(measurements$subject,3,pad="0")
measurements_all<-measurements[,-c(2:14,32:42,44)]
measurements_all<- measurements_all %>%
  distinct()
measurements_all$subject<-str_pad(measurements_all$subject,3,pad="0")
measurements_all<-subset(measurements_all, !is.na(nose_gap_area))
measurements_all$subject<-as.character(measurements_all$subject)
all_measures<-left_join(measurements_all, stats_FFE, by="subject")
write.csv(all_measures,"Craniometrics Paper/raw_scaled_measures_clusters_full.csv")

#all_measure<-as.data.frame(all_measures)
colnames(all_measure)[11]<-"Menton_sellion_length"
msl_mean<-all_measure %>% group_by(k_cluster) %>% 
  summarise(meanMSL = mean(Menton_sellion_length))

colnames(all_measures)[c(8,9,14)]<-c("ear.length.raw","lip.length.raw","bigonial.breadth.raw")
colnames(all_measures)[17]<-"head.length.raw"

ba_el_mean<-all_measures %>% group_by(k_cluster) %>%
  summarise(mean_el = mean(ear.length.raw),
            mean_ll = mean(lip.length.raw),
            mean_bb = mean(bigonial.breadth.raw),
            mean_hl = mean(head.length.raw))



select_var<-c("Bizygomatic Breadth","Nose Length","nose_gap_area","Neck Circumference","Ear Breadth")
select_var<-c("subject","Condition", select_var)
#measurements<-measurements[,c(1,2,15:31,43)]
measurements<-measurements[,c(select_var)]

# measure_clust<-left_join(measurements, FFE_subset, by = c("subject" = "Subject"))
# measure_clust<- measure_clust[,-2] %>%
#   distinct()
measurements$subject<-as.character(measurements$subject)

measure_clust<-left_join(measurements, stats_FFE[,-c(3:6,9:13,15,17,18,20)], by = "subject")
 measure_clust<- measure_clust[,-2] %>%
  distinct()  


measure_clust<-left_join(measure_clust, FFE_subset[,c(1,4,5,19)], by = c("subject"= "Subject"))
measure_clust<-subset(measure_clust, !is.na(measure_clust$k_clust))
colnames(measure_clust)[1:16]<-c("Subject","Bizygomatic.Breadth","Nose.Length","nose_gap_area","Neck.Circumference","Ear.Breadth",
                           "Sex","Neck.Circumference_scale","Ear.Breadth_scale","Nose.Length_scale", "Bizygomatic.Breadth_scale",
                           "nose_gap_area_scaled","Age_scale","Height_scale","Weight_scale", "BMI_scale")
write.csv(measure_clust,"Craniometrics Paper/raw_scaled_measures_clusters_select.csv")

measure_clust_avg<- measure_clust %>%
  group_by(k_clust) %>%
  summarise(Bizygomatic.Breadth = mean(Bizygomatic.Breadth, na.rm=TRUE),
            Nose.Length = mean(Nose.Length, na.rm = TRUE),
            nose_gap_area = mean(nose_gap_area, na.rm = TRUE),
            Neck.Circumference  = mean(Neck.Circumference, na.rm = TRUE),
            Ear.Breadth = mean(Ear.Breadth, na.rm = TRUE))

lapply(measure_clust[which(measure_clust$k_clust=="4"),3:7],sd, na.rm = TRUE)


#select_stats<-subset(select_stats, !(subject %in% outliers))
select_stats<-subset(select_stats, !is.na(nose_gap_area))
stats_clust<-left_join(select_stats, FFE_subset, by = c("subject" = "Subject"))

write.csv(stats_clust, "Craniometrics Paper/measurement_clusters_scaled.csv")

###################################################################################
## raw measurements significance test by cluster
###################################################################################
colnames(measure_clust)[2:6]<-c("Bizygomatic.Breadth.raw","Nose.Length.raw","nose_gap_area.raw",
                                "Neck.Circumference.raw","Ear.Breadth.raw")

model1<-lm(Bizygomatic.Breadth.raw~k_cluster , data = measure_clust)
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ k_cluster,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

BB<-ggplot(data = mod_means, aes(x = k_cluster, y = emmean, color = k_cluster)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.5) +
  geom_point(size = 4) +
  annotate(geom="text", x="1", y=142, label = "Bizygomatic Breadth", color = "black", size =6, fontface = "bold")+
  geom_text(aes(label = gsub(" ","", .group)),
            position = position_nudge(x = 0.2), color = "black", size = 8) +
  labs(#caption = "Means followed by a common letter are \nnot significantly different according to the Tukey-test",
       color = "Cluster", y = "") +
  scale_color_manual(values = c("#FFAE00" , "#0C7D74", "#0F353D","#D55E00")) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 1.5),
        panel.grid.major.y = element_line(color = "gray80"),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.text = element_text(size=18, face = "bold"),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(colour="black"),
        legend.position = "top")

model1<-lm(Nose.Length.raw~k_cluster , data = measure_clust)
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ k_cluster,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

NL<-ggplot(data = mod_means, aes(x = k_cluster, y = emmean, color = k_cluster)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.5) +
  geom_point(size = 4) +
  annotate(geom="text", x="1", y=58, label = "Nose Length", color = "black", size =6, fontface = "bold")+
  geom_text(aes(label = gsub(" ","", .group)),
            position = position_nudge(x = 0.2), color = "black", size = 8) +
  labs(#caption = "Means followed by a common letter are \nnot significantly different according to the Tukey-test",
    color = "Cluster", y = "") +
  scale_color_manual(values = c("#FFAE00" , "#0C7D74", "#0F353D","#D55E00")) + 
  ylim(c(45,60)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 1.5),
        panel.grid.major.y = element_line(color = "gray80"),
        plot.margin = unit(c(0,0,0,0.38),"cm"),
        axis.text.y = element_text(size=18, face = "bold"),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(colour="black"),
        legend.position = "top")

model1<-lm(nose_gap_area.raw~k_cluster , data = measure_clust)
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ k_cluster,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

NGA<-ggplot(data = mod_means, aes(x = k_cluster, y = emmean, color = k_cluster)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.5) +
  geom_point(size = 4) +
  annotate(geom="text", x="1", y=348, label = "Nose Gap Area", color = "black", size =6, fontface = "bold")+
  geom_text(aes(label = gsub(" ","", .group)),
            position = position_nudge(x = 0.2), color = "black", size = 8) +
  labs(#caption = "Means followed by a common letter are \nnot significantly different according to the Tukey-test",
    color = "Cluster", y = "") +
  scale_color_manual(values = c("#FFAE00" , "#0C7D74", "#0F353D","#D55E00")) +
  theme(axis.title = element_blank(),
        panel.grid.major.y = element_line(color = "gray80"),
        plot.margin = unit(c(0,0,0.2,0),"cm"),
        axis.text = element_text(size=18, face = "bold"),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(colour="black"),
        legend.position = "top")

model1<-lm(Neck.Circumference.raw~k_cluster , data = measure_clust)
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ k_cluster,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

NC<-ggplot(data = mod_means, aes(x = k_cluster, y = emmean, color = k_cluster)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.5) +
  geom_point(size = 4) +
  annotate(geom="text", x="1", y=415, label = "Neck Circumference", color = "black", size =6, fontface = "bold")+
  geom_text(aes(label = gsub(" ","", .group)),
            position = position_nudge(x = 0.2), color = "black", size = 8) +
  labs(#caption = "Means followed by a common letter are \nnot significantly different according to the Tukey-test",
    color = "Cluster", y = "") +
  scale_color_manual(values = c("#FFAE00" , "#0C7D74", "#0F353D","#D55E00")) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 1.5),
        panel.grid.major.y = element_line(color = "gray80"),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.text = element_text(size=18, face = "bold"),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(colour="black"),
        legend.position = "top")

model1<-lm(Ear.Breadth.raw~k_cluster , data = measure_clust)
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ k_cluster,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

EB<-ggplot(data = mod_means, aes(x = k_cluster, y = emmean, color = k_cluster)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.5) +
  geom_point(size = 4) +
  annotate(geom="text", x="1", y=38, label = "Ear Breadth", color = "black", size =6, fontface = "bold")+
  geom_text(aes(label = gsub(" ","", .group)),
            position = position_nudge(x = 0.2), color = "black", size = 8) +
  labs(#caption = "Means followed by a common letter are \nnot significantly different according to the Tukey-test",
    color = "Cluster", y = "") +
  scale_color_manual(values = c("#FFAE00" , "#0C7D74", "#0F353D","#D55E00")) +
  ylim(c(25,40)) +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 1.5),
        panel.grid.major.y = element_line(color = "gray80"),
        plot.margin = unit(c(0,0,0,0.38),"cm"),
        axis.text = element_text(size=18, face = "bold"),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(colour="black"),
        legend.position = "top")

ggpubr::ggarrange(BB,NC,NL,EB,NGA, ncol =1, common.legend = TRUE, legend = "top")

####################################################################################
## Test for multivariate normality
####################################################################################
library(energy)
mvnorm.etest(select_var_stats[,c(3:7)], R = 100)
############################################################################################
## Self Measurements
## Rerun above clustering, but with variables that are easier to take individually
############################################################################################
# select_var<-c("Bitragion.Chin.Arc","Bitragion.Coronal.Arc","nose_gap_area","Neck.Circumference")
# select_var_stats<-scale_stats_subset[,c("subject", "Sex", select_var)]
# 
# ## k-means
# #k=4
# set.seed(91)
# km.out<-kmeans(select_var_stats[,-c(1,2)], 4, nstart=20)
# km.out
# 
# plot(select_var_stats, col = (km.out$cluster +1), pch = 20, cex = 2)
# 
# #k=5
# # set.seed(91)
# # km.out<-kmeans(select_var_stats[,-c(1,2)], 5, nstart=20)
# # km.out
# # 
# # plot(select_var_stats, col = (km.out$cluster +1), pch = 20, cex = 2)
# 
# FFE_subset<-subset(FFE, FFE$Subject %in% select_var_stats$subject)
# FFE_subset$k_clust<-as.factor(km.out$cluster)
# 
# ggplot(FFE_subset, aes(x = KN95, y = KN95_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ggplot(FFE_subset, aes(x = KF94, y = KF94_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ggplot(FFE_subset, aes(x = Surgical, y = Surgical_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ggplot(FFE_subset, aes(x = MKF94, y = MKF94_clip, color = k_clust)) +
#   geom_point(size = 3)
# ## make scatter plot of baseline and clip with clusters as colors.
# ############################### PCA ##################################
# pr.out<-prcomp(select_var_stats[,-c(1,2)])
# biplot(pr.out)
# 
# pv_dat<-as.data.frame(pr.out$x)
# 
# pr.var<-(pr.out$sdev)
# eigenval<-pr.var^2
# pve<-pr.var/sum(pr.var)
# 
# plot(pve, xlab="Principal Component", ylab="Proportion of Variance Explained", type = "b")
# plot(cumsum(pve), xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", type="b")
# cumsum(pve)
# 
# pv_dat$k_clust<-as.factor(km.out$cluster)
# 
# ggplot(pv_dat, aes(x = PC1, y=PC2, color = k_clust)) +
#   geom_point(size = 3)
# #########################################################################
# ## scatterplots of classes
# ## FFE baseline vs FFE with clip
# #########################################################################
# ggplot(FFE_subset, aes(x = KN95, y = KN95_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ggplot(FFE_subset, aes(x = KF94, y = KF94_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ggplot(FFE_subset, aes(x = Surgical, y = Surgical_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ggplot(FFE_subset, aes(x = MKF94, y = MKF94_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ###########################################################################
# # Correlation with FFE
# ###########################################################################
# all_stats_FFE<-left_join(select_stats, FFE, by = c("subject" = "Subject"))
# all_stats_FFE<-all_stats_FFE[,-c(20:28)]
# 
# cormatrix<-cor_mat(all_stats_FFE[,-1], method = "pearson")
# cor_pvals<-cor_gather(cormatrix)
# cor_pvals<-as.data.frame(cor_pvals)
# masks<-colnames(all_stats_FFE)[20:32]
# cor_pvals<-subset(cor_pvals, var2 %in% masks)
# cor_pvals<-subset(cor_pvals, !(var1 %in% masks))
# #write.csv(cor_pvals,"Craniometrics Paper/cor_matrix_data.csv")
# 
# 
# ##############################################################################
# ## Average Overall FFE By Cluster
# ##############################################################################
# avg_FFE<-FFE_subset %>% group_by(k_clust) %>%
#   summarise(N95 = mean(N95), 
#             KN95 = mean(KN95),
#             KN95_clip = mean(KN95_clip),
#             Surgical = mean(Surgical),
#             Surgical_clip = mean(Surgical_clip),
#             KF94 = mean(KF94),
#             KF94_clip  = mean(KF94_clip),
#             MKF94 = mean(MKF94, na.rm = TRUE),
#             MKF94_clip = mean(MKF94_clip, na.rm = TRUE))
# 
# #############################################################################
# ## Average Raw Measurements by Cluster
# #############################################################################
# measurements<-fread("PythonOutput/total_table_8-14-2023.csv")
# measurements<-as.data.frame(measurements)
# measurements$subject<-str_pad(measurements$subject,3,pad="0")
# 
# select_var<-c("Bitragion Chin Arc","Bitragion Coronal Arc","nose_gap_area","Neck Circumference", "subject")
# 
# #measurements<-measurements[,c(1,2,15:31,43)]
# measurements<-measurements[,c(select_var)]
# 
# measure_clust<-left_join(measurements, FFE_subset, by = c("subject" = "Subject"))
# measure_clust<-subset(measure_clust, !is.na(measure_clust$k_clust))
# colnames(measure_clust)[c(1:5)]<-c("Bitragion.Chin.Arc","Bitragion.Coronal.Arc","nose_gap_area","Neck.Circumference","Subject")
# 
# measure_clust_avg<- measure_clust %>%
#   group_by(k_clust) %>%
#   summarise(Bitragion.Chin.Arc = mean(Bitragion.Chin.Arc, na.rm=TRUE),
#             Bitragion.Coronal.Arc = mean(Bitragion.Coronal.Arc, na.rm = TRUE),
#             nose_gap_area = mean(nose_gap_area, na.rm = TRUE),
#             Neck.Circumference  = mean(Neck.Circumference, na.rm = TRUE))
# sd(measure_clust$Bitragion.Chin.Arc[which(measure_clust$k_clust == "1")])
# sd(measure_clust$Bitragion.Chin.Arc[which(measure_clust$k_clust == "2")])
# sd(measure_clust$Bitragion.Chin.Arc[which(measure_clust$k_clust == "3")])
# sd(measure_clust$Bitragion.Chin.Arc[which(measure_clust$k_clust == "4")])
# sd(measure_clust$Bitragion.Coronal.Arc[which(measure_clust$k_clust == "1")])
# sd(measure_clust$Bitragion.Coronal.Arc[which(measure_clust$k_clust == "2")])
# sd(measure_clust$Bitragion.Coronal.Arc[which(measure_clust$k_clust == "3")])
# sd(measure_clust$Bitragion.Coronal.Arc[which(measure_clust$k_clust == "4")])
# sd(measure_clust$nose_gap_area[which(measure_clust$k_clust == "1")])
# sd(measure_clust$nose_gap_area[which(measure_clust$k_clust == "2")])
# sd(measure_clust$nose_gap_area[which(measure_clust$k_clust == "3")])
# sd(measure_clust$nose_gap_area[which(measure_clust$k_clust == "4")])
# sd(measure_clust$Neck.Circumference[which(measure_clust$k_clust == "1")])
# sd(measure_clust$Neck.Circumference[which(measure_clust$k_clust == "2")])
# sd(measure_clust$Neck.Circumference[which(measure_clust$k_clust == "3")])
# sd(measure_clust$Neck.Circumference[which(measure_clust$k_clust == "4")])
# 
# 
# select_stats<-subset(select_stats, !(subject %in% outliers))
# select_stats<-subset(select_stats, !is.na(nose_gap_area))
# stats_clust<-left_join(select_stats, FFE_subset, by = c("subject" = "Subject"))
# 
# 
# 
# 
# #########################################################
# ## New Select Variables
# ########################################################
# select_var<-c("Bitragion.Coronal.Arc","nose_gap_area","Neck.Circumference", "Ear.Breadth", "ratio_sel.lip")
# select_var_stats<-scale_stats_subset[,c("subject", "Sex", select_var)]
# 
# ## k-means
# #k=4
# set.seed(91)
# km.out<-kmeans(select_var_stats[,-c(1,2)], 4, nstart=20)
# km.out
# 
# plot(select_var_stats, col = (km.out$cluster +1), pch = 20, cex = 2)
# 
# #k=5
# # set.seed(91)
# # km.out<-kmeans(select_var_stats[,-c(1,2)], 5, nstart=20)
# # km.out
# # 
# # plot(select_var_stats, col = (km.out$cluster +1), pch = 20, cex = 2)
# 
# FFE_subset<-subset(FFE, FFE$Subject %in% select_var_stats$subject)
# FFE_subset$k_clust<-as.factor(km.out$cluster)
# 
# ggplot(FFE_subset, aes(x = KN95, y = KN95_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ggplot(FFE_subset, aes(x = KF94, y = KF94_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ggplot(FFE_subset, aes(x = Surgical, y = Surgical_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ggplot(FFE_subset, aes(x = MKF94, y = MKF94_clip, color = k_clust)) +
#   geom_point(size = 3)
# ## make scatter plot of baseline and clip with clusters as colors.
# ############################### PCA ##################################
# pr.out<-prcomp(select_var_stats[,-c(1,2)])
# biplot(pr.out)
# 
# pv_dat<-as.data.frame(pr.out$x)
# 
# pr.var<-(pr.out$sdev)
# eigenval<-pr.var^2
# pve<-pr.var/sum(pr.var)
# 
# plot(pve, xlab="Principal Component", ylab="Proportion of Variance Explained", type = "b")
# plot(cumsum(pve), xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", type="b")
# cumsum(pve)
# 
# pv_dat$k_clust<-as.factor(km.out$cluster)
# 
# ggplot(pv_dat, aes(x = PC1, y=PC2, color = k_clust)) +
#   geom_point(size = 3)
# #########################################################################
# ## scatterplots of classes
# ## FFE baseline vs FFE with clip
# #########################################################################
# ggplot(FFE_subset, aes(x = KN95, y = KN95_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ggplot(FFE_subset, aes(x = KF94, y = KF94_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ggplot(FFE_subset, aes(x = Surgical, y = Surgical_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ggplot(FFE_subset, aes(x = MKF94, y = MKF94_clip, color = k_clust)) +
#   geom_point(size = 3)
# 
# ###########################################################################
# # Correlation with FFE
# ###########################################################################
# all_stats_FFE<-left_join(select_stats, FFE, by = c("subject" = "Subject"))
# all_stats_FFE<-all_stats_FFE[,-c(20:28)]
# 
# cormatrix<-cor_mat(all_stats_FFE[,-1], method = "pearson")
# cor_pvals<-cor_gather(cormatrix)
# cor_pvals<-as.data.frame(cor_pvals)
# masks<-colnames(all_stats_FFE)[20:32]
# cor_pvals<-subset(cor_pvals, var2 %in% masks)
# cor_pvals<-subset(cor_pvals, !(var1 %in% masks))
# #write.csv(cor_pvals,"Craniometrics Paper/cor_matrix_data.csv")
# 
# 
# ##############################################################################
# ## Average Overall FFE By Cluster
# ##############################################################################
# avg_FFE<-FFE_subset %>% group_by(k_clust) %>%
#   summarise(N95 = mean(N95), 
#             KN95 = mean(KN95),
#             KN95_clip = mean(KN95_clip),
#             Surgical = mean(Surgical),
#             Surgical_clip = mean(Surgical_clip),
#             KF94 = mean(KF94),
#             KF94_clip  = mean(KF94_clip),
#             MKF94 = mean(MKF94, na.rm = TRUE),
#             MKF94_clip = mean(MKF94_clip, na.rm = TRUE))
# 
# #############################################################################
# ## Average Raw Measurements by Cluster
# #############################################################################
# measurements<-fread("PythonOutput/total_table_8-14-2023.csv")
# measurements<-as.data.frame(measurements)
# measurements$subject<-str_pad(measurements$subject,3,pad="0")
# 
# #select_var<-c( "subject","Bitragion Coronal Arc","nose_gap_area","Neck Circumference", "Ear Breadth", "ratio_sel-lip")
# 
# #measurements<-measurements[,c(1,2,15:31,43)]
# #measurements<-measurements[,c(select_var)]
# 
# measure_clust<-left_join(measurements, FFE_subset, by = c("subject" = "Subject"))
# #measure_clust<-subset(measure_clust, !is.na(measure_clust$k_clust))
# #colnames(measure_clust)[c(1:6)]<-c("Subject","Bitragion.Coronal.Arc","nose_gap_area","Neck.Circumference", "Ear.Breadth", "ratio_sel.lip")
# #measure_clust<-na.omit(measure_clust)
# 
# 
# measure_clust_avg<- measure_clust %>%
#   group_by(k_clust) %>%
#   summarise(Ear.Breadth = mean(Ear.Breadth, na.rm=TRUE),
#             Bitragion.Coronal.Arc = mean(Bitragion.Coronal.Arc, na.rm = TRUE),
#             nose_gap_area = mean(nose_gap_area, na.rm = TRUE),
#             Neck.Circumference  = mean(Neck.Circumference, na.rm = TRUE),
#             ratio_sel.lip  = mean(ratio_sel.lip, na.rm = TRUE))
# sd(measure_clust$Bitragion.Chin.Arc[which(measure_clust$k_clust == "1")])
# sd(measure_clust$Bitragion.Chin.Arc[which(measure_clust$k_clust == "2")])
# sd(measure_clust$Bitragion.Chin.Arc[which(measure_clust$k_clust == "3")])
# sd(measure_clust$Bitragion.Chin.Arc[which(measure_clust$k_clust == "4")])
# sd(measure_clust$Bitragion.Coronal.Arc[which(measure_clust$k_clust == "1")])
# sd(measure_clust$Bitragion.Coronal.Arc[which(measure_clust$k_clust == "2")])
# sd(measure_clust$Bitragion.Coronal.Arc[which(measure_clust$k_clust == "3")])
# sd(measure_clust$Bitragion.Coronal.Arc[which(measure_clust$k_clust == "4")])
# sd(measure_clust$nose_gap_area[which(measure_clust$k_clust == "1")])
# sd(measure_clust$nose_gap_area[which(measure_clust$k_clust == "2")])
# sd(measure_clust$nose_gap_area[which(measure_clust$k_clust == "3")])
# sd(measure_clust$nose_gap_area[which(measure_clust$k_clust == "4")])
# sd(measure_clust$Neck.Circumference[which(measure_clust$k_clust == "1")])
# sd(measure_clust$Neck.Circumference[which(measure_clust$k_clust == "2")])
# sd(measure_clust$Neck.Circumference[which(measure_clust$k_clust == "3")])
# sd(measure_clust$Neck.Circumference[which(measure_clust$k_clust == "4")])
# 
# 
# select_stats<-subset(select_stats, !(subject %in% outliers))
# select_stats<-subset(select_stats, !is.na(nose_gap_area))
# stats_clust<-left_join(select_stats, FFE_subset, by = c("subject" = "Subject"))


stats_FFE[,c(7:9,12:16,19,38)] %>% group_by(k_cluster) %>% 
  summarise(neck.circum= mean(Neck.Circumference),
            ear.breadth = mean(Ear.Breadth),
            ear.length = mean(Ear.Length),
            msl = mean(Menton.Sellion.Length),
            nose.breadth = mean(Nose.Breadth),
            nose.length = mean(Nose.Length),
            bigonial.breadth = mean(Bigonial.Breadth),
            bizyg.breadth = mean(Bizygomatic.Breadth),
            nga = mean(nose_gap_area))
