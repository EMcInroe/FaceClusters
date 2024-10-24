###############################################################################################
## FACE FIT 
## Cleaner Version of PCA_explor.R
## Date Original Created: December 13, 2023
## Melissa McInroe
## Updated February 2, 2024 
##
## This script continues from  "L:\PRIV\EPHD_CRB\FACEFIT\DATA\Craniometrics Paper\PCA_cluster_selectvars_v2.R"
## to create plots and summmaries of the data.  
## 99 Participants
## Data Input
## Data Output
## Update includes altering plots to fit journal submission
## Updated version changes cluster numbers (1-4) to shape assignments (D, P, R, T)
## Updated version changes asthetics of FFE boxplots by cluster.
################################################################################################
## Load Packages

library(data.table)
library(ggplot2)
library(corrplot)
library(stringr)
library(dplyr)
library(tidyr)
library(cluster)
library(leaps)
library(rstatix)
library(RColorBrewer)
library(MASS)
library(tree)
library(ggpattern)
library(grid)
library(gridExtra)


###############################################################################
## Load Data
## Demographics
## Craniometric Measurements
## Mask Performance (FFE)
###############################################################################
setwd("L:/PRIV/EPHD_CRB/FACEFIT/DATA")

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

# scale demographics
scale_stats<-data.frame(select_stats[,c(1,21)], apply(select_stats[,-c(1,21)], 2, scale))

scale_stats<-subset(scale_stats, !is.na(nose_gap_area))

#scale_stats_subset<-subset(scale_stats, !(subject %in% outliers))
FFE$Subject<-str_sub(FFE$Subject, -3,-1)
stats_FFE<-left_join(scale_stats, FFE, by=c("subject"="Subject"))
stats_FFE<-stats_FFE[,-c(25:28)]

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

#k=4
set.seed(91)
km.out<-kmeans(select_var_stats[,-c(1,2,8,9)], 4, nstart=20)
km.out

col_pal<-c("#D55E00","#FFAE00","#0C7D74","#0F353D")  
plot(select_var_stats[,-c(8,9)], col = col_pal[km.out$cluster], pch = 20, cex = 2)


FFE_subset<-subset(FFE, FFE$Subject %in% select_var_stats$subject)
FFE_subset$k_clust<-as.factor(km.out$cluster)
stats_FFE$k_cluster<-as.factor(km.out$cluster)


############################################################################
## Box and Whisker plots
## separated by cluster and condition (unmodified, modified/with clip)
## clusters as shapes and colors

## with pattern and with fill/no fill
############################################################################
# ### Pattern boxes
# FFE_long<-pivot_longer(FFE_subset, col = 7:14, names_to="Condition", values_to="FFE")
# FFE_long$k_clust<-factor(FFE_long$k_clust, labels = c("D","P","R","T"))
# 
# #KN95 
# FFE_long_KN95<-subset(FFE_long, Condition == "KN95" | Condition == "KN95_clip")
# 
# KN95<-ggplot(FFE_long_KN95, aes(x = k_clust, y = FFE, fill = interaction(k_clust,Condition), dodge = Condition)) +                     #shape = k_clust
#   geom_boxplot_pattern(pattern = c("none","circle","none","circle","none","circle","none","circle"),
#                        pattern_fill = "white",
#                        pattern_density = 0.35) +
#   scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D","#D55E00","#FFAE00","#0C7D74","#0F353D"), 
#                     labels = c("D Unmodified","P Unmodified","R Unmodified","T Unmodified", 
#                                "d Modified","P Modified","R Modified","T Modified"),
#                     name = "") +                                      #"gray10","gray70"
#                      #guide = guide_legend(override.aes = list(size = .5))) +
#   annotate(geom="text", x="P", y=30, label = "KN95", color = "black", size =5, fontface = "bold", hjust=0.05)+
#   geom_hline(yintercept = 80, linetype = 2) +
#   geom_hline(yintercept = 65, linetype = 2) +
#   scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
#      theme(panel.grid.major.x = element_blank(),
#            panel.grid.minor.y = element_blank(),
#            legend.position = "none",                                                                   #"bottom"
#            #legend.text = element_text(size = 14),
#            #legend.title = element_text(size = 14),
#            axis.title.y = element_text(size = 18, face = "bold"),
#            axis.text.y = element_text(size = 18, face = "bold"),
#            axis.title.x = element_blank(),
#            axis.text.x = element_blank(),
#            axis.ticks.x = element_blank(),
#            plot.margin = unit(c(0,0,0,0),"cm"),
#            panel.background = element_blank(),
#            panel.border = element_rect(color = "black", fill = NA)) +
#   labs(y = "FFE (%)")
# 
# #Surgical
# FFE_long_surgical<-subset(FFE_long, Condition == "Surgical" | Condition == "Surgical_clip")
# 
# surgical<-ggplot(FFE_long_surgical, aes(x = k_clust, y = FFE, fill = interaction(k_clust,Condition), dodge = Condition)) +                     #shape = k_clust
#   geom_boxplot_pattern(pattern = c("none","circle","none","circle","none","circle","none","circle"),
#                        pattern_fill = "white",
#                        pattern_density = 0.35) +
#   scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D","#D55E00","#FFAE00","#0C7D74","#0F353D"), 
#                     labels = c("D Unmodified","P Unmodified","R Unmodified","T Unmodified", 
#                                "d Modified","P Modified","R Modified","T Modified"),
#                     name = "") +    
#   # scale_fill_manual(values = c("gray10","gray70"), labels = c("Unmodified", "Modified"),
#   #                   guide = guide_legend(override.aes = list(size = .5))) +
#   annotate(geom="text", x="P", y=30, label = "Surgical", color = "black", size =5, fontface = "bold", hjust = 0.2)+
#   geom_hline(yintercept = 80, linetype = 2) +
#   geom_hline(yintercept = 65, linetype = 2) +
#   scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         legend.position = "none",                                                                   #"bottom"
#         #legend.text = element_text(size = 14),
#         #legend.title = element_text(size = 14),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.line.y.right = element_line("black"),
#         panel.background = element_blank(),
#         plot.margin = unit(c(0,0,0,0),"cm"),
#         panel.border = element_rect(color = "black", fill = NA))
# 
# #KF94
# FFE_long_KF94<-subset(FFE_long, Condition == "KF94" | Condition == "KF94_clip")
# 
# KF94<-ggplot(FFE_long_KF94, aes(x = k_clust, y = FFE, fill = interaction(k_clust,Condition), dodge = Condition)) +                     #shape = k_clust
#   geom_boxplot_pattern(pattern = c("none","circle","none","circle","none","circle","none","circle"),
#                        pattern_fill = "white",
#                        pattern_density = 0.35) +
#   scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D","#D55E00","#FFAE00","#0C7D74","#0F353D"), 
#                     labels = c("D Unmodified","P Unmodified","R Unmodified","T Unmodified", 
#                                "d Modified","P Modified","R Modified","T Modified"),
#                     name = "") +  
#   annotate(geom="text", x="P", y=30, label = "KF94", color = "black", size =5, fontface = "bold", hjust=0.08)+
#   geom_hline(yintercept = 80, linetype = 2) +
#   geom_hline(yintercept = 65, linetype = 2) +
#   scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
#   scale_x_discrete(labels = c("D"="D\u2666","P"="P\u2BC2","R"="R\u25AA","T"="T\u2BC6")) +
#   #scale_x_discrete((labels = c("1 \u2666","2 \u2B1F","3 \u25AA","4 \u2BC6"))) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         axis.text.x = element_text(size=18, face = "bold", colour = c("#D55E00","#FFAE00","#0C7D74","#0F353D")),
#         axis.text.y = element_text(size = 18, face = "bold"),
#         axis.title.y = element_text(size = 18, face = "bold"),
#         axis.title.x = element_blank(),
#         legend.position = "none",                                                                   #"bottom"
#         #legend.text = element_text(size = 14),
#         #legend.title = element_text(size = 14),
#         panel.background = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA),
#         plot.margin = unit(c(0,0,0,0),"cm")) +
#   labs(y = "FFE (%)", x = "")
# 
# #MKF94
# FFE_long_MKF94<-subset(FFE_long, Condition == "MKF94" | Condition == "MKF94_clip")
# FFE_long_MKF94<-subset(FFE_long_MKF94, !is.na(FFE))
# 
# MKF94<-ggplot(FFE_long_MKF94, aes(x = k_clust, y = FFE, fill = interaction(k_clust,Condition), dodge = Condition)) +                     #shape = k_clust
#   geom_boxplot_pattern(pattern = c("none","circle","none","circle","none","circle","none","circle"),
#                        pattern_fill = "white",
#                        pattern_density = 0.35) +
#   scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D","#D55E00","#FFAE00","#0C7D74","#0F353D"), 
#                     labels = c("D Unmodified","P Unmodified","R Unmodified","T Unmodified", 
#                                "d Modified","P Modified","R Modified","T Modified"),
#                     name = "") +  
#   geom_hline(yintercept = 80, linetype = 2) +
#   geom_hline(yintercept = 65, linetype = 2) +
#   scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
#   annotate(geom="text", x="P", y=30, label = "MKF94", color = "black", size =5, fontface = "bold", hjust = 0.15)+
#   scale_x_discrete(labels = c("D"="D\u2666","P"="P\u2BC2","R"="R\u25AA","T"="T\u2BC6")) +
#   theme(legend.position = "none",                                                                   #"bottom"
#         #legend.text = element_text(size = 14),
#         #legend.title = element_text(size = 14),
#         axis.text.x = element_text(size = 18, face = "bold", colour = c("#D55E00","#FFAE00","#0C7D74","#0F353D")),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.line.y.right = element_line("black"),
#         panel.background = element_blank(),
#         plot.margin = unit(c(0,0,0,0),"cm"),
#         panel.border = element_rect(color = "black", fill = NA))
# 
# ggpubr::ggarrange(KN95, surgical, KF94, MKF94, nrow =2, ncol = 2, widths = c(1,.85),heights = c(1.8,2))


### Outline Boxes
FFE_long_KN95<-subset(FFE_long, Condition == "KN95" | Condition == "KN95_clip")

KN95<-ggplot(FFE_long_KN95, aes(x = k_clust, y = FFE, fill = interaction(k_clust,Condition), color = interaction(k_clust,Condition), dodge = Condition)) +                     #shape = k_clust
  geom_boxplot(position = position_dodge(width = .6),width = .4, size = 1.5) +
  scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D","white","white","white","white"), 
                    labels = c("D Unmodified","P Unmodified","R Unmodified","T Unmodified", 
                               "d Modified","P Modified","R Modified","T Modified"),
                    name = "") + 
  scale_color_manual(values = c("black","black","black","black","#D55E00","#FFAE00","#0C7D74","#0F353D")) +
  #guide = guide_legend(override.aes = list(size = .5))) +
  annotate(geom="text", x="P", y=30, label = "KN95", color = "black", size =5, fontface = "bold", hjust=0.05)+
  geom_hline(yintercept = 80, linetype = 2) +
  geom_hline(yintercept = 65, linetype = 2) +
  scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",                                                                   #"bottom"
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
  labs(y = "FFE (%)")

#Surgical
FFE_long_surgical<-subset(FFE_long, Condition == "Surgical" | Condition == "Surgical_clip")

surgical<-ggplot(FFE_long_surgical, aes(x = k_clust, y = FFE, fill = interaction(k_clust,Condition), color = interaction(k_clust,Condition), dodge = Condition)) +                     #shape = k_clust
  geom_boxplot(position = position_dodge(width = .6),width = .4, size = 1.5) +
  scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D","white","white","white","white"), 
                    labels = c("D Unmodified","P Unmodified","R Unmodified","T Unmodified", 
                               "d Modified","P Modified","R Modified","T Modified"),
                    name = "") + 
  scale_color_manual(values = c("black","black","black","black","#D55E00","#FFAE00","#0C7D74","#0F353D")) +    
  # scale_fill_manual(values = c("gray10","gray70"), labels = c("Unmodified", "Modified"),
  #                   guide = guide_legend(override.aes = list(size = .5))) +
  annotate(geom="text", x="P", y=30, label = "Surgical", color = "black", size =5, fontface = "bold", hjust = 0.2)+
  geom_hline(yintercept = 80, linetype = 2) +
  geom_hline(yintercept = 65, linetype = 2) +
  scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "none",                                                                   #"bottom"
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y.right = element_line("black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.border = element_rect(color = "black", fill = NA))

#KF94
FFE_long_KF94<-subset(FFE_long, Condition == "KF94" | Condition == "KF94_clip")

KF94<-ggplot(FFE_long_KF94, aes(x = k_clust, y = FFE, fill = interaction(k_clust,Condition), color = interaction(k_clust,Condition), dodge = Condition)) +                     #shape = k_clust
  geom_boxplot(position = position_dodge(width = .6),width = .4, size = 1.5)  +
  scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D","white","white","white","white"), 
                    labels = c("D Unmodified","P Unmodified","R Unmodified","T Unmodified", 
                               "d Modified","P Modified","R Modified","T Modified"),
                    name = "") + 
  scale_color_manual(values = c("black","black","black","black","#D55E00","#FFAE00","#0C7D74","#0F353D")) +  
  annotate(geom="text", x="P", y=30, label = "KF94", color = "black", size =5, fontface = "bold", hjust=0.08)+
  geom_hline(yintercept = 80, linetype = 2) +
  geom_hline(yintercept = 65, linetype = 2) +
  scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
  scale_x_discrete(labels = c("D"="D\u2666","P"="P\u2BC2","R"="R\u25AA","T"="T\u2BC6")) +
  #scale_x_discrete((labels = c("1 \u2666","2 \u2B1F","3 \u25AA","4 \u2BC6"))) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(size=18, face = "bold", colour = c("#D55E00","#FFAE00","#0C7D74","#0F353D")),
        axis.text.y = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "none",                                                                   #"bottom"
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14),
        panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        plot.margin = unit(c(0,0,0,0),"cm")) +
  labs(y = "FFE (%)", x = "")

#MKF94
FFE_long_MKF94<-subset(FFE_long, Condition == "MKF94" | Condition == "MKF94_clip")
FFE_long_MKF94<-subset(FFE_long_MKF94, !is.na(FFE))

MKF94<-ggplot(FFE_long_MKF94, aes(x = k_clust, y = FFE, fill = interaction(k_clust,Condition), color = interaction(k_clust,Condition), dodge = Condition)) +                     #shape = k_clust
  geom_boxplot(position = position_dodge(width = .6),width = .4, size = 1.5)  +
  scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D","white","white","white","white"), 
                    labels = c("D Unmodified","P Unmodified","R Unmodified","T Unmodified", 
                               "d Modified","P Modified","R Modified","T Modified"),
                    name = "") + 
  scale_color_manual(values = c("black","black","black","black","#D55E00","#FFAE00","#0C7D74","#0F353D")) +  
  geom_hline(yintercept = 80, linetype = 2) +
  geom_hline(yintercept = 65, linetype = 2) +
  scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
  annotate(geom="text", x="P", y=30, label = "MKF94", color = "black", size =5, fontface = "bold", hjust = 0.15)+
  scale_x_discrete(labels = c("D"="D\u2666","P"="P\u2BC2","R"="R\u25AA","T"="T\u2BC6")) +
  theme(legend.position = "none",                                                                   #"bottom"
        #legend.text = element_text(size = 14),
        #legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 18, face = "bold", colour = c("#D55E00","#FFAE00","#0C7D74","#0F353D")),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y.right = element_line("black"),
        panel.background = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"),
        panel.border = element_rect(color = "black", fill = NA))

ggpubr::ggarrange(KN95, surgical, KF94, MKF94, nrow =2, ncol = 2, widths = c(1,.85),heights = c(1.8,2))

## black and white legend

blk_white<-ggplot(FFE_long_MKF94, aes(x = k_clust, y = FFE, fill = Condition)) +                     
  geom_boxplot(position = position_dodge(width = .6),width = .4)  +
  scale_fill_manual(values = c("black","white"), 
                    labels = c("Unmodified","Modified"),
                    name = "") +
  theme(legend.text = element_text(size = 22),
        legend.title = element_text(size = 22),
        legend.key.size = unit(3, 'cm'),
        legend.key = element_rect(colour = "transparent", fill = "white"),)

legend<-cowplot::get_legend(blk_white)
grid.newpage()
grid.draw(legend)


# ### Outline Boxes Alt
# FFE_long_KN95<-subset(FFE_long, Condition == "KN95" | Condition == "KN95_clip")
# 
# KN95<-ggplot(FFE_long_KN95, aes(x = k_clust, y = FFE, fill = interaction(k_clust,Condition), color = interaction(k_clust,Condition), dodge = Condition)) +                     #shape = k_clust
#   geom_boxplot(width = .5) +
#   scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D","gray","gray","gray","gray"), 
#                     labels = c("D Unmodified","P Unmodified","R Unmodified","T Unmodified", 
#                                "d Modified","P Modified","R Modified","T Modified"),
#                     name = "") + 
#   scale_color_manual(values = c("black","black","black","black","#D55E00","#FFAE00","#0C7D74","#0F353D")) +
#   #guide = guide_legend(override.aes = list(size = .5))) +
#   annotate(geom="text", x="P", y=30, label = "KN95", color = "black", size =5, fontface = "bold", hjust=0.05)+
#   geom_hline(yintercept = 80, linetype = 2) +
#   geom_hline(yintercept = 65, linetype = 2) +
#   scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         legend.position = "none",                                                                   #"bottom"
#         #legend.text = element_text(size = 14),
#         #legend.title = element_text(size = 14),
#         axis.title.y = element_text(size = 18, face = "bold"),
#         axis.text.y = element_text(size = 18, face = "bold"),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         plot.margin = unit(c(0,0,0,0),"cm"),
#         panel.background = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA)) +
#   labs(y = "FFE (%)")
# 
# #Surgical
# FFE_long_surgical<-subset(FFE_long, Condition == "Surgical" | Condition == "Surgical_clip")
# 
# surgical<-ggplot(FFE_long_surgical, aes(x = k_clust, y = FFE, fill = interaction(k_clust,Condition), color = interaction(k_clust,Condition), dodge = Condition)) +                     #shape = k_clust
#   geom_boxplot(width = .5) +
#   scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D","gray","gray","gray","gray"), 
#                     labels = c("D Unmodified","P Unmodified","R Unmodified","T Unmodified", 
#                                "d Modified","P Modified","R Modified","T Modified"),
#                     name = "") + 
#   scale_color_manual(values = c("black","black","black","black","#D55E00","#FFAE00","#0C7D74","#0F353D")) +    
#   # scale_fill_manual(values = c("gray10","gray70"), labels = c("Unmodified", "Modified"),
#   #                   guide = guide_legend(override.aes = list(size = .5))) +
#   annotate(geom="text", x="P", y=30, label = "Surgical", color = "black", size =5, fontface = "bold", hjust = 0.2)+
#   geom_hline(yintercept = 80, linetype = 2) +
#   geom_hline(yintercept = 65, linetype = 2) +
#   scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         legend.position = "none",                                                                   #"bottom"
#         #legend.text = element_text(size = 14),
#         #legend.title = element_text(size = 14),
#         axis.title = element_blank(),
#         axis.text = element_blank(),
#         axis.ticks = element_blank(),
#         axis.line.y.right = element_line("black"),
#         panel.background = element_blank(),
#         plot.margin = unit(c(0,0,0,0),"cm"),
#         panel.border = element_rect(color = "black", fill = NA))
# 
# #KF94
# FFE_long_KF94<-subset(FFE_long, Condition == "KF94" | Condition == "KF94_clip")
# 
# KF94<-ggplot(FFE_long_KF94, aes(x = k_clust, y = FFE, fill = interaction(k_clust,Condition), color = interaction(k_clust,Condition), dodge = Condition)) +                     #shape = k_clust
#   geom_boxplot(width = .5) +
#   scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D","gray","gray","gray","gray"), 
#                     labels = c("D Unmodified","P Unmodified","R Unmodified","T Unmodified", 
#                                "d Modified","P Modified","R Modified","T Modified"),
#                     name = "") + 
#   scale_color_manual(values = c("black","black","black","black","#D55E00","#FFAE00","#0C7D74","#0F353D")) +  
#   annotate(geom="text", x="P", y=30, label = "KF94", color = "black", size =5, fontface = "bold", hjust=0.08)+
#   geom_hline(yintercept = 80, linetype = 2) +
#   geom_hline(yintercept = 65, linetype = 2) +
#   scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
#   scale_x_discrete(labels = c("D"="D\u2666","P"="P\u2BC2","R"="R\u25AA","T"="T\u2BC6")) +
#   #scale_x_discrete((labels = c("1 \u2666","2 \u2B1F","3 \u25AA","4 \u2BC6"))) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.y = element_blank(),
#         axis.text.x = element_text(size=18, face = "bold", colour = c("#D55E00","#FFAE00","#0C7D74","#0F353D")),
#         axis.text.y = element_text(size = 18, face = "bold"),
#         axis.title.y = element_text(size = 18, face = "bold"),
#         axis.title.x = element_blank(),
#         legend.position = "none",                                                                   #"bottom"
#         #legend.text = element_text(size = 14),
#         #legend.title = element_text(size = 14),
#         panel.background = element_blank(),
#         panel.border = element_rect(color = "black", fill = NA),
#         plot.margin = unit(c(0,0,0,0),"cm")) +
#   labs(y = "FFE (%)", x = "")
# 
# #MKF94
# FFE_long_MKF94<-subset(FFE_long, Condition == "MKF94" | Condition == "MKF94_clip")
# FFE_long_MKF94<-subset(FFE_long_MKF94, !is.na(FFE))
# 
# MKF94<-ggplot(FFE_long_MKF94, aes(x = k_clust, y = FFE, fill = interaction(k_clust,Condition), color = interaction(k_clust,Condition), dodge = Condition)) +                     #shape = k_clust
#   geom_boxplot(width = .5) +
#   scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D","gray","gray","gray","gray"), 
#                     labels = c("D Unmodified","P Unmodified","R Unmodified","T Unmodified", 
#                                "d Modified","P Modified","R Modified","T Modified"),
#                     name = "") + 
#   scale_color_manual(values = c("black","black","black","black","#D55E00","#FFAE00","#0C7D74","#0F353D")) +  
#   geom_hline(yintercept = 80, linetype = 2) +
#   geom_hline(yintercept = 65, linetype = 2) +
#   scale_y_continuous(limits = c(25,105), breaks =c(20,30,40,50,60,70,80,90,100)) +
#   annotate(geom="text", x="P", y=30, label = "MKF94", color = "black", size =5, fontface = "bold", hjust = 0.15)+
#   scale_x_discrete(labels = c("D"="D\u2666","P"="P\u2BC2","R"="R\u25AA","T"="T\u2BC6")) +
#   theme(legend.position = "none",                                                                   #"bottom"
#         #legend.text = element_text(size = 14),
#         #legend.title = element_text(size = 14),
#         axis.text.x = element_text(size = 18, face = "bold", colour = c("#D55E00","#FFAE00","#0C7D74","#0F353D")),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.line.y.right = element_line("black"),
#         panel.background = element_blank(),
#         plot.margin = unit(c(0,0,0,0),"cm"),
#         panel.border = element_rect(color = "black", fill = NA))
# 
# ggpubr::ggarrange(KN95, surgical, KF94, MKF94, nrow =2, ncol = 2, widths = c(1,.85),heights = c(1.8,2))
#################################################################################
## tukey of FFE by cluster
#################################################################################
FFE_long$Condition<-factor(FFE_long$Condition, levels = c("KN95","Surgical","KF94","MKF94","KN95_clip",
                                                          "Surgical_clip","KF94_clip","MKF94_clip"),
                           labels =c("KN95","Surgical","KF94","MKF94","Mod-KN95",
                                     "Mod-Surgical","Mod-KF94","Mod-MKF94"))

## cluster 1
model1<-lm(FFE~Condition , data = FFE_long, subset = c(k_clust == "D"))
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ Condition,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

clust_1<-ggplot(data = mod_means, aes(x = Condition, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, linewidth = 1.5, color = "#D55E00") +
  geom_point(size = 4, color = "#D55E00") +
  annotate(geom="text", x="Surgical", y=93, label = "Cluster D", color = "black", size =8, fontface = "bold")+
  annotate(geom="text", x ="Surgical", y = 93, label = "\u2666", color = "#D55E00", size = 8, hjust = -4.5)+
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = mod_means$lower.CL[which(mod_means$Condition == "KF94")], ymax = mod_means$upper.CL[which(mod_means$Condition == "KF94")], 
           alpha = 0.3, fill = "gray90") +
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = mod_means$lower.CL[which(mod_means$Condition == "Surgical")], ymax = 60, 
           alpha = 0.3, fill = "gray70") +
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = mod_means$lower.CL[which(mod_means$Condition == "Mod-Surgical")], ymax = mod_means$upper.CL[which(mod_means$Condition == "Mod-Surgical")], 
           alpha = 0.3, fill = "gray50") +
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = mod_means$lower.CL[which(mod_means$Condition == "Mod-KN95")], ymax = mod_means$upper.CL[which(mod_means$Condition == "Mod-KN95")], 
           alpha = 0.3, fill = "gray10") +
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", ymin = 72, ymax = 79, alpha = 0.3, fill = "gray30") +
  annotate(geom="text", x="MKF94", y=44, label = "a", color = "black", size =6, fontface = "bold", hjust = -3)+
  annotate(geom="text", x="MKF94", y=55, label = "b", color = "black", size =6, fontface = "bold", hjust = -3)+
  annotate(geom="text", x="MKF94", y=67, label = "c", color = "black", size =6, fontface = "bold", hjust = -3)+
  annotate(geom="text", x="MKF94", y=76, label = "d", color = "black", size =6, fontface = "bold", hjust = -3)+
  annotate(geom="text", x="MKF94", y=84, label = "e", color = "black", size =6, fontface = "bold", hjust = -3)+
  #annotate("rect", xmin = "Mod-KN95", xmax = "Mod-MKF94", ymin = 35, ymax = 100, alpha = 0.2, fill = "gray") +
  #geom_text(aes(label = gsub(" ","", .group)),
   #         position = position_nudge(x = -0.3), color = "black", size = 8) +
  ylim(c(35,100)) +
  labs(y = "FFE") +
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.x = element_line(color = "black", linewidth = 1.5),
        panel.grid.major.y = element_line(color = "gray80"),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.text = element_text(size=18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(colour="black"))

## cluster 2
model1<-lm(FFE~Condition , data = FFE_long, subset = c(k_clust == "P"))
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ Condition,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

clust_2<-ggplot(data = mod_means, aes(x = Condition, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, linewidth = 1.5, color = "#FFAE00") +
  geom_point(size = 4, color = "#FFAE00") +
  annotate(geom="text", x="Surgical", y=93, label = "Cluster P", color = "black", size =8, fontface = "bold")+
  annotate(geom="text", x ="Surgical", y = 93, label = "\u2605", color = "#FFAE00", size = 11, hjust = -2.5, vjust = .3)+
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = mod_means$lower.CL[which(mod_means$Condition == "KF94")], ymax = 72, 
           alpha = 0.3, fill = "gray80") +
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = 73, ymax = 84, 
           alpha = 0.3, fill = "gray50") +
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = mod_means$lower.CL[which(mod_means$Condition == "Mod-KN95")], ymax = mod_means$upper.CL[which(mod_means$Condition == "Mod-KN95")], 
           alpha = 0.3, fill = "gray20") +
  annotate(geom="text", x="MKF94", y=67, label = "a", color = "black", size =6, fontface = "bold", hjust = -3)+
  annotate(geom="text", x="MKF94", y=79, label = "b", color = "black", size =6, fontface = "bold", hjust = -3)+
  annotate(geom="text", x="MKF94", y=90, label = "c", color = "black", size =6, fontface = "bold", hjust = -3)+
  # annotate("rect", xmin = "Mod-KN95", xmax = "Mod-MKF94", ymin = 35, ymax = 100, alpha = 0.2, fill = "gray") +
  # geom_text(aes(label = gsub(" ","", .group)),
  #           position = position_nudge(x = -0.3), color = "black", size = 8) +
  labs(y = "FFE") +
  ylim(c(35,100)) +
  theme(axis.title = element_blank(),
        #axis.ticks.x = element_blank(),
        #axis.text.x = element_blank(),
        # = element_line(color = "black", linewidth = 1.5),
        panel.grid.major.y = element_line(color = "gray80"),
        plot.margin = unit(c(0,0,0,0),"cm"),
        axis.text.x = element_text(size=18, face = "bold", angle = 90),
        axis.text.y = element_text(size=18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"),
        panel.background = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.line = element_line(colour="black"))

## Cluster 3
model1<-lm(FFE~Condition , data = FFE_long, subset = c(k_clust == "R"))
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ Condition,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

clust_3<-ggplot(data = mod_means, aes(x = Condition, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, linewidth = 1.5, color = "#0C7D74") +
  geom_point(size = 4, color = "#0C7D74") +
  annotate(geom="text", x="Surgical", y=93, label = "Cluster R", color = "black", size =8, fontface = "bold")+
  annotate(geom="text", x ="Surgical", y = 93, label = "\u25AA", color = "#0C7D74", size = 14, hjust = -4)+
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = mod_means$lower.CL[which(mod_means$Condition == "Surgical")], ymax = mod_means$upper.CL[which(mod_means$Condition == "Surgical")], 
           alpha = 0.3, fill = "gray80") +
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = 64, ymax = 66.5, 
           alpha = 0.3, fill = "gray60") +
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = 67, ymax = 67.5, 
           alpha = 0.3, fill = "gray40") +
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94",
           ymin = mod_means$lower.CL[which(mod_means$Condition == "Mod-KF94")], ymax = mod_means$upper.CL[which(mod_means$Condition == "Mod-KN95")], 
           alpha = 0.3, fill = "gray20") +
  annotate(geom="text", x="MKF94", y=60, label = "a", color = "black", size =6, fontface = "bold", hjust = -3)+
  annotate(geom="text", x="MKF94", y=65.25, label = "b", color = "black", size =6, fontface = "bold", hjust = -3)+
  annotate(geom="text", x="MKF94", y=67.45, label = "c", color = "black", size =6, fontface = "bold", hjust = -3)+
  annotate(geom="text", x="MKF94", y=72, label = "d", color = "black", size =6, fontface = "bold", hjust = -3)+
  # annotate("rect", xmin = "Mod-KN95", xmax = "Mod-MKF94", ymin = 35, ymax = 100, alpha = 0.2, fill = "gray") +
  # geom_text(aes(label = gsub(" ","", .group)),
  #           position = position_nudge(x = -0.4), color = "black", size = 8) +
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
model1<-lm(FFE~Condition , data = FFE_long, subset = c(k_clust == "T"))
summary(model1)
TukeyHSD(aov(model1))
mod_mean_contr<-emmeans::emmeans(object = model1,
                                 pairwise ~ Condition,
                                 adjust = "tukey")
mod_means<-multcomp::cld(object = mod_mean_contr$emmeans,
                         Letters = letters)

clust_4<-ggplot(data = mod_means, aes(x = Condition, y = emmean)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.3, linewidth = 1.5, color = "#0F353D") +
  geom_point(size = 4, color = "#0F353D") +
  annotate(geom="text", x="Surgical", y=93, label = "Cluster T", color = "black", size =8, fontface = "bold")+
  annotate(geom="text", x ="Surgical", y = 93, label = "\u2BC6", color = "#0F353D", size = 10, hjust = -3)+
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = mod_means$lower.CL[which(mod_means$Condition == "KF94")], ymax = mod_means$upper.CL[which(mod_means$Condition == "Surgical")], 
           alpha = 0.3, fill = "gray80") +
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = mod_means$lower.CL[which(mod_means$Condition == "Mod-Surgical")], ymax = mod_means$upper.CL[which(mod_means$Condition == "Mod-KF94")],
           alpha = 0.3, fill = "gray50") +
  annotate("rect", xmin = "KN95", xmax = "Mod-MKF94", 
           ymin = mod_means$lower.CL[which(mod_means$Condition == "Mod-KN95")], ymax = mod_means$upper.CL[which(mod_means$Condition == "Mod-KN95")], 
           alpha = 0.3, fill = "gray20") +
  annotate(geom="text", x="MKF94", y=46, label = "a", color = "black", size =6, fontface = "bold", hjust = -3)+
  annotate(geom="text", x="MKF94", y=65, label = "b", color = "black", size =6, fontface = "bold", hjust = -3)+
  annotate(geom="text", x="MKF94", y=80, label = "c", color = "black", size =6, fontface = "bold", hjust = -3)+
  # annotate("rect", xmin = "Mod-KN95", xmax = "Mod-MKF94", ymin = 35, ymax = 100, alpha = 0.2, fill = "gray") +
  # geom_text(aes(label = gsub(" ","", .group)),
  #           position = position_nudge(x = -0.2), color = "black", size = 8) +
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

ggpubr::ggarrange(clust_1,clust_4,clust_2,clust_3, ncol =2, nrow = 2, heights = c(1.5,2), widths = c(1.75,1.5))

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




#############################################################################
## Average Raw Measurements by Cluster
#############################################################################
measurements<-fread("PythonOutput/total_table_8-14-2023.csv")
measurements<-as.data.frame(measurements)
#measurements<-na.omit(measurements)
summary(measurements)

measurements$subject<-str_pad(measurements$subject,3,pad="0")
measurements_all<-measurements[,-c(2:14,32:42,44)]
measurements_all<- measurements_all %>%
  distinct()
measurements_all$subject<-str_pad(measurements_all$subject,3,pad="0")
measurements_all<-subset(measurements_all, !is.na(nose_gap_area))
measurements_all$subject<-as.character(measurements_all$subject)
all_measures<-left_join(measurements_all, stats_FFE, by="subject")
#write.csv(all_measures,"Craniometrics Paper/raw_scaled_measures_clusters_full.csv")


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

measure_clust<-subset(measure_clust, !is.na(measure_clust$k_clust))
colnames(measure_clust)[1:16]<-c("Subject","Bizygomatic.Breadth","Nose.Length","nose_gap_area","Neck.Circumference","Ear.Breadth",
                           "Sex","Neck.Circumference_scale","Ear.Breadth_scale","Nose.Length_scale", "Bizygomatic.Breadth_scale",
                           "nose_gap_area_scaled","Age_scale","Height_scale","Weight_scale", "BMI_scale")
#write.csv(measure_clust,"Craniometrics Paper/raw_scaled_measures_clusters_select.csv")

measure_clust$k_cluster<-factor(measure_clust$k_cluster, labels = c("D","P","R","T"))
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
  geom_point(size = 3) +
  annotate(geom="text", x="D", y=142, label = "Bizygomatic Breadth", color = "black", size =6, fontface = "bold")+
  annotate("rect", xmin = "D", xmax = "T", 
           ymin = mod_means$lower.CL[which(mod_means$k_cluster == "T")], ymax = mod_means$upper.CL[which(mod_means$k_cluster == "D")], 
           alpha = 0.3, fill = "gray70") +
  annotate("rect", xmin = "D", xmax = "T", 
           ymin = mod_means$lower.CL[which(mod_means$k_cluster == "P")], ymax = mod_means$upper.CL[which(mod_means$k_cluster == "R")],
           alpha = 0.3, fill = "gray30") +
  annotate(geom="text", x="P", y=130, label = "a", color = "black", size =6, fontface = "bold", hjust = -5)+
  annotate(geom="text", x="P", y=142, label = "b", color = "black", size =6, fontface = "bold", hjust = -5)+
  # geom_text(aes(label = gsub(" ","", .group)),
  #           position = position_nudge(x = 0.2), color = "black", size = 8) +
  labs(#caption = "Means followed by a common letter are \nnot significantly different according to the Tukey-test",
       color = "Cluster", y = "") +
  scale_color_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D")) +
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
  geom_point(size = 3) +
  annotate(geom="text", x="D", y=58, label = "Nose Length", color = "black", size =6, fontface = "bold")+
  annotate("rect", xmin = "D", xmax = "T", 
           ymin = mod_means$lower.CL[which(mod_means$k_cluster == "D")], ymax = mod_means$upper.CL[which(mod_means$k_cluster == "D")], 
           alpha = 0.3, fill = "gray80") +
  annotate("rect", xmin = "D", xmax = "T", 
           ymin = mod_means$lower.CL[which(mod_means$k_cluster == "P")], ymax = mod_means$upper.CL[which(mod_means$k_cluster == "P")],
           alpha = 0.3, fill = "gray50") +
  annotate("rect", xmin = "D", xmax = "T", 
           ymin = mod_means$lower.CL[which(mod_means$k_cluster == "T")], ymax = mod_means$upper.CL[which(mod_means$k_cluster == "R")],
           alpha = 0.3, fill = "gray20") +
  annotate(geom="text", x="P", y=48, label = "a", color = "black", size =6, fontface = "bold", hjust = -5)+
  annotate(geom="text", x="P", y=50.5, label = "b", color = "black", size =6, fontface = "bold", hjust = -5)+
  annotate(geom="text", x="P", y=54, label = "c", color = "black", size =6, fontface = "bold", hjust = -5)+
  # geom_text(aes(label = gsub(" ","", .group)),
  #           position = position_nudge(x = 0.2), color = "black", size = 8) +
  labs(#caption = "Means followed by a common letter are \nnot significantly different according to the Tukey-test",
    color = "Cluster", y = "") +
  scale_color_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D")) +
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
  geom_point(size = 3) +
  annotate(geom="text", x="D", y=348, label = "Nose Gap Area", color = "black", size =6, fontface = "bold")+
  annotate("rect", xmin = "D", xmax = "T", 
           ymin = mod_means$lower.CL[which(mod_means$k_cluster == "P")], ymax = mod_means$upper.CL[which(mod_means$k_cluster == "T")], 
           alpha = 0.3, fill = "gray70") +
  annotate("rect", xmin = "D", xmax = "T", 
           ymin = mod_means$lower.CL[which(mod_means$k_cluster == "R")], ymax = mod_means$upper.CL[which(mod_means$k_cluster == "R")],
           alpha = 0.3, fill = "gray30") +
  annotate(geom="text", x="P", y=220, label = "a", color = "black", size =6, fontface = "bold", hjust = -5)+
  annotate(geom="text", x="P", y=340, label = "b", color = "black", size =6, fontface = "bold", hjust = -5)+
  # geom_text(aes(label = gsub(" ","", .group)),
  #           position = position_nudge(x = 0.2), color = "black", size = 8) +
  labs(#caption = "Means followed by a common letter are \nnot significantly different according to the Tukey-test",
    color = "Cluster", y = "") +
  scale_color_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D")) +
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
  geom_point(size = 3) +
  annotate(geom="text", x="D", y=415, label = "Neck Circumference", color = "black", size =6, fontface = "bold")+
  annotate("rect", xmin = "D", xmax = "T", 
           ymin = mod_means$lower.CL[which(mod_means$k_cluster == "D")], ymax = mod_means$upper.CL[which(mod_means$k_cluster == "T")], 
           alpha = 0.3, fill = "gray80") +
  annotate("rect", xmin = "D", xmax = "T", 
           ymin = mod_means$lower.CL[which(mod_means$k_cluster == "P")], ymax = mod_means$upper.CL[which(mod_means$k_cluster == "P")],
           alpha = 0.3, fill = "gray50") +
  annotate("rect", xmin = "D", xmax = "T", 
           ymin = mod_means$lower.CL[which(mod_means$k_cluster == "R")], ymax = mod_means$upper.CL[which(mod_means$k_cluster == "R")],
           alpha = 0.3, fill = "gray20") +
  annotate(geom="text", x="P", y=365, label = "a", color = "black", size =6, fontface = "bold", hjust = -5)+
  annotate(geom="text", x="P", y=395, label = "b", color = "black", size =6, fontface = "bold", hjust = -5)+
  annotate(geom="text", x="P", y=412, label = "c", color = "black", size =6, fontface = "bold", hjust = -5)+
  # geom_text(aes(label = gsub(" ","", .group)),
  #           position = position_nudge(x = 0.2), color = "black", size = 8) +
  labs(#caption = "Means followed by a common letter are \nnot significantly different according to the Tukey-test",
    color = "Cluster", y = "") +
  scale_color_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D")) +
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

EB<-ggplot(data = mod_means, aes(x = k_cluster, y = emmean, color = k_cluster,)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, linewidth = 1.5) +
  geom_point(size = 3) +
  annotate(geom="text", x="D", y=38, label = "Ear Breadth", color = "black", size =6, fontface = "bold")+
  annotate("rect", xmin = "D", xmax = "T", 
           ymin = mod_means$lower.CL[which(mod_means$k_cluster == "T")], ymax = mod_means$upper.CL[which(mod_means$k_cluster == "T")], 
           alpha = 0.3, fill = "gray70") +
  annotate("rect", xmin = "D", xmax = "T", 
           ymin = mod_means$lower.CL[which(mod_means$k_cluster == "R")], ymax = mod_means$upper.CL[which(mod_means$k_cluster == "P")],
           alpha = 0.3, fill = "gray30") +
  annotate(geom="text", x="P", y=32, label = "a", color = "black", size =6, fontface = "bold", hjust = -5)+
  annotate(geom="text", x="P", y=36.5, label = "b", color = "black", size =6, fontface = "bold", hjust = -5)+
  # geom_text(aes(label = gsub(" ","", .group)),
  #           position = position_nudge(x = 0.2), color = "black", size = 8) +
  labs(#caption = "Means followed by a common letter are \nnot significantly different according to the Tukey-test",
    color = "Cluster", y = "") +
  scale_color_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D")) +
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


#######################################################################################
### LDA star plot
#######################################################################################
ggplotLDAPrep <- function(x){
  if (!is.null(Terms <- x$terms)) {
    data <- model.frame(x)
    X <- model.matrix(delete.response(Terms), data)
    g <- model.response(data)
    xint <- match("(Intercept)", colnames(X), nomatch = 0L)
    if (xint > 0L) 
      X <- X[, -xint, drop = FALSE]
  }
  means <- colMeans(x$means)
  X <- scale(X, center = means, scale = FALSE) %*% x$scaling
  rtrn <- as.data.frame(cbind(X,labels=as.character(g)))
  rtrn <- data.frame(X,labels=as.character(g))
  return(rtrn)
}


set.seed(5)
lda.fit<-lda(k_cluster~., data = measure_clust[,-c(1,7:25)])
lda.fit
fitGraph<-ggplotLDAPrep(lda.fit)
fitGraph_sum<-fitGraph %>%
  group_by(labels) %>%
  summarise(LD1 = mean(LD1), LD2 = mean(LD2))

fitGraph_all<-left_join(fitGraph,fitGraph_sum, by="labels")
fit_clust1<-subset(fitGraph[,-3], labels == "D")
fit_clust1<-rbind(fit_clust1, fitGraph_sum[1,] %>% slice(rep(1:n(), times = 24)))
fit_clust1$grp<-as.factor(rep(1:24, times = 1))

fit_clust2<-subset(fitGraph[,-3], labels == "P")
fit_clust2<-rbind(fit_clust2, fitGraph_sum[2,] %>% slice(rep(1:n(), times = 20)))
fit_clust2$grp<-as.factor(rep(1:20, times = 1))

fit_clust3<-subset(fitGraph[,-3], labels == "R")
fit_clust3<-rbind(fit_clust3, fitGraph_sum[3,] %>% slice(rep(1:n(), times = 34)))
fit_clust3$grp<-as.factor(rep(1:34, times = 1))

fit_clust4<-subset(fitGraph[,-3], labels == "T")
fit_clust4<-rbind(fit_clust4, fitGraph_sum[4,] %>% slice(rep(1:n(), times = 21)))
fit_clust4$grp<-as.factor(rep(1:21, times = 1))

colnames(fitGraph_sum)[1]<-"Cluster"
                                                  #shapes "\u2666","\u2605","\u25AA","\u2BC6")
ggplot(fit_clust1)+
  geom_point(aes(LD1,LD2, group = grp),color = "#D55E00", size = 2) +
  geom_line(aes(LD1,LD2, group = grp),color = "#D55E00") +
  geom_point(data = fit_clust2, aes(LD1,LD2, group = grp), color = "#FFAE00", size = 2) +
  geom_line(data = fit_clust2, aes(LD1,LD2, group = grp),color = "#FFAE00") +
  geom_point(data = fit_clust3, aes(LD1,LD2, group = grp), color = "#0C7D74", size = 2) +
  geom_line(data = fit_clust3, aes(LD1,LD2, group = grp),color = "#0C7D74") +
  geom_point(data = fit_clust4, aes(LD1,LD2, group = grp), color = "#0F353D", size = 2) +
  geom_line(data = fit_clust4, aes(LD1,LD2, group = grp),color = "#0F353D") +
  geom_point(data = fitGraph_sum, mapping = aes(x = LD1, y = LD2, color = Cluster), size = 4) +
  scale_color_manual(values = col_pal) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.x = element_line(color = "gray80"),
        panel.grid.major.y = element_line(color = "gray80"),
        axis.text = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"))



#####################################################################################
## Decision Tree
## Tree should use measurements to lead to cluster
#####################################################################################

#library(datasets) 
library(caTools)
library(party)
library(magrittr)

tree_mod<-ctree(k_cluster~.,data = measure_clust[,-c(1,7:25)], controls = ctree_control(mincriterion = .8))
plot(tree_mod)


## Decision Tree 
set.seed(2)

tree.clust<-tree(k_cluster~., measure_clust[,-c(1,7:25)])
tree.clust

plot(tree.clust)
text(tree.clust, pretty = 5)


library(rpart)
library(rpart.plot)
rpart_tree<-rpart(k_cluster~., measure_clust[,-c(1,7:25)], control = rpart.control(minsplit = 14))
rpart.plot(rpart_tree)

##########################################################################################
## Standardized means bar plot grouped by craniometric and cluster
##########################################################################################

vars_long<-pivot_longer(stats_FFE, cols = c(7,8,14,16,19), names_to="cran_metric", values_to="measure")
vars_long<-as.data.frame(vars_long[,c(1,29:31)])
vars_long$k_cluster<-factor(vars_long$k_cluster, labels=c("D","P","R","T"))

vars_sum<-vars_long %>% group_by(k_cluster, cran_metric) %>%
  summarise(mean = mean(measure, na.rm=T))

library(tidytext)

ggplot(vars_sum, aes(y = mean, x = reorder_within(k_cluster, by =  mean, within = cran_metric), fill = k_cluster)) +
  geom_bar(position = "dodge", stat = "summary", fun = "mean", width = .8) +
  scale_fill_manual(values = c("#D55E00","#FFAE00","#0C7D74","#0F353D"), ) +
  scale_x_reordered() +
  ylim(-1.1,1.1) +
  facet_wrap(~cran_metric, scales = "free") +
  #scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,by= .2)) +
  labs(y = "Average", x = "", fill="Cluster") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.x = element_line(color = "gray80"),
        #panel.grid.major.y = element_line(color = "gray80"),
        axis.text = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"),
        strip.text = element_text(size = 18),
        aspect.ratio = 2/3) +
  coord_flip() 
  





all_measures[,c(12,52)] %>% group_by(k_cluster) %>% 
  summarise(mean= mean(`Nose Breadth`))