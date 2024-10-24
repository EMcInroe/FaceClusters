### Predictive Modeling for Craniofacial measurements
### In this script we try multiple methods to predict the clustering created in 
## the script PCA_cluster_selectvars.R
## Only the select variables are used for the predictions
## Models are validated using cross validations and a confusion matrix 
## Created October 23, 2023
## By: Melissa McInroe
#-----------------------------------------------------------------------------------

library(data.table)
library(tree)
library(class)
library(MASS)
library(randomForest)
library(gbm)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(corrplot)
library(rstatix)

## load data
setwd("L:/PRIV/EPHD_CRB/FACEFIT/DATA")
measure_clust<-fread("Craniometrics Paper/Old and Misc/raw_scaled_measures_clusters_select.csv", header = TRUE, drop = 1)

measure_clust$k_clust<-as.factor(measure_clust$k_clust)

# measure_clust_long<-pivot_longer(measure_clust, cols = c(2:4,6), names_to = "aspects", values_to = "measurements")

## Plots of clusters
col_pal<-c("#FFAE00","#0F353D", "#0C7D74", "#D55E00")  
#col_pal<-c("#78160C","#164C45", "#508C9B","#C4830A")        #C4830A


## Scatterplots of craniometric variables by cluster
# plot1<-ggplot(measure_clust, aes(x = Bizygomatic.Breadth, y = Neck.Circumference, color = k_clust)) +
#   geom_point(size = 3) +
#   scale_color_manual(values = col_pal) +
#   labs(x = "Bizygomatic Breadth", y = "Neck Circumference", color = "Cluster") +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 
# plot2<-ggplot(measure_clust, aes(x = Ear.Breadth, y = Neck.Circumference, color = k_clust)) +
#   geom_point(size = 3) +
#   scale_color_manual(values = col_pal) +
#   labs(x = "Ear.Breadth", y = "Neck Circumference", color = "Cluster") +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 
# plot3<-ggplot(measure_clust, aes(x = Nose.Length, y = Neck.Circumference, color = k_clust)) +
#   geom_point(size = 3) +
#   scale_color_manual(values = col_pal) +
#   labs(x = "Nose Length", y = "Neck Circumference", color = "Cluster") +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 
# plot4<-ggplot(measure_clust, aes(x = nose_gap_area, y = Neck.Circumference, color = k_clust)) +
#   geom_point(size = 3) +
#   scale_color_manual(values = col_pal) +
#   labs(x = "Nose Gap Area", y = "Neck Circumference", color = "Cluster") +
#   theme(axis.line = element_line(colour = "black"),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill=NA, linewidth=1))
# ggarrange(plot1, plot2,plot3,plot4, common.legend = TRUE)

###Correlation
cormatrix<-cor_mat(measure_clust[,-c(1,7:16,24:29,32)], method = "pearson")
cor_pvals<-cor_gather(cormatrix)
cor_pvals<-as.data.frame(cor_pvals)

corrplot(cor(measure_clust[,-c(1,7:16,24:29,32)], method = "pearson"))


## Decision Tree 
set.seed(2)
train<-sample(1:nrow(measure_clust), 82)
measure_clust_test<-measure_clust[-train,]
tree.clust<-tree(k_clust~.-Subject, measure_clust, subset = train)
tree.pred<-predict(tree.clust, measure_clust_test, type = "class")
table(tree.pred, measure_clust_test$k_clust)
6/10 * 100

train<-sample(1:nrow(measure_clust), 79)
measure_clust_test<-measure_clust[-train,]
tree.clust<-tree(k_clust~.-Subject, measure_clust, subset = train)
tree.pred<-predict(tree.clust, measure_clust_test, type = "class")
table(tree.pred, measure_clust_test$k_clust)
(4+5+4)/20*100

plot(tree.clust)
text(tree.clust, pretty = 5)

cv.cluster<-cv.tree(tree.clust, FUN = prune.misclass)
cv.cluster
plot(cv.cluster$size, cv.cluster$dev, type = "b")
plot(cv.cluster$k, cv.cluster$dev, type = "b")
prune.cluster<- prune.misclass(tree.clust, best = 7)
plot(prune.cluster)
text(prune.cluster, pretty =0)
tree.pred<-predict(prune.cluster, measure_clust_test, type = "class")
table(tree.pred, measure_clust_test$k_clust)
## 7 nodes is no different from 9 nodes

###########################################################
## Linear Discriminant Analysis
###########################################################
set.seed(5)
train<-sample(1:nrow(measure_clust), 79)
measure_clust_test<-measure_clust[-train,]
lda.fit<-lda(k_clust~., data = measure_clust[,-c(1,7:32)], subset = train)
lda.fit
plot(lda.fit)
lda.pred<-predict(lda.fit, measure_clust_test[,-c(1,7:33)])
lda.class<-lda.pred$class
table(lda.class, measure_clust_test$k_clust)
# (3+5+7+4)/20*100             #90% accurate
# mean(lda.class == measure_clust_test$k_clust)



## cross_validated
set.seed(10)
lda.fit<-lda(k_clust~., data = measure_clust[,-c(1,7:32)], CV=TRUE)
lda.fit
measure_clust$class<-lda.fit$class
table(measure_clust$k_clust, measure_clust$class)
measure_clust$class<-NULL

posterior<-as.data.frame(lda.fit$posterior)
posterior$class<-lda.fit$class
posterior$cluster<-measure_clust$k_clust
posterior$misclass<-ifelse(posterior$class == posterior$cluster, "NO", "YES")

subset(posterior, misclass == "YES")

#find way to plot with colors of clusters and centroids
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
lda.fit<-lda(k_cluster~., data = measure_clust[,-c(1,7:32)])
lda.fit
fitGraph<-ggplotLDAPrep(lda.fit)
fitGraph_sum<-fitGraph %>%
  group_by(labels) %>%
  summarise(LD1 = mean(LD1), LD2 = mean(LD2))

fitGraph_all<-left_join(fitGraph,fitGraph_sum, by="labels")
fit_clust1<-subset(fitGraph[,-3], labels == "1")
fit_clust1<-rbind(fit_clust1, fitGraph_sum[1,] %>% slice(rep(1:n(), times = 24)))
fit_clust1$grp<-as.factor(rep(1:24, times = 2))

fit_clust2<-subset(fitGraph[,-3], labels == "2")
fit_clust2<-rbind(fit_clust2, fitGraph_sum[2,] %>% slice(rep(1:n(), times = 20)))
fit_clust2$grp<-as.factor(rep(1:20, times = 2))

fit_clust3<-subset(fitGraph[,-3], labels == "3")
fit_clust3<-rbind(fit_clust3, fitGraph_sum[3,] %>% slice(rep(1:n(), times = 34)))
fit_clust3$grp<-as.factor(rep(1:34, times = 2))

fit_clust4<-subset(fitGraph[,-3], labels == "4")
fit_clust4<-rbind(fit_clust4, fitGraph_sum[4,] %>% slice(rep(1:n(), times = 21)))
fit_clust4$grp<-as.factor(rep(1:21, times = 2))

colnames(fitGraph_sum)[1]<-"Cluster"

ggplot(fit_clust1)+
  geom_point(aes(LD1,LD2, group = grp),color = "#78160C") +
  geom_line(aes(LD1,LD2, group = grp),color = "#78160C") +
  geom_point(data = fit_clust2, aes(LD1,LD2, group = grp), color = "#164C45") +
  geom_line(data = fit_clust2, aes(LD1,LD2, group = grp),color = "#164C45") +
  geom_point(data = fit_clust3, aes(LD1,LD2, group = grp), color = "#508C9B") +
  geom_line(data = fit_clust3, aes(LD1,LD2, group = grp),color = "#508C9B") +
geom_point(data = fit_clust4, aes(LD1,LD2, group = grp), color = "#C4830A") +
  geom_line(data = fit_clust4, aes(LD1,LD2, group = grp),color = "#C4830A") +
  geom_point(data = fitGraph_sum, mapping = aes(x = LD1, y = LD2, shape = Cluster, color = Cluster), size = 4) +
  scale_shape_manual(values = c(15,16,17,18)) +
  scale_color_manual(values = col_pal) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        panel.grid.major.x = element_line(color = "gray80"),
        panel.grid.major.y = element_line(color = "gray80"),
        axis.text = element_text(size = 18, face = "bold"),
        axis.title = element_text(size = 18, face = "bold"))
  

##############################################################
## Random Forest
#############################################################
set.seed(3)
# boost.clust<-gbm(k_clust~.-Subject, data = measure_clust[train,],
#                  distribution = "gaussian", n.trees = 5000, interaction.dept = 4)
# yhat.boost<-predict(boost.clust, newdata = measure_clust_test, n.trees = 5000, type = "response")

## Random Forest
cluster.rf<-randomForest(k_clust~.-Subject, data = measure_clust, subset = train, mtry = 2, importance = TRUE)
yhat.rf<-predict(cluster.rf, newdata = measure_clust_test)
table(yhat.rf, measure_clust_test$k_clust)
(1+4+7+4)/20*100      # 80% accurate
 


#####################################################################################################
## Clustering with different measurements 
## Used to create a prediction model for shiny app
## Compare new clusters to old clusters 
#####################################################################################################
# measurements<-fread("PythonOutput/total_table_8-14-2023.csv")
# measurements<-as.data.frame(measurements)
# #measurements<-na.omit(measurements)
# summary(measurements)
# 
# measurements$subject<-str_pad(measurements$subject,3,pad="0")
# measurements_all<-measurements[,-c(2:14,32:42,44)]
# measurements_all<- measurements_all %>%
#   distinct()
# measurements_all$subject<-str_pad(measurements_all$subject,3,pad="0")
# measurements_all<-subset(measurements_all, !is.na(nose_gap_area))
# measurements_all$subject<-as.character(measurements_all$subject)
# all_measures<-left_join(measurements_all, stats_FFE, by="subject")
# 
# select_var<-c("Bizygomatic Breadth","Nose Length","nose_gap_area","Neck Circumference","Ear Breadth")
# select_var<-c("subject","Condition", select_var)
# #measurements<-measurements[,c(1,2,15:31,43)]
# measurements<-measurements[,c(select_var)]
# 
# # measure_clust<-left_join(measurements, FFE_subset, by = c("subject" = "Subject"))
# # measure_clust<- measure_clust[,-2] %>%
# #   distinct()
# measurements$subject<-as.character(measurements$subject)
# 
# measure_clust<-left_join(measurements, stats_FFE[,-c(3:6,9:13,15,17,18,20)], by = "subject")
# measure_clust<- measure_clust[,-2] %>%
#   distinct()  
# 
# measure_clust<-subset(measure_clust, !is.na(measure_clust$k_clust))
# colnames(measure_clust)[1:16]<-c("Subject","Bizygomatic.Breadth","Nose.Length","nose_gap_area","Neck.Circumference","Ear.Breadth",
#                                  "Sex","Neck.Circumference_scale","Ear.Breadth_scale","Nose.Length_scale", "Bizygomatic.Breadth_scale",
#                                  "nose_gap_area_scaled","Age_scale","Height_scale","Weight_scale", "BMI_scale")
# #write.csv(measure_clust,"Craniometrics Paper/raw_scaled_measures_clusters_select.csv")
# 
# measure_clust$k_cluster<-factor(measure_clust$k_cluster, labels = c("D","P","R","T"))