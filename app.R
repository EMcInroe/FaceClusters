#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(tidyverse)
library(ggthemes)
library(purrr)
library(DT)
library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(leaps)
library(rstatix)
library(RColorBrewer)
library(MASS)
options(tigris_class = "sf")
options(tigris_use_cache = TRUE)

setwd("")
measurements<-fread("PythonOutput/total_table_8-14-2023.csv")
#measurements<-fread("www/total_table_8-14-2023.csv")
names(measurements)
# [1] "subject"               "Condition"             "RH"                    "temp"                  "Bending Mean"          "Bending SD"            "Reading Mean"          "Reading SD"            "LR Mean"              
# [10] "LR SD"                 "UD Mean"               "UD SD"                 "Overall Mean"          "Overall SD"            "Bitragion Chin Arc"    "Bitragion Coronal Arc" "Bitragion Frontal Arc" "Head Circumference"   
# [19] "Neck Circumference"    "Ear Breadth"           "Ear Length"            "Lip Length"            "Upper Facial Breadth"  "Menton-Sellion Length" "Nose Breadth"          "Nose Length"           "Bigonial Breadth"     
# [28] "Bizygomatic Breadth"   "Head Breadth"          "Head Length"           "nose_gap_area"         "nose_gap_curve_length" "Bizyg_1_11"            "MenSel_6_89"           "EarLength_2_9"         "LipLength_59_65"      
# [37] "sel-earb"              "nose-earb"             "chin-earb"             "lip-earb"              "sphenomaxillary_angle" "ratio_sel-nose_tip"    "ratio_sel-lip"         "ratio_sel-chin"       
measure_clust<-fread("Craniometrics Paper/Old and Misc/raw_scaled_measures_clusters_select.csv", header = TRUE, drop = 1)
#measure_clust<-fread("/www/raw_scaled_measures_clusters_select.csv")
names(measure_clust)
# [1] "Subject"                   "Bizygomatic.Breadth"       "Nose.Length"               "nose_gap_area"             "Neck.Circumference"        "Ear.Breadth"               "Sex"                      
# [8] "Neck.Circumference_scale"  "Ear.Breadth_scale"         "Nose.Length_scale"         "Bizygomatic.Breadth_scale" "nose_gap_area_scaled"      "Age_scale"                 "Height_scale"             
# [15] "Weight_scale"              "BMI_scale"                 "N95"                       "KN95"                      "KN95_clip"                 "Surgical"                  "Surgical_clip"            
# [22] "KF94"                      "KF94_clip"                 "MKF94"                     "MKF94_clip"                "KN95_diff"                 "Surgical_diff"             "KF94_diff"                
# [29] "MKF94_diff"                "k_cluster"                 "Age"                       "BMI"                       "k_clust"       

col_pal<-c("#D55E00","#FFAE00","#0C7D74","#0F353D") 
measure_clust$k_cluster<-factor(measure_clust$k_cluster, labels = c("D","P","R","T"))
measure_clust<-left_join(measure_clust,measurements[,c(1,32)], by = c("Subject"="subject"))
measure_clust<-distinct(measure_clust)
###########################################################
## Linear Discriminant Analysis
###########################################################
set.seed(5)
lda.fit<-lda(k_clust~., data = measure_clust[,-c(1,4,7:32)])
lda.fit
# plot(lda.fit)
# lda.pred<-predict(lda.fit, measure_clust_test[,-c(1,7:33)])
# lda.class<-lda.pred$class
# table(lda.class, measure_clust_test$k_clust)

# sample<-data.frame(matrix(c(130.3,48.9,203,388,34.7),nrow = 1, ncol=5))
# colnames(sample)<-c("Bizygomatic.Breadth","Nose.Length","nose_gap_area","Neck.Circumference","Ear.Breadth")
# lda.pred_samp<-predict(lda.fit,sample)
# lda.pred_samp$class


lda_fun<-function(bb,nl,nac,nc,eb, model = lda.fit) {
  eb_mm<-eb*10
  nl_mm<-nl*10
  bb_dub<-bb*2*10
  nac_dub<-nac*2*10
  nc_dub<-nc*2*10
  sample<-data.frame(matrix(c(bb_dub,nl_mm,nac_dub,nc_dub,eb_mm), nrow = 1, ncol=5))
  colnames(sample)<-c("Bizygomatic.Breadth","Nose.Length","nose_gap_curve_length","Neck.Circumference","Ear.Breadth")
  lda.pred_samp<-predict(model,sample)
  return(lda.pred_samp$class)
}


#######################################################################################################################
############################################# Set Up Dashboard ########################################################
#######################################################################################################################

ui<-dashboardPage(
  dashboardHeader(title = "Face Fit 2.0"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Overview", tabName = "overview", icon = icon("home",lib = "glyphicon")),
      menuItem("Measurement Input", tabName = "fittest", icon = icon("stats", lib = "glyphicon")),
      menuItem("Group Description", tabName = "results", icon = icon("book", lib = "glyphicon"))
    )
  ),
  dashboardBody(
    tags$head(tags$style(HTML(".skin-blue .main-header .logo {background-color:#60978B;}
                              .skin-blue .main-header .navbar {background-color:#D1CF7C;}
                              .skin-blue .main-sidebar {background-color;#7E8F57;}"))),
    tabItems(
      ################################################# Overview ###########################################################
      tabItem(tabName = "overview",
              fluidPage(
                tags$head(
                  tags$style(HTML("
            code {
                display:block;
                padding:9.5px;
                margin:0 0 10px;
                margin-top:10px;
                font-size:13px;
                line-height:20px;
                word-break:break-all;
                word-wrap:break-word;
                white-space:pre-wrap;
                background-color:#F5F5F5;
                border:1px solid rgba(0,0,0,0.15);
                border-radius:4px; 
                font-family:monospace;
            }"))),
                titlePanel("What is Face Fit"),
                mainPanel(
                  p("Face Fit is a study conducted at the EPA's Human Studies Facility in Chapel Hill to determine how well disposable masks perform under different conditions. The study measured 
                  filtering efficiency for a diverse sample of the population while they performed a series of simple movements. Facial measurements were also taken for each particpant as well 
                  as demographic information. The results were then compliled and used to create groups based on the facial measurements and mask performance. We can now use these groupings to
                  predict the most effective mask for others."),
                  p("By entering a few simple measurements, you will be able to know the", span("most efficient mask", style = "color:blue"),
                    "for your face."),
                  fluidRow(column(6,imageOutput("star_plot")),
                  column(6,box(p(strong("Proceed to the next page to enter your measurements.")))))
                )
              )),
      ########################################## Input Measurements ######################################################
      tabItem(tabName = "fittest",
              fluidPage(
                titlePanel("Craniofacial Measurements"),
                fluidRow(
                  tabBox(tags$head(tags$style(HTML("#tabBox{height:90vh !important;}"))),
                         width = 12,
                    tabPanel(
                      fluidRow(
                      column(6,imageOutput("measure_image",inline = T)),
                      column(6,imageOutput("measure_image2",inline =T))
                    ))),
                ),
                h2("Enter Your Measurements Here"),
                fluidRow(
                  column(4,
                         box(numericInput("bb","Face Width",0))),
                  column(4,
                         box(numericInput("nl","Nose Length",0))),
                  column(4,
                         box(numericInput("nga","Nose Arc",0)))
                ),
                fluidRow(
                  column(4,
                         box(numericInput("nc","Neck Circumference",0))),
                  column(4, 
                         box(numericInput("eb","Ear Width",0))))
              ),
              fluidRow(
                box(title = "Your Group is:",
                    span(textOutput("lda_group_results"), style = "color:blue;font-size:30px"))
              )
              ),
################################## Group Description #######################################
        tabItem(
          tabName = "results",
          fluidPage(
            titlePanel("Group Description"),
            fluidRow(
              box(span(textOutput("lda_group_results2"), style = "color:blue;font-size:40px"))
            ),
            fluidRow(
              box(imageOutput("clust_image", inline = T))
            ),
            fluidRow(
              span(htmlOutput("text1"), style = "font-size:20px")
            ),
            fluidRow(
              span(htmlOutput("text2"), style = "font-size:20px")
            ),
            fluidRow(
              span(htmlOutput("text3"), style = "font-size:20px")
            ),
            fluidRow(
              span(htmlOutput("text4"), style = "font-size:20px")
            )
          )
        )

  )))


####################################################################################################################
################################################# Server ###########################################################
####################################################################################################################
server<-function(input,output) {
  output$lda_group_results<-renderText({
    if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==1) {
      "D - Diamond"
    } else if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==2) {
      "P - Pentagon"
    } else if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==3) {
      "R - Rectangle"
    } else {
      "T - Triangle"
    }
  })
  
  output$measure_image <- renderImage({
    filename <- "./www/pg1directions.png"
    list(src = filename,
         width = "100%")
  }, deleteFile = FALSE)
  
  output$measure_image2 <- renderImage({
    filename <- "./www/pg2directions.png"
    list(src = filename,
         width = "100%")
  }, deleteFile = FALSE)
  
  output$lda_group_results2<- renderText({
    if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==1) {
      "Diamond:"
    } else if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==2) {
      "Pentagon:"
    } else if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==3) {
      "Rectangle:"
    } else {
      "Triangle:"
    }
  })
  
  output$clust_image<-renderImage({
    if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==1) {
      filename<-"./www/diamond_clust.png"
    } else if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==2) {
      filename<-"./www/pentagon_clust.png"
    } else if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==3) {
      filename<-"./www/rectangle_clust.png"
    } else {
      filename<-"./www/triangle_clust.png"
    }
    list(src = filename,
         width = "100%")
  }, deleteFile = FALSE)
  
  output$text1 <-renderUI({
    if (lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==3){
      HTML(paste0("<b>","Should you consider wearing a clip: ", "</b>", "MAYBE"))
    } else {
      HTML(paste0("<b>","Should you consider wearing a clip: ", "</b>", "YES"))
    }
  })
  
  output$text2 <-renderUI({
    if (lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==2){
      HTML(paste0("<b>","What should you wear if you want to get above 80% filtration from your disposable mask?","</b><li>", 
                  "N95, KN95 (clipped), Large KF94 (clipped), or Medium KF94 (clipped)", "</li>"))
    } else if (lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==3){
      HTML(paste0("</b>","What should you wear if you want to get above 80% filtration from your disposable mask?","</b><li>", 
                  "N95", "</li>"))
    } else {
      HTML(paste0("<b>","What should you wear if you want to get above 80% filtration from your disposable mask?","</b><li>", 
                  "N95, KN95 (clipped), or Medium KF94 (clipped)", "</li>"))
    }
  })
  
  output$text3 <- renderUI({
    if (lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==1){
      HTML(paste0("<b>","What should you wear if you want to get above 60% filtration from your disposable mask?","</b><li>", 
                  "N95, KN95 (clipped), Medium KF94 (clipped), Surgical (clipped), KF94 (clipped)", "</li>"))
    } else if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==2) {
      HTML(paste0("<b>","What should you wear if you want to get above 60% filtration from your disposable mask?","</b><li>", 
                  "N95, KN95 (clipped or unclipped), Surgical (clipped or unclipped), Large KF94 (clipped or unclipped), Medium KF94 (clipped or unclipped", "</li>"))
    } else if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==3) {
      HTML(paste0("<b>","What should you wear if you want to get above 60% filtration from your disposable mask?","</b><li>", 
                  "N95, KN95 (clipped or unclipped), Surgical (clipped), Large KF94 (clipped), Medium KF94 (clipped or unclipped)", "</li>"))
    } else {
      HTML(paste0("<b>","What should you wear if you want to get above 60% filtration from your disposable mask?","</b><li>", 
                  "N95, KN95 (clipped), MKF94 (clipped), Large KF94 (clipped)", "</li>"))
      }
  })
  
  output$text4 <- renderUI ({
    if (lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==1){
      HTML(paste0("<b>","You are likely to substantially improve the protection you receive from any surgical/procedure, KN95, or KF94 mask by wearing a clip. 
           Using a clip increases your chances of having a fitted filtration efficiency of above 80%.","</b>"))
    } else if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==2) {
      HTML(paste0("<b>","You are likely to substantially improve the protection you receive from any surgical/procedure, KN95, or KF94 mask by wearing a clip. 
           Using a clip increases your chances of having a fitted filtration efficiency of above 80%. You are likely to receive 60% protection without.","</b>"))
    } else if(lda_fun(input$bb, input$nl, input$nga, input$nc, input$eb)==3) {
      HTML(paste0("<b>","You are not likely to substantially improve the protection you receive from any surgical/procedure, KN95, or KF94 mask. 
           There is a chance that using a clip will hurt your mask performance. The N95 is your best chance of getting 80% or higher protection.","</b>"))
    } else {
      HTML(paste0("<b>","You are likely to substantially improve the protection you receive from any surgical/procedure, KN95, or KF94 mask by wearing a clip. You are unlikely to 
           receive 60% or higher protection in a surgical mask whether clipped or unclipped. Your best chance to get above 60% protection is by wearing a N95, 
           KN95 (clipped), or KF94 (clipped). ","</b>"))
    }
  })

  
  
  output$star_plot<-renderImage({
    filename <- "./www/star_plot.png"
    list(src = filename,
         width = "80%")
  }, deleteFile = FALSE)

}

#--------------------------------------------------------------------------------------------------------------------

shinyApp(ui=ui, server = server)
