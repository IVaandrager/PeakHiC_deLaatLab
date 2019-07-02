####################################################
###########create a model to predict TPM############
####################################################

library(caret)
library(mlbench)
library(caretEnsemble)
library(MLmetrics)

set.seed(1243)
#statObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/statObjTSS.rds")
#statObj <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/avg_chip/statObjTSS.rds")
statObj2 <- readRDS("/delaat/group/iwan/Hap1_peakHiC/statistics/statObj/FullStatObj3.rds")
statObj = statObj2[statObj2$vp.TSS==1 | statObj2$anch.TSS==1,]

expression.cat <- function(x){
  if(x<10){
    return('low')
  }
  if(x<30){
    return('medium')
  }else{
    return('high')
  }
}


######### for only VP gene prediction in WT ##################
#statObj_vp.Hap1 <- data.frame(HiC.reads=statObj$Hap1, loopdist=statObj$loopDist, statObj[,28:41], statObj[,46:56], 
#                              TPM=statObj$vp.TPM.Hap1)
statObj <- statObj[,-c(1,2,4:12,14,17:19,21:27,81:84,85:94,96:109)]
statObj_vp.Hap1 <- statObj[!is.na(statObj$avg.Hap1.TPM),]
names(statObj_vp.Hap1)[which(names(statObj_vp.Hap1)=="avg.Hap1.TPM")] <- "TPM"
statObj_vp.Hap1$TPM <- factor((lapply(statObj_vp.Hap1$TPM, expression.cat)), levels=c("low","medium","high"))

df <- statObj_vp.Hap1

###general models###
control <- trainControl(method="repeatedcv", number=10, repeats=3)
metric <- "Accuracy"
# Bagged CART
fit.treebag <- train(TPM~., data=df, method="treebag", metric=metric, trControl=control,na.action=na.exclude)
# Random Forest
fit.rf <- train(TPM~., data=df, method="rf", metric=metric, trControl=control,na.action=na.exclude)
# summarize results
bagging_results <- resamples(list(treebag=fit.treebag, rf=fit.rf))
summary(bagging_results)
dotplot(bagging_results)


#alternative with CV
inTrain <- createDataPartition(
  y = df$TPM,
  ## the outcome data are needed
  p = .75,
  ## The percentage of data in the
  ## training set
  list = FALSE
)

training <- df[ inTrain,]
testing  <- df[-inTrain,]

ctrl <- trainControl(
  method = "repeatedcv", 
  repeats = 3,
  classProbs = TRUE, 
  summaryFunction = multiClassSummary
)


plsFit <- train(
  TPM ~ .,
  data = training,
  method = "pls",
  preProc = c("center", "scale"),
  tuneLength = 15,
  trControl = ctrl,
  metric = "ROC",
  na.action=na.exclude
)

plsFit
ggplot(plsFit)

plsClasses <- predict(plsFit, newdata = testing)
str(plsClasses)
#>  Factor w/ 2 levels "M","R": 2 1 1 2 1 2 2 2 2 2 ...
plsProbs <- predict(plsFit, newdata = testing, type = "prob")
head(plsProbs)

confusionMatrix(data = plsClasses, testing$TPM)


logit <- glm(TPM~., data = training, family = 'binomial')
summary(logit)
   
# Deviance Residuals: 
#   Min       1Q   Median       3Q      Max  
# -3.3619  -1.0911   0.6282   0.9465   2.1308  
# 
# Coefficients: (1 not defined because of singularities)
# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                 -6.319e+00  1.344e+00  -4.703 2.57e-06 ***
#   Hap1                         6.463e-03  1.665e-03   3.882 0.000104 ***
#   maxV4Cscore                 -3.146e-03  1.083e-01  -0.029 0.976816    
# delta                       -2.396e+00  6.669e-01  -3.593 0.000327 ***
#   ratio                        5.195e+00  1.231e+00   4.220 2.44e-05 ***
#   loopDist                     3.726e-07  1.618e-07   2.303 0.021269 *  
#   vp.TSS                      -7.365e-01  9.346e-02  -7.880 3.27e-15 ***
#   anch.TSS                    -4.950e-01  1.973e-01  -2.509 0.012116 *  
#   vp.CTCF_WT                   4.731e-01  4.230e-01   1.118 0.263371    
# anch.CTCF_WT                 5.897e-01  4.253e-01   1.387 0.165592    
# vp.SMC1_WT                  -4.271e-02  1.111e-01  -0.385 0.700572    
# anch.SMC1_WT                -1.452e-01  1.108e-01  -1.310 0.190315    
# vp.H3K27ac                   6.976e-02  5.462e-02   1.277 0.201499    
# anch.H3K27ac                 1.401e-01  5.552e-02   2.524 0.011607 *  
#   vp.H3K4me3                  -1.585e-01  5.344e-02  -2.967 0.003008 ** 
#   anch.H3K4me3                -1.692e-02  5.857e-02  -0.289 0.772620    
# vp.H3K56ac                   1.966e-02  6.001e-02   0.328 0.743184    
# anch.H3K56ac                -1.247e-02  6.231e-02  -0.200 0.841442    
# vp.RNAPol                   -3.892e-02  4.289e-02  -0.907 0.364162    
# anch.RNAPol                  6.090e-02  4.029e-02   1.511 0.130673    
# vp.CTCFinGene                6.724e-01  2.380e-01   2.825 0.004729 ** 
#   anch.CTCFinGene              2.604e-01  9.268e-02   2.809 0.004964 ** 
#   ctcf.ctfc                   -6.828e-01  4.342e-01  -1.573 0.115820    
# prom.enh                     1.698e-01  8.399e-02   2.022 0.043201 *  
#   prom.ctcf                   -1.018e-01  8.339e-02  -1.221 0.222187    
# active.gene                  9.200e-01  9.332e-02   9.859  < 2e-16 ***
#   vp.Hap1_CTCF_chip            5.977e-03  9.614e-03   0.622 0.534108    
# anch.Hap1_CTCF_chip          1.404e-02  9.987e-03   1.406 0.159783    
# vp.Hap1_H3k27ac_chip        -1.458e-02  7.801e-03  -1.869 0.061680 .  
# anch.Hap1_H3k27ac_chip      -1.241e-02  7.610e-03  -1.631 0.102921    
# vp.Hap1_SMC1_chip           -1.683e-02  2.863e-02  -0.588 0.556588    
# anch.Hap1_SMC1_chip         -3.857e-02  2.836e-02  -1.360 0.173748    
# vp.Wapl3_3_CTCF_chip        -7.137e-03  1.193e-02  -0.598 0.549697    
# anch.Wapl3_3_CTCF_chip      -1.800e-02  1.235e-02  -1.457 0.145068    
# vp.Wapl3_3_SMC1_chip         5.235e-04  1.471e-02   0.036 0.971606    
# anch.Wapl3_3_SMC1_chip       7.782e-03  1.488e-02   0.523 0.600904    
# vp.Hap1_ATAC_chip            6.340e-03  9.360e-03   0.677 0.498158    
# anch.Hap1_ATAC_chip          6.040e-03  9.281e-03   0.651 0.515155    
# vp.Hap1_H2AK119ub1_chip     -4.756e-02  4.049e-02  -1.174 0.240249    
# anch.Hap1_H2AK119ub1_chip   -8.465e-02  4.086e-02  -2.072 0.038279 *  
#   vp.Hap1_H3K27me3_chip       -9.774e-02  3.655e-02  -2.674 0.007496 ** 
#   anch.Hap1_H3K27me3_chip     -6.167e-02  3.604e-02  -1.711 0.087071 .  
# max.Hap1_CTCF_chip          -1.785e-02  1.068e-02  -1.672 0.094469 .  
# max.Hap1_H3k27ac_chip        9.974e-03  8.938e-03   1.116 0.264460    
# max.Hap1_SMC1_chip           3.472e-02  3.131e-02   1.109 0.267501    
# max.Wapl3_3_CTCF_chip        2.216e-02  1.332e-02   1.665 0.095990 .  
# max.Wapl3_3_SMC1_chip       -7.256e-03  1.616e-02  -0.449 0.653405    
# max.Hap1_ATAC_chip           7.493e-04  1.067e-02   0.070 0.944027    
# max.Hap1_H2AK119ub1_chip     6.810e-02  4.592e-02   1.483 0.138098    
# max.Hap1_H3K27me3_chip       7.909e-02  4.141e-02   1.910 0.056150 .  
# vp.LFC_CTCF_chip             1.220e-01  9.020e-02   1.352 0.176311    
# anch.LFC_CTCF_chip           1.192e-01  9.464e-02   1.260 0.207666    
# max.LFC_CTCF_chip           -3.095e-01  1.295e-01  -2.390 0.016831 *  
#   vp.LFC_SMC1_chip             4.335e-02  9.807e-02   0.442 0.658461    
# anch.LFC_SMC1_chip          -1.856e-02  1.013e-01  -0.183 0.854590    
# max.LFC_SMC1_chip           -8.182e-03  1.297e-01  -0.063 0.949687    
# loopdepth                    9.961e-03  7.670e-03   1.299 0.194066    
# TADid                        8.221e-05  3.295e-05   2.495 0.012603 *  
#   TADstretch                  -6.924e-02  8.014e-02  -0.864 0.387570    
# max.CTCF_WT                 -4.846e-01  4.317e-01  -1.123 0.261596    
# max.SMC1_WT                 -4.878e-02  1.228e-01  -0.397 0.691266    
# max.H3K27ac                  5.151e-02  5.926e-02   0.869 0.384682    
# max.H3K4me3                  5.556e-02  6.264e-02   0.887 0.375082    
# max.H3K56ac                 -5.312e-03  6.479e-02  -0.082 0.934659    
# max.RNAPol                   6.809e-02  4.645e-02   1.466 0.142697    
# max.CTCFinGene                      NA         NA      NA       NA    
# window.Hap1_CTCF_chip        1.421e-01  1.138e+00   0.125 0.900674    
# window.Hap1_H3k27ac_chip     9.196e-01  1.146e-01   8.023 1.03e-15 ***
#   window.Hap1_SMC1_chip        3.146e+00  7.335e-01   4.289 1.79e-05 ***
#   window.Wapl3_3_CTCF_chip     3.476e-01  1.301e+00   0.267 0.789357    
# window.Wapl3_3_SMC1_chip    -2.490e+00  5.623e-01  -4.429 9.46e-06 ***
#   window.Hap1_ATAC_chip       -4.324e-01  2.475e-01  -1.747 0.080614 .  
# window.Hap1_H2AK119ub1_chip -5.910e-01  1.371e-01  -4.311 1.62e-05 ***
#   window.Hap1_H3K27me3_chip    3.039e-01  1.048e-01   2.900 0.003732 ** 
#   window.LFC_CTCF_chip        -4.789e-01  7.778e-01  -0.616 0.538057    
# window.LFC_SMC1_chip         2.551e+00  4.340e-01   5.878 4.15e-09 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
# Null deviance: 11051.4  on 8177  degrees of freedom
# Residual deviance:  9746.5  on 8103  degrees of freedom
# (2956 observations deleted due to missingness)
# AIC: 9896.5
# 
# Number of Fisher Scoring iterations: 4
