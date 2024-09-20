#devtools::install_github("stephenrho/pminternal")

library(pminternal)
library(glmnet)
library(caret)

setwd("/Users/fabioaffaticati/Desktop/Work/activ_covid-and-omics")

data <- read.csv("pminternal.csv")

data <- data[, which(colnames(data) != "Subjectnr")]
data <- data[, which(colnames(data) != "Group")]
data$labels <- ifelse(data$labels == 'C 1', 1, 0)


dmy <- dummyVars(" ~ .", data = data)
data <- data.frame(predict(dmy, newdata = data))

#data$Genderm<-as.factor(data$Genderm)
#data$Genderf<-as.factor(data$Genderf)
#data$labels<-factor(data$labels)

lasso_fun <- function(data, ...){
  
  y <- data$labels
  x <- as.matrix(data[, which(colnames(data) != "labels")])
  
  glmnet(x=x, y=y, family=binomial(link='logit'))
}

lasso_predict <- function(model, data, ...){
  
  y <- data$labels
  x <- as.matrix(data[, which(colnames(data) != "labels")])
  
  predict(model, newx = x, type = "response")[,1]
}


(val <- validate(data = data, outcome = "labels", 
                 model_fun = lasso_fun, pred_fun = lasso_predict, 
                 method = "boot_simple", B = 10))









stepglm <- function(data, ...){
  
  m <- cv.glmnet(x=x, y=y, family=binomial(link='logit'))
}

steppred <- function(model, data, ...){
  predict(model, newdata = data, type = "response")
}

validate(data = gusto, outcome = "y", model_fun = stepglm, 
         pred_fun = steppred, method = "boot_simple", B = 10)



prediction_stability(val, smooth_bounds = TRUE)

mape_stability(val)

calibration_stability(val)
