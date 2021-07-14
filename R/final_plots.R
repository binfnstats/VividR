resultsData = bind_rows(data4, data2, data3)

data4 = timeValsLM[,-1] %>% melt()
colnames(data4) = c("Fold", "Run", "Value")
data4$Metric = "Time"
data4$Method = "LIMMA"

data2 = accValsLM[,-1] %>% melt()
colnames(data2) = c("Fold", "Run", "Value")
data2$Metric = "ACC"
data2$Method = "LIMMA"

data3 = aucValsLM[,-1] %>% melt()
colnames(data3) = c("Fold", "Run", "Value")
data3$Metric = "AUC"
data3$Method = "LIMMA"


######################
accValsRFE = matrix(data = 0,
                  nrow = k,
                  ncol = repB)
aucValsRFE = matrix(data = 0,
                  nrow = k,
                  ncol = repB)
timeValsRFE = matrix(data = 0,
                   nrow = k,
                   ncol = repB)


kFolds = readRDS("final_folds.rds")

B = 100
k = 10
repB = 5
setSeed = 1234567
nCores = parallel::detectCores() - 1
compareMethod = "BIC"

for(i in c(1,2,4,5)){
  kFold = kFolds[[i]]
  features2 = list()
  
  for(j in 1:k){
    trainTest = unlist(kFolds[[i]][-j])
    xTrain = x[trainTest,]
    yTrain = y[trainTest]
    xTest = x[-trainTest,]
    yTest = y[-trainTest]
    
    preProcValues = caret::preProcess(xTrain, method = c("center", "scale"))
    xTrain = stats::predict(preProcValues, xTrain)
    xTest = stats::predict(preProcValues, xTest)
    
    tic()
    RFE = sigFeature::sigFeature(X = xTrain,
                                 Y = yTrain)
    time = toc()
    saveRDS(RFE,paste0("RFE_fold_",i,"_",j,".rds",collapse=""))
    newX = xTrain[,RFE > (p - topN)]
    RFEGLM = glmnet::cv.glmnet(newX, 
                               yTrain, 
                               alpha = 0, 
                               family = "binomial")
    newXTest = xTest[,RFE > (p - topN)]
    RFEPred = predict(RFEGLM,
                      s = "lambda.1se",
                      newx = newXTest)
    fitted = c(exp(RFEPred)/(1+exp(RFEPred)))
    
    
    modelPred = 1*(fitted > 0.5)
    binary = yTest
    accValsRFE[j,i] = mean(1*(modelPred == binary))
    roc = pROC::roc(binary, fitted)
    aucValsRFE[j,i] = roc$auc
    timeValsRFE[j, i] = time$toc - time$tic
  }
}

######################

# PLot




borutaData = readRDS('borutaResults.rds')
VIVIDData = readRDS('VIVIDResults.rds')
LIMMAData = readRDS('LIMMAResults.rds')

LIMMAData %>% 
  group_by(Metric, Method) %>% 
  summarise(mean_value = mean(Value))

finalData = bind_rows(VIVIDData, borutaData, LIMMAData)

ggplot(finalData, aes(x = Method, y = Value, colour = Method)) + 
  geom_boxplot() +
  facet_wrap(~Metric, scale = "free") + 
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle=90, hjust=1)) + 
  ylab(NULL)

