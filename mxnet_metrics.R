mx.metric.mlogloss <- mx.metric.custom("mlogloss", function(label, pred){
  label <- as.factor(label)
  #label_mat <- matrix(0, nrow = length(label), ncol = length(levels(label)))
  label_mat <- matrix(0, nrow = length(label), ncol = length(pred) / length(label) )
  
  sample_levels <- as.integer(label)
  for (i in 1:length(label)) label_mat[i, sample_levels[i]] <- 1
  label <- label_mat
  
  label <- as.vector(t(label))
  eps <- 1e-15
  nrow <- nrow(pred)
  batch <- ncol(pred)
  pred <- pmax(pmin(pred, 1 - eps), eps)
  MultiLogLoss <- (-1/batch) * sum(label * log(pred))
  return(MultiLogLoss)
})


mx.metric.logloss <- mx.metric.custom("logloss", function(label, pred){
  eps <- 1e-15
  batch <- length(pred)
  pred <- pmax(pmin(pred, 1 - eps), eps)
  LogLoss <- (-1/batch) * sum(label * log(pred))
  return(LogLoss)
})


mx.metric.ExpMAE <- mx.metric.custom("ExpMAE", function(label, pred){
  label <- exp(label)
  pred<- exp(pred)
  eps <- 1e-15
  batch <- ncol(pred)
  ExpMAE <- (1/batch) * sum(abs(label-pred))
  return(ExpMAE)
})
