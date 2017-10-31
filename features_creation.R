library(data.table)
library(dplyr)

# knn creator: function to return target value smoothed over k nearest neighbors
knn_creator <- function(features, train_flag, target_tot, k) {
  
  knn_train_tot <- FNN::knn.reg(train = features[train_flag, ], test = data_tot_matrix[train_flag, ], y = target_tot, k = k+1)$pred
  knn_train_tot <- (knn_train_tot*(k+1)-target_tot)/k
  
  knn_holdout <- FNN::knn.reg(train = features[train_flag, ], test = features[!train_flag, ], y = target_tot, k = k)$pred
  return(c(knn_train_tot, knn_holdout))
}


# cross creator: returns interactions between nbins of
cross_creator <- function(var1, var2, bins, sparse = T) {
  
  var1_bin  <- cut(var1, breaks = sort(unique(quantile(var1, probs = 0:bins/bins))), include.lowest = T, labels = F)
  var2_bin<- cut(var2, breaks = sort(unique(quantile(var2, probs = 0:bins/bins))), include.lowest = T, labels = F)
  
  var_cross <- paste(var1_bin, var2_bin, sep="_")
  
  if (sparse) { 
    return(sparse.model.matrix(~.-1, data.frame(var_cross)))
  } else {
    return(model.matrix(~.-1, data.frame(var_cross)))
  }
}


# cat_estimates: returns the mean target associated with given each cat level
cat_estimates <- function(var, target, train_flag, bayes = F, weight = NULL) {
  
  if (is.null(weight)){
    data = data.frame(var = as.factor(var), target = target, train_flag=train_flag, weight = 1)
  } else {
    data = data.frame(var = as.factor(var), target = target, train_flag=train_flag, weight = weight)
  }
  
  # bayes
  # VHM: variance of the hypothetical means
  data_VHM <- data %>% 
    filter(train_flag) %>% 
    dplyr::group_by(var) %>% 
    summarise(target = weighted.mean(target, weight),
              weight = sum(weight)) %>% 
    ungroup() %>% 
    summarise(VHM = sum(weight/sum(weight) * (target - weighted.mean(target, weight)) ^ 2))
  VHM <- data_VHM$VHM
  
  # EVPV: expected value of the process variance
  data_EVPV <- data %>% 
    filter(train_flag) %>% 
    dplyr::group_by(var) %>% 
    summarise(EVPV = sum(weight/sum(weight) * (target - weighted.mean(target, weight)) ^ 2),
              target = weighted.mean(target, weight),
              weight = sum(weight)) %>% 
    ungroup() %>% 
    mutate(target_mean = weighted.mean(target, weight),
           EVPV = weighted.mean(EVPV, weight) / weight,
           Z = VHM / (VHM + EVPV),
           target_smooth = Z * target + (1 - Z) * target_mean) 
  
  target_smooth <- data_EVPV$target_smooth
  names(target_smooth) <- data_EVPV$var
  
  return(target_smooth[data$var])
}

# cat_estimates: returns the mean target associated with given each cat level
cat_bayes_estimates <- function(var, target, train_flag, bayes = T, weight = NULL) {
  
  if (is.null(weight)){
    data = data.frame(var = as.character(var), x = target, train_flag=train_flag, w = 1, stringsAsFactors = F)
  } else {
    data = data.frame(var = as.character(var), x = target, train_flag=train_flag, w = weight, stringsAsFactors = F)
  }

  if(bayes) {
    # VHM: variance of the hypothetical means
    data_cred_train <- data %>% 
      filter(train_flag) %>% 
      # mutate(w = w / sum(w)) %>% 
      dplyr::group_by(var) %>% 
      mutate(x_cat = sum(w*x),
             x2_cat = sum(w*x^2),
             w_cat = sum(w),
             mean_cat = x_cat/w_cat,
             PV_cat = x2_cat/w_cat - mean_cat^2,
             w_cat_adj = w_cat-w,
             mean_cat_adj = (x_cat - w*x)/(w_cat_adj),
             PV_cat_adj = (x2_cat - w*x^2)/(w_cat_adj) - mean_cat_adj^2) %>% 
      ungroup() %>% 
      mutate(w_tot =sum(w),
             mean_cat_tot = sum(w*mean_cat)/w_tot,
             mean_tot_adj = (sum(w*x) - x*w)/(w_tot-w),
             mean_cat2_tot = sum(w*mean_cat^2)/w_tot,
             VHM = mean_cat2_tot - mean_cat_tot^2,
             VHM_adj = (w_tot*mean_cat2_tot - w_cat*mean_cat^2 + w_cat_adj*mean_cat_adj^2)/(w_tot - w) - 
               ((w_tot*mean_cat_tot - w_cat*mean_cat + w_cat_adj*mean_cat_adj)/(w_tot - w))^2,
             EVPV = sum(w*PV_cat)/w_tot, 
             EVPV_adj = (w_tot*EVPV - w_cat*PV_cat + w_cat_adj * PV_cat_adj)/(w_tot-w)/w_cat_adj,
             Z = VHM_adj / (VHM_adj + EVPV_adj),
             x_smooth = Z * mean_cat_adj + (1- Z) * mean_tot_adj,
             x_smooth = ifelse(is.na(x_smooth), mean_tot_adj, x_smooth))
    

    
  } else {
    data_cred_train <- data %>% 
      filter(train_flag) %>% 
      dplyr::group_by(var) %>% 
      mutate(x_cat = sum(w*x),
             w_cat = sum(w),
             mean_cat = x_cat/w_cat,
             w_cat_adj = w_cat-w,
             x_smooth = (x_cat - w*x)/(w_cat_adj))
  }
  
  # attach metric on test dataset
  cat_smooth_key <- data_cred_train %>% 
    group_by(var) %>% summarise(x_smooth = mean(x_smooth)) %>% 
    as.data.table()
  setkeyv(cat_smooth_key, cols = "var")
  
  data_cred_test <- data %>% 
    filter(!train_flag) %>% 
    as.data.table()
  
  data_cred_test <- cat_smooth_key[data_cred_test]
  
  cat_smooth <- numeric(nrow(data))
  cat_smooth[train_flag] <- data_cred_train$x_smooth
  cat_smooth[!train_flag] <- data_cred_test$x_smooth
  cat_smooth[is.na(cat_smooth)] <- mean(cat_smooth, na.rm=T)
  
  return(cat_smooth)
}
