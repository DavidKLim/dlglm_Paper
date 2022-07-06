run_mice = function(dir_name, sim.params, miss.params, X, Y, mask_x, mask_y, g, data_types_x, niws = 500){
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(gtable)
  library(mice)
  
  N=sim.params$N; P=sim.params$P; data_types_x=sim.params$data_types_x; family=sim.params$family
  
  link=if(family=="Gaussian"){"identity"}else if(family=="Multinomial"){"mlogit"}else if(family=="Poisson"){"log"}
  
  pi = miss.params$pi
  mechanism=miss.params$mechanism
  sim_index=sim.params$sim_index
  miss_pct_features = miss.params$miss_pct_features
  
  data_type_x = if(all(data_types_x==data_types_x[1])){data_types_x[1]}else{"mixed"}
  
  # dir_name = sprintf("Results_X%s_Y%s/%s%s/miss_%s/phi%d/sim%d", data_type_x,data_type_y,prefix, dataset, case, miss.params$phi0, sim_index)   # to save interim results
  # ifelse(!dir.exists(sprintf("%s/Diagnostics",dir_name)), dir.create(sprintf("%s/Diagnostics",dir_name)), F)
  
  data.fname = sprintf("%s/data_%s_%d_%d.RData", dir_name, mechanism, miss_pct_features, pi*100)
  print(paste("Data file: ", data.fname))
  if(!file.exists(data.fname)){
    stop("Data file does not exist..")
  } else{
    load(data.fname)
  }
  
  if(family=="Multinomial"){
    levels_Y = levels(factor(Y))   # convert into numeric.
    Y = as.numeric(factor(Y))
    Y = Y-min(Y)
  }else{
    levels_Y = NA
  }   # set to start from 0, not 1
  
  data_types_x_0 = data_types_x
  # mask_x = (res$mask_x)^2; mask_y = (res$mask_y)^2
  if(sum(data_types_x=="cat") == 0){
    X_aug = X
    mask_x_aug = mask_x
  } else{
    # reorder to real&count covariates first, then augment cat dummy vars
    X_aug = X[,!(data_types_x %in% c("cat"))]
    mask_x_aug = mask_x[,!(data_types_x %in% c("cat"))]
    
    ## onehot encode categorical variables
    # X_cats = X[,data_types_x=="cat"]
    Cs = rep(0,sum(data_types_x=="cat"))
    # X_cats_onehot = matrix(nrow=N,ncol=0)
    cat_ids = which(data_types_x=="cat")
    for(i in 1:length(cat_ids)){
      X_cat = as.numeric(as.factor(X[,cat_ids[i]]))-1
      Cs[i] = length(unique(X_cat))
      X_cat_onehot = matrix(ncol = Cs[i], nrow=length(X_cat))
      for(ii in 1:Cs[i]){
        X_cat_onehot[,ii] = (X_cat==ii-1)^2
      }
      # X_cats_onehot = cbind(X_cats_onehot, X_cat_onehot)
      X_aug = cbind(X_aug, X_cat_onehot)
      mask_x_aug = cbind(mask_x_aug, matrix(mask_x[,cat_ids[i]], nrow=N, ncol=Cs[i]))
    }
    
    
    ## column bind real/count and one-hot encoded cat vars
    data_types_x = c( data_types_x[!(data_types_x %in% c("cat"))], rep("cat",sum(Cs)) )
  }
  
  Xs = split(data.frame(X), g)        # split by $train, $test, and $valid
  Xs_aug = split(data.frame(X_aug), g)        # split by $train, $test, and $valid
  Ys = split(data.frame(Y), g)        # split by $train, $test, and $valid
  Rxs = split(data.frame(mask_x), g)
  Rys = split(data.frame(mask_y), g)
  
  
  X_MICE = X; X_MICE[mask_x==0]=NA
  Xs_MICE = split(data.frame(X_MICE), g)
  Y_MICE = Y; Y_MICE[mask_y==0]=NA
  Ys_MICE = split(data.frame(Y_MICE), g)
  
  dat_MICE = data.frame(cbind(Xs_MICE$train,Ys_MICE$train)); colnames(dat_MICE)[ncol(dat_MICE)] = "y"
  # res_MICE = mice::mice(dat_MICE, m=res$train_params$L * res$train_params$M)
  res_MICE = mice::mice(dat_MICE, m=niws)
  dat_MICE_test = data.frame(cbind(Xs_MICE$test,Ys_MICE$test)); colnames(dat_MICE_test)[ncol(dat_MICE_test)] = "y"
  res_MICE_test = mice::mice(dat_MICE_test, m=niws)
  
  xyhat_mice = matrix(0,nrow=nrow(Xs$test),ncol=ncol(Xs$test))
  for(i in 1:niws){
    xyhat_mice = xyhat_mice + mice::complete(res_MICE_test,i)/niws
  }  # mean across imputations
  xhat_mice = xyhat_mice[,-ncol(xyhat_mice)]
  
  res = list(xhat_mice=xhat_mice, xyhat_mice=xyhat_mice,
             res_MICE=res_MICE, res_MICE_test=res_MICE_test,
             levels_Y=levels_Y)
  
  return(res)
}