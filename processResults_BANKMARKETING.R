unbalanced = F
# unbalanced = T
pref = if(!unbalanced){""}else{"unbalanced/"}

dir_name = "Results_BANKMARKETING/miss_x/phiNA/sim1"
load(sprintf("%s/data_NA_NA_NA.RData",dir_name))
load(sprintf("%s/%sIgnorable/res_dlglm.RData",dir_name,pref))
ires=res
load(sprintf("%s/%sIgnorable/NA_NA_NA/with_normalization/opt_train.out",dir_name,pref))
ires_train = res_train
load(sprintf("%s/%sres_dlglm.RData",dir_name,pref))
load(sprintf("%s/%sNA_NA_NA/with_normalization/opt_train.out",dir_name,pref))
load(sprintf("%s/res_mice.RData",dir_name))

X_aug = X
mask_x_aug = mask_x
Y = as.numeric(as.factor(Y))
Y = Y-min(Y)
Y_aug = Y
mask_y_aug = mask_y

Xs = split(data.frame(X_aug), g)        # split by $train, $test, and $valid
Ys = split(data.frame(Y_aug), g)        # split by $train, $test, and $valid
Rxs = split(data.frame(mask_x_aug), g)
Rys = split(data.frame(mask_y_aug), g)

metrics = data.frame(matrix(nrow=4,ncol=11))
rownames(metrics) = c("dlglm","idlglm","mean","mice")
colnames(metrics) = c("ARI","PA","AUC","sens","PPV","FPR", "kappa", "brier","brierScaled","F1","MCC")

compute_sens_PP = function(truth,pred, prs){
  sens = sum(truth==1 & pred==1)/sum(truth==1)
  PP = sum(truth==1 & pred==1)/sum(pred==1)
  FPR = sum(truth==0 & pred==1)/sum(truth==0)
  kappa = measures::KAPPA(truth,pred)
  brier = measures::Brier(prs, truth, 0, 1)
  brierScaled = measures::BrierScaled(prs, truth, 0, 1)
  f1 = measures::F1(truth, pred, 1)
  mcc = yardstick::mcc_vec(truth=factor(truth), estimate=factor(pred))
  return(list(sens=sens, PP=PP, FPR=FPR, kappa=kappa,
              brier=brier, brierScaled=brierScaled, f1=f1, mcc=mcc))
}
#### dlglm
library(mclust)
library(nnet)
library(yardstick)
prs = res$results$all_params$y$probs; prs_train = res_train$all_params$y$probs
# cutoffs = seq(0.01,1,0.01)
# PPs_train = matrix(nrow=length(cutoffs),ncol=8)  # 8 metrics
# colnames(PPs_train) = c("sens","PP","FPR","kappa","brier","brierScaled","f1","MCC")
# rownames(PPs_train) = cutoffs
# for(c in 1:length(cutoffs)){
#   pred_classes_train = apply(prs_train,1,function(x){(x[2]>cutoffs[c])^2})
#   PPs_train[c,] = unlist(compute_sens_PP(Ys$train$Y_aug, pred_classes_train, prs_train))
# }
opt_cutoff = 0.2
# pred_classes = apply(prs,1, which.max)-1   # 0.5 cutoff
pred_classes = apply(prs,1, function(x){(x[2]>opt_cutoff)^2})
table(pred_classes,Ys$test$Y_aug)
ARI = adjustedRandIndex(pred_classes,Ys$test$Y_aug)
PA = mean(pred_classes==Ys$test$Y_aug)
library(pROC)
library(measures)
AUC = pROC::auc(Ys$test$Y_aug,prs[,2])
PP = compute_sens_PP(Ys$test$Y_aug, pred_classes, prs)
metrics[1,] = c(ARI, PA, AUC, unlist(PP))

#### idlglm
iprs = ires$results$all_params$y$probs; iprs_train = ires_train$all_params$y$probs
# PPs_train = matrix(nrow=length(cutoffs),ncol=8)  # 8 metrics
# colnames(PPs_train) = c("sens","PP","FPR","kappa","brier","brierScaled","f1","MCC")
# rownames(PPs_train) = cutoffs
# for(c in 1:length(cutoffs)){
#   pred_classes_train = apply(iprs_train,1,function(x){(x[2]>cutoffs[c])^2})
#   PPs_train[c,] = unlist(compute_sens_PP(Ys$train$Y_aug, pred_classes_train, iprs_train))
# }
# opt_cutoff = 0.2
opt_cutoff = 0.27
# ipred_classes = apply(iprs,1, which.max)-1
ipred_classes = apply(iprs,1, function(x){(x[2]>opt_cutoff)^2})
table(ipred_classes,Ys$test$Y_aug)
ARI = adjustedRandIndex(ipred_classes,Ys$test$Y_aug)
PA = mean(ipred_classes==Ys$test$Y_aug)
AUC = pROC::auc(Ys$test$Y_aug, iprs[,2])
PP = compute_sens_PP(Ys$test$Y_aug, ipred_classes, iprs)
metrics[2,] = c(ARI, PA, AUC, unlist(PP))

#### MEAN
X_mean = X_aug; Y_mean = Y_aug
for(i in 1:ncol(X_mean)){
  # print(mean(X[mask_x[,i]==1,i]))
  if(is.factor(X_mean[,i])){
    X_mean[mask_x[,i]==0,i] = names(table(X_mean[,i]))[which.max(table(X_mean[,i]))]
  }else{
    X_mean[mask_x[,i]==0,i] = mean(X[mask_x[,i]==1,i])
  }
}
Y_mean[mask_y==0] = mean(Y[mask_y==1])
Xs_mean = split(data.frame(X_mean), g)
Ys_mean = split(data.frame(Y_mean), g)
xhat_mean = Xs_mean$test; yhat_mean = Ys_mean$test

dat = cbind( Xs_mean$train ,Ys_mean$train )
fit_mean = multinom(as.factor(Y_mean) ~ 0+., data=dat)

# pred_probs_train_mean = predict(fit_mean, type="probs", cbind(Xs_mean$train,row.names = NULL))
# PPs_train = matrix(nrow=length(cutoffs),ncol=8)  # 8 metrics
# colnames(PPs_train) = c("sens","PP","FPR","kappa","brier","brierScaled","f1","MCC")
# rownames(PPs_train) = cutoffs
# for(c in 1:length(cutoffs)){
#   pred_classes_train = (pred_probs_train_mean> cutoffs[c])^2 #apply(iprs_train,1,function(x){(x[2]>cutoffs[c])^2})
#   PPs_train[c,] = unlist(compute_sens_PP(Ys$train$Y_aug, pred_classes_train, pred_probs_train_mean))
# }

# opt_cutoff = 0.2
opt_cutoff = 0.21

pred_probs_mean = predict(fit_mean, newdata=cbind(Xs_mean$test,row.names = NULL), type="probs")
# pred_classes_mean = predict(fit_mean, newdata=cbind(Xs_mean$test,row.names = NULL))  # 0.5 cutoff
pred_classes_mean = (pred_probs_mean > opt_cutoff)^2

table(pred_classes_mean,Ys$test$Y_aug)
ARI = adjustedRandIndex(pred_classes_mean,Ys$test$Y_aug)
PA = mean(pred_classes_mean==Ys$test$Y_aug)
AUC = pROC::auc(Ys$test$Y_aug, pred_probs_mean)
PP = compute_sens_PP(Ys$test$Y_aug, pred_classes_mean, pred_probs_mean)
metrics[3,] = c(ARI, PA, AUC, unlist(PP))


#### MICE
library(mice)
### USE MICE ON TRAIN SET --> FIT GLM --> PREDICT ON PRED SET USING POOLED ESTIMATES

fits_MICE = list()
for(i in 1:res_mice$res_MICE$m){
  fits_MICE[[i]] = multinom(y ~ 0+., data=complete(res_mice$res_MICE,i))
  fits_MICE[[i]]$data = NULL; fits_MICE[[i]]$model=NULL
}
fit_MICE = pool(fits_MICE)
dummy_fit_MICE = fits_MICE[[1]]; rm(fits_MICE)
dummy_fit_MICE$coefficients = fit_MICE$pooled$estimate

Xs_test_mice = complete(res_mice$res_MICE_test,1)
for(j in 2:res_mice$res_MICE_test$m){
  Xs_test_mice = Xs_test_mice + complete(res_mice$res_MICE_test,j)
}
Xs_test_mice = Xs_test_mice/res_mice$res_MICE_test$m          # average of multiply-imputed test sets

Xs_train_mice = complete(res_mice$res_MICE,1)
for(j in 2:res_mice$res_MICE$m){
  Xs_train_mice = Xs_train_mice + complete(res_mice$res_MICE,j)
}
Xs_train_mice = Xs_train_mice/res_mice$res_MICE_test$m

# pred_probs_train_mice = predict(dummy_fit_MICE, type="probs", newdata=Xs_train_mice)
# PPs_train = matrix(nrow=length(cutoffs),ncol=8)  # 8 metrics
# colnames(PPs_train) = c("sens","PP","FPR","kappa","brier","brierScaled","f1","MCC")
# rownames(PPs_train) = cutoffs
# for(c in 1:length(cutoffs)){
#   pred_classes_train = (pred_probs_train_mice> cutoffs[c])^2 #apply(iprs_train,1,function(x){(x[2]>cutoffs[c])^2})
#   PPs_train[c,] = unlist(compute_sens_PP(Ys$train$Y_aug, pred_classes_train, pred_probs_train_mice))
# }
# opt_cutoff = 0.2
opt_cutoff = 0.21

pred_probs_mice = predict(dummy_fit_MICE, newdata = Xs_test_mice, type="probs")    # 0.5 cutoff
# pred_classes_mice = predict(dummy_fit_MICE, newdata = Xs_test_mice)
pred_classes_mice = (pred_probs_mice > opt_cutoff)^2

table(pred_classes_mice,Ys$test$Y_aug)
ARI = adjustedRandIndex(pred_classes_mice,Ys$test$Y_aug)
PA = mean(pred_classes_mice==Ys$test$Y_aug)
AUC = pROC::auc(Ys$test$Y_aug, pred_probs_mice)
PP = compute_sens_PP(Ys$test$Y_aug, pred_classes_mice, pred_probs_mice)
metrics[4,] = c(ARI, PA, AUC, unlist(PP))



metrics
# ### USE MICE ON TEST SET (PERFECT PERFORMANCFE: this would be wrong since we're feeding correct y)
# pred_classes_mice = complete(res_mice$res_MICE_test,1)$y
# table(pred_classes_mice,Ys$test$Y_aug)
# mean(pred_classes_mice==Ys$test$Y_aug)
# pROC::auc(Ys$test$Y_aug, pred_classes_mice)

save(metrics, file=sprintf("%s/%smetrics.out",dir_name,pref))


## looking at the coef est's
w = res$results$w   # 0 HLs, 2 x 16.
betas = w[2,]-w[1,]
names(betas) = names(X)
betas  # print coef ests for 0 HL res
