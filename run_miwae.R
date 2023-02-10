source_python("MIWAE.py")    # MIWAE class
source_python("notMIWAE.py") # notMIWAE class
source_python("trainer.py")  # "train" function
source_python("utils.py")    # "imputationRMSE" and "not_imputationRMSE" functions
np = import("numpy")
tf = import("tensorflow")

# ### TESTING
# 
# ## things needed for both
# X = np$array(matrix(rnorm(1000,0,1),nrow=1000,ncol=10))
# Xnan = X; Xnan[1:5,1:5] = NA  # need NA missing version
# X0 = Xnan; X0[is.na(Xnan)]=0 # zero preimputed version
# Xval = np$array(matrix(rnorm(1000,0,1),nrow=100,ncol=10))
# Xvalnan = Xval; Xvalnan[1:5,1:5] = NA
# Xval0 = Xvalnan; Xval0[is.na(Xvalnan)]=0
# 
# R = np$array(!np$isnan(Xnan), dtype=np$float)


# hyperparams = list(sigma="elu", bss=c(1000L), lrs=c(0.01,0.001), impute_bs = 1000L, arch="IWAE",
#                    niws=5L, n_imps = 500L, n_epochss=2002L, n_hidden_layers = c(0L,1L,2L), n_hidden_layers_y = c(0L), n_hidden_layers_r = c(0L,1L),
#                    h=c(128L,64L), h_y=NULL, h_r=c(16L,32L),
#                    dim_zs = c(as.integer(floor(ncol(X)/12)),as.integer(floor(ncol(X)/4)), as.integer(floor(ncol(X)/2)), as.integer(floor(3*ncol(X)/4))),
#                    L1_weights = 0)
### for MIWAE  # learning rate is fixed at 0.001 in their code. also, no option to change # HLs
run_miwae = function(dir_name, X, Y, mask_x, mask_y, g, hyperparams, niws = 500L){
  X = as.matrix(X)
  Xnan = X; Xnan[mask_x==0] = NA
  X0 = X; X0[mask_x==0] = 0
  
  Xs = split(data.frame(X), g)
  X0s =  split(data.frame(X0), g)
  Xnans =  split(data.frame(Xnan), g)
  Rxs = split(data.frame(mask_x), g)
  R = mask_x
  ## unsupervised method: doesn't utilize Y or Ry
  # Ys = split(data.frame(Y), g)
  # Rys = split(data.frame(mask_y), g)
  
  ## hyperparameters to tune
  bs = hyperparams$bss
  d = hyperparams$dim_zs
  h = hyperparams$h
  lrs = hyperparams$lrs
  n_epochs = hyperparams$n_epochss
  
  opt_val_loss = Inf
  for(i in 1:length(bs)){for(j in 1:length(d)){for(k in 1:length(h)){for(l in 1:length(lrs)){
    mod_name = sprintf("%s/miwae-Z%d_h%d_bs%d_lr%f",dir_name, d[j], h[k], bs[i], lrs[l])
    miwae = MIWAE(X=np$array(as.matrix(Xnans$train)), Xval=np$array(as.matrix(Xnans$valid)),
                  n_latent=d[j], n_samples=5L, n_hidden=h[k], lr=lrs[l], name=mod_name)
    train(miwae, batch_size=bs[i], max_iter=n_epochs)
    
    val_loss = miwae$val_batch()
    
    if(val_loss < opt_val_loss){
      opt_d = d[j]; opt_h = h[k]; opt_bs = bs[i]; opt_lr = lrs[l]
      opt_mod_name = mod_name
      miwae$save(sprintf("%s/saved_model",mod_name))
      # tf$keras$models$save_model(sprintf("%s/saved_model.pth",mod_name))
      opt_val_loss = val_loss
    }
  }}}}
  
  print("Finished training")
  miwae_test = MIWAE(X=np$array(as.matrix(Xnans$test)), Xval=np$array(as.matrix(Xnans$valid)),
                n_latent=opt_d, n_samples=5L, n_hidden=opt_h, lr=opt_lr, name=opt_mod_name)
  print(opt_mod_name)
  miwae_test$load(sprintf("%s/saved_model",opt_mod_name))
  
  print("Loaded optimal model")
  ## on test set
  res_MIWAE = imputationRMSE(miwae_test,np$array(X),np$array(X0),np$array(Xnan),np$array(R),niws)  # put in Xs$test here, with Rx$test
  print("Computed importance-weighted xhat")
  res_MIWAE[[2]]  # imputations
  res_MIWAE[[2]]-X
  Xhat = res_MIWAE[[2]]; Xhat[R==1] = X[R==1]
  
  return(Xhat)
}

### for nOT MIWAE
run_notmiwae = function(dir_name, X, Y, mask_x, mask_y, g, mprocess="linear", hyperparams, niws = 500L){
  X = as.matrix(X)
  Xnan = X; Xnan[mask_x==0] = NA
  X0 = X; X0[mask_x==0] = 0
  
  Xs = split(data.frame(X), g)
  X0s =  split(data.frame(X0), g)
  Xnans =  split(data.frame(Xnan), g)
  Rxs = split(data.frame(mask_x), g)
  R = mask_x
  ## unsupervised method: doesn't utilize Y or Ry
  # Ys = split(data.frame(Y), g)
  # Rys = split(data.frame(mask_y), g)
  
  ## hyperparameters to tune
  bs = hyperparams$bss
  d = hyperparams$dim_zs
  h = hyperparams$h
  lrs = hyperparams$lrs
  n_epochs = hyperparams$n_epochss
  
  opt_val_loss = Inf
  for(i in 1:length(bs)){for(j in 1:length(d)){for(k in 1:length(h)){for(l in 1:length(lrs)){
    mod_name = sprintf("%s/notmiwae-Z%d_h%d_bs%d_lr%f", dir_name, d[j], h[k], bs[i], lrs[l])
    notmiwae = notMIWAE(X=np$array(as.matrix(Xnans$train)), Xval=np$array(as.matrix(Xnans$valid)),
                        n_latent=d[j], n_samples=5L, n_hidden=h[k], missing_process=mprocess, lr=lrs[l], name=mod_name)
    train(notmiwae, batch_size=bs[i], max_iter=n_epochs)
    val_loss = notmiwae$val_batch()
    if(val_loss < opt_val_loss){
      opt_d = d[j]; opt_h = h[k]; opt_bs = bs[i]; opt_lr = lrs[l]
      opt_mod_name = mod_name
      notmiwae$save(sprintf("%s/saved_model",mod_name))
      # notmiwae$save_model(sprintf("%s/saved_model",mod_name))
      opt_val_loss = val_loss
    }
  }}}}
  print("Finished training")
  
  notmiwae_test = notMIWAE(X=np$array(as.matrix(Xnans$test)), Xval=np$array(as.matrix(Xnans$valid)),
                        n_latent=opt_d, n_samples=5L, n_hidden=opt_h, missing_process=mprocess, lr=opt_lr, name=opt_mod_name)
  print(opt_mod_name)
  notmiwae_test$load(sprintf("%s/saved_model",opt_mod_name))
  print("Loaded optimal model")
  
  res_notMIWAE = not_imputationRMSE(notmiwae_test,np$array(X),np$array(X0),np$array(Xnan),np$array(R),niws)
  print("Computed importance-weighted xhat")
  res_notMIWAE[[2]]  # imputations
  res_notMIWAE[[2]]-X   # not 0 for nonmissing entries: these are samples for observed entries too.
  Xhat = res_notMIWAE[[2]]; Xhat[R==1] = X[R==1]
  return(Xhat)
}