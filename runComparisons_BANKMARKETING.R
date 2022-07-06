library("dlglm")
source("run_mice.R")

learn_r = T
covars = "ally"
early_stop=T

dataset = "BANKMARKETING"

trace=F; draw_miss=T   # default.

unbalanced = F
normalize=T

methods=c("dlglm","idlglm","mice"); niws=50 ## 500?

dir_name = sprintf("Results_BANKMARKETING/miss_x/phiNA/sim1")

if(unbalanced){
  fname = sprintf("%s/unbalanced/res_dlglm.RData",dir_name)
  ifname = sprintf("%s/unbalanced/Ignorable/res_dlglm.RData",dir_name)
  ifelse(!dir.exists(sprintf("%s/unbalanced",dir_name)), dir.create(sprintf("%s/unbalanced",dir_name), recursive=T), F)
  ifelse(!dir.exists(sprintf("%s/unbalanced/Ignorable",dir_name)), dir.create(sprintf("%s/unbalanced/Ignorable",dir_name), recursive=T), F)
  dir_name0 = sprintf("%s/unbalanced/NA_NA_NA", dir_name)
  idir_name0 = sprintf("%s/unbalanced/Ignorable/NA_NA_NA", dir_name)
}else{
  fname = sprintf("%s/res_dlglm.RData",dir_name)
  ifname = sprintf("%s/Ignorable/res_dlglm.RData",dir_name)
  dir_name0 = sprintf("%s/NA_NA_NA", dir_name)
  idir_name0 = sprintf("%s/Ignorable/NA_NA_NA", dir_name)
}
fname_mice = sprintf("%s/res_mice.RData",dir_name)

mechanism=NA; miss_pct_features=NA; pi=NA
load( sprintf("%s/data_%s_%d_%d.RData", dir_name, mechanism, miss_pct_features, pi*100) )  # loads "X","Y","mask_x","mask_y","g"

P = ncol(X); N= nrow(X)
# data_types_x = c("real",rep("cat",7),rep("real",8))     # try this too?
data_types_x = rep("real",P)
family = "Multinomial"; link = "mlogit"
sim.params=list(N=N, P=P, data_types_x=data_types_x, family=family, sim_index=1);  miss.params = list(pi=NA, mechanism=NA, miss_pct_features=NA)



covars_r_y = if(grepl("y",covars)){ 1L }else{ 0L }
covars_r_x = if(grepl("miss",covars)){ (colMeans(mask_x)!=1)^2    # SELECT COVARIATES FOR MNAR
} else if(grepl("obs", covars)){ (colMeans(mask_x)==1)^2
} else if(grepl("all", covars)){ covars_r_x = rep(1L,P) }

## Physionet/UCI hyperparams
hyperparams = list(sigma="elu", bss=c(1000L), lrs=c(0.01,0.001), impute_bs = 1000L, arch="IWAE",
                   niws=5L, n_imps = 500L, n_epochss=2002L,
                   # n_hidden_layers = c(0L,1L,2L), n_hidden_layers_y = c(0L), n_hidden_layers_r = c(0L,1L),  # simulations
                   # n_hidden_layers = c(1L,2L,4L), n_hidden_layers_y = c(0L), n_hidden_layers_r = c(0L,1L,2L,4L),  # Physionet
                   n_hidden_layers = c(1L,2L,4L), n_hidden_layers_y = c(0L,1L,2L), n_hidden_layers_r = c(0L,1L,2L),  # Physionet
                   h=c(128L,64L), h_y=NULL, h_r=c(16L,32L),
                   # dim_zs = as.integer(floor(c(ncol(X)/4, ncol(X)/2, 3*ncol(X)/4))), # simulations
                   dim_zs = as.integer(floor(  c(8L, ncol(X)/4, ncol(X)/2, 3*ncol(X)/4)  )),  # Physionet
                   L1_weights = 0)

if(!file.exists(fname) & "dlglm" %in% methods){
  ifelse(!dir.exists(dir_name0), dir.create(dir_name0, recursive=T), F)
  res = dlglm(dir_name0, X, Y, mask_x, mask_y, g,
              covars_r_x, covars_r_y, learn_r, data_types_x, Ignorable=F,
              family, link, normalize, early_stop, trace, draw_miss, unbalanced, hyperparams=hyperparams)
  save(res, file=fname)
  rm(res)
}

if(!file.exists(ifname) & "idlglm" %in% methods){
  ifelse(!dir.exists(idir_name0), dir.create(idir_name0, recursive=T), F)

  res = dlglm(idir_name0, X, Y, mask_x, mask_y, g,
              covars_r_x, covars_r_y, learn_r, data_types_x, Ignorable=T,
              family, link, normalize, early_stop, trace, draw_miss, unbalanced, hyperparams=hyperparams)
  save(res, file=ifname)
  rm(res)
}

if(!file.exists(fname_mice) & "mice" %in% methods){

  # sim.params = list(N=N, D=D, P=P, data_types=data_types, family=family, sim_index=sim_indexes[s], ratios=c(train=.8,valid=.1,test=.1),
  #                   beta=beta, C=C, Cy=Cy, NL_x=NL_x, NL_y=NL_y)
  # miss.params=list(scheme="UV", mechanism=mechanisms[m], NL_r=NL_r, pi=pi, phi0=phi0, miss_pct_features=miss_pct_features, miss_cols=miss_cols, ref_cols=ref_cols)
  res_mice = run_mice(dir_name, sim.params, miss.params, X, Y, mask_x, mask_y, g, data_types_x, family, data_type_y_pref, niws = niws)
  save(res_mice, file=fname_mice)
  rm(res_mice)
}




# sbatch -p gpu --gres=gpu:1 --qos=gpu_access -N 1 -n 1 --mem=24g -t 3- -o /pine/scr/d/e/deelim/dump/dlglm.out -J sim_dlglm --wrap='R CMD BATCH run_dlglm.R'

