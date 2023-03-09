
library("dlglm")
source("run_mice.R")
source("run_miwae.R")
Ns=c(1e4,1e5); Ps=c(25,50); Ds=c(2,8)
family="Multinomial"; C=2
link = if(family=="Gaussian"){ "identity"
} else if(family=="Multinomial"){ "mlogit"
} else if(family=="Poisson"){ "log" }

beta=0.25
case = "x"; pi=0.3; phi0=5; miss_pct_features = 50

learn_r = T#; covars_r_x = rep(1L,P); covars_r_y = 1L  # include all
covars = "ally"   # "obs", "missy", "obsy", "all", "ally"
normalize=F; mechanisms=c("MNAR","MAR","MCAR"); sim_indexes = 1:30
early_stop=T
trace=F; draw_miss=T   # default.
methods=c("dlglm","idlglm","mice","miwae","notmiwae"); niws=50
init_r = "default"   # or "alt"
dlglm_pref = if(init_r=="default"){""}else if(init_r=="alt"){"alt_init/"}

## run linear miss model simulations:
NL_x=F; NL_y=F; NL_r=F
for(ns in 1:length(Ns)){for(ps in 1:length(Ps)){for(ds in 1:length(Ds)){
  data_types_x = c( rep("real",Ps[ps]) )
  data_type_x_pref = if(all(data_types_x == data_types_x[1])){ data_types_x[1] } else{ "mixed" }
  data_type_y_pref = if(family=="Gaussian"){"real"}else if(family=="Multinomial"){"cat"}else if(family=="Poisson"){"cts"}

  prefix = sprintf("Xr%dct%dcat%d_beta%f_pi%d/",
                   sum(data_types_x=="real"), sum(data_types_x=="count"), sum(data_types_x=="cat"),
                   beta, pi*100)
  dataset = sprintf("SIM_N%d_P%d_D%d", Ns[ns], Ps[ps],Ds[ds])
  for(m in 1:length(mechanisms)){for(s in 1:length(sim_indexes)){

    iNL_x = if(NL_x){"NL"}else{""}
    iNL_y = if(NL_y){"NL"}else{""}
    iNL_r = if(NL_r){"NL_"}else{""}
    dir_name = sprintf("Results_%sX%s_%sY%s/%s%s/miss_%s%s/phi%d/sim%d",
                       iNL_x, data_type_x_pref, iNL_y, data_type_y_pref, prefix, dataset, iNL_r, case, phi0, sim_indexes[s])

    fname = sprintf("%s/%sres_dlglm_%s_%d_%d.RData",dir_name,dlglm_pref,mechanisms[m],miss_pct_features,pi*100)
    ifname = sprintf("%s/Ignorable/res_dlglm_%s_%d_%d.RData",dir_name,mechanisms[m],miss_pct_features,pi*100)
    fname_mice = sprintf("%s/res_mice_%s_%d_%d.RData",dir_name,mechanisms[m],miss_pct_features,pi*100)
    fname_miwae = sprintf("%s/res_miwae_%s_%d_%d.RData",dir_name,mechanisms[m],miss_pct_features,pi*100)
    fname_notmiwae = sprintf("%s/res_notmiwae_%s_%d_%d.RData",dir_name,mechanisms[m],miss_pct_features,pi*100)

    load( sprintf("%s/data_%s_%d_%d.RData", dir_name, mechanisms[m], miss_pct_features, pi*100) )  # loads "X","Y","mask_x","mask_y","g"
    load( sprintf("%s/params_%s_%d_%d.RData", dir_name, mechanisms[m], miss_pct_features, pi*100))  # sim.params, miss.params, sim.data, sim.mask

    dir_name0 = sprintf("%s/%s%s_%d_%d", dir_name, dlglm_pref, mechanisms[m], miss_pct_features, pi*100)
    idir_name0 = sprintf("%s/Ignorable/%s_%d_%d", dir_name, mechanisms[m], miss_pct_features, pi*100)

    covars_r_y = if(grepl("y",covars)){ 1L }else{ 0L }
    covars_r_x = if(grepl("miss",covars)){ (colMeans(mask_x)!=1)^2    # SELECT COVARIATES FOR MNAR
    } else if(grepl("obs", covars)){ (colMeans(mask_x)==1)^2
    } else if(grepl("all", covars)){ covars_r_x = rep(1L,P) }

    hyperparams = list(sigma="elu", bss=c(1000L), lrs=c(0.01,0.001), impute_bs = 1000L, arch="IWAE",
                       niws=5L, n_imps = 500L, n_epochss=2002L, n_hidden_layers = c(0L,1L,2L), n_hidden_layers_y = c(0L), n_hidden_layers_r = c(0L,1L),
                       h=c(128L,64L), h_y=NULL, h_r=c(16L,32L),
                       dim_zs = c(as.integer(floor(ncol(X)/12)),as.integer(floor(ncol(X)/4)), as.integer(floor(ncol(X)/2)), as.integer(floor(3*ncol(X)/4))),
                       L1_weights = 0)

    if(!file.exists(fname) & "dlglm" %in% methods){
      ifelse(!dir.exists(dir_name0), dir.create(dir_name0, recursive=T), F)
      res = dlglm(dir_name0, X, Y, mask_x, mask_y, g,
                  covars_r_x, covars_r_y, learn_r, data_types_x, Ignorable=F,
                  family, link, normalize, early_stop, trace, draw_miss, init_r=init_r, hyperparams=hyperparams)
      save(res, file=fname)
      rm(res)
    }

    if(!file.exists(ifname) & "idlglm" %in% methods){
      ifelse(!dir.exists(idir_name0), dir.create(idir_name0, recursive=T), F)

      res = dlglm(idir_name0, X, Y, mask_x, mask_y, g,
                  covars_r_x, covars_r_y, learn_r, data_types_x, Ignorable=T,
                  family, link, normalize, early_stop, trace, draw_miss, hyperparams=hyperparams)
      save(res, file=ifname)
      rm(res)
    }

    if(!file.exists(fname_mice) & "mice" %in% methods){
      res_mice = run_mice(dir_name, sim.params, miss.params, X, Y, mask_x, mask_y, g, data_types_x, niws = niws)
      save(res_mice, file=fname_mice)
      rm(res_mice)
    }

    if(!file.exists(fname_miwae) & "miwae" %in% methods){
      res_mice = run_miwae(dir_name, sim.params, miss.params, X, Y, mask_x, mask_y, g, data_types_x, niws = niws)
      save(res_mice, file=fname_mice)
      rm(res_mice)
    }

    if(!file.exists(fname_notmiwae) & "notmiwae" %in% methods){
      res_mice = run_notmiwae(dir_name, sim.params, miss.params, X, Y, mask_x, mask_y, g, data_types_x, niws = niws)
      save(res_mice, file=fname_mice)
      rm(res_mice)
    }
  }}
}}}


## run nonlinear miss model simulations:
NL_x=F; NL_y=F; NL_r=T
for(ns in 1:length(Ns)){for(ps in 1:length(Ps)){for(ds in 1:length(Ds)){
  data_types_x = c( rep("real",Ps[ps]) )
  data_type_x_pref = if(all(data_types_x == data_types_x[1])){ data_types_x[1] } else{ "mixed" }
  data_type_y_pref = if(family=="Gaussian"){"real"}else if(family=="Multinomial"){"cat"}else if(family=="Poisson"){"cts"}

  prefix = sprintf("Xr%dct%dcat%d_beta%f_pi%d/",
                   sum(data_types_x=="real"), sum(data_types_x=="count"), sum(data_types_x=="cat"),
                   beta, pi*100)
  dataset = sprintf("SIM_N%d_P%d_D%d", Ns[ns], Ps[ps],Ds[ds])
  for(m in 1:length(mechanisms)){for(s in 1:length(sim_indexes)){

    iNL_x = if(NL_x){"NL"}else{""}
    iNL_y = if(NL_y){"NL"}else{""}
    iNL_r = if(NL_r){"NL_"}else{""}
    dir_name = sprintf("Results_%sX%s_%sY%s/%s%s/miss_%s%s/phi%d/sim%d",
                       iNL_x, data_type_x_pref, iNL_y, data_type_y_pref, prefix, dataset, iNL_r, case, phi0, sim_indexes[s])

    fname = sprintf("%s/%sres_dlglm_%s_%d_%d.RData",dir_name,dlglm_pref,mechanisms[m],miss_pct_features,pi*100)
    ifname = sprintf("%s/Ignorable/res_dlglm_%s_%d_%d.RData",dir_name,mechanisms[m],miss_pct_features,pi*100)
    fname_mice = sprintf("%s/res_mice_%s_%d_%d.RData",dir_name,mechanisms[m],miss_pct_features,pi*100)
    fname_miwae = sprintf("%s/res_miwae_%s_%d_%d.RData",dir_name,mechanisms[m],miss_pct_features,pi*100)
    fname_notmiwae = sprintf("%s/res_notmiwae_%s_%d_%d.RData",dir_name,mechanisms[m],miss_pct_features,pi*100)

    load( sprintf("%s/data_%s_%d_%d.RData", dir_name, mechanisms[m], miss_pct_features, pi*100) )  # loads "X","Y","mask_x","mask_y","g"
    load( sprintf("%s/params_%s_%d_%d.RData", dir_name, mechanisms[m], miss_pct_features, pi*100))  # sim.params, miss.params, sim.data, sim.mask

    dir_name0 = sprintf("%s/%s%s_%d_%d", dir_name, dlglm_pref, mechanisms[m], miss_pct_features, pi*100)
    idir_name0 = sprintf("%s/Ignorable/%s_%d_%d", dir_name, mechanisms[m], miss_pct_features, pi*100)

    covars_r_y = if(grepl("y",covars)){ 1L }else{ 0L }
    covars_r_x = if(grepl("miss",covars)){ (colMeans(mask_x)!=1)^2    # SELECT COVARIATES FOR MNAR
    } else if(grepl("obs", covars)){ (colMeans(mask_x)==1)^2
    } else if(grepl("all", covars)){ covars_r_x = rep(1L,P) }

    hyperparams = list(sigma="elu", bss=c(1000L), lrs=c(0.01,0.001), impute_bs = 1000L, arch="IWAE",
                       niws=5L, n_imps = 500L, n_epochss=2002L, n_hidden_layers = c(0L,1L,2L), n_hidden_layers_y = c(0L), n_hidden_layers_r = c(0L,1L),
                       h=c(128L,64L), h_y=NULL, h_r=c(16L,32L),
                       dim_zs = c(as.integer(floor(ncol(X)/12)),as.integer(floor(ncol(X)/4)), as.integer(floor(ncol(X)/2)), as.integer(floor(3*ncol(X)/4))),
                       L1_weights = 0)

    if(!file.exists(fname) & "dlglm" %in% methods){
      ifelse(!dir.exists(dir_name0), dir.create(dir_name0, recursive=T), F)
      res = dlglm(dir_name0, X, Y, mask_x, mask_y, g,
                  covars_r_x, covars_r_y, learn_r, data_types_x, Ignorable=F,
                  family, link, normalize, early_stop, trace, draw_miss, init_r=init_r, hyperparams=hyperparams)
      save(res, file=fname)
      rm(res)
    }

    if(!file.exists(ifname) & "idlglm" %in% methods){
      ifelse(!dir.exists(idir_name0), dir.create(idir_name0, recursive=T), F)

      res = dlglm(idir_name0, X, Y, mask_x, mask_y, g,
                  covars_r_x, covars_r_y, learn_r, data_types_x, Ignorable=T,
                  family, link, normalize, early_stop, trace, draw_miss, hyperparams=hyperparams)
      save(res, file=ifname)
      rm(res)
    }

    if(!file.exists(fname_mice) & "mice" %in% methods){
      res_mice = run_mice(dir_name, sim.params, miss.params, X, Y, mask_x, mask_y, g, data_types_x, niws = niws)
      save(res_mice, file=fname_mice)
      rm(res_mice)
    }

    if(!file.exists(fname_miwae) & "miwae" %in% methods){
      res_mice = run_miwae(dir_name, sim.params, miss.params, X, Y, mask_x, mask_y, g, data_types_x, niws = niws)
      save(res_mice, file=fname_mice)
      rm(res_mice)
    }

    if(!file.exists(fname_notmiwae) & "notmiwae" %in% methods){
      res_mice = run_notmiwae(dir_name, sim.params, miss.params, X, Y, mask_x, mask_y, g, data_types_x, niws = niws)
      save(res_mice, file=fname_mice)
      rm(res_mice)
    }
  }}
}}}


# sbatch -p gpu --gres=gpu:1 --qos=gpu_access -N 1 -n 1 --mem=24g -t 3- -o /pine/scr/d/e/deelim/dump/dlglm.out -J sim_dlglm --wrap='R CMD BATCH run_dlglm.R'

