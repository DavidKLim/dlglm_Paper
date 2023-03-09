library("dlglm")
source("run_mice.R")

phi0=5; pi=0.3; sim_indexes=1
miss_pct_features = 50
mechanisms=c("MNAR","MAR","MCAR")
case="x"
miss_cols=NULL; ref_cols=NULL

datasets=c("DRYBEAN","LETTER","SHUTTLE","RED","WHITE","BANKNOTE","SPAM")  # something doesn't match up with RED. 6 cats in test, 5 in train?

NL_r = F

learn_r = T
covars = "ally"
normalize=T
early_stop=T

trace=F; draw_miss=T   # default.
methods=c("dlglm","idlglm","mice")

init_r = "default"   # or "alt"
dlglm_pref = if(init_r=="default"){""}else if(init_r=="alt"){"alt_init/"}



for(d in 1:length(datasets)){
  dataset = datasets[d]
  print(dataset)
  #### DEFINE data_types_x, data_type_y, family, link, Cy, C, by datasets[d]
  # datasets=c("SPAM","IRIS","ADULT","WINE","RED","WHITE","YEAST","CONCRETE","BANKNOTE")


  for(m in 1:length(mechanisms)){for(s in 1:length(sim_indexes)){

    iNL_r = if(NL_r){"NL_"}else{""}
    dir_name = sprintf("Results_%s/miss_%s%s/phi%d/sim%d", dataset, iNL_r, case, phi0, sim_indexes[s])

    fname = sprintf("%s/%sres_dlglm_%s_%d_%d.RData",dir_name,dlglm_pref,mechanisms[m],miss_pct_features,pi*100)
    ifname = sprintf("%s/Ignorable/res_dlglm_%s_%d_%d.RData",dir_name,mechanisms[m],miss_pct_features,pi*100)
    fname_mice = sprintf("%s/res_mice_%s_%d_%d.RData",dir_name,mechanisms[m],miss_pct_features,pi*100)

    load( sprintf("%s/data_%s_%d_%d.RData", dir_name, mechanisms[m], miss_pct_features, pi*100) )  # loads "X","Y","mask_x","mask_y","g"
    load( sprintf("%s/params_%s_%d_%d.RData", dir_name, mechanisms[m], miss_pct_features, pi*100))  # sim.params, miss.params, sim.data, sim.mask

    if(dataset %in% c("SPAM","IRIS","WINE","RED","WHITE","YEAST","BANKNOTE",
                      "LETTER","SHUTTLE", "DRYBEAN","CAREVALUATION","RAISIN")){
      C=NULL
      Cy = length(unique(Y))
      data_type_y = "cat"
      family="Multinomial"; link="mlogit"
    } else if(dataset%in%c("CONCRETE","ABALONE")){
      C=NULL
      Cy = NULL
      data_type_y = "real"
      family="Gaussian"; link="identity"
    }
    #### UCI hyperparams
    hyperparams = list(sigma="elu", bss=c(1000L), lrs=c(0.01,0.001), impute_bs = 1000L, arch="IWAE",
                       niws=5L, n_imps = 500L, n_epochss=2002L,
                       # n_hidden_layers = c(0L,1L,2L), n_hidden_layers_y = c(0L), n_hidden_layers_r = c(0L,1L),  # simulations
                       n_hidden_layers = c(1L,2L,4L), n_hidden_layers_y = c(0L,1L,2L), n_hidden_layers_r = c(0L,1L),  # Physionet
                       h=c(128L,64L), h_y=NULL, h_r=c(16L,32L),
                       # dim_zs = as.integer(floor(c(ncol(X)/4, ncol(X)/2, 3*ncol(X)/4))), # simulations
                       dim_zs = as.integer(floor(  c(8L, ncol(X)/4, ncol(X)/2, 3*ncol(X)/4)  )),  # Physionet
                       L1_weights = 0)


    if(dataset=="CAREVALUATION"){
      data_types_x=rep("cat", ncol(X))
    }else if(dataset=="ABALONE"){
      data_types_x=c("cat",rep("real", ncol(X)-1))
    }else{ data_types_x=rep("real", ncol(X)) }

    N=nrow(X); P=ncol(X); D=NA


    dir_name0 = sprintf("%s/%s%s_%d_%d", dir_name, dlglm_pref, mechanisms[m], miss_pct_features, pi*100)
    idir_name0 = sprintf("%s/Ignorable/%s_%d_%d", dir_name, mechanisms[m], miss_pct_features, pi*100)

    covars_r_y = if(grepl("y",covars)){ 1L }else{ 0L }
    covars_r_x = if(grepl("miss",covars)){ (colMeans(mask_x)!=1)^2    # SELECT COVARIATES FOR MNAR
    } else if(grepl("obs", covars)){ (colMeans(mask_x)==1)^2
    } else if(grepl("all", covars)){ covars_r_x = rep(1L,P) }


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
      res_mice = run_mice(dir_name, sim.params, miss.params, X, Y, mask_x, mask_y, g, data_types_x, family, niws = niws)
      save(res_mice, file=fname_mice)
      rm(res_mice)
    }
    # if(!file.exists(fname_miwae) & "miwae" %in% methods){
    #   res_mice = run_miwae(dir_name, sim.params, miss.params, X, Y, mask_x, mask_y, g, data_types_x, niws = niws)
    #   save(res_mice, file=fname_mice)
    #   rm(res_mice)
    # }
    # if(!file.exists(fname_notmiwae) & "notmiwae" %in% methods){
    #   res_mice = run_notmiwae(dir_name, sim.params, miss.params, X, Y, mask_x, mask_y, g, data_types_x, niws = niws)
    #   save(res_mice, file=fname_mice)
    #   rm(res_mice)
    # }
  }}
}




# sbatch -p gpu --gres=gpu:1 --qos=gpu_access -N 1 -n 1 --mem=24g -t 3- -o /pine/scr/d/e/deelim/dump/dlglm.out -J sim_dlglm --wrap='R CMD BATCH run_dlglm.R'

