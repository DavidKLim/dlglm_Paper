# source('processResults.R')
# source('processResults_script.R')
library(dlglm)
run_processScript = function(beta=0.25, pi=0.3, miss_pct_features = 50,
                             Ns, Ps, Ds, phi0=5, case="x", data_type_y="real", Cy=NULL,
                             sim_indexes=1:5, mechanisms=c("MCAR","MAR","MNAR"),
                             NL_x=F, NL_y=F, NL_r=F, methods=c("mean","mice","idlglm","dlglm"),
                             init_r="default"){
  family = if(data_type_y=="real"){"Gaussian"}else if(data_type_y=="cat"){"Multinomial"}else if(data_type_y=="cts"){"Poisson"}
  dlglm_pref = if(init_r=="default"){""}else if(init_r=="alt"){"alt_init/"}



  for(i in 1:length(Ns)){for(j in 1:length(Ps)){for(k in 1:length(Ds)){
  N=Ns[i]; P=Ps[j]; D=Ds[k]
  data_types_x = rep("real",P); C=NULL
  data_type_x = if(all(data_types_x==data_types_x[1])){data_types_x[1]}else{"mixed"}
  P_real=sum(data_types_x=="real"); P_count=sum(data_types_x=="count"); P_cat=sum(data_types_x=="cat")
  prefix=sprintf("Xr%dct%dcat%d_beta%f_pi%d/",P_real,P_count,P_cat,beta,pi*100)

  iNL_x = if(NL_x){"NL"}else{""}
  iNL_y = if(NL_y){"NL"}else{""}
  iNL_r = if(NL_r){"NL_"}else{""}
  dataset = sprintf("SIM_N%d_P%d_D%d", N, P, D)
  df = data.frame()
  df_ratio = data.frame()

  dir_name = sprintf("Results_%sX%s_%sY%s/%s%s/miss_%s%s/phi%d/%s",
                     iNL_x,data_type_x,iNL_y,data_type_y,
                     prefix, dataset,
                     iNL_r, case, phi0, dlglm_pref)
  ifelse(!dir.exists(dir_name), dir.create(dir_name), F)
  for(m in 1:length(mechanisms)){

    fname = sprintf("%slist_res_%s.out",
                    dir_name, mechanisms[m])
    if(file.exists(fname)){
      load(fname)
    } else{
      list_res = list()
      for(s in 1:length(sim_indexes)){
        sim_index=sim_indexes[s]; mechanism=mechanisms[m]
        sim.params = list(N=N, P=P, D=D, data_types=rep(data_type_x,P),
                          family=family,
                          sim_index=sim_index, ratios=c(train=.8,valid=.1,test=.1),
                          beta=beta, C=C, Cy=Cy, NL_x=NL_x, NL_y=NL_y)
        miss.params = list(scheme="UV", mechanism=mechanism, pi=pi, phi0=phi0, miss_pct_features=miss_pct_features,NL_r=F)

        list_res[[s]] = dlglm::processResults(prefix=prefix,data.file.name = data.file.name, mask.file.name=mask.file.name,
                                       sim.params = sim.params,
                                       miss.params = miss.params,
                                       case=case,
                                       data_types_x=data_types_x, data_type_y = data_type_y, methods=methods, init_r=init_r)
        ## list_res[[1]]$imp$L1s
        ## list_res[[1]]$inf$df
        ## list_res[[1]]$pred$complete
        ## list_res[[1]]$pred$imputed
      }
      save(list_res, file=fname)
    }

    #### PLOT IMPUTATION barplots with error bars
    df1 = lapply(list_res,function(x) x$imp$L1s)
    df1 = lapply(df1, function(x) x[methods])
    # List of 5
    # $ :List of 5
    # ..$ zero  : num [1:7381] 1.131 0.833 0.774 1.805 0.2 ...
    # ..$ mean  : num [1:7381] 1.325 1.623 1.682 0.651 2.256 ...
    # ..$ dlglm : num [1:7381] 0.0365 0.3191 0.0306 0.6849 0.51 ...
    # ..$ idlglm: num [1:7381] 1.297 1.676 1.445 0.619 1.554 ...
    # ..$ mice  : num [1:7381] 1.1 1.663 1.233 0.692 1.621 ...
    # $ :List of 5
    # ..$ zero  : num [1:7405] 0.309 1.238 1.271 1.809 1.111 ...
    # ..$ mean  : num [1:7405] 2.164 1.235 1.202 0.664 1.362 ...
    # ..$ dlglm : num [1:7405] 0.791 0.271 0.192 0.786 0.192 ...
    # ..$ idlglm: num [1:7405] 2.086 0.933 1.012 0.756 1.206 ...
    # ..$ mice  : num [1:7405] 2.42 1.054 0.857 0.733 1.281 ...
    # $ :List of 5
    # ..$ zero  : num [1:7576] 0.509 1.539 1.596 1.85 0.711 ...
    # ..$ mean  : num [1:7576] 1.954 0.924 0.867 0.613 1.752 ...
    # ..$ dlglm : num [1:7576] 0.73 0.486 0.397 0.792 0.225 ...
    # ..$ idlglm: num [1:7576] 1.754 0.474 0.694 0.219 1.227 ...
    # ..$ mice  : num [1:7576] 1.553 0.718 0.684 0.098 1.338 ...
    # $ :List of 5
    # ..$ zero  : num [1:7480] 0.731 1.769 1.315 0.973 0.838 ...
    # ..$ mean  : num [1:7480] 1.716 0.678 1.131 1.474 1.609 ...
    # ..$ dlglm : num [1:7480] 0.408 0.734 0.142 0.113 0.45 ...
    # ..$ idlglm: num [1:7480] 1.801 0.598 1.326 1.319 0.886 ...
    # ..$ mice  : num [1:7480] 1.552 0.342 1.236 1.361 1.206 ...
    # $ :List of 5
    # ..$ zero  : num [1:7556] 1.36 1.131 0.817 1.471 0.615 ...
    # ..$ mean  : num [1:7556] 1.104 1.333 1.647 0.994 1.849 ...
    # ..$ dlglm : num [1:7556] 0.0587 0.2456 0.0523 0.2801 0.0509 ...
    # ..$ idlglm: num [1:7556] 1.201 0.746 0.906 0.651 0.908 ...
    # ..$ mice  : num [1:7556] 0.955 0.777 0.92 0.825 0.981 ...
    avg_L1 = sapply(df1, function(x) sapply(x,mean))
    ## compute SD across sim_indexes?

    #### PLOT INFERENCE: mean/sd PB, RB
    df2 = lapply(list_res,function(x) x$inf$df)
    # List of 5
    # $ : num [1:51, 1:15] 0.07566 -0.03848 0.04127 -0.00362 -0.04606 ...
    # ..- attr(*, "dimnames")=List of 2
    # .. ..$ : chr [1:51] "1" "2" "3" "4" ...
    # .. ..$ : chr [1:15] "RB_zero" "RB_mean" "RB_dlglm" "RB_idlglm" ...
    # $ : num [1:51, 1:15] -0.0481 0.0823 0.0208 -0.1189 0.0805 ...
    # ..- attr(*, "dimnames")=List of 2
    # .. ..$ : chr [1:51] "1" "2" "3" "4" ...
    # .. ..$ : chr [1:15] "RB_zero" "RB_mean" "RB_dlglm" "RB_idlglm" ...
    # $ : num [1:51, 1:15] 0.03542 0.05888 0.10666 0.00946 -0.06864 ...
    # ..- attr(*, "dimnames")=List of 2
    # .. ..$ : chr [1:51] "1" "2" "3" "4" ...
    # .. ..$ : chr [1:15] "RB_zero" "RB_mean" "RB_dlglm" "RB_idlglm" ...
    # $ : num [1:51, 1:15] 0.00735 -0.00848 -0.0431 0.0197 0.09996 ...
    # ..- attr(*, "dimnames")=List of 2
    # .. ..$ : chr [1:51] "1" "2" "3" "4" ...
    # .. ..$ : chr [1:15] "RB_zero" "RB_mean" "RB_dlglm" "RB_idlglm" ...
    # $ : num [1:51, 1:15] -0.1323 0.1096 0.0875 -0.0212 -0.1344 ...
    # ..- attr(*, "dimnames")=List of 2
    # .. ..$ : chr [1:51] "1" "2" "3" "4" ...
    # .. ..$ : chr [1:15] "RB_zero" "RB_mean" "RB_dlglm" "RB_idlglm" ...

    df2_avg = sapply(df2,function(x)x[nrow(df2[[1]]),])

    metric = "PB"  # PB, RB, or SE. PB is most informative
    avg_metric = df2_avg[grepl(metric,rownames(df2_avg)),]
    avg_metric = avg_metric[match(methods,sapply(strsplit(rownames(avg_metric),"_"),function(x)x[2])), ]
    ## compute SD across sim_indexes?

    #### PLOT PREDICTION: barplots with error bars (complete and imputed)
    ref_method = "mean"
    id_ref = which(methods == ref_method)          # reference: dlglm. Show all metrics divided by metric of dlglm (relative)
    if(family=="Multinomial"){
      df3 = lapply(list_res,function(x) unlist(lapply(x$pred$complete,mean)));     df3 = lapply(df3, function(x) x[methods])
      df4 = lapply(list_res,function(x) unlist(lapply(x$pred$imputed,mean)));     df4 = lapply(df4, function(x) x[methods])
      # df3 = lapply(list_res,function(x) unlist(lapply(x$pred$complete,function(y){mean(-log(y))})));     df3 = lapply(df3, function(x) x[methods])
      # df4 = lapply(list_res,function(x) unlist(lapply(x$pred$imputed,function(y){mean(-log(y))})));     df4 = lapply(df4, function(x) x[methods])
      df5 = lapply(list_res,function(x) unlist(x$pred$AUC_complete));     df5 = lapply(df5, function(x) x[methods])
      df6 = lapply(list_res,function(x) unlist(x$pred$AUC_imputed));     df6 = lapply(df6, function(x) x[methods])
      df7 = lapply(list_res,function(x) unlist(x$pred$Briers1));     df7 = lapply(df7, function(x) x[methods])
      df8 = lapply(list_res,function(x) unlist(x$pred$Briers2));     df8 = lapply(df8, function(x) x[methods])


      # [[1]]
      # zero      mean      dlglm     idlglm       mice
      # [1,] 0.06591256 0.1039643 0.02416338 0.09320132 0.09634262
      #
      # [[2]]
      # zero       mean     dlglm     idlglm       mice
      # [1,] 0.03251658 0.03979739 0.0200064 0.03396019 0.03456542
      #
      # [[3]]
      # zero      mean      dlglm     idlglm      mice
      # [1,] 0.07217582 0.1061206 0.02238534 0.09978853 0.1056182
      #
      # [[4]]
      # zero       mean      dlglm     idlglm       mice
      # [1,] 0.06602664 0.08124393 0.03557768 0.09308788 0.07734531
      #
      # [[5]]
      # zero       mean     dlglm     idlglm       mice
      # [1,] 0.04061242 0.04114375 0.0225036 0.03072247 0.03074503
      avg_predC_L1 = sapply(df3,function(x) x )
      avg_predI_L1 = sapply(df4,function(x) x )
      avg_AUCC = sapply(df5,function(x) x )
      avg_AUCI = sapply(df6,function(x) x )
      avg_BrierC = sapply(df7,function(x) x )
      avg_BrierI = sapply(df8,function(x) x )

      # add_df = cbind(apply(avg_L1,1,mean), apply(avg_L1,1,sd),
      #                apply(avg_metric,1,mean), apply(avg_metric,1,sd),
      #                apply(avg_predC_L1,1,mean), apply(avg_predC_L1,1,sd),
      #                apply(avg_predI_L1,1,mean), apply(avg_predI_L1,1,sd))             # avg_L1, avg_metric, avg_pred1_L1, avg_pred2_L1
      # colnames(add_df) = c("mean_L1","sd_L1",
      #                      sprintf("mean_%s",metric), sprintf("sd_%s",metric),
      #                      "mean_predC_L1","sd_predC_L1","mean_predI_L1","sd_predI_L1")

      add_df = cbind(apply(avg_L1,1,mean), apply(avg_L1,1,sd),
                     apply(avg_metric,1,mean), apply(avg_metric,1,sd),
                     apply(avg_predC_L1,1,mean), apply(avg_predC_L1,1,sd),
                     apply(avg_predI_L1,1,mean), apply(avg_predI_L1,1,sd),
                     apply(avg_AUCC,1,mean), apply(avg_AUCC,1,sd),
                     apply(avg_AUCI,1,mean), apply(avg_AUCI,1,sd),
                     apply(avg_BrierC, 1, mean), apply(avg_BrierC, 1, sd),
                     apply(avg_BrierI, 1, mean), apply(avg_BrierI, 1, sd))             # avg_L1, avg_metric, avg_pred1_L1, avg_pred2_L1
      colnames(add_df) = c("mean_L1","sd_L1",
                           sprintf("mean_%s",metric), sprintf("sd_%s",metric),
                           "mean_predC_L1","sd_predC_L1","mean_predI_L1","sd_predI_L1",
                           "mean_AUCC","sd_AUCC","mean_AUCI","sd_AUCI", "mean_BrierC", "sd_BrierC","mean_BrierI","sd_BrierI")


      add_df_ratio = cbind(apply(avg_L1[-id_ref,]/avg_L1[id_ref,],1,mean), apply(avg_L1[-id_ref,]/avg_L1[id_ref,],1,sd),
                           apply(avg_metric[-id_ref,]/avg_metric[id_ref,],1,mean), apply(avg_metric[-id_ref,]/avg_metric[id_ref,],1,sd),
                           apply(avg_predC_L1[-id_ref,]/avg_predC_L1[id_ref,],1,mean), apply(avg_predC_L1[-id_ref,]/avg_predC_L1[id_ref,],1,sd),
                           apply(avg_predI_L1[-id_ref,]/avg_predI_L1[id_ref,],1,mean), apply(avg_predI_L1[-id_ref,]/avg_predI_L1[id_ref,],1,sd),
                           apply(avg_AUCC[-id_ref,]/avg_AUCC[id_ref,],1,mean), apply(avg_AUCC[-id_ref,]/avg_AUCC[id_ref,],1,sd),
                           apply(avg_AUCI[-id_ref,]/avg_AUCI[id_ref,],1,mean), apply(avg_AUCI[-id_ref,]/avg_AUCI[id_ref,],1,sd),
                           apply(avg_BrierC[-id_ref,]/avg_BrierC[id_ref,],1,mean), apply(avg_BrierC[-id_ref,]/avg_BrierC[id_ref,],1,sd),
                           apply(avg_BrierI[-id_ref,]/avg_BrierI[id_ref,],1,mean), apply(avg_BrierI[-id_ref,]/avg_BrierI[id_ref,],1,sd))
      colnames(add_df_ratio) = c("mean_L1_ratio","sd_L1_ratio",
                                 sprintf("mean_%s_ratio",metric), sprintf("sd_%s_ratio",metric),
                                 "mean_predC_L1_ratio","sd_predC_L1_ratio","mean_predI_L1_ratio","sd_predI_L1_ratio",
                                 "mean_AUCC_ratio","sd_AUCC_ratio","mean_AUCI_ratio","sd_AUCI_ratio","mean_BrierC_ratio", "sd_BrierC_ratio",
                                 "mean_BrierI_ratio", "sd_BrierI_ratio")

    } else if(family=="Gaussian"){
      df3 = lapply(list_res,function(x) unlist(lapply(x$pred$complete,mean)));     df3 = lapply(df3, function(x) x[methods])
      df4 = lapply(list_res,function(x) unlist(lapply(x$pred$imputed,mean)));     df4 = lapply(df4, function(x) x[methods])


      avg_predC_L1 = sapply(df3,function(x) x )
      avg_predI_L1 = sapply(df4,function(x) x )
      add_df = cbind(apply(avg_L1,1,mean), apply(avg_L1,1,sd),
                     apply(avg_metric,1,mean), apply(avg_metric,1,sd),
                     apply(avg_predC_L1,1,mean), apply(avg_predC_L1,1,sd),
                     apply(avg_predI_L1,1,mean), apply(avg_predI_L1,1,sd))             # avg_L1, avg_metric, avg_pred1_L1, avg_pred2_L1
      colnames(add_df) = c("mean_L1","sd_L1",
                           sprintf("mean_%s",metric), sprintf("sd_%s",metric),
                           "mean_predC_L1","sd_predC_L1","mean_predI_L1","sd_predI_L1")


      add_df_ratio = cbind(apply(avg_L1[-id_ref,]/avg_L1[id_ref,],1,mean), apply(avg_L1[-id_ref,]/avg_L1[id_ref,],1,sd),
                           apply(avg_metric[-id_ref,]/avg_metric[id_ref,],1,mean), apply(avg_metric[-id_ref,]/avg_metric[id_ref,],1,sd),
                           apply(avg_predC_L1[-id_ref,]/avg_predC_L1[id_ref,],1,mean), apply(avg_predC_L1[-id_ref,]/avg_predC_L1[id_ref,],1,sd),
                           apply(avg_predI_L1[-id_ref,]/avg_predI_L1[id_ref,],1,mean), apply(avg_predI_L1[-id_ref,]/avg_predI_L1[id_ref,],1,sd))
      colnames(add_df_ratio) = c("mean_L1_ratio","sd_L1_ratio",
                                 sprintf("mean_%s_ratio",metric), sprintf("sd_%s_ratio",metric),
                                 "mean_predC_L1_ratio","sd_predC_L1_ratio","mean_predI_L1_ratio","sd_predI_L1_ratio")


    }
    add_df = data.frame(add_df)
    add_df$mechanism=mechanisms[m]
    add_df$method = factor(rownames(add_df), levels=methods)
    rownames(add_df) = NULL
    add_df_ratio = data.frame(add_df_ratio)
    add_df_ratio$mechanism=mechanisms[m]
    add_df_ratio$method = factor(rownames(add_df_ratio), levels=methods[-id_ref])
    rownames(add_df_ratio) = NULL
    df = rbind(df, add_df)
    df_ratio = rbind(df_ratio, add_df_ratio)
  }

  df$mechanism = factor(df$mechanism, levels=c("MCAR","MAR","MNAR"))
  df_ratio$mechanism = factor(df_ratio$mechanism, levels=c("MCAR","MAR","MNAR"))

  df = df[df$method %in% methods,]
  df_ratio = df_ratio[df_ratio$method %in% methods,]
  # perf_metrics = c("L1","PB","predC_L1","predI_L1")
  if(family=="Multinomial"){
    perf_metrics = c("L1","PB","predC_L1","predI_L1","AUCC","AUCI", "BrierC","BrierI")
  }else if(family=="Gaussian"){
    perf_metrics = c("L1","PB","predC_L1","predI_L1")
  }

  for(l in 1:length(perf_metrics)){
    print(perf_metrics[l])
    p=ggplot(df,aes(x=mechanism, y=eval(parse(text=paste("mean_",perf_metrics[l],sep=""))), fill=method, color=method))+
      geom_bar(stat="identity",position=position_dodge(.9),alpha=0.4,color="black")+#ylim(c(0,3))+#ylim(c(0,0.5))+
      geom_errorbar(aes(ymin=eval(parse(text=paste("mean_",perf_metrics[l],sep="")))-eval(parse(text=paste("sd_",perf_metrics[l],sep=""))),
                        ymax=eval(parse(text=paste("mean_",perf_metrics[l],sep="")))+eval(parse(text=paste("sd_",perf_metrics[l],sep="")))),
                    width=.2, position=position_dodge(.9),color="black") +
      labs(title=bquote( list("n =" ~ .(N) ~ ", p =" ~ .(P) ~ ", d =" ~ .(D), ~ mu[phi]==.(phi0)) ), y = perf_metrics[l], x="Mechanism") +
      theme(legend.title = element_blank(), text=element_text(size = 20))#,axis.text.x = element_text(colour = c(rep("black",nlevels(bar_df$mechanism)-1),"red")))

    ggsave(filename=sprintf("%s%s.png",
                            dir_name, perf_metrics[l]),
           plot=p, width = 12, height=7, units="in")

    p=ggplot(df_ratio,aes(x=mechanism, y=eval(parse(text=paste("mean_",perf_metrics[l],"_ratio",sep=""))), fill=method, color=method))+
      geom_bar(stat="identity",position=position_dodge(.9),alpha=0.4,color="black")+#ylim(c(0,3))+#ylim(c(0,0.5))+
      geom_errorbar(aes(ymin=eval(parse(text=paste("mean_",perf_metrics[l],"_ratio",sep="")))-eval(parse(text=paste("sd_",perf_metrics[l],"_ratio",sep=""))),
                        ymax=eval(parse(text=paste("mean_",perf_metrics[l],"_ratio",sep="")))+eval(parse(text=paste("sd_",perf_metrics[l],"_ratio",sep="")))),
                    width=.2, position=position_dodge(.9),color="black") +
      labs(title=bquote( list("n =" ~ .(N) ~ ", p =" ~ .(P) ~ ", d =" ~ .(D), ~ mu[phi]==.(phi0)) ), y = sprintf("%s_ratio (ref: %s)",perf_metrics[l], ref_method), x="Mechanism") +
      theme(legend.title = element_blank(), text=element_text(size = 20))#,axis.text.x = element_text(colour = c(rep("black",nlevels(bar_df$mechanism)-1),"red")))

    ggsave(filename=sprintf("%s%s_ratio.png",
                            dir_name, perf_metrics[l]),
           plot=p, width = 12, height=7, units="in")
  }
}}}
}

beta=0.25; pi=0.3; miss_pct_features = 50
# Ns=c(1e4); Ps=c(100,50,25); Ds=c(2,8)        # 1e4, 50 done. N=1e5 and P=25 needs to be done
# Ns=c(1e5); Ps=c(100,50,25); Ds=c(2,8)        # 1e4, 50 done. N=1e5 and P=25 needs to be done
phi0 = 5
data.file.name = NULL; mask.file.name=NULL
case="x"

# data_type_y = "real"; Cy=NULL   # real, cat, cts
# data_type_y = "cat"; Cy=2   # real, cat, cts

sim_indexes = 1:5; mechanisms=c("MNAR","MAR","MCAR") #; Ignorables=c(F,T)
methods=c("mean","mice","idlglm","dlglm","miwae","notmiwae")

init_r = "default"   # or "alt"

# process linear miss model sims
NL_x=F; NL_y=F; NL_r=F
run_processScript(beta, pi, miss_pct_features,
                  Ns = c(1e4,1e5), Ps=c(50, 25), Ds=c(2, 8), phi0, case, data_type_y="cat", Cy=2,
                  sim_indexes, mechanisms,
                  NL_x, NL_y, NL_r, methods, init_r)

# process nonlinear miss model sims
NL_x=F; NL_y=F; NL_r=T
run_processScript(beta, pi, miss_pct_features,
                  Ns = c(1e4,1e5), Ps=c(50, 25), Ds=c(2, 8), phi0, case, data_type_y="cat", Cy=2,
                  sim_indexes, mechanisms,
                  NL_x, NL_y, NL_r, methods, init_r)
