
library(dlglm)

run_processScript = function(dataset="SIM", pi=0.3, miss_pct_features = 50,
                             phi0=5, case="x", sim_indexes=1, mechanisms=c("MCAR","MAR","MNAR"),
                             methods=c("mean","mice","idlglm","dlglm"),
                             init_r="default"){

  dlglm_pref = if(init_r=="default"){""}else if(init_r=="alt"){"alt_init/"}

  dir_name = sprintf("Results_%s/miss_%s/phi%d%s",
                     dataset, case, phi0, dlglm_pref)
  ifelse(!dir.exists(dir_name), dir.create(dir_name), F)

  prefix=""
  df = data.frame()
  df_ratio = data.frame()


  for(m in 1:length(mechanisms)){
    data_fname = sprintf("%s/sim1/data_%s_%d_%d.RData",dir_name,mechanisms[m],miss_pct_features,pi*100)
    load(data_fname)

    if(dataset=="CAREVALUATION"){
      data_types_x=rep("cat", ncol(X))
    }else if(dataset=="ABALONE"){
      data_types_x=c("cat",rep("real", ncol(X)-1))
    }else{ data_types_x=rep("real", ncol(X)) }

    N=nrow(X); P=ncol(X); D=NULL
    P_real=sum(data_types_x=="real"); P_count=sum(data_types_x=="count"); P_cat=sum(data_types_x=="cat")

    if(dataset %in% c("SPAM","IRIS","WINE","RED","WHITE","YEAST","BANKNOTE", "BREAST",
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

    family = if(data_type_y=="real"){"Gaussian"}else if(data_type_y=="cat"){"Multinomial"}else if(data_type_y=="cts"){"Poisson"}

    fname = sprintf("%s/list_res_%s.out",
                    dir_name, mechanisms[m])
    if(file.exists(fname)){
      load(fname)
    } else{
      list_res = list()
      for(s in 1:length(sim_indexes)){
        sim_index=sim_indexes[s]; mechanism=mechanisms[m]
        sim.params = list(N=N, P=P, D=D, data_types=data_types_x,
                          family=family,
                          sim_index=sim_index, ratios=c(train=.8,valid=.1,test=.1),
                          C=C, Cy=Cy, NL_x=F, NL_y=F)
        miss.params = list(scheme="UV", mechanism=mechanism, pi=pi, phi0=phi0, miss_pct_features=miss_pct_features,NL_r=NL_r)

        list_res[[s]] = dlglm::processResults(dataset=dataset,prefix=prefix,data.file.name = data.file.name, mask.file.name=mask.file.name,
                                       sim.params = sim.params,
                                       miss.params = miss.params,
                                       case=case,
                                       data_types_x=data_types_x, data_type_y = data_type_y, methods=methods, init_r=init_r, normalize=T)
        ## list_res[[1]]$imp$L1s
        ## list_res[[1]]$inf$df
        ## list_res[[1]]$pred$complete
        ## list_res[[1]]$pred$imputed
        gc()
      }
      save(list_res, file=fname)
    }

    #### PLOT IMPUTATION barplots with error bars
    df1 = lapply(list_res,function(x) x$imp$L1s)
    df1 = lapply(df1, function(x) x[methods])
    avg_L1 = sapply(df1, function(x) sapply(x,mean))
    ## compute SD across sim_indexes?

    #### PLOT PREDICTION: barplots with error bars (complete and imputed)
    ref_method = "mean"
    id_ref = which(methods == ref_method)          # reference: dlglm. Show all metrics divided by metric of dlglm (relative)
    if(family=="Multinomial"){
      # df3 = lapply(list_res,function(x) unlist(lapply(x$pred$complete,mean)));     df3 = lapply(df3, function(x) x[methods])
      # df4 = lapply(list_res,function(x) unlist(lapply(x$pred$imputed,mean)));     df4 = lapply(df4, function(x) x[methods])
      df5 = lapply(list_res,function(x) unlist(x$pred$AUC_complete));     df5 = lapply(df5, function(x) x[methods])
      df6 = lapply(list_res,function(x) unlist(x$pred$AUC_imputed));     df6 = lapply(df6, function(x) x[methods])
      df7 = lapply(list_res,function(x) unlist(x$pred$Briers1));     df7 = lapply(df7, function(x) x[methods])
      df8 = lapply(list_res,function(x) unlist(x$pred$Briers2));     df8 = lapply(df8, function(x) x[methods])
      df9 = lapply(list_res,function(x) unlist(x$pred$ARI_complete));     df9 = lapply(df9, function(x) x[methods])
      df10 = lapply(list_res,function(x) unlist(x$pred$ARI_imputed));     df10 = lapply(df10, function(x) x[methods])

      df11 = lapply(list_res,function(x) unlist(x$pred$TPR_complete));     df11 = lapply(df11, function(x) x[methods])
      df12 = lapply(list_res,function(x) unlist(x$pred$TPR_imputed));     df12 = lapply(df12, function(x) x[methods])
      df13 = lapply(list_res,function(x) unlist(x$pred$FPR_complete));     df13 = lapply(df13, function(x) x[methods])
      df14 = lapply(list_res,function(x) unlist(x$pred$FPR_imputed));     df14 = lapply(df14, function(x) x[methods])

      df15 = lapply(list_res,function(x) unlist(x$pred$F1_complete));     df15 = lapply(df15, function(x) x[methods])
      df16 = lapply(list_res,function(x) unlist(x$pred$F1_imputed));     df16 = lapply(df16, function(x) x[methods])

      df17 = lapply(list_res,function(x) unlist(x$pred$kappa1));     df17 = lapply(df17, function(x) x[methods])
      df18 = lapply(list_res,function(x) unlist(x$pred$kappa2));     df18 = lapply(df18, function(x) x[methods])

      df19 = lapply(list_res,function(x) unlist(x$pred$mcc1));     df19 = lapply(df19, function(x) x[methods])
      df20 = lapply(list_res,function(x) unlist(x$pred$mcc2));     df20 = lapply(df20, function(x) x[methods])

      df21 = lapply(list_res,function(x) unlist(x$pred$PPV_complete));     df21 = lapply(df21, function(x) x[methods])
      df22 = lapply(list_res,function(x) unlist(x$pred$PPV_imputed));     df22 = lapply(df22, function(x) x[methods])

      # avg_predC_L1 = sapply(df3,function(x) x )
      # avg_predI_L1 = sapply(df4,function(x) x )
      avg_AUCC = sapply(df5,function(x) x )
      avg_AUCI = sapply(df6,function(x) x )
      avg_BrierC = sapply(df7,function(x) x )
      avg_BrierI = sapply(df8,function(x) x )
      avg_ARIC = sapply(df9, function(x) x)
      avg_ARII = sapply(df10, function(x) x)

      avg_TPRC = sapply(df11, function(x) x)
      avg_TPRI = sapply(df12, function(x) x)
      avg_FPRC = sapply(df13, function(x) x)
      avg_FPRI = sapply(df14, function(x) x)
      avg_F1C = sapply(df15, function(x) x)
      avg_F1I = sapply(df16, function(x) x)
      avg_kappaC = sapply(df17, function(x) x)
      avg_kappaI = sapply(df18, function(x) x)

      avg_MCCC = sapply(df19, function(x) x)
      avg_MCCI = sapply(df20, function(x) x)

      avg_PPVC = sapply(df21, function(x) x)
      avg_PPVI = sapply(df22, function(x) x)


      add_df = cbind(apply(avg_L1,1,mean), #apply(avg_L1,1,sd),
                     # apply(avg_predC_L1,1,mean), apply(avg_predC_L1,1,sd),
                     # apply(avg_predI_L1,1,mean), apply(avg_predI_L1,1,sd),
                     apply(avg_AUCC,1,mean), #apply(avg_AUCC,1,sd),
                     apply(avg_AUCI,1,mean), #apply(avg_AUCI,1,sd),
                     apply(avg_BrierC, 1, mean), #apply(avg_BrierC, 1, sd),
                     apply(avg_BrierI, 1, mean), #apply(avg_BrierI, 1, sd),
                     apply(avg_ARIC,1,mean), #apply(avg_ARIC,1,sd),
                     apply(avg_ARII,1,mean),#, apply(avg_ARII,1,sd)
                     apply(avg_TPRC,1,mean), apply(avg_TPRI,1,mean), apply(avg_FPRC,1,mean), apply(avg_FPRI,1,mean),
                     apply(avg_F1C,1,mean), apply(avg_F1I,1,mean),apply(avg_kappaC,1,mean), apply(avg_kappaI,1,mean),
                     apply(avg_MCCC,1,mean), apply(avg_MCCI,1,mean), apply(avg_PPVC,1,mean), apply(avg_PPVI,1,mean)
      )             # avg_L1, avg_metric, avg_pred1_L1, avg_pred2_L1
      colnames(add_df) = c("mean_L1",#"sd_L1",
                           # "mean_predC_L1","sd_predC_L1",
                           # "mean_predI_L1","sd_predI_L1",
                           "mean_AUCC",#"sd_AUCC",
                           "mean_AUCI",#"sd_AUCI",
                           "mean_BrierC",#"sd_BrierC",
                           "mean_BrierI",#"sd_BrierI"
                           "mean_ARIC",#"sd_ARIC",
                           "mean_ARII",#"sd_ARII"
                           "mean_TPRC","mean_TPRI","mean_FPRC","mean_FPRI","mean_F1C","mean_F1I","mean_kappaC","mean_kappaI",
                           "mean_MCCC","mean_MCCI","mean_PPVC","mean_PPVI"
      )


      add_df_ratio = cbind(avg_L1[-id_ref,]/avg_L1[id_ref,],
                           # apply(avg_predC_L1[-id_ref,]/avg_predC_L1[id_ref,],1,mean), apply(avg_predC_L1[-id_ref,]/avg_predC_L1[id_ref,],1,sd),
                           # apply(avg_predI_L1[-id_ref,]/avg_predI_L1[id_ref,],1,mean), apply(avg_predI_L1[-id_ref,]/avg_predI_L1[id_ref,],1,sd),
                           avg_AUCC[-id_ref,]/avg_AUCC[id_ref,],
                           avg_AUCI[-id_ref,]/avg_AUCI[id_ref,],
                           avg_BrierC[-id_ref,]/avg_BrierC[id_ref,],
                           avg_BrierI[-id_ref,]/avg_BrierI[id_ref,],
                           avg_ARIC[-id_ref,]/avg_ARIC[id_ref,],
                           avg_ARII[-id_ref,]/avg_ARII[id_ref,],
                           avg_TPRC[-id_ref,]/avg_TPRC[id_ref,], avg_TPRI[-id_ref,]/avg_TPRI[id_ref,],
                           avg_FPRC[-id_ref,]/avg_FPRC[id_ref,], avg_FPRI[-id_ref,]/avg_FPRI[id_ref,],
                           avg_F1C[-id_ref,]/avg_F1C[id_ref,], avg_F1I[-id_ref,]/avg_F1I[id_ref,],
                           avg_kappaC[-id_ref,]/avg_kappaC[id_ref,], avg_kappaI[-id_ref,]/avg_kappaI[id_ref,],
                           avg_MCCC[-id_ref,]/avg_MCCC[id_ref,], avg_MCCI[-id_ref,]/avg_MCCI[id_ref,],
                           avg_PPVC[-id_ref,]/avg_PPVC[id_ref,], avg_PPVI[-id_ref,]/avg_PPVI[id_ref,])
      colnames(add_df_ratio) = c("mean_L1_ratio",
                                 # "mean_predC_L1_ratio","mean_predI_L1_ratio",
                                 "mean_AUCC_ratio","mean_AUCI_ratio","mean_BrierC_ratio",
                                 "mean_BrierI_ratio", "mean_ARIC_ratio","mean_ARII_ratio",
                                 "mean_TPRC_ratio","mean_TPRI_ratio","mean_FPRC_ratio","mean_FPRI_ratio",
                                 "mean_F1C_ratio","mean_F1I_ratio","mean_kappaC_ratio","mean_kappaI_ratio",
                                 "mean_MCCC_ratio","mean_MCCI_ratio","mean_PPVC_ratio","mean_PPVI_ratio")

    } else if(family=="Gaussian"){
      df3 = lapply(list_res,function(x) unlist(lapply(x$pred$complete,mean)));     df3 = lapply(df3, function(x) x[methods])
      df4 = lapply(list_res,function(x) unlist(lapply(x$pred$imputed,mean)));     df4 = lapply(df4, function(x) x[methods])


      avg_predC_L1 = sapply(df3,function(x) x )
      avg_predI_L1 = sapply(df4,function(x) x )
      add_df = cbind(apply(avg_L1,1,mean),# apply(avg_L1,1,sd),
                     apply(avg_predC_L1,1,mean),# apply(avg_predC_L1,1,sd),
                     apply(avg_predI_L1,1,mean)#, apply(avg_predI_L1,1,sd)
      )             # avg_L1, avg_metric, avg_pred1_L1, avg_pred2_L1
      colnames(add_df) = c("mean_L1",#"sd_L1",
                           "mean_predC_L1",#"sd_predC_L1",
                           "mean_predI_L1"#,"sd_predI_L1"
      )


      add_df_ratio = cbind(avg_L1[-id_ref,]/avg_L1[id_ref,],
                           avg_predC_L1[-id_ref,]/avg_predC_L1[id_ref,],
                           avg_predI_L1[-id_ref,]/avg_predI_L1[id_ref,])
      colnames(add_df_ratio) = c("mean_L1_ratio","mean_predC_L1_ratio","mean_predI_L1_ratio")


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
    # perf_metrics = c("L1","predC_L1","predI_L1","AUCC","AUCI", "BrierC","BrierI")
    perf_metrics = c("L1","AUCC","AUCI", "BrierC","BrierI", "ARIC","ARII", "TPRC","TPRI",
                     "FPRC","FPRI","F1C","F1I","kappaC","kappaI","MCCC","MCCI","PPVC","PPVI")
    # if(Cy==2){ perf_metrics = c("L1","AUCC","AUCI", "BrierC","BrierI", "ARIC","ARII")
    # } else{ perf_metrics = c("L1","BrierC","BrierI", "ARIC","ARII") }
  }else if(family=="Gaussian"){
    perf_metrics = c("L1","predC_L1","predI_L1")
  }

  for(l in 1:length(perf_metrics)){
    # library(scales)
    ## AUC: (0,1), Brier: (0, 1), F1: (0, 1), ARI: (-0.1, 1), TPR: (0, 1), FPR: (0, 1), kappa: (-inf, 1)
    print(perf_metrics[l])

    p=ggplot(df,aes(x=mechanism, y=eval(parse(text=paste("mean_",perf_metrics[l],sep=""))), fill=method, color=method))+
      geom_bar(stat="identity",position="dodge",alpha=0.4,color="black")+
      # geom_errorbar(aes(ymin=eval(parse(text=paste("mean_",perf_metrics[l],sep="")))-eval(parse(text=paste("sd_",perf_metrics[l],sep=""))),
      #                   ymax=eval(parse(text=paste("mean_",perf_metrics[l],sep="")))+eval(parse(text=paste("sd_",perf_metrics[l],sep="")))),
      #               width=.2, position=position_dodge(.9),color="black") +
      labs(title=bquote( list("n =" ~ .(N) ~ ", p =" ~ .(P) ~ "," ~ mu[phi]==.(phi0)) ), y = perf_metrics[l], x="Mechanism") +
      theme(legend.title = element_blank(), text=element_text(size = 20))#,axis.text.x = element_text(colour = c(rep("black",nlevels(bar_df$mechanism)-1),"red")))
    if(perf_metrics[l] %in% c("AUCC","AUCI", "ARIC","ARII", "TPRC","TPRI","FPRC","FPRI","F1C","F1I","MCCC","MCCI","PPVC","PPVI")){
      p = p + scale_y_continuous(limits= c(0, 1), oob=scales::rescale_none)
    }
    ggsave(filename=sprintf("%s/%s.png",
                            dir_name, perf_metrics[l]),
           plot=p, width = 12, height=7, units="in")

    p=ggplot(df_ratio,aes(x=mechanism, y=eval(parse(text=paste("mean_",perf_metrics[l],"_ratio",sep=""))), fill=method, color=method))+
      geom_bar(stat="identity",position=position_dodge(.9),alpha=0.4,color="black")+#ylim(c(0,3))+#ylim(c(0,0.5))+
      # geom_errorbar(aes(ymin=eval(parse(text=paste("mean_",perf_metrics[l],"_ratio",sep="")))-eval(parse(text=paste("sd_",perf_metrics[l],"_ratio",sep=""))),
      #                   ymax=eval(parse(text=paste("mean_",perf_metrics[l],"_ratio",sep="")))+eval(parse(text=paste("sd_",perf_metrics[l],"_ratio",sep="")))),
      #               width=.2, position=position_dodge(.9),color="black") +
      labs(title=bquote( list("n =" ~ .(N) ~ ", p =" ~ .(P) ~ "," ~ mu[phi]==.(phi0)) ), y = sprintf("%s_ratio (ref: %s)",perf_metrics[l], ref_method), x="Mechanism") +
      theme(legend.title = element_blank(), text=element_text(size = 20))#,axis.text.x = element_text(colour = c(rep("black",nlevels(bar_df$mechanism)-1),"red")))

    ggsave(filename=sprintf("%s/%s_ratio.png",
                            dir_name, perf_metrics[l]),
           plot=p, width = 12, height=7, units="in")
  }
  save(df,file=sprintf("%s/df.out", dir_name, perf_metrics[l]))
  save(df_ratio,file= sprintf("%s/df_ratio.out", dir_name, perf_metrics[l]))
}

pi=0.3; miss_pct_features = 50
phi0 = 5
data.file.name = NULL; mask.file.name=NULL
case="x"

sim_indexes = 1; mechanisms=c("MNAR","MAR","MCAR")
NL_x=F; NL_y=F; NL_r=F
methods=c("mean","mice","idlglm","dlglm")

init_r = "default"   # or "alt"
datasets=c("DRYBEAN","LETTER","SHUTTLE","RED","WHITE","BANKNOTE","SPAM")
for(d in 1:length(datasets)){
  run_processScript(datasets[d], pi, miss_pct_features,
                    phi0, case,
                    sim_indexes, mechanisms,
                    methods, init_r)
}
