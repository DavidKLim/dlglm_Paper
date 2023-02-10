
library(dlglm)
phi0=5; pi=0.3; sim_index=1
miss_pct_features = 50
mechanisms=c("MCAR","MAR","MNAR")
case="x"
miss_cols=NULL; ref_cols=NULL


datasets=c("RED","WHITE","BANKNOTE","CAREVALUATION","SPAM")
datasets=c("DRYBEAN","LETTER","SHUTTLE")
beta=NA

# family="Gaussian"; Cy=NULL
family="Multinomial"; Cy=NA       ## 3 classes for cat vars
N=NA; D=NA; P=NA
C=NULL       ## 3 classes for cat vars
data_types=NA
NL_x = F; NL_y = F; NL_r = F

# res = list()
# index = 1



for(d in 1:length(datasets)){
for(m in 1:length(mechanisms)){
  for(i in sim_index){
    ## by default: save datasets as file.
    # res[[index]] = dlglm::prepareData(dataset = datasets[d], sim.params = list(N=N, D=D, P=P, data_types=data_types, family=family, sim_index=i, ratios=c(train=.8,valid=.1,test=.1),
    #                                                                            beta=beta, C=C, Cy=Cy, NL_x=NL_x, NL_y=NL_y),
    #                                   miss.params=list(scheme="UV", mechanism=mechanisms[m], NL_r=NL_r, pi=pi, phi0=phi0, miss_pct_features=miss_pct_features,
    #                                                    miss_cols=miss_cols, ref_cols=ref_cols),
    #                                   case=case)
    dlglm::prepareData(dataset = datasets[d], sim.params = list(N=N, D=D, P=P, data_types=data_types, family=family, sim_index=i, ratios=c(train=.8,valid=.1,test=.1),
                                                                               beta=beta, C=C, Cy=Cy, NL_x=NL_x, NL_y=NL_y),
                                      miss.params=list(scheme="UV", mechanism=mechanisms[m], NL_r=NL_r, pi=pi, phi0=phi0, miss_pct_features=miss_pct_features,
                                                       miss_cols=miss_cols, ref_cols=ref_cols),
                                      case=case)
    # index = index+1
  }
}
}

# X as miss model

datasets="SIM"
beta=0.25

Ns=c(10000,100000); Ps=c(25, 50); Ds=c(2,8)
family="Multinomial"; Cy=2       ## 3 classes for cat vars
NL_x = F; NL_y = F; NL_r = F

for(d in 1:length(datasets)){
  for(ns in 1:length(Ns)){for(ps in 1:length(Ps)){ for(ds in 1:length(Ds)){
    P_real=Ps[ps]; P_count=0; P_cat=0; C=NULL       ## 3 classes for cat vars
    data_types = c( rep("real",P_real), rep("count",P_count), rep("cat",P_cat) )
      for(m in 1:length(mechanisms)){
        for(i in sim_index){
          dlglm::prepareData(dataset = datasets[d], sim.params = list(N=Ns[ns], D=Ds[ds], P=Ps[ps], data_types=data_types, family=family, sim_index=i, ratios=c(train=.8,valid=.1,test=.1),
                                                                      beta=beta, C=C, Cy=Cy, NL_x=NL_x, NL_y=NL_y),
                             miss.params=list(scheme="UV", mechanism=mechanisms[m], NL_r=NL_r, pi=pi, phi0=phi0, miss_pct_features=miss_pct_features,
                                              miss_cols=miss_cols, ref_cols=ref_cols),
                             case=case)
        }
      }
  }}}}

## log X miss model

datasets="SIM"
beta=0.25

Ns=c(10000,100000); Ps=c(25, 50); Ds=c(2,8)
family="Multinomial"; Cy=2       ## 3 classes for cat vars
NL_x = F; NL_y = F; NL_r = T

for(d in 1:length(datasets)){
  for(ns in 1:length(Ns)){for(ps in 1:length(Ps)){ for(ds in 1:length(Ds)){
    P_real=Ps[ps]; P_count=0; P_cat=0; C=NULL       ## 3 classes for cat vars
    data_types = c( rep("real",P_real), rep("count",P_count), rep("cat",P_cat) )
    for(m in 1:length(mechanisms)){
      for(i in sim_index){
        dlglm::prepareData(dataset = datasets[d], sim.params = list(N=Ns[ns], D=Ds[ds], P=Ps[ps], data_types=data_types, family=family, sim_index=i, ratios=c(train=.8,valid=.1,test=.1),
                                                                    beta=beta, C=C, Cy=Cy, NL_x=NL_x, NL_y=NL_y),
                           miss.params=list(scheme="UV", mechanism=mechanisms[m], NL_r=NL_r, pi=pi, phi0=phi0, miss_pct_features=miss_pct_features,
                                            miss_cols=miss_cols, ref_cols=ref_cols),
                           case=case)
      }
    }
  }}}}




phi0=NA; pi=NA; sim_index=1
miss_pct_features = NA
mechanisms=NA
case="x"
miss_cols=NULL; ref_cols=NULL


datasets="BANKMARKETING"
beta=NA

# family="Gaussian"; Cy=NULL
family="Multinomial"; Cy=NA       ## 3 classes for cat vars
N=NA; D=NA; P=NA
C=NULL       ## 3 classes for cat vars
data_types=NA
NL_x = F; NL_y = F; NL_r = F

dlglm::prepareData(dataset = datasets[d], sim.params = list(N=N, D=D, P=P, data_types=data_types, family=family, sim_index=i, ratios=c(train=.8,valid=.1,test=.1),
                                                            beta=beta, C=C, Cy=Cy, NL_x=NL_x, NL_y=NL_y),
                   miss.params=list(scheme="UV", mechanism=mechanisms[m], NL_r=NL_r, pi=pi, phi0=phi0, miss_pct_features=miss_pct_features,
                                    miss_cols=miss_cols, ref_cols=ref_cols),
                   case=case)
