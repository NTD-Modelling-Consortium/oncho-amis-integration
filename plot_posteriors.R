
library(truncnorm)
library(invgamma)
library(ggplot2)
library(gridExtra)
library(sf)
library(magrittr)
library(dplyr)
library(viridis)
library(vioplot)

setwd("~/Documents/oncho-endgame-apr25/")
load("Maps/ALL_prevalence_map_multipletimespoints.rds")
load("Maps/iu_task_lookup.rds")

plotTraj=T

map_all_mtp[[1]]$TaskID = as.integer(map_all_mtp[[1]]$TaskID)
map_all_mtp[[2]]$TaskID = as.integer(map_all_mtp[[2]]$TaskID)
map_all_mtp[[3]]$TaskID = as.integer(map_all_mtp[[3]]$TaskID)

#just loop over batches rerun
failed_ids = c() # get these from multipletimepoints_projections_inputs.R
table_iu_idx=read.csv(paste0("outputs/table_iu_idx.csv")) 
ids=setdiff(1:max(table_iu_idx$TaskID),failed_ids)

# Load parameters subsampled on cluster 
load("outputs/InputPars_MTP_allIUs.rds") # generated in multipletimepoints_projections_inputs.R
pars_all = sampled_params_all %>%
  select(-c(location))

# Calculate all ESS
ess_all=c() 
ius_with_sigma0.025 = c()
for (id in ids){
  
  load(paste0("outputs/output_",id,".Rdata")) # this is just to get ESS initial run with sigma=0.0025
  
  # colnames of weight_matrix should reflect IUID
  iu_list = colnames(output$weight_matrix)
  ess_iu = rep(NA,length(iu_list))
  
  # find ESS with ESS<200
  ess = output$ess
  iu_names_lt200 = iu_list[ess<200]
  iu_names_ge200 = iu_list[!ess<200]
  
  # keep track of IUs that should use sigma0.025 results (those that achieved ESS < 200 when sigma=0.0025)
  ius_with_sigma0.025 = c(ius_with_sigma0.025,iu_names_lt200)
  
  load(paste0("outputs/output_",id,".Rdata")) # loads amis_output
  ess_iu[!ess<200] = output$ess[!ess<200]
  
  if(file.exists(paste0("outputs/output_",id,"_sigma0.025.Rdata"))){
    load(paste0("outputs/output_",id,"_sigma0.025.Rdata")) # loads amis_output
    ess_iu[ess<200] = output$ess[ess<200]
  }
  
  ess_all = data.frame(rbind(ess_all,cbind(IUID=iu_list,ESS=ess_iu)))
}
ess_all = ess_all %>%
  mutate(ESS = as.numeric(ESS)) %>%
  arrange(ESS) 
 
# ess < 200 even with sigma0.025 
ess_lt200 = ess_all %>%
  filter(ESS<200)
write.csv(ess_lt200, file="outputs/plots/ess_lt200.csv",row.names = F)

dinvgammatrunc <- function(x, shape, rate, lower, upper, log=F) {
  res=ifelse(x < lower | x > upper, 0, dinvgamma(x, shape, rate) / diff(pinvgamma(c(lower,upper), shape, rate)))
  if (log==TRUE){
    res=log(res)
  }
  return(res)
}
rinvgammatrunc <- function(n,shape, rate, lower, upper){
  samples=rep(NA,n)
  for (i in 1:n){
    samples[i]=rinvgamma(1, shape, rate)
    while (samples[i] < lower | samples[i] > upper) {samples[i]=rinvgamma(1, shape, rate)}
  }
  return(samples)
}

qinvgammatrunc <- function(p, shape, rate, lower, upper) {
  res = rinvgammatrunc(100000,shape,rate,lower,upper)
  quant = quantile(res,probs=p)
  #res= diff(pinvgamma(c(lower,p), shape, rate)) / diff(pinvgamma(c(lower,upper), shape, rate))
  return(quant)
}

dat <- pars_all 

# Priors
k_alpha = 10
k_beta = 3
k_min = 0.01
k_max = 4
abr_mean = log(1000)
abr_sd = log(5000)-log(1000)
abr_min = log(0.01)
abr_max = log(200000)

priorkmed <- qinvgammatrunc(0.5, shape=k_alpha, rate=k_beta, lower=k_min, upper=k_max) 
priorklwr <- qinvgammatrunc(0.025, shape=k_alpha, rate=k_beta, lower=k_min, upper=k_max) 
priorkupr <- qinvgammatrunc(0.975, shape=k_alpha, rate=k_beta, lower=k_min, upper=k_max) 

priorABRmed <- exp(qnorm(0.5, abr_mean, abr_sd))
priorABRlwr <- exp(qnorm(0.025,  abr_mean, abr_sd))
priorABRupr <- exp(qnorm(0.975,  abr_mean, abr_sd))


mnABR <- aggregate(exp(bite_rate_per_person_per_year)~IUID, data=dat, FUN="median")
lwrABR <- aggregate(exp(bite_rate_per_person_per_year)~IUID, data=dat, FUN=quantile, probs=0.05)
uprABR <- aggregate(exp(bite_rate_per_person_per_year)~IUID, data=dat, FUN=quantile, probs=0.95)

ABR <- (cbind(mnABR, lwrABR[,2], uprABR[,2]))
colnames(ABR) <- c("IUID", "ABR", "ABRlwr", "ABRupr")

mnK <- aggregate(individual_exposure~IUID, data=dat, FUN="median")
lwrK <- aggregate(individual_exposure~IUID, data=dat, FUN=quantile, probs=0.05)
uprK <- aggregate(individual_exposure~IUID, data=dat, FUN=quantile, probs=0.95)

k <- (cbind(mnK, lwrK[,2], uprK[,2]))
colnames(k) <- c("IUID", "k", "klwr", "kupr")

mnprev_t1 <- aggregate(prev_t1~IUID, data=dat, FUN="median")
lwrprev_t1 <- aggregate(prev_t1~IUID, data=dat, FUN=quantile, probs=0.05)
uprprev_t1 <- aggregate(prev_t1~IUID, data=dat, FUN=quantile, probs=0.95)

prev_t1 <- (cbind(mnprev_t1, lwrprev_t1[,2], uprprev_t1[,2]))
colnames(prev_t1) <- c("IUID", "prev_t1", "prev_t1lwr", "prev_t1upr")

mnprev_t2 <- aggregate(prev_t2~IUID, data=dat, FUN="median")
lwrprev_t2 <- aggregate(prev_t2~IUID, data=dat, FUN=quantile, probs=0.05)
uprprev_t2 <- aggregate(prev_t2~IUID, data=dat, FUN=quantile, probs=0.95)

prev_t2 <- (cbind(mnprev_t2, lwrprev_t2[,2], uprprev_t2[,2]))
colnames(prev_t2) <- c("IUID", "prev_t2", "prev_t2lwr", "prev_t2upr")

mnprev_t3 <- aggregate(prev_t3~IUID, data=dat, FUN="median")
lwrprev_t3 <- aggregate(prev_t3~IUID, data=dat, FUN=quantile, probs=0.05)
uprprev_t3 <- aggregate(prev_t3~IUID, data=dat, FUN=quantile, probs=0.95)

prev_t3 <- (cbind(mnprev_t3, lwrprev_t3[,2], uprprev_t3[,2]))
colnames(prev_t3) <- c("IUID", "prev_t3", "prev_t3lwr", "prev_t3upr")

df <- cbind(ABR, k[,2:4], prev_t1[,2:4], prev_t2[,2:4], prev_t3[,2:4]) 
df <- df[order(df$ABR),]
df$IUN <- seq(1,nrow(df))

#trearment naive IUs
tn_ius = iu_task_lookup$IUID[which(iu_task_lookup$TaskID %in% c(675,720,721))]
df$treatment_naive = factor(df$IUID %in% tn_ius)


pABR <- ggplot(data=df, aes(y=IUN, x=ABR))+
  geom_segment(aes(x=ABRlwr, xend=ABRupr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point(aes(col=treatment_naive),alpha=0.5) +
  xlab("ABR") +
  theme(axis.text.y = element_blank()) +
  theme_bw() +
  scale_y_continuous(name="IU")+
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::number_format(accuracy = 1),
    limits = c(0.01,170000))

#labels = scales::trans_format("log10", scales::math_format(10^.x)))

pk <- ggplot(data=df, aes(y=IUN, x=k))+
  geom_segment(aes(x=klwr, xend=kupr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point(aes(col=treatment_naive),alpha=0.5) +
  xlab("Exposure heterogeneity") +
  theme(axis.text.y = element_blank()) +
  theme_bw() +
  scale_y_continuous(name="IU")+
  scale_x_continuous(trans="log10",
                     breaks = scales::trans_breaks("log10", function(x) 10^x),
                     labels = scales::number_format(),
                     limits = c(0.01,4))


pprev_t1 <- ggplot(data=df, aes(y=IUN, x=prev_t1))+
  geom_segment(aes(x=prev_t1lwr, xend=prev_t1upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point(aes(col=treatment_naive),alpha=0.5) +
  xlab("Prevalence at baseline (1975) ") +
  theme(axis.text.y = element_blank()) +
  theme_bw() +
  scale_y_continuous(name="IU") +
  xlim(0,1)
#scale_x_log10(
#  breaks = scales::trans_breaks("log10", function(x) 10^x),
#  labels = scales::number_format())

pprev_t2 <- ggplot(data=df, aes(y=IUN, x=prev_t2))+
  geom_segment(aes(x=prev_t2lwr, xend=prev_t2upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point(aes(col=treatment_naive),alpha=0.5) +
  xlab("Prevalence in 2000") +
  theme(axis.text.y = element_blank()) +
  theme_bw() +
  scale_y_continuous(name="IU") +
  xlim(0,1)

pprev_t3 <- ggplot(data=df, aes(y=IUN, x=prev_t3))+
  geom_segment(aes(x=prev_t3lwr, xend=prev_t3upr, y=IUN, yend=IUN), col="grey", linewidth=0.1)+
  geom_point(aes(col=treatment_naive),alpha=0.5) +
  xlab("Prevalence in 2018") +
  theme(axis.text.y = element_blank()) +
  theme_bw() +
  scale_y_continuous(name="IU") +
  xlim(0,1)


pdf("outputs/plots/posteriors_mtp.pdf" , width=14, height=10)
grid.arrange(
  pABR + 
    annotate(geom="segment", x=priorABRlwr, xend=priorABRupr, y=mean(df$IUN), yend=mean(df$IUN), col="red") +
    annotate(geom="point", x=priorABRmed, y=mean(df$IUN), col="red"),
  
  pk +   annotate(geom="segment", x=priorklwr, xend=priorkupr, y=mean(df$IUN), yend=mean(df$IUN), col="red") +
    annotate(geom="point", x=priorkmed, y=mean(df$IUN), col="red"),
  
  ggplot() + theme_void(),
  pprev_t1,
  pprev_t2,
  pprev_t3,
  
  
  nrow=2)
dev.off()

pdf("outputs/plots/posteriors_mtp_parameters.pdf" , width=10, height=6)
grid.arrange(
  pABR + 
    annotate(geom="segment", x=priorABRlwr, xend=priorABRupr, y=mean(df$IUN), yend=mean(df$IUN), col="red") +
    annotate(geom="point", x=priorABRmed, y=mean(df$IUN), col="red"),
  
  pk +   annotate(geom="segment", x=priorklwr, xend=priorkupr, y=mean(df$IUN), yend=mean(df$IUN), col="red") +
    annotate(geom="point", x=priorkmed, y=mean(df$IUN), col="red"),
  nrow=1)
dev.off()


pABRk <- ggplot(data=df, aes(y=k, x=ABR))+
  geom_point(aes(col=treatment_naive),alpha=0.5) +
  theme_bw() +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::number_format(accuracy = 1))

pkprev_t1 = ggplot(data=df, aes(y=k, x=prev_t1))+
  geom_point(aes(col=treatment_naive),alpha=0.5) +
  xlab("Prevalence at baseline") +
  theme_bw() 


pkprev_t2 = ggplot(data=df, aes(y=k, x=prev_t2))+
  geom_point(aes(col=treatment_naive),alpha=0.5) +
  xlab("Prevalence at 2000") +
  theme_bw() 


pkprev_t3 = ggplot(data=df, aes(y=k, x=prev_t3))+
  geom_point(aes(col=treatment_naive),alpha=0.5) +
  xlab("Prevalence at 2018") +
  theme_bw() 

pdf("outputs/plots/k_plots_mtp.pdf" , width=10, height=10)
grid.arrange(
  pABRk,
  pkprev_t1,
  pkprev_t2,
  pkprev_t3,
  nrow=2)
dev.off()


# Maps in Evandro's style
th <- theme(plot.title = element_text(size=12, hjust=0.5))
lwd_borders <- 0.1
colour <- "D"

# Map median posterior prevalence, abr and k
shape <- read_sf(dsn = "../ESPEN_IU_2021/", layer = "ESPEN_IU_2021")

#Prevalence
prev_t1_medians = mnprev_t1
colnames(prev_t1_medians) = c("IUID","prev_1975")
espen_prev_t1 = shape %>%
  dplyr::select(IU_CODE,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  mutate(IUID = paste0(ADMIN0ISO3,substr(IU_CODE,(nchar(IU_CODE)-4),nchar(IU_CODE)))) %>%
  left_join(prev_t1_medians, by="IUID")
map_prev_t1 = st_as_sf(espen_prev_t1)

prev_t1_map=ggplot(map_prev_t1) +
  geom_sf(aes(fill=prev_1975), lwd=lwd_borders)+
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80", limits=c(0,1)) +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Prev. 1975")) +
  ggtitle("Prevalence at baseline (1975)") + th

prev_t2_medians = mnprev_t2
colnames(prev_t2_medians) = c("IUID","prev_2000")
espen_prev_t2 = shape %>%
  dplyr::select(IU_CODE,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  mutate(IUID = paste0(ADMIN0ISO3,substr(IU_CODE,(nchar(IU_CODE)-4),nchar(IU_CODE)))) %>%
  left_join(prev_t2_medians, by="IUID")
map_prev_t2 = st_as_sf(espen_prev_t2)

prev_t2_map=ggplot(map_prev_t2) +
  geom_sf(aes(fill=prev_2000), lwd=lwd_borders)+
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80", limits=c(0,1)) +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Prev. 2000")) +
  ggtitle("Prevalence in 2000") + th

prev_t3_medians = mnprev_t3
colnames(prev_t3_medians) = c("IUID","prev_2018")
espen_prev_t3 = shape %>%
  dplyr::select(IU_CODE,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  mutate(IUID = paste0(ADMIN0ISO3,substr(IU_CODE,(nchar(IU_CODE)-4),nchar(IU_CODE)))) %>%
  left_join(prev_t3_medians, by="IUID")
map_prev_t3 = st_as_sf(espen_prev_t3)

prev_t3_map=ggplot(map_prev_t3) +
  geom_sf(aes(fill=prev_2018), lwd=lwd_borders)+
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80", limits=c(0,1)) +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="Prev. 2018")) +
  ggtitle("Prevalence in 2018") + th



#ABR
ABR_medians = mnABR
colnames(ABR_medians) = c("IUID","ABR")
espen_ABR = shape %>%
  dplyr::select(IU_CODE,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  mutate(IUID = paste0(ADMIN0ISO3,substr(IU_CODE,(nchar(IU_CODE)-4),nchar(IU_CODE)))) %>%
  left_join(ABR_medians, by="IUID")
map_ABR = st_as_sf(espen_ABR)

ABR_map=ggplot(map_ABR) +
  geom_sf(aes(fill=ABR), lwd=lwd_borders)+
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80", limits=c(0,3100)) +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="ABR")) +
  ggtitle("ABR") + th



#k
K_medians = mnK
colnames(K_medians) = c("IUID","Exposure")
espen_K = shape %>%
  dplyr::select(IU_CODE,ADMIN0ISO3,Shape_Leng,Shape_Area,geometry) %>%
  mutate(IUID = paste0(ADMIN0ISO3,substr(IU_CODE,(nchar(IU_CODE)-4),nchar(IU_CODE)))) %>%
  left_join(K_medians, by="IUID")
map_K = st_as_sf(espen_K)

K_map=ggplot(map_K) +
  geom_sf(aes(fill=Exposure), lwd=lwd_borders)+
  scale_fill_viridis(option=colour, direction = 1, na.value="gray80", limits=c(0,0.5)) +
  theme_bw() + theme(legend.position="right") +
  guides(fill=guide_colorbar(title="k")) +
  ggtitle("k") + th
  

pdf("outputs/plots/maps_all_mtp.pdf" , width=14, height=10)
grid.arrange(
  ABR_map,
  K_map,
  ggplot() + theme_void(),
  prev_t1_map,
  prev_t2_map,
  prev_t3_map,
  nrow=2)
dev.off()

if(plotTraj == TRUE){
  for(id in ids){
    
    mda_history= read.csv(paste0("oncho-venv/EPIONCHO-IBM/model_output/InputMDA_MTP_",id,".csv"))
    vc_history= read.csv(paste0("oncho-venv/EPIONCHO-IBM/model_output/InputVC_MTP_",id,".csv")) %>%
      filter(vector_control==1)
    
    # IUIDs 
    wh=which(map_all_mtp[[1]]$TaskID==id)
    iu_list=map_all_mtp[[1]]$IUID[wh]
    
    prevalence_map = lapply(1:length(map_all_mtp), function(t) {
      map_t = map_all_mtp[[t]]
      wh<-which(map_t$TaskID==id)
      map_t<-map_t[wh,]
      rownames(map_t)<-iu_list
      map_t = as.matrix(map_t[,(ncol(map_t)-1000+1):ncol(map_t),drop=FALSE])
      return(list(data=map_t))
    })
    
    # plot settings 
    L <- nrow(prevalence_map[[1]]$data)
    panel_ncols <- 4
    panel_nrows <- ceiling(L/panel_ncols)
    map_years = c(1975,2000,2018)
    
    # draws from prior
    if(id %in% c(1,675,676)){ # batches with no MDA are large
      pdf(paste0("outputs/plots/trajectories_prior_",id,".pdf") ,width=10,height=panel_nrows*2.5)
      
    } else {
      png(paste0("outputs/plots/trajectories_prior_",id,".png") ,width=10,height=panel_nrows*2.5,units="in",res=900)
    }
    par(mfrow=c(panel_nrows,panel_ncols))
    for (j in 1:length(iu_list)){
      iu = iu_list[j]
      
      if(iu %in% ius_with_sigma0.025){
        load(paste0("outputs/output_",id,"_sigma0.025.Rdata")) # loads amis_output
        load(paste0("outputs/trajectories_",id,"_sigma0.025.Rdata"))
      } else {
        load(paste0("outputs/output_",id,".Rdata")) # loads amis_output
        load(paste0("outputs/trajectories_",id,".Rdata"))
      }
      
      if(output$ess[j]>200){
        plot(y=trajectories[1,],x=c(1975:2018),type="l",ylim=c(0,1),
             main=paste0(iu, "\n ESS ",round(output$ess[j],digits=2)), xlab="year", ylab="prevalence")
        
        for(i in seq(1,500,10)){
          lines(y=trajectories[i,],x=1975:2018)
        }
        
        for (year_ind in 1:length(prevalence_map)){
          if (length(which(!is.na(prevalence_map[[year_ind]]$data[j,]))) > 0){
            vioplot(prevalence_map[[year_ind]]$data[j,],add=T,col="red",at=map_years[year_ind],wex=2)
          }
        }
        
        for (mda_ind in 1:nrow(mda_history)){
          year = as.numeric(mda_history$Year[mda_ind])
          abline(v=year,col="green")
        }
      }
    }
    dev.off()
    
    # draws from posterior"
    if(id %in% c(1,675,676)){
      pdf(paste0("outputs/plots/trajectories_posterior_",id,".pdf") ,width=10,height=panel_nrows*2.5)
      
    } else {
      png(paste0("outputs/plots/trajectories_posterior_",id,".png") ,width=10,height=panel_nrows*2.5,units="in",res=900)
    }
    
    par(mfrow=c(panel_nrows,panel_ncols))
    for (j in 1:length(iu_list)){
      iu = iu_list[j]
      
      if(iu %in% ius_with_sigma0.025){
        load(paste0("outputs/output_",id,"_sigma0.025.Rdata")) # loads amis_output
        load(paste0("outputs/trajectories_",id,"_sigma0.025.Rdata"))
      } else {
        load(paste0("outputs/output_",id,".Rdata")) # loads amis_output
        load(paste0("outputs/trajectories_",id,".Rdata"))
      }
      
      if(output$ess[j]>200){
        sampled_params_iu =  pars_all %>%
          filter(IUID == iu)
        plot(y=trajectories[sampled_params_iu$seed[1],],x=1975:2018,type="l",ylim=c(0,1),
             main=paste0(iu, "\n ESS ",round(output$ess[j],digits=2)), xlab="year", ylab="prevalence")
        
        for(i in 2:nrow(sampled_params_iu)){
          lines(y=trajectories[sampled_params_iu$seed[i],],x=1975:2018)
        }
        
        for (year_ind in 1:length(prevalence_map)){
          if (length(which(!is.na(prevalence_map[[year_ind]]$data[j,]))) > 0){
            vioplot(prevalence_map[[year_ind]]$data[j,],add=T,col="red",at=map_years[year_ind],wex=2)
          }
        }
        
        for (mda_ind in 1:nrow(mda_history)){
          year = as.numeric(mda_history$Year[mda_ind])
          abline(v=year,col="green")
        }
        
        for (vc_ind in 1:nrow(vc_history)){
          year = as.numeric(vc_history$Year[vc_ind])
          abline(v=year,col="darkgreen")
        }
      }
    }
    dev.off() 
    
    
  }
}

# just plot IUs that didn't reach ESS=200 even with sigma=0.025
if(plotTraj == TRUE){

  # plot settings 
  L <- nrow(ess_lt200)
  panel_ncols <- 4
  panel_nrows <- ceiling(L/panel_ncols)
  map_years = c(1975,2000,2018)
  
  pdf(paste0("outputs/plots/trajectories_posterior_lowESS.pdf") ,width=10,height=panel_nrows*2.5)
  
  par(mfrow=c(panel_nrows,panel_ncols))
  
  for(iu in ess_lt200$IUID){
    
    id = table_iu_idx$TaskID[which(table_iu_idx$IUID == iu)]
    
    load(paste0("outputs/output_",id,"_sigma0.025.Rdata")) # loads amis_output
    load(paste0("outputs/trajectories_",id,"_sigma0.025.Rdata"))
    mda_history= read.csv(paste0("oncho-venv/EPIONCHO-IBM/model_output/InputMDA_MTP_",id,".csv"))
    vc_history= read.csv(paste0("oncho-venv/EPIONCHO-IBM/model_output/InputVC_MTP_",id,".csv")) %>%
      filter(vector_control==1)
    
    # colnames of weight_matrix should reflect IUID
    iu_list = colnames(output$weight_matrix)
    
    prevalence_map = lapply(1:length(map_all_mtp), function(t) {
      map_t = map_all_mtp[[t]]
      wh<-which(map_t$TaskID==id)
      map_t<-map_t[wh,]
      rownames(map_t)<-iu_list
      map_t = as.matrix(map_t[,(ncol(map_t)-1000+1):ncol(map_t),drop=FALSE])
      return(list(data=map_t))
    })

    # draws from posterior"
    j=which(iu_list==iu)
    sampled_params_iu =  pars_all %>%
      filter(IUID == iu)
    plot(y=trajectories[sampled_params_iu$seed[1],],x=1975:2018,type="l",ylim=c(0,1),
         main=paste0(iu, "\n ESS ",round(output$ess[j],digits=2)), xlab="year", ylab="prevalence")
    
    for(i in 2:nrow(sampled_params_iu)){
      lines(y=trajectories[sampled_params_iu$seed[i],],x=1975:2018)
    }
    
    for (year_ind in 1:length(prevalence_map)){
      if (length(which(!is.na(prevalence_map[[year_ind]]$data[j,]))) > 0){
        vioplot(prevalence_map[[year_ind]]$data[j,],add=T,col="red",at=map_years[year_ind],wex=2)
      }
    }
    
    for (mda_ind in 1:nrow(mda_history)){
      year = as.numeric(mda_history$Year[mda_ind])
      abline(v=year,col="green")
    }
    
    for (vc_ind in 1:nrow(vc_history)){
      year = as.numeric(vc_history$Year[vc_ind])
      abline(v=year,col="darkgreen")
    }
  }
  dev.off() 
}

