

#### univariable analysis
run_univar_analysis <- function(data, variables){
  
  precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))

  
  var_names <- variables
  models<-list()
  list<-list()
  
  for(i in 1:length(var_names)){
    
    formula <- logRt ~ 1 + 
        f(State, model="iid", hyper=precision.prior)+
        f(YEAR_NO, model="rw1", replicate=ISLAND_NO, hyper=precision.prior)+
        data[,var_names[i]] 
    
    
    models[[i]] <- inla(formula = formula,
                        data = data, family="gaussian",
                        control.family=list(link="identity"),
                        control.compute = list(waic = TRUE, dic=TRUE),
                        control.inla=list(int.strategy = "eb"), verbose=F)
    
    waic<- models[[i]]$waic$waic
    dic <- models[[i]]$dic$dic
    
    beta<-models[[i]]$summary.fixed[2,]  
    
    
    list[[i]]<-data.frame(var = var_names[i],
                          waic = waic,
                          dic = dic,
                          beta_mean = c(beta$mean),
                          beta_low = c(beta$`0.025quant`),
                          beta_upp = c(beta$`0.975quant`))
    
    
    models[[i]]<-NA
    print(paste0("finished model ", i, " out of ", length(variables)))
    
  }
  
    result<- bind_rows(list)
    return(result)
}

plot_univariable_beta <- function(univar_results_table, choice = "all"){
  
  if(choice!="all" & choice !="only_significant"){
    warning("choice must be all or only_significant")
  }
  
  res<-univar_results_table
  res$group <- c(rep(c("hols","rain", "max_hum","mean_hum", "min_hum",
                                      "max_temp", "mean_temp", "min_temp"), times=15)) 
  
  
  res$group2<- c(rep(0:14, each=8))
  
  res$Direction <-ifelse(res$beta_mean>0,"Positive", "Negative")
  
  if(choice=="only_significant"){
    res<-res[res$beta_low < 0 & res$beta_upp <0 | res$beta_low > 0 & res$beta_upp >0,]
  }  
  
    p<-res%>% 
      mutate(group=sub("hols", "Proportion of school holiday days", group),
             group=sub("mean_temp", "Mean temperature", group),
             group=sub("min_temp", "Minimum temperature", group),
             group=sub("max_temp", "Maximum temperature", group),
             group=sub("mean_hum", "Mean humidity", group),
             group=sub("min_hum", "Minimum humidity", group),
             group=sub("max_hum", "Maximum humidity", group),
             group=sub("rain", "Cumulative rainfall", group))%>%
      ggplot()+ facet_wrap("group",scales="free")+
      geom_pointrange(aes(x=group2, y=beta_mean, col=Direction,
                          ymin=beta_low, 
                          ymax=beta_upp), 
                      size=0.3,position=position_dodge(width=0.2)) + 
      geom_hline(yintercept=0,linetype="dashed")+
      labs(x = "Lag", y = "Beta value", col=NULL) +
      coord_flip() + #makes horizontal
      theme_classic() +
      theme(legend.position="bottom",
            legend.text=element_text(size=15),
            strip.background = element_blank())
    
  return(p)
  
}

best_univar <- function(univar_results_table){
  univar_results_table %>%filter(beta_low < 0 & beta_upp <0 | beta_low > 0 & beta_upp >0)->sig
  sig$var[sig$waic==min(sig$waic)] ->var
  waic<-min(sig$waic)
  dic<-min(sig$dic)
  return(c(var,waic,dic))
}

collinarity_check<-function(data,all_variables, variable_chosen){
  
  m.cor = abs(cor(data[,all_variables],use="complete.obs",method="pearson"))
  var<-variable_chosen
  sub1 <- m.cor[,colnames(m.cor)==var] #subset just the chosen variable
  test_df1 <- data.frame(names = all_variables,
                         values = sub1)
  remove<-test_df1[which(test_df1$values>0.6 | test_df1$values < -0.6),1] #names of the variables colinear not to be included in next round
  
  return(remove)
}


run_multivar_analysis <- function(data, variables_to_try, variables_chosen_already = NA, stage = NA){
  
  precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
  
  var_names <- variables_to_try
  
  models<-list()
  list<-list()
  
  for(i in 1:length(var_names)){
    
    if(length(variables_chosen_already)==1){
    formula <- logRt ~ 1  + data[,variables_chosen_already] +
      f(State, model="iid", hyper=precision.prior)+
      f(YEAR_NO, model="rw1", replicate=ISLAND_NO, hyper=precision.prior)+
      data[,var_names[i]]
    }
    if(length(variables_chosen_already)==2){
      formula <- logRt ~ 1  + data[,variables_chosen_already[1]] + data[,variables_chosen_already[2]] +
        f(State, model="iid", hyper=precision.prior)+
        f(YEAR_NO, model="rw1", replicate=ISLAND_NO, hyper=precision.prior)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==3){
      formula <- logRt ~ 1  + data[,variables_chosen_already[1]] +
                                        data[,variables_chosen_already[2]] +
                                        data[,variables_chosen_already[3]] +
        f(State, model="iid", hyper=precision.prior)+
        f(YEAR_NO, model="rw1", replicate=ISLAND_NO, hyper=precision.prior)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==4){
      formula <- logRt ~ 1  + data[,variables_chosen_already[1]] +
                                        data[,variables_chosen_already[2]] +
                                        data[,variables_chosen_already[3]] +
                                        data[,variables_chosen_already[4]] +
        f(State, model="iid", hyper=precision.prior)+
        f(YEAR_NO, model="rw1", replicate=ISLAND_NO, hyper=precision.prior)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==5){
      formula <- logRt ~ 1  + data[,variables_chosen_already[1]] +
                                        data[,variables_chosen_already[2]] +
                                        data[,variables_chosen_already[3]] +
                                        data[,variables_chosen_already[4]] +
                                        data[,variables_chosen_already[5]] +
        f(State, model="iid", hyper=precision.prior)+
        f(YEAR_NO, model="rw1", replicate=ISLAND_NO, hyper=precision.prior)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==6){
      formula <- logRt ~ 1  + data[,variables_chosen_already[1]] +
                                        data[,variables_chosen_already[2]] +
                                        data[,variables_chosen_already[3]] +
                                        data[,variables_chosen_already[4]] +
                                        data[,variables_chosen_already[5]] +
                                        data[,variables_chosen_already[6]] +
        f(State, model="iid", hyper=precision.prior)+
        f(YEAR_NO, model="rw1", replicate=ISLAND_NO, hyper=precision.prior)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==7){
      formula <- logRt ~ 1 + data[,variables_chosen_already[1]] +
                                        data[,variables_chosen_already[2]] +
                                        data[,variables_chosen_already[3]] +
                                        data[,variables_chosen_already[4]] +
                                        data[,variables_chosen_already[5]] +
                                        data[,variables_chosen_already[6]] +
                                        data[,variables_chosen_already[7]] +
        f(State, model="iid", hyper=precision.prior)+
        f(YEAR_NO, model="rw1", replicate=ISLAND_NO, hyper=precision.prior)+
        data[,var_names[i]]
    }
    if(length(variables_chosen_already)==8){
      formula <- logRt ~ 1 + data[,variables_chosen_already[1]] +
                                        data[,variables_chosen_already[2]] +
                                        data[,variables_chosen_already[3]] +
                                        data[,variables_chosen_already[4]] +
                                        data[,variables_chosen_already[5]] +
                                        data[,variables_chosen_already[6]] +
                                        data[,variables_chosen_already[7]] +
                                        data[,variables_chosen_already[8]] +
        f(State, model="iid", hyper=precision.prior)+
        f(YEAR_NO, model="rw1", replicate=ISLAND_NO, hyper=precision.prior)+
        data[,var_names[i]]
    }
    
    models[[i]] <- inla(formula = formula,
                        data = data, family="gaussian",
                        control.family=list(link="identity"),
                        control.compute = list(waic = TRUE, dic=TRUE),
                        control.inla=list(int.strategy = "eb"), verbose=F)
    
    waic<- models[[i]]$waic$waic
    dic<- models[[i]]$dic$dic
    beta<-models[[i]]$summary.fixed[2+stage,]
    
    list[[i]]<-data.frame(var = var_names[i],
                          waic = waic,
                          dic = dic,
                          beta_mean = c(beta$mean),
                          beta_low = c(beta$`0.025quant`),
                          beta_upp = c(beta$`0.975quant`))
    
    
    models[[i]]<-NA
    print(paste0("finished model ", i, " out of ", length(var_names), " in stage ", stage))
    
  }
  
  result <- bind_rows(list)
  
  # choose the best fitting
  result %>%filter(beta_low < 0 & beta_upp <0 | beta_low > 0 & beta_upp >0)->sig
  sig$var[sig$waic==min(sig$waic)] ->var_waic
  sig$var[sig$dic==min(sig$dic)] ->var_dic
  waic<-sig$waic[sig$var==var_waic]
  dic<-sig$dic[sig$var==var_dic]
  
  return(c(var_waic,waic,var_dic, dic))
}


#### run final multivariable model
run_final_model <- function(data){
  
  ## prior for RE in INLA models
  precision.prior <- list(prec = list(prior = "pc.prec", param = c(0.5, 0.01)))
  
  formula <- logRt ~ 1 + mean_temp + prop_hols + max_hum + min_hum14 + ev71 + cum_rf +
    f(State, model="iid", hyper=precision.prior)+
    f(YEAR_NO, model="rw1", replicate=ISLAND_NO, hyper=precision.prior)

  model <- inla(formula = formula,
                data = data, family="gaussian",
                control.family=list(link="identity"),
                control.compute = list(waic = TRUE, dic=TRUE, config=T),
                control.inla=list(int.strategy = "eb"), verbose=F)

return(model)
}

#### sampling
sampling <- function(model, data, nsamp=1000){
  
samples=inla.posterior.sample(nsamp,model)
f2 <- matrix(NA, nrow(data), nsamp)
xx.s <- inla.posterior.sample.eval(function(...) c(Predictor[1:nrow(data)]), samples)
for(n in 1:nsamp){
  f2[,n] = (xx.s[, n])
  }

data$est<- apply(f2, 1, median)
data$est_low<-apply(f2, 1, quantile, probs = c(0.025))
data$est_upp<- apply(f2, 1, quantile, probs = c(0.975))

return(data)
}


#### plotting final model fit
plot_final_fit <- function(data_sampled, choice){

if(choice!="Peninsular" & choice !="East"){
  warning("choice must be Peninsular or East")
}
  
if(choice=="East"){
p<-ggplot(data_sampled[data_sampled$Island!="West",])+
  facet_grid("State",scales="free")+
  geom_line(aes(x=date,y=exp(logRt),group=EPI_PERIOD_NO,col="Calculated"))+
  labs(x="Date",y="Rt",col=NULL)+
  geom_errorbar(data=data_sampled[data_sampled$Island!="West",],
                aes(x=date,ymin=exp(est_low),ymax=exp(est_upp),group=EPI_PERIOD_NO),col="pink")+
  geom_line(data=data_sampled[data_sampled$Island!="West",],
            aes(x=date,y=exp(est),group=EPI_PERIOD_NO,col="Estimated"))+
  geom_hline(aes(yintercept=1),linetype="dashed")+
  theme_bw()+
  theme(legend.position="bottom")+
  scale_color_manual(values = c("grey","red"))+
  guides(col=guide_legend(override.aes=list(linewidth=5)))
return(p)
}
if(choice=="Peninsular"){
p1<-ggplot(data_sampled[data_sampled$State=="JOHOR"|data_sampled$State=="MELAKA"|data_sampled$State=="NEGERI SEMBILAN",])+
  facet_grid("State",scales="free")+
  geom_line(aes(x=date,y=exp(logRt),group=EPI_PERIOD_NO,col="Calculated"))+
  labs(x="Date",y="Rt",col=NULL)+
  geom_errorbar(data=data_sampled[data_sampled$State=="JOHOR"|data_sampled$State=="MELAKA"|data_sampled$State=="NEGERI SEMBILAN",],
                aes(x=date,ymin=exp(est_low),ymax=exp(est_upp),group=EPI_PERIOD_NO),col="pink")+
  geom_line(data=data_sampled[data_sampled$State=="JOHOR"|data_sampled$State=="MELAKA"|data_sampled$State=="NEGERI SEMBILAN",],
            aes(x=date,y=exp(est),group=EPI_PERIOD_NO,col="Estimated"))+
  geom_hline(aes(yintercept=1),linetype="dashed")+
  theme_bw()+
  theme(legend.position="bottom")+
  scale_color_manual(values = c("grey","red"))+
  guides(col=guide_legend(override.aes=list(linewidth=5)))

p2<-ggplot(data_sampled[data_sampled$State=="PERAK"|data_sampled$State=="PULAU PINANG"|data_sampled$State=="SELANGOR",])+
  facet_grid("State",scales="free")+
  geom_line(aes(x=date,y=exp(logRt),group=EPI_PERIOD_NO,col="Calculated"))+
  labs(x="Date",y="Rt",col=NULL)+
  geom_errorbar(data=data_sampled[data_sampled$State=="PERAK"|data_sampled$State=="PULAU PINANG"|data_sampled$State=="SELANGOR",],
                aes(x=date,ymin=exp(est_low),ymax=exp(est_upp),group=EPI_PERIOD_NO),col="pink")+
  geom_line(data=data_sampled[data_sampled$State=="PERAK"|data_sampled$State=="PULAU PINANG"|data_sampled$State=="SELANGOR",],
            aes(x=date,y=exp(est),group=EPI_PERIOD_NO,col="Estimated"))+
  geom_hline(aes(yintercept=1),linetype="dashed")+
  theme_bw()+
  theme(legend.position="bottom")+
  scale_color_manual(values = c("grey","red"))+
  guides(col=guide_legend(override.aes=list(linewidth=5)))

p3<-ggplot(data_sampled[data_sampled$State=="KUALA LUMPUR"|data_sampled$State=="KELANTAN"|data_sampled$State=="KEDAH",])+
  facet_grid("State",scales="free")+
  geom_line(aes(x=date,y=exp(logRt),group=EPI_PERIOD_NO,col="Calculated"))+
  labs(x="Date",y="Rt",col=NULL)+
  geom_errorbar(data=data_sampled[data_sampled$State=="KUALA LUMPUR"|data_sampled$State=="KELANTAN"|data_sampled$State=="KEDAH",],
                aes(x=date,ymin=exp(est_low),ymax=exp(est_upp),group=EPI_PERIOD_NO),col="pink")+
  geom_line(data=data_sampled[data_sampled$State=="KUALA LUMPUR"|data_sampled$State=="KELANTAN"|data_sampled$State=="KEDAH",],
            aes(x=date,y=exp(est),group=EPI_PERIOD_NO,col="Estimated"))+
  geom_hline(aes(yintercept=1),linetype="dashed")+
  theme_bw()+
  theme(legend.position="bottom")+
  scale_color_manual(values = c("grey","red"))+
  guides(col=guide_legend(override.aes=list(linewidth=5)))

return(list(p1,p2,p3))
}
}


#### plotting final model random effects
plot_final_RE <- function(model_RE, choice, shapefile=NA){
  
  if(choice!="State_line"&choice!="State_map" & choice !="Year"){
    warning("choice must be State_line or State_map or Year")
  }
  
  if(choice=="State_line"){
  random<-model_RE$State
  colnames(random)<-c("ID","mean","sd","low","med","upp","mode","kld")
  A<-ggplot(random)+geom_point(aes(x=ID,y=mean))+geom_errorbar(aes(x=ID, ymin=low, ymax=upp))+theme_bw()+
  labs(y="random intercept on state",x="state")
  
  return(A)
  }
  if(choice=="State_map"){
    random<-model_RE$State
    colnames(random)<-c("ID","mean","sd","low","med","upp","mode","kld")
    
    shapefile$NAME_1 <- toupper(shapefile$NAME_1)
    shapefile <-sf::st_as_sf(shapefile)
    
    shapefile%>%left_join(random, by=c("NAME_1"="ID"))->test
    
    myCols <- adjustcolor(colorRampPalette(brewer.pal(n=9, 'Reds'))(20), .85)
    A<-ggplot(data=test)+
      geom_sf(aes(fill=mean),col="black") +
      scale_fill_gradient2(na.value = "grey", 
                           low = ("blue"),
                           #mid = "white",
                           high = ("red"), 
                           #midpoint=0
      )+
      theme_bw()+
      labs(fill="State-level intercept",x=NULL, y=NULL)+
      theme(legend.position = "bottom",
            legend.key.width=unit(2, "cm"),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            strip.background = element_blank())+
      guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5),
             size = guide_legend(title.position="top", title.hjust = 0.5))
    
    return(A)
  }
  if(choice=="Year"){
    random<-model_RE$YEAR_NO
    colnames(random)<-c("ID","mean","sd","low","med","upp","mode","kld")
    random$Island<- rep(c("East","Peninsular"), each=8)
    
    A<-ggplot(random)+geom_line(aes(x=ID,y=mean,col=Island))+
      geom_ribbon(aes(x=ID, ymin=low, ymax=upp,fill=Island),alpha=0.5)+theme_bw()+
      scale_x_continuous(breaks=c(1:8),labels=c(seq(from=min(data$Year),to=max(data$Year),by=1)))+
      labs(y="random intercept on year",x="year", fill=NULL, col=NULL)+theme(legend.position = c(0.1,0.8))
    
    return(A) 
    
  }
}


#### marginal r-squared
rsquared <- function(model, data){
  y_fixed <- model$summary.linear.predictor$mean
  
  # Calculate the observed values variance and predicted values variance
  var_obs <- var(data$logRt)
  var_fixed <- var(y_fixed)
  
  # Compute the marginal R-squared
  r_squared <- var_fixed / var_obs
  return(r_squared)
}

#### Examination of the trend
# Detrending the observed time-series by subtracting the predicted values. The detrended series should appear stationary (i.e. with no clear trend)
plot_detrended_ts <- function(choice, data){
  
  if(choice!="Peninsular" & choice !="East"){
    warning("choice must be Peninsular or East")
  }
  
  if(choice=="East"){
    p<-ggplot(data[data$Island!="West",])+
      facet_grid("State",scales="free")+
      geom_bar(aes(x=date, y = logRt - pred_mean),stat = 'identity')+
      labs(x="",y="Residuals",col=NULL)+
      geom_hline(aes(yintercept=0),linetype="dashed")+
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      theme_bw()+
      guides(col=guide_legend(override.aes=list(linewidth=5)))
    return(p)
  }
  if(choice=="Peninsular"){
    p1<-ggplot(data[data$State=="JOHOR"|data$State=="MELAKA"|data$State=="NEGERI SEMBILAN",])+
      facet_grid("State",scales="free")+
      geom_bar(aes(x=date, y = logRt - pred_mean),stat = 'identity')+
      labs(x="",y="Residuals",col=NULL)+
      geom_hline(aes(yintercept=0),linetype="dashed")+
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      theme_bw()+
      guides(col=guide_legend(override.aes=list(linewidth=5)))
    
    p2<-ggplot(data[data$State=="PERAK"|data$State=="PULAU PINANG"|data$State=="SELANGOR",])+
      facet_grid("State",scales="free")+
      geom_bar(aes(x=date, y = logRt - pred_mean),stat = 'identity')+
      labs(x="",y="Residuals",col=NULL)+
      geom_hline(aes(yintercept=0),linetype="dashed")+
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      theme_bw()+
      guides(col=guide_legend(override.aes=list(linewidth=5)))
    
    p3<-ggplot(data[data$State=="KUALA LUMPUR"|data$State=="KELANTAN"|data$State=="KEDAH",])+
      facet_grid("State",scales="free")+
      geom_bar(aes(x=date, y = logRt - pred_mean),stat = 'identity')+
      labs(x="",y="Residuals",col=NULL)+
      geom_hline(aes(yintercept=0),linetype="dashed")+
      scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
      theme_bw()+
      guides(col=guide_legend(override.aes=list(linewidth=5)))
    
    return(list(p1,p2,p3))
  }
}
