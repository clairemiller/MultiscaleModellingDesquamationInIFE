source("scripts/paper_plot_theme.R")
results_dir <- paste0(results_parent_dir,"HeterogeneousInhibitor/")


# Steady state height ------------------------------------------------------
# Process data
d <- applyFunctionToResultsDirectories(getHeights,maindir=results_dir)
d <- bind_rows(d)
d$n_mut <- as.numeric(gsub(".*NAbnormal|pct","",d$setup))
d$n_healthy <- 100-d$n_mut
d$time <- d$time-min(d$time)

d_ave <- summarise(group_by(filter(d,time >= 30),time,n_healthy),h=mean(top_mean)-4)
d_ave$time <- d_ave$time-min(d_ave$time)
d_ss <- summarise(group_by(d_ave,n_healthy),ss=mean(h), ss_min=min(h), ss_max = max(h))

# Calculate the fit (quadratic)
d_ss_fit <- lm(ss~poly(n_healthy,2),data=d_ss)
d_ss_predict <- data.frame(n_healthy=seq(0,100,1))
d_ss_predict$ss = predict(d_ss_fit,newdata=d_ss_predict)

# Steady state summary plot
p_ss <- ggplot(d_ss,aes(x=n_healthy,y=ss)) +
  geom_point(size=10) + geom_errorbar(aes(ymin=ss_min,ymax=ss_max),width=5,size=2) +
  geom_line(data=d_ss_predict,linetype=2,size=2) +
  labs(x="Prop. healthy SC [%]",y=TeX("Corneum $\\tau_{ss}$ \\[CD\\]"))
p_ss <- add_paper_theme(p_ss) 
output_pdf_plot(p_ss,"heterogeneous_ss.pdf")


# Reactants ---------------------------------------------------------------
d_r <- getReactants(paste0(results_dir,"NAbnormal050/Seed00/results_from_time_360"))
d_r <- filter(d_r,t==max(t))
d_r <- tidyr::pivot_longer(d_r,cols=c("s","e","i"),names_to="reactant")
d_r$reactant <- factor(d_r$reactant,levels=c("s","e","i")) # Need to ensure consistent order

# Main plot
p_r <- ggplot(d_r) +
  geom_point(aes(x=z,y=value,colour=reactant,alpha=1),size=6) +
  geom_vline(xintercept=d_ss$ss[3] + 4, size=2,linetype=2) +
  labs(x="z [CD]",y=NULL,colour=NULL) +
  scale_alpha(range=c(0.2,0.2)) + guides(alpha=F) +
  ylim(0,1)
# Sub plot of i
p_sub <- ggplot(filter(d_r,reactant=="i")) + 
  geom_point(aes(x=z,y=value,colour=reactant),show.legend=F,alpha=0.2) +
  labs(y=NULL,x=NULL) + scale_colour_manual(values=brewer_colours[3]) +
  scale_y_continuous(breaks = c(0,1e-5), labels=c("0",TeX("10^{-5}")))
# Add paper theme
p_r <- add_paper_theme_colourplot(p_r,panel.grid.major.x = element_line(colour = "lightgrey"))
p_sub <- add_paper_theme(p_sub,plot.background= element_rect(fill="white",colour="black"))
vp <- viewport(width=0.3,height=0.3,x=0.72,y=0.68)
# Manual output to add the extra subplot
pdf("heterogeneous_reactants.pdf",width=15,height=10)
print(p_r)
print(p_sub,vp=vp)
dev.off()


# Velocity ----------------------------------------------------------------
# Process data
d_v <- applyFunctionToResultsDirectories(getZVelocity,maindir=results_dir)
d_v <- bind_rows(d_v,.id="setup")
d_v$n_healthy <- 100-as.numeric(gsub(".*NAbnormal|Seed.*","",d_v$setup))

d_v <- filter(d_v,time >= (min(time)+30))
d_v$time <- d_v$time-min(d_v$time)

mean_v <- summarise(group_by(d_v,n_healthy),mean_v=mean(v),min_v=min(v),max_v=max(v))

# Plot
p_v <- ggplot() +
  geom_errorbar(aes(x=as.factor(n_healthy),ymin=min_v,ymax=max_v), 
                linetype=1, width=0.5, size=2, data=mean_v) +
  geom_boxplot(aes(x=as.factor(n_healthy),y=v),coef=10, size=2, data=d_v) +
  labs(x="Prop. healthy SC [%]",y=TeX("z velocity \\[CD.hr^{-1}\\]"))
p_v <- add_paper_theme(p_v)
output_pdf_plot(p_v,"heterogeneous_velocity.pdf")


# Turnover time -----------------------------------------------------------
d_tt <- applyFunctionToResultsDirectories(getTurnoverTime,maindir=results_dir)
d_tt <- bind_rows(d_tt)
d_tt <- filter(d_tt,time >= (min(time) + 30))
d_tt$n_healthy <- 100-as.numeric(gsub(".*NAbnormal|/Seed.*","",d_tt$setup))

d_tt_ave <- summarise(group_by(d_tt,n_healthy),mean_tt=mean(age),max_tt=max(age),min_tt=min(age),median_tt=median(age))

p_tt <- ggplot(d_tt,aes(x=as.factor(n_healthy))) +
  geom_errorbar(aes(ymin=min_tt,ymax=max_tt), data=d_tt_ave, linetype=1, width=0.5, size=2) +
  stat_boxplot(aes(y=age),geom='boxplot',coef=20, size=2) +
  labs(x="Prop. healthy SC [%]",y="Turnover time [days]")
p_tt <- add_paper_theme(p_tt)
output_pdf_plot(p_tt,"heterogeneous_turnovertime.pdf")
