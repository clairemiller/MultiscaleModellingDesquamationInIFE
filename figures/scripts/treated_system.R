source("scripts/paper_plot_theme.R")
results_dir <- paste0(results_parent_dir,"Treated/")


# Steady state height ------------------------------------------------------
# Process data
d <- applyFunctionToResultsDirectories(getHeights,maindir=results_dir)
d <- rbindlist(d)
d$n_mut <- as.numeric(gsub(".*NDiseased|pct","",d$setup))
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
output_pdf_plot(p_ss,"treated_ss.pdf")


# Reactants ---------------------------------------------------------------
d_r <- getReactants(paste0(results_dir,"NDiseased050/Seed00/results_from_time_360"))
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
pdf("treated_reactants.pdf",width=15,height=10)
print(p_r)
print(p_sub,vp=vp)
dev.off()


# Velocity ----------------------------------------------------------------
# Process data
d_v <- applyFunctionToResultsDirectories(getZVelocity,maindir=results_dir)
d_v <- rbindlist(d_v,idcol="setup")
d_v$n_healthy <- 100-as.numeric(gsub(".*NDiseased|Seed.*","",d_v$setup))

d_v <- filter(d_v,time >= (min(time)+30))
d_v$time <- d_v$time-min(d_v$time)

ts_v <- summarise(group_by(d_v,n_healthy,time),mean_v=mean(v), min_v = min(v), max_v=max(v))
mean_v <- summarise(group_by(d_v,n_healthy),mean_v=mean(v),min_v=min(v),max_v=max(v))
var_v =  (max(mean_v$mean_v)-min(mean_v$mean_v))/mean(mean_v$mean_v)

# Plot
p_v <- ggplot(mean_v,aes(x=n_healthy,y=mean_v)) +
  geom_point(size=10) + 
  geom_errorbar(aes(ymin=min_v,ymax=max_v),width=5,size=2) +
  labs(x="Prop. healthy SC [%]",y=TeX("z velocity \\[CD.hr^{-1}\\]"))
p_v <- add_paper_theme(p_v)
output_pdf_plot(p_v,"treated_velocity.pdf")


# Turnover time -----------------------------------------------------------
d_tt <- applyFunctionToResultsDirectories(getTurnoverTime,maindir=results_dir)
d_tt <- rbindlist(d_tt)
d_tt <- filter(d_tt,time >= (min(time) + 30))
d_tt$n_healthy <- 100-as.numeric(gsub(".*NDiseased|/Seed.*","",d_tt$setup))

d_tt_ave <- summarise(group_by(d_tt,n_healthy),mean_tt=mean(age),max_tt=max(age),min_tt=min(age),median_tt=median(age))

# Against height
d_tt_ss <- inner_join(d_tt_ave,d_ss,by="n_healthy")
fit <- lm(median_tt~ss, data=d_tt_ss)
p_tt_ss <- ggplot(d_tt_ss,aes(x=ss,y=median_tt, colour=n_healthy)) + 
  geom_point(size=10) +
  geom_smooth(method="lm",se=F, colour="black", linetype=2, size=2) +
  labs(y="Turnover time [days]", x=TeX("Corneum $\\tau_{ss}$ \\[CD\\]"), colour="Prop. healthy\nSC [%]") +
  xlim(7.7,11) + ylim(9,12.3)
p_tt_ss <- add_paper_theme(p_tt_ss,
                           panel.grid.major.x = element_line(colour = "lightgrey"),
                           legend.position=c(0.83,0.3))
output_pdf_plot(p_tt_ss,"treated_heightvturnover.pdf") 