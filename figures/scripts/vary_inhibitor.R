source("scripts/paper_plot_theme.R")

results_dir <- paste0(results_parent_dir,"DiseasedTissue/")

addITScale <- function(p,levels_prop,...) {
  p + scale_color_gradient(breaks=levels_prop, 
                         labels = TeX(paste0("$",format(levels_prop,digits=2),"e_T$")),...)
}

# Steady state height plot ------------------------------------------------------
# Process data
d <- applyFunctionToResultsDirectories(getHeights, maindir=results_dir)
d <- rbindlist(d)
d$iT_pct <- as.numeric(gsub(".*InhibitorLevel|pct","",d$setup))

d <- mutate(filter(d,time >= 30+min(time)), time = time-min(time), iT = iT_pct/100)
d_ss_perseed <- summarise(group_by(d,time,iT), h_t=mean(top_mean)-4)
d_ss <- summarise(group_by(d_ss_perseed,iT),ss=mean(h_t), ss_max = max(h_t), ss_min=min(h_t))
d_ss <- mutate(d_ss,pct_change = (d_ss$ss[5]-ss)/d_ss$ss[5])

# Calculate the fit
fit <- lm(h_t~iT,data=d_ss_perseed)
d_ss$predictlm <- predict(fit,newdata=d_ss)

# Plot
p_ss <- ggplot(d_ss,aes(x=iT,y=ss)) +
  geom_point(size=10) + 
  geom_errorbar(aes(ymin=ss_min,ymax=ss_max),width=0.1,size=2) +
  geom_line(aes(y=predictlm),linetype=2,size=2) +
  labs(x=TeX("$\\frac{i_T}{e_T}$ \\[nM\\]"),y=TeX("Corneum $\\tau_{ss}$ \\[CD\\]"))
p_ss <- add_paper_theme(p_ss) 
output_pdf_plot(p_ss,"varyinhibitor_ss.pdf",pwidth=12)


# Reactants ---------------------------------------------------------------
# Process data
d_r <- applyFunctionToResultsDirectories(getReactants,maindir=results_dir,pattern="Seed00")
d_r <- rbindlist(d_r)
d_r$iT <- as.numeric(gsub(".*InhibitorLevel|pct/Seed.*","",d_r$id))/100

# Determine total amount of free enzyme in the system
filter(d_r,z > 4 & z < d_ss$ss[d_ss$iT==0]) %>%
  filter(t > 30*24 + min(t)) %>%
  group_by(t) %>%
  summarise(total_e = sum(e)) %>%
  summarise(mean_e = mean(total_e))
# mean_e = 1527.

# Plot substrate (corneodesmosome)
d_r <- filter(d_r,t==max(t))
d_ss_hl <- filter(d_ss,iT %in% 0:1)
p_s <- ggplot(d_r) +
  geom_point(aes(x=z,y=s,colour=iT),size=4,alpha=0.2) +
  geom_vline(aes(xintercept=ss+4,colour=iT),data=d_ss_hl,linetype=2,size=2) +
  labs(x="z [CD]",y="s",colour=TeX("$i_T$ \\[nM\\]")) +
  ylim(0,1)
p_s <- addITScale(p_s,unique(d_r$iT),high="#FF1E20",low="#4D0609")
p_s <- add_paper_theme(p_s,legend.position=c(0.85,0.75)) + 
  theme(panel.grid.major.x = element_line(colour = "lightgrey"))
output_pdf_plot(p_s,"varyinhibitor_reactants_s.pdf")

# Plot enzyme (KLK)
p_e <- ggplot(d_r) +
  geom_point(aes(x=z,y=e,colour=iT),size=4,alpha=0.2) +
  geom_vline(aes(xintercept=ss+4,colour=iT),data=d_ss_hl,linetype=2,size=2) +
  labs(x="z [CD]",y="e",colour=TeX("$i_T$ \\[nM\\]")) +
  ylim(0,1)
p_e <- addITScale(p_e,rev(unique(d_r$iT)))
p_e <- add_paper_theme(p_e,legend.position=c(0.85,0.3)) +
  theme(panel.grid.major.x = element_line(colour = "lightgrey"))
output_pdf_plot(p_e,"varyinhibitor_reactants_e.pdf")



# Turnover time -----------------------------------------------------------
# Process data
d_tt <- applyFunctionToResultsDirectories(getTurnoverTime,maindir=results_dir)
d_tt <- rbindlist(d_tt)
d_tt$iT <- as.numeric(gsub(".*InhibitorLevel|pct/Seed.*","",d_tt$setup))/100
d_tt <- filter(d_tt,time >= (min(time) + 30))
d_tt_ave <- summarise(group_by(d_tt,iT),mean_tt=mean(age),max_tt=max(age),min_tt=min(age),median_tt=median(age))

#change_tt <- filter(d_tt_ave,iT %in% c(0,1))
#change_tt <- 1-change_tt$median_tt[1]/change_tt$median_tt[2]

# Box plot
p_tt_violin <- ggplot(d_tt) +
  stat_boxplot(aes(x=as.factor(iT),y=age),geom="boxplot",coef=100) +
  #geom_point(aes(x=as.factor(iT),y=median_tt),data=d_tt_ave,size=10) +
  labs(x=TeX("$\\frac{i_T}{e_T}$ \\[nM\\]"),y="Turnover time [days]") +
  ylim(5,23)
p_tt_violin <- add_paper_theme(p_tt_violin)
output_pdf_plot(p_tt_violin,"varyinhibitor_turnovertime.pdf")

# Against height
d_tt_ss <- inner_join(d_tt_ave,d_ss,by="iT")
fit <- lm(median_tt~ss,data=d_tt_ss)
p_tt_ss <- ggplot(d_tt_ss,aes(x=ss,y=median_tt,colour=iT)) + 
  geom_point(size=10) +
  geom_smooth(method="lm",se=F, colour="black", linetype=2, size=2) +
  labs(y="Turnover time [days]", x=TeX("Corneum $\\tau_{ss}$ \\[CD\\]"), colour=TeX("$i_T$ \\[nM\\]"))
p_tt_ss <- addITScale(p_tt_ss,unique(d_tt_ss$iT))
p_tt_ss <- add_paper_theme(p_tt_ss,panel.grid.major.x = element_line(colour = "lightgrey"),legend.position=c(0.85,0.3))
output_pdf_plot(p_tt_ss,"varyinhibitor_heightvturnover.pdf")  


# Cell velocity -----------------------------------------------------------
# Process data
d_v <- applyFunctionToResultsDirectories(getZVelocity,maindir=results_dir)
d_v <- rbindlist(d_v,idcol="setup")
d_v$iT <- as.numeric(gsub(".*InhibitorLevel|pctSeed.*","",d_v$setup))/100
d_v <- filter(d_v,time >= (min(time)+30))
d_v$time <- d_v$time-min(d_v$time)

ts_v <- summarise(group_by(d_v,iT,time),mean_v=mean(v), min_v = min(v), max_v=max(v))
mean_v <- summarise(group_by(d_v,iT),mean_v=mean(v),min_v=min(v),max_v=max(v))
var_v =  (max(mean_v$mean_v)-min(mean_v$mean_v))/mean(mean_v$mean_v)

# Plot
p_v <- ggplot(mean_v,aes(x=iT,y=mean_v)) +
  geom_point(size=10) + 
  geom_errorbar(aes(ymin=min_v,ymax=max_v),width=0.05,size=2) +
  labs(x=TeX("\\frac{i_T}{e_T} \\[nM\\]"),y=TeX("z velocity \\[CD.hr^{-1}\\]"))
p_v <- add_paper_theme(p_v) 
output_pdf_plot(p_v,"varyinhibitor_velocity.pdf")


