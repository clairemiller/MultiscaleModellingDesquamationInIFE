source("scripts/paper_plot_theme.R")
results_dir <- paste0(results_parent_dir,"TwoProliferativePopulations/HarmonicMean/")
healthy_dir <- paste0(results_parent_dir,"Treated/NDiseased000")


# Steady state height plot ------------------------------------------------------
# Get data
d <- applyFunctionToResultsDirectories(getHeights, maindir=results_dir)
d <- rbindlist(d)
d <- mutate(d,deltaT = as.numeric(gsub(".*DeltaT|hr_.*","",setup)))
# Add the deltaT=0 data
d0 <- applyFunctionToResultsDirectories(getHeights,healthy_dir)
d0 <- rbindlist(d0)
d0 <- mutate(d0, deltaT=0)
d <- rbind(d,d0)

d <- mutate(filter(d,time >= 30+min(time)), time = time-min(time))
d_ss_perseed <- summarise(group_by(d,time,deltaT), h_t=mean(top_mean)-4)
d_ss <- summarise(group_by(d_ss_perseed,deltaT),ss=mean(h_t), ss_max = max(h_t), ss_min=min(h_t))

# Plot
p_ss <- ggplot(d_ss,aes(x=deltaT,y=ss)) +
  geom_point(size=10) + 
  geom_errorbar(aes(ymin=ss_min,ymax=ss_max),size=2,width=0.5) +
  labs(x=TeX("$\\Delta t$ \\[hr\\]"),y=TeX("Corneum $\\tau_{ss}$ \\[CD\\]")) +
  ylim(9.75,12)
p_ss <- add_paper_theme(p_ss) 
output_pdf_plot(p_ss,"twoprolifpops_ss.pdf",pwidth=12)


# Cell velocity -----------------------------------------------------------
# Process data
d_v <- applyFunctionToResultsDirectories(getZVelocity,maindir=results_dir)
d_v <- rbindlist(d_v,idcol="setup")
d_v$deltaT <- as.numeric(gsub(".*DeltaT|hr_.*","",d_v$setup))
# Add in the comparison
d0 <- applyFunctionToResultsDirectories(getZVelocity,maindir=healthy_dir)
d0 <- rbindlist(d0,idcol="setup")
d0$deltaT <- 0
d_v <- rbind(d_v,d0)

d_v <- filter(d_v,time >= (min(time)+30))
d_v$time <- d_v$time-min(d_v$time)

ts_v <- summarise(group_by(d_v,deltaT,time),mean_v=mean(v), min_v = min(v), max_v=max(v))
mean_v <- summarise(group_by(d_v,deltaT),mean_v=mean(v),min_v=min(v),max_v=max(v))

# Plot
p_v <- ggplot(mean_v,aes(x=deltaT,y=mean_v)) +
  geom_point(size=10) + 
  geom_errorbar(aes(ymin=min_v,ymax=max_v),width=0.3,size=2) +
  labs(x=TeX("$\\Delta t$ \\[hr\\]"),y=TeX("z velocity \\[CD.hr^{-1}\\]"))
p_v <- add_paper_theme(p_v)
output_pdf_plot(p_v,"twoprolifpop_velocity.pdf",pwidth=12)


