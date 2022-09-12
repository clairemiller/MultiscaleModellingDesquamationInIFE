source("scripts/paper_plot_theme.R")
results_dir <- paste0(results_parent_dir,"Treated/NDiseased000")

# Process height data ------------------------------------------------------
d <- applyFunctionToResultsDirectories(getHeights, maindir=results_dir)
d <- rbindlist(d)
d$time <- d$time-min(d$time)
d <- filter(d,time >= 30)
d$time <- d$time-min(d$time)
d_ss <- mean(d$top_mean)

d_ave <- summarise(group_by(d,time),h=mean(top_mean),h_min=min(top_mean),h_max=max(top_mean))

# Height Plot --------------------------------------------------------------
p_healthy <- ggplot(d_ave) +
  geom_line(aes(x=time,y=h,colour="Mean tissue"),size=2) +
  geom_ribbon(aes(x=time,ymin=h_min,ymax=h_max),alpha=0.4,fill=brewer_colours[3]) +
  geom_line(aes(x=time,y=h-4,colour="Mean corneum"),size=2) +
  geom_ribbon(aes(x=time,ymin=h_min-4,ymax=h_max-4),alpha=0.4,fill=brewer_colours[2]) +
  geom_hline(aes(yintercept=d_ss,colour="Steady state est."),linetype=2,size=2) +
  geom_hline(aes(yintercept=d_ss-4,colour="Steady state est."),linetype=2,size=2) +
  labs(x="Time [days]",y="Thickness [CD]",colour=NULL) +
  scale_colour_manual(values=c(brewer_colours[c(2,3)],"black")) +
  ylim(0,ceiling(max(d_ave$h_max))) 
p_healthy <- add_paper_theme(p_healthy, legend.position=c(0.8,0.2)) 
output_pdf_plot(p_healthy,"healthy_height.pdf")


# Cell velocities ---------------------------------------------------------
v <- applyFunctionToResultsDirectories(getZVelocity,maindir=results_dir)
v <- rbindlist(v)

v <- filter(v,time >= (30 + min(time)))
v$time <- v$time - min(v$time)

mean_v <- summarise(group_by(v,time),mean_v = mean(v),max_v=max(v),min_v=min(v))

p_v <- ggplot(mean_v) +
  geom_line(aes(x=time,y=mean_v,colour="Daily mean"),size=2) +
  geom_ribbon(aes(x=time,ymin=min_v,ymax=max_v),fill="black",linetype=2,alpha=0.1) + 
  labs(x="Time [days]",y=TeX("z velocity $\\[CD.hr^{-1}\\]$"),colour=NULL) +
  scale_colour_manual(values="black")
p_v <- add_paper_theme(p_v,legend.position=c(0.85,0.13))
output_pdf_plot(p_v,"healthy_velocity.pdf")

v_corneum <- applyFunctionToResultsDirectories(getCorneumVelocity,maindir = results_dir)
v_corneum <- rbindlist(v_corneum)
v_corneum <- filter(v_corneum,time >= (30 + min(time)))
v_corneum <- summarise(group_by(v_corneum,mut),mean_v=mean(v))

# Check the velocity profiles from a couple of sims
# to justify constant velocity in SC assumption
getDailySCVelocity <- function(filepath) {
  v_ex <- readChasteResultsFile("nodevelocities.dat",filepath,
                        columns=c("id","x","y","z","vx","vy","vz")) %>%
    mutate(time = time/24, vx = vx*120, vy = vy*120, vz = vz*120) %>%    
    filter(time > (30 + min(time))) %>%
    filter(id > 100) %>%
    filter(z > 4, z < d_ss) %>%
    group_by(id) %>%
    summarise(n_dp = n(),
              max_v = max(diff(z)),
              min_v = min(diff(z)),
              mean_v = mean(diff(z))) %>%
    mutate(range_v = max_v-min_v,
           rel_range = (max_v-min_v)/mean_v,
           rel_var = abs((max_v+min_v)/(2*mean_v)-1)) %>%
    filter(n_dp > 3)
  return(v_ex)
}
v_ex <- applyFunctionToResultsDirectories(getDailySCVelocity,results_dir)
v_ex <- bind_rows(v_ex, .id="Seed")
  
summarise(v_ex, 
            mean_relvar = mean(rel_var),
            ci95_relvar = quantile(rel_var,p=0.95)
            )
quantile(v_ex$rel_var,p=0.95)

# Process reactant data ---------------------------------------------------
filepath <- paste0(results_dir,"/Seed00/results_from_time_360")
d_e <- readCellData(filepath,"KLKLevel")
d_i <- readCellData(filepath,"LEKTILevel")
d_s <- readCellData(filepath,"AdhesiveProteinLevel")
d_ages <- readChasteResultsFile(filename="cellages.dat",filepath=filepath,columns = c("id","x","y","z","age"))

# Determine total amount of free enzyme in the system
filter(d_e,z > 4 & z < d_ss) %>%
  filter(time > 30*24 + min(time)) %>%
  group_by(time) %>%
  summarise(total_e = sum(KLKLevel)) %>%
  summarise(mean_e = mean(total_e))
# mean_e = 804.

# Get top of tissue cells
d_mut <- readChasteResultsFile("results.vizmutationstates",filepath,columns = "mut")
d_top <- d_s[d_mut$mut == 6,]

# Rearrange for plotting
x <- data.frame(t=d_e$time,z=d_e$z,s=d_s$AdhesiveProteinLevel,e=d_e$KLKLevel,i=d_i$LEKTILevel,age=d_ages$age)
x <- filter(x,t == max(t))
x <- tidyr::pivot_longer(x,cols=c("s","e","i"))

# Read in the Matlab data
filename <- "/Users/millerc2/MATLAB/projects/project3/data/multicellularmodel_solutions.csv"
d_matlab <- read.csv(filename)
d <- mutate(d, z=mean(mean_v$mean_v)*t_hr + 4)
d <- tidyr::pivot_longer(d[,c("z","s","e")],cols = c("s","e"))
d <- rbind(data.frame(z=0,name=c("s","e"),value=c(1,0)),d)

# investigate the time step
d_matlab <- mutate(d_matlab, dt = c(NA,diff(t_hr)))
ggplot(d_matlab) + geom_point(aes(x=t_hr,y=dt)) + scale_y_log10()
quantile(d_matlab$dt[-1])

# Calculate the variation
d_s$z_bin <- round(d_s$z/0.1)*0.1
var_s <- summarise(group_by(filter(d_s,time==max(time) & z > 4),z_bin),s_max = max(AdhesiveProteinLevel),s_min=min(AdhesiveProteinLevel))
var_s$var <- var_s$s_max-var_s$s_min
mean(var_s$var)


# Reactants plot ----------------------------------------------------------
p_reactants <- ggplot(x) + 
  geom_point(aes(x=z,y=value,colour=name,alpha=1),size=6) +
  geom_line(aes(x=z,y=value,group=name),data=d,size=4,linetype=2) +
  geom_vline(xintercept=d_ss,linetype=3,size=2) + 
  labs(x="z [CD]",y=NULL,colour=NULL) +
  scale_alpha(range=c(0.2,0.2)) + guides(alpha=F)
# Smaller plot of i
p_sub <- ggplot(filter(x,name=="i")) + 
  geom_point(aes(x=z,y=value,colour=name),show.legend=F) +
  labs(y=NULL,x=NULL) + scale_colour_manual(values=brewer_colours[2]) +
  scale_y_continuous(breaks = c(0,1e-5), labels=c("0",TeX("10^{-5}")))
# Add paper theme
p_reactants <- add_paper_theme_colourplot(p_reactants,panel.grid.major.x = element_line(colour = "lightgrey"))
p_sub <- add_paper_theme(p_sub,plot.background= element_rect(fill="white",colour="black"))
vp <- viewport(width=0.35,height=0.35,x=0.65,y=0.65)
# Manual output to add the extra subplot
pdf("healthy_reactants.pdf",width=15,height=10)
  print(p_reactants)
  print(p_sub,vp=vp)
dev.off()



# Turnover time -----------------------------------------------
d_tt <- applyFunctionToResultsDirectories(getTurnoverTime,maindir=results_dir)
d_tt <- rbindlist(d_tt)
d_tt <- filter(d_tt,time >= (min(time) + 30))
median(d_tt$age)*24

p_tt <- ggplot(d_tt) + 
  geom_histogram(aes(x=age,y=stat(density)),fill="lightgrey",colour="black",binwidth=0.5) +
  labs(x="Cell age [days]",y="Proportion of cells")
p_tt <- add_paper_theme(p_tt)
output_pdf_plot(p_tt,"healthy_turnovertime.pdf")
