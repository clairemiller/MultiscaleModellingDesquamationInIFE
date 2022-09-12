source("scripts/paper_plot_theme.R")

# Read in data and process ------------------------------------------------
filename <- "~/matlab/data/movingcell_solutions.csv"
d <- read.csv(filename)
d$t_day <- d$t_hr/24
d <- filter(d,klk==5)

d <- filter(d,t_day <= 20)
d <- pivot_longer(d,cols = c("s","e","cs","i","ci"), names_to="reactant", names_transform = list(reactant = as.factor))
levels(d$reactant) <- gsub("c([a-z])","c[\\1]",levels(d$reactant))



# Enzyme solution calculated from literature ----------------------------------
d_exp <- filter(d,s0==1e-5 & n_iT == 1 & eT == 1e-6)
pA <- ggplot(d_exp) +
    geom_line(aes(x=t_day,y=value,colour=reactant), size=4) +
    scale_color_manual(values=brewer_colours,labels = parse(text=levels(d$reactant))) +
    labs(x="t [day]",y=NULL,colour=NULL)
pA <- add_paper_theme(pA,legend.position=c(0.8,0.5),panel.grid.major.x = element_line(colour = "lightgrey"))

output_pdf_plot(pA,"movingcell_actual.pdf", pwidth=12)
  
# Effective enzyme solution ----------------------------------
pE <- ggplot(filter(d,s0==1e-5 & n_iT == 1 & eT == 1e-10)) +
  geom_line(aes(x=t_day,y=value,colour=reactant), size=4) +
  scale_color_manual(values=brewer_colours,labels = parse(text=levels(d$reactant))) +
  labs(x="t [day]",y=NULL,colour=NULL)
pE <- add_paper_theme(pE,legend.position=c(0.9,0.6),panel.grid.major.x = element_line(colour = "lightgrey"))

output_pdf_plot(pE,"movingcell_effective.pdf",pwidth=12)
