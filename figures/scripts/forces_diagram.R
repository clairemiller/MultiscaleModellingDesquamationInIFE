source("scripts/paper_plot_theme.R")

# Adhesion
calculatePalssonAdhesionForce <- function(s = seq(0,1.5,0.01), alpha = 0.2)
{
  lambda = 7.0
  c1 = 1/sqrt(2*lambda)
  s = s/0.5;
  Fij = alpha*( (s+c1)*exp(-lambda*(s+c1)^2) - c1*exp(-lambda*c1*c1)*exp(-lambda*s*s) )
  return(Fij)
}
da <- data.frame(s = seq(0.0,1.0,0.01))
da$f = -1.0*calculatePalssonAdhesionForce(da$s,alpha=1/0.002669)
da$t = "Adhesion"

# Exponential repulsion
dr <- data.frame(s = seq(-0.5,0.0,0.01))
dr$f <- 15* log(1+dr$s)
dr$t <- "Repulsion"

# Combine
d <- rbind(da,dr)


p <- ggplot(d) + 
  geom_hline(yintercept=0) +
  geom_line(aes(x=s,y=-f,colour=t),size=2) + 
  #coord_cartesian(ylim=c(-0.001,0.005)) +
  labs(x=TeX("$r_{ij}$ \\[CD\\]"),y=TeX("$F_{ij}\\;\\[\\mu N\\]$"), colour="Force" )
print(p)
p <- add_paper_theme_colourplot(p,legend.position=c(0.8,0.75))


# Note: need to remove line from plot function to run this code
output_pdf_plot(p, "force.pdf")
