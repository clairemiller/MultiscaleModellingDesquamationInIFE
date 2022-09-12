library(latex2exp)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(grid)
library(tidyr)


# Plotting functions/data
#---------------------------------------------------

brewer_palette <- "Set1"
brewer_colours <- brewer.pal(9,brewer_palette)
brewer_colours <- brewer_colours[-6] # Removing yellow

add_paper_theme <- function(p_in,...) {
  p_out <- p_in + theme(
    panel.background = element_rect(colour = "black",fill="white"),
    panel.grid.major.y = element_line(colour = "lightgrey"),
    legend.key.size = unit(30,"pt"),
    legend.text=element_text(size=rel(3), margin = margin(t=10, r=0, b = 10, l=0), hjust=0),
    legend.title = element_text(face="bold", margin = margin(b=20)),
    legend.margin=margin(t=20,l=20,r=20,b=20),
    legend.background = element_rect(colour="darkgrey",fill="white"),
    legend.key = element_rect(colour=NA,fill=NA),
    plot.margin = unit(c(30,30,10,10),"pt"),
    title=element_text(size=rel(3)),
    plot.title = element_text(hjust=0.5),
    axis.text=element_text(size=rel(3)),
    strip.text=element_text(size=rel(4)),
    strip.text.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)),
    strip.background = element_rect(fill=NA),
    panel.spacing = unit(30,"pt"),
    axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0),size=rel(1.5)),
    axis.title.x = element_text(margin = margin(t = 20, r = 0, b = 0, l = 0),size=rel(1.5)),
    ...)
  return(p_out)
}

add_paper_theme_colourplot <- function(p_in,...)
{
  p_out <- add_paper_theme(p_in,...)
  p_out <- p_out + scale_color_manual(values=brewer_colours)
  return(p_out)
}

add_paper_theme_fillplot <- function(p_in,...)
{
  p_out <- add_paper_theme(p_in,...)
  p_out <- p_out + scale_fill_manual(values=brewer_colours)
  return(p_out)
}

add_paper_theme_gradientcolour <- function(p_in,...)
{
  p_out <- add_paper_theme(p_in,...)
  return(p_out)
}

output_pdf_plot <- function(p, filename, pwidth=15,pheight=10) {
  pdf(filename,width=pwidth,height=pheight)
  tryCatch(
    {
      print(p)
      dev.off()
    },
    error=function(e) {
      dev.off()
      print(e)
    }
  )
}



# Data read functions
#-----------------------------------------
results_parent_dir <- "~/Results/MultiscaleModellingDesquamationInIFE/"

getHeights <- function(filepath)
{
  print(filepath)
  heights <- readChasteResultsFile("tissueheight.dat",filepath=filepath,
                                   columns = c("wild_mut","wild_mean","wild_max","top_mut","top_mean","top_max"))
  heights$time <- heights$time/24
  
  heights$setup <- gsub(results_parent_dir,"",filepath)
  heights$seed <- gsub(".*Seed|/results_from_time.*","",heights$setup)
  heights$setup <- gsub("/Seed.*","",heights$setup)
  
  return(heights)
}

getReactants <- function(filepath) {
  print(filepath)
  d_e <- readCellData(filepath,"KLKLevel")
  d_i <- readCellData(filepath,"LEKTILevel")
  d_s <- readCellData(filepath,"AdhesiveProteinLevel")
  #d_ages <- readChasteResultsFile(filename="cellages.dat",filepath=filepath,columns = c("id","x","y","z","age"))
  
  x <- data.frame(t=d_e$time,z=d_e$z,s=d_s$AdhesiveProteinLevel,e=d_e$KLKLevel,i=d_i$LEKTILevel)#,age=d_ages$age)
  #x <- filter(x,t == max(t))
  #x <- melt(x,id.vars=c("t","z"))#,"age"))
  
  x$seed <- gsub(".*Seed|/results_from_time_.*","",filepath)
  x$id <- filepath
  
  return(x)
}

getTurnoverTime <- function(filepath) {
  print(filepath)
  x <- readChasteResultsFile("cellageatdeath.dat",filepath,columns = c("id","age"))
  x$seed <- gsub(".*Seed|/results_from.*","",filepath)
  x$setup <- filepath
  x$age <- x$age/24
  x$time <- x$time/24
  return(x)
}

getZVelocity <- function(filepath) {
  print(filepath)
  d <- readChasteResultsFile("nodevelocities.dat",filepath,columns=c("id","x","y","z","vx","vy","vz"))
  celltypes <- readChasteResultsFile("results.vizcelltypes",filepath)
  cellmut <- readChasteResultsFile("results.vizmutationstates",filepath)
  i_diff <- which(celltypes$v != 0 & cellmut$v != 6)
  d <- d[i_diff,]
  output <- as.data.frame(summarise(group_by(d,time),v=mean(vz)))
  output$seed <- gsub(".*Seed|/results_from_time_.*","",filepath)
  #output$id <- filepath
  output$time <- output$time/24
  # Adjust for the time step (this is multiplied in the output cpp so need to now divide)
  output$v <- output$v*120.0;
  return(output)
}

getCorneumVelocity <- function(filepath) {
  print(filepath)
  d <- readChasteResultsFile("nodevelocities.dat",filepath,columns=c("id","x","y","z","vx","vy","vz"))
  d$mut <- readChasteResultsFile("results.vizmutationstates",filepath)$v
  d <- filter(d,z >= 4)
  output <- as.data.frame(summarise(group_by(d,time,mut),v=mean(vz)))
  output$seed <- gsub(".*Seed|/results_from_time_.*","",filepath)
  #output$id <- filepath
  output$time <- output$time/24
  # Adjust for the time step (this is multiplied in the output cpp so need to now divide)
  output$v <- output$v*120.0;
  return(output) 
}

applyFunctionToResultsDirectories <- function(fun,maindir="./",pattern="",...)
{
  subdirs <- grep(pattern,list.dirs(path=maindir),value=T)
  subdirs <- as.list(grep("results_from_time_",subdirs,value=T))
  names(subdirs) <- gsub("\\.|/|results_from_time_[0-9]+","",subdirs)
  return( lapply(subdirs,fun,...) )
}
