# This R script is to make graphs of the number of TE insertion sites per 10kb (window size 9) in each chromosome.
# The input file is "_somatic_ins_rmevade_10kb_sliding_size9.csv" from TED-seq analysis.

rm(list=ls())
library(ggplot2)
library(dplyr)
library(readr)
setwd("") # set working directory
files  <- list.files(pattern = "_somatic_ins_rmevade_10kb_sliding_size9.csv")  
for (file.name in files) 
{
  X1 <- read.table(file.name, header = FALSE, sep = "\t", col.names = c("c1","c2","c3","c4","c5","c6","c7")) 
  View(X1)
  ff.base <- sub('.csv', "", file.name)    
  
  #to set the scale of y-axis 
  if(max(X1$c4)<=10)
  {
    main_scale = 10
    minor_scale = 5
    limity = 10
  }
  else if ((max(X1$c4)>10 & max(X1$c4)<=100))
  {
    main_scale = trunc(max(X1$c4)/10)*10
    minor_scale = trunc(max(X1$c4)/10)*10/2
    limity = max(X1$c4)
  }
  else if ((max(X1$c4)>100 & max(X1$c4)<=1000))
  {
    main_scale = trunc(max(X1$c4)/100)*100
    minor_scale = trunc(max(X1$c4)/100)*100/2
    limity = max(X1$c4)
  }
  else if ((max(X1$c4)>1000 & max(X1$c4)<=10000))
  {
    main_scale = trunc(max(X1$c4)/1000)*1000
    minor_scale = trunc(max(X1$c4)/1000)*1000/2
    limity = max(X1$c4)
  }
  else
  {
    main_scale = "out of scale"
    minor_scale= "out of scale"
    limity = "out of scale"
  }
  
  print(file.name)
  print(main_scale)
  print(minor_scale)
  print(limity)
  
################Chr1###############################
df = filter(.data = X1, c1 %in% c("Chr1")) 
ggplot(df, aes(x=c2, y=c4)) +  
  geom_bar(stat = "identity", fill = 'blue') +
  xlab("") +
  ylab("") +
  scale_x_continuous(breaks = seq(0,30000000,5000000), 
                     labels = NULL,
                     limits = c(-10000, 32990000),　
expand= c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(0, main_scale, minor_scale),
                     labels = NULL,
                     limits = c(0, limity),
                     expand= c(0.03,0.03)) +
  theme(text = element_text(size = 0),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(size = 0.2),
        axis.line = element_line(size = 0.2),
        panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin= unit(c(2, 1, 0, 0), "pt"))
ggsave(paste0(ff.base, "_Chr1.pdf"), device = "pdf", width=2, height=0.3) 

################Chr2###############################
df = filter(.data = X1, c1 %in% c("Chr2")) 
ggplot(df, aes(x=c2, y=c4)) + 
  geom_bar(stat = "identity", fill = 'blue') +
  xlab("") +
  ylab("") +
  scale_x_continuous(breaks = seq(0,30000000,5000000), 
                     labels = NULL,
                     limits = c(-10000, 32990000), 　
                     expand= c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(0, main_scale, minor_scale),
                     labels = NULL,
                     limits = c(0, limity),
                     expand= c(0.03,0.03)) +
  theme(text = element_text(size = 0),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(size = 0.2),
        axis.line = element_line(size = 0.2),
        panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin= unit(c(2, 1, 0, 0), "pt"))
ggsave(paste0(ff.base, "_Chr2.pdf"), device = "pdf", width=2, height=0.3) 

################Chr3###############################
df = filter(.data = X1, c1 %in% c("Chr3")) 
ggplot(df, aes(x=c2, y=c4)) + 
  geom_bar(stat = "identity", fill = 'blue') +
  xlab("") +
  ylab("") +
  scale_x_continuous(breaks = seq(0,30000000,5000000), 
                     labels = NULL,
                     limits = c(-10000, 32990000),　
                     expand= c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(0, main_scale, minor_scale),
                     labels = NULL,
                     limits = c(0, limity),
                     expand= c(0.03,0.03)) +
  theme(text = element_text(size = 0),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(size = 0.2),
        axis.line = element_line(size = 0.2),
        panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin= unit(c(2, 1, 0, 0), "pt"))
ggsave(paste0(ff.base, "_Chr3.pdf"), device = "pdf", width=2, height=0.3) 

################Chr4###############################
df = filter(.data = X1, c1 %in% c("Chr4")) 
ggplot(df, aes(x=c2, y=c4)) + 
  geom_bar(stat = "identity", fill = 'blue') +
  xlab("") +
  ylab("") +
  scale_x_continuous(breaks = seq(0,30000000,5000000), 
                     labels = NULL,
                     limits = c(-10000, 32990000),　
                     expand= c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(0, main_scale, minor_scale),
                     labels = NULL,
                     limits = c(0, limity),
                     expand= c(0.03,0.03)) +
  theme(text = element_text(size = 0),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(size = 0.2),
        axis.line = element_line(size = 0.2),
        panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin= unit(c(2, 1, 0, 0), "pt"))
ggsave(paste0(ff.base, "_Chr4.pdf"), device = "pdf", width=2, height=0.3) 

################Chr5###############################
df = filter(.data = X1, c1 %in% c("Chr5")) 
ggplot(df, aes(x=c2, y=c4)) + 
  geom_bar(stat = "identity", fill = 'blue') +
  xlab("") +
  ylab("") +
  scale_x_continuous(breaks = seq(0,30000000,5000000), 
                     labels = NULL,
                     limits = c(-10000, 32990000),　
                     expand= c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(0, main_scale, minor_scale),
                     labels = NULL,
                     limits = c(0, limity),
                     expand= c(0.03,0.03)) +
  theme(text = element_text(size = 0),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(size = 0.2),
        axis.line = element_line(size = 0.2),
        panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid = element_blank(),
        plot.margin= unit(c(2, 1, 0, 0), "pt"))
ggsave(paste0(ff.base, "_Chr5.pdf"), device = "pdf", width=2, height=0.3) 

}
