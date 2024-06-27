# This R script plots the number of TE insertions per 10 kb in each CEN region. 
# The input file is "_somatic_ins_10 kb_windows.csv" from TED-seq analysis. 
rm(list=ls())
library(ggplot2)
library(dplyr)
library(readr)
setwd("") # set working directory
files  <- list.files(pattern = "_somatic_ins_10kb_windows.csv")   
for (file.name in files) 
{
  X1 <- read.table(file.name, header = FALSE, sep = "\t", col.names = c("c1","c2","c3","c4","c5","c6","c7"))  # quote = '"',ファイルをデータファイルに読み込む
  View(X1)
  ff.base <- sub('.csv', "", file.name)     
  
  if((max(X1$c4)>=0 & max(X1$c4)<=10))
  {
    main_scale = 10
    minor_scale = 5
  }
  else if ((max(X1$c4)>10 & max(X1$c4)<=100))
  {
    main_scale = trunc(max(X1$c4)/10)*10
    minor_scale = trunc(max(X1$c4)/10)*10/2
  }
  else if ((max(X1$c4)>100 & max(X1$c4)<=1000))
  {
    main_scale = trunc(max(X1$c4)/100)*100
    minor_scale = trunc(max(X1$c4)/100)*100/2
  }
  else if ((max(X1$c4)>1000 & max(X1$c4)<=10000))
  {
    main_scale = trunc(max(X1$c4)/1000)*1000
    minor_scale = trunc(max(X1$c4)/1000)*1000/2
  }
  else
  {
    main_scale = "out of scale"
    minor_scale= "out of scale"
  }
  print(file.name)
  print(main_scale)
  print(minor_scale)
  
################CEN1###############################
df = filter(.data = X1, c1 %in% c("Chr1") & (c2>13200000 & c2<=19200000)) #CEN1
#df = filter(.data = X1, c1 %in% c("Chr2") & (c2>1935000 & c2<=7935000)) #CEN2
#df = filter(.data = X1, c1 %in% c("Chr3") & (c2>11665000 & c2<=17665000)) #CEN3
#df = filter(.data = X1, c1 %in% c("Chr4") & (c2>2590000 & c2<=8590000)) #CEN4
#df = filter(.data = X1, c1 %in% c("Chr5") & (c2>10170000 & c2<=16170000)) #CEN5
ggplot(df, aes(x=c2, y=c4)) + 
  geom_bar(stat = "identity", fill = 'blue') +
  xlab("") + ylab("") +
  scale_x_continuous(breaks = seq(14000000,19000000,1000000), #CEN1
  #scale_x_continuous(breaks = seq(2000000,7000000,1000000), #CEN2
  #scale_x_continuous(breaks = seq(12000000,17000000,1000000), #CEN3
  #scale_x_continuous(breaks = seq(3000000,8000000,1000000), #CEN4
  #scale_x_continuous(breaks = seq(11000000,16000000,1000000), #CEN5
                     labels = NULL,
                     expand= c(0,0)) +
  scale_y_continuous(breaks = seq(0, main_scale, minor_scale),  
                     labels = NULL,
                     limits = c(0, max(X1$c4))) +
  theme(text = element_text(size = 0),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(size = 0.2),
        panel.border = element_rect(fill = NA, size=0.3),
        panel.background = element_blank(),
        panel.grid = element_blank())

ggsave(paste0(ff.base, "_CEN1.pdf"), device = "pdf", width=1.25, height=0.89) 

###########CEN2######################
#df = filter(.data = X1, c1 %in% c("Chr1") & (c2>13200000 & c2<=19200000)) #CEN1
df = filter(.data = X1, c1 %in% c("Chr2") & (c2>1935000 & c2<=7935000)) #CEN2
#df = filter(.data = X1, c1 %in% c("Chr3") & (c2>11665000 & c2<=17665000)) #CEN3
#df = filter(.data = X1, c1 %in% c("Chr4") & (c2>2590000 & c2<=8590000)) #CEN4
#df = filter(.data = X1, c1 %in% c("Chr5") & (c2>10170000 & c2<=16170000)) #CEN5
ggplot(df, aes(x=c2, y=c4)) + 
  geom_bar(stat = "identity", fill = 'blue') +
  xlab("") + ylab("") +
  #scale_x_continuous(breaks = seq(14000000,19000000,1000000), #CEN1
  scale_x_continuous(breaks = seq(2000000,7000000,1000000), #CEN2
  #scale_x_continuous(breaks = seq(12000000,17000000,1000000), #CEN3
  #scale_x_continuous(breaks = seq(3000000,8000000,1000000), #CEN4
  #scale_x_continuous(breaks = seq(11000000,16000000,1000000), #CEN5
                     labels = NULL,
                     expand= c(0,0)) +
  scale_y_continuous(breaks = seq(0, main_scale, minor_scale),  
                     labels = NULL,
                     limits = c(0, max(X1$c4))) +
  theme(text = element_text(size = 0),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(size = 0.2),
        panel.border = element_rect(fill = NA, size=0.3),
        panel.background = element_blank(),
        panel.grid = element_blank())

ggsave(paste0(ff.base, "_CEN2.pdf"), device = "pdf", width=1.25, height=0.89) 

###########CEN3######################
#df = filter(.data = X1, c1 %in% c("Chr1") & (c2>13200000 & c2<=19200000)) #CEN1
#df = filter(.data = X1, c1 %in% c("Chr2") & (c2>1935000 & c2<=7935000)) #CEN2
df = filter(.data = X1, c1 %in% c("Chr3") & (c2>11665000 & c2<=17665000)) #CEN3
#df = filter(.data = X1, c1 %in% c("Chr4") & (c2>2590000 & c2<=8590000)) #CEN4
#df = filter(.data = X1, c1 %in% c("Chr5") & (c2>10170000 & c2<=16170000)) #CEN5
ggplot(df, aes(x=c2, y=c4)) + 
  geom_bar(stat = "identity", fill = 'blue') +
  xlab("") + ylab("") +
  #scale_x_continuous(breaks = seq(14000000,19000000,1000000), #CEN1
  #scale_x_continuous(breaks = seq(2000000,7000000,1000000), #CEN2
  scale_x_continuous(breaks = seq(12000000,17000000,1000000), #CEN3
  #scale_x_continuous(breaks = seq(3000000,8000000,1000000), #CEN4
  #scale_x_continuous(breaks = seq(11000000,16000000,1000000), #CEN5
                     labels = NULL,
                     expand= c(0,0)) +
  scale_y_continuous(breaks = seq(0, main_scale, minor_scale),  
                     labels = NULL,
                     limits = c(0, max(X1$c4))) +
  theme(text = element_text(size = 0),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(size = 0.2),
        panel.border = element_rect(fill = NA, size=0.3),
        panel.background = element_blank(),
        panel.grid = element_blank())
       
ggsave(paste0(ff.base, "_CEN3.pdf"), device = "pdf", width=1.25, height=0.89) 

###########CEN4######################
#df = filter(.data = X1, c1 %in% c("Chr1") & (c2>13200000 & c2<=19200000)) #CEN1
#df = filter(.data = X1, c1 %in% c("Chr2") & (c2>1935000 & c2<=7935000)) #CEN2
#df = filter(.data = X1, c1 %in% c("Chr3") & (c2>11665000 & c2<=17665000)) #CEN3
df = filter(.data = X1, c1 %in% c("Chr4") & (c2>2590000 & c2<=8590000)) #CEN4
#df = filter(.data = X1, c1 %in% c("Chr5") & (c2>10170000 & c2<=16170000)) #CEN5
ggplot(df, aes(x=c2, y=c4)) + 
  geom_bar(stat = "identity", fill = 'blue') +
  xlab("") + ylab("") +
  #scale_x_continuous(breaks = seq(14000000,19000000,1000000), #CEN1
  #scale_x_continuous(breaks = seq(2000000,7000000,1000000), #CEN2
  #scale_x_continuous(breaks = seq(12000000,17000000,1000000), #CEN3
  scale_x_continuous(breaks = seq(3000000,8000000,1000000), #CEN4
  #scale_x_continuous(breaks = seq(11000000,16000000,1000000), #CEN5
                     labels = NULL,
                     expand= c(0,0)) +
  scale_y_continuous(breaks = seq(0, main_scale, minor_scale),  
                     labels = NULL,
                     limits = c(0, max(X1$c4))) +
  theme(text = element_text(size = 0),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(size = 0.2),
        panel.border = element_rect(fill = NA, size=0.3),
        panel.background = element_blank(),
        panel.grid = element_blank())
       
ggsave(paste0(ff.base, "_CEN4.pdf"), device = "pdf", width=1.25, height=0.89) 

###########CEN5######################
#df = filter(.data = X1, c1 %in% c("Chr1") & (c2>13200000 & c2<=19200000)) #CEN1
#df = filter(.data = X1, c1 %in% c("Chr2") & (c2>1935000 & c2<=7935000)) #CEN2
#df = filter(.data = X1, c1 %in% c("Chr3") & (c2>11665000 & c2<=17665000)) #CEN3
#df = filter(.data = X1, c1 %in% c("Chr4") & (c2>2590000 & c2<=8590000)) #CEN4
df = filter(.data = X1, c1 %in% c("Chr5") & (c2>10170000 & c2<=16170000)) #CEN5
ggplot(df, aes(x=c2, y=c4)) + 
  geom_bar(stat = "identity", fill = 'blue') +
  xlab("") + ylab("") +
  #scale_x_continuous(breaks = seq(14000000,19000000,1000000), #CEN1
  #scale_x_continuous(breaks = seq(2000000,7000000,1000000), #CEN2
  #scale_x_continuous(breaks = seq(12000000,17000000,1000000), #CEN3
  #scale_x_continuous(breaks = seq(3000000,8000000,1000000), #CEN4
  scale_x_continuous(breaks = seq(11000000,16000000,1000000), #CEN5
                     labels = NULL,
                     expand= c(0,0)) +
  scale_y_continuous(breaks = seq(0, main_scale, minor_scale), 
                     labels = NULL,
                     limits = c(0, max(X1$c4))) +
  theme(text = element_text(size = 0),
        axis.ticks.length = unit(0.5, "mm"),
        axis.ticks = element_line(size = 0.2),
        panel.border = element_rect(fill = NA, size=0.3),
        panel.background = element_blank(),
        panel.grid = element_blank())

ggsave(paste0(ff.base, "_CEN5.pdf"), device = "pdf", width=1.25, height=0.89) 

}
