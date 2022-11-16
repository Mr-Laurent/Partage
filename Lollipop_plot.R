setwd("~/Lollipop/")

df<-read.csv(file="cohorteBatch2.csv",sep=";",header=T)

library(ggplot2)
library(scales)
theme_set(theme_classic())

# Reorder manually by batch :
df$CODE<- factor(df$CODE, levels = rev(c("CT3","CR4","HR4","CT4","CR5","HR5","CT5","CR6","HR6")))
# Adding the batch information
# Each samples have 3 time points, except the first sample from Run10 (4 time points) :
df$Run<-c(rep("Run8",3),rep("Run9",3),rep("Run10",4),rep(c(rep("Run8",3),rep("Run9",3),rep("Run10",3)),2) )
df$Run<- factor(df$Run, levels = c("Run8","Run9","Run10"))


pdf("Lollipop_plot_Batch2.pdf")

ggplot(df, aes(x=CODE, y=delai))+ 
  geom_segment(aes(x=CODE,xend=CODE,y=0,yend=delai,colour=Run),size=2,)+
  geom_point(aes(x=CODE, y=delai, colour=Run), size=5)+
  geom_point(data=df,aes(CODE,temps.rejet,color="Rejection"),shape=18,size=8)+
  scale_color_manual(values = c("Rejection"="#D8A131","Run8"="#64BED4","Run9"="#428CBE","Run10"="#305DAD"))+
  coord_flip()+ labs(title="Blood samples through time")+
  theme(axis.text=element_text(size=17,family="sans"),
        axis.title=element_text(size=14, family="sans"),
        plot.title = element_text(size=20,hjust = 0.5),
        legend.title = element_blank())+xlab(label = '') +
  ylab(label = 'post-transplantation time (in days)')+ ylim(-1, max(df$delai))
  
dev.off()
