library(plyr)
library(dplyr)
library(ggplot2)
library(Seurat)
library(Matrix)
library(ggalluvial)

obj<-readRDS(file="myobject.rds")


df1<-data.frame(names=colnames(obj),data.frame(obj@meta.data$temps1))
colnames(df1)<-c("cells","obj_temps1")
df2<-data.frame(names=colnames(obj),data.frame(obj@meta.data$temps2))
colnames(df2)<-c("cells","obj_temps2")

dfall<-cbind(df1,df2$obj_temps2)
colnames(dfall)<-c("cells","obj_temps1","obj_temps2")

#Transform data.frame to a readable table for ggplot : use melt()
repartition<-melt(table(dfall$obj_pred,dfall$obj_pred2))
colnames(repartition)<-c("Temps1","Temps2","Freq")

# Taken from https://stackoverflow.com/questions/71661440/making-an-alluvial-sankey-diagram-using-the-first-axis-as-the-fill

pdf("Sankey_alluvial_ChangeTemps.pdf")
ggplot(data = repartition,
       aes(y = Freq, axis1 = Temps1, 
           axis2 = Temps2)) +
  geom_alluvium(aes(fill = Temps1),
                width = 1/12) +
  geom_stratum(width = 1/12, fill ="grey", color = "black") +
  geom_text(x = 0.95, stat = "stratum", 
            aes(label = Temps1),
            color = 'black', hjust = 1) +
  geom_text(x = 2.05, stat = "stratum", 
            aes(label = Temps2),
            color = 'black', hjust = 0) +
  scale_x_discrete(limits = c("Before", "After"), expand = c(.2, .2)) +
  ggtitle("Changes in cell proportions through time") +
  theme_void() +
  theme(legend.position = 'none')
dev.off()

