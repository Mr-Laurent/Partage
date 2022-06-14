library(ggplot2)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(gplots)
library(circlize)
library(pals)
library(plyr)


arg1<-"name_experiment"
# Number of modules you want:
nclust<-"100"


ds=obj@assays$RNA@counts
# Remove rows with genes not expressed  (>0 to not have bug for breaks, >1 to not have bug for z
ds<-ds[which(Matrix::rowSums(ds)>1),]
ds_mean<-Matrix::rowMeans(ds)
ds_var<-apply(ds, 1, var)

varmean_df=data.frame(m=ds_mean,v=ds_var,gene=rownames(ds))
rownames(varmean_df)=rownames(ds)

#Compute loess curve to get the variable genes
x=log10(varmean_df$m)
breaks=seq(min(x),max(x),.2)
lv=log2(varmean_df$v/varmean_df$m)
z=sapply(split(lv,cut(x,breaks)),min,na.rm=T)
maskinf=is.infinite(z)
z=z[!maskinf]
b=breaks[-length(breaks)]
b=b[!maskinf]
lo=loess(z~b)

# Selection of the threshold (Adapt with the plot, or the number of genes you want)
inVarMean_MeanThresh=-1
inVarMean_VarmeanThresh=1

#Plot of the loess curve for variable gene selection 
plot(log10(varmean_df$m),log2(varmean_df$v/varmean_df$m),xlab="Log10(mean)",ylab="log2(var/mean)",panel.first=grid())
x1=seq(-1,3.5,l=100)
lline=predict(lo,newdata =x1)
lines(x1,lline+as.numeric(inVarMean_VarmeanThresh),col=2)
abline(v=inVarMean_MeanThresh,col=2)
lline2=predict(lo,newdata =log10(varmean_df$m))


## Selection of the genes 

# Use the same threshold but this time to subset the genes to keep from the matrix
geneModuleMask<-log10(varmean_df$m)>as.numeric(inVarMean_MeanThresh)&log2(varmean_df$v/varmean_df$m)>lline2+as.numeric(inVarMean_VarmeanThresh)
#See how many genes (TRUE) are kept for the module selection
table(geneModuleMask) 

effivar<-rownames(ds)[which(geneModuleMask==T)]
# Just for information, see the number of reads minimum in the object (to be sure we're not taking badly captured genes in our variable genes)
min(rowSums(ds[which(rownames(ds)%in%effivar),])) 

# computing gene-gene correlations on a reduced matrix with only the selected genes, and k-mean clustering to select the wanted number of clusters
ds_select<-ds[effivar,]
res2<-cor(as.matrix(t(ds_select)), method = c("pearson"))
res2_dist=dist(res2,method = 'euclidean')
res2_clust=hclust(res2_dist,method = 'complete')
res2_cut<-cutree(res2_clust, k = nclust)

# Dataframe of the associated module for each gene
dfcut<-as.data.frame(names(res2_cut))
dfcut$related_genes_module<-res2_cut
colnames(dfcut)<-c("X","related_genes_module")
rownames(dfcut)<-rownames(dfcut)


## Writting a csv file with genelist by module:

newmod<-c()
for(i in 1:nclust){ newmod<-append(newmod, paste0("My_Mod_",i))}

for( i in levels(as.factor(dfcut$related_genes_module))){
  assign(paste0("My_Mod_",i),dfcut$X[which(dfcut$related_genes_module==i)])
}

allmod<-c()
alllist<-list()
for(i in 1:nclust){
  submod<-data.frame(Mod=paste0(newmod[i]),Gen=paste0(eval(as.name(newmod[i]))))
  allmod<-rbind(allmod,submod)
  alllist[i]<-list(paste0( paste0(eval(as.name(newmod[i])),collapse=",") ) )
}

write.table(unlist(alllist),file=paste0("Modulelist_",arg1,"_",nclust,"mod_mTh",gsub("-","m",as.character(inVarMean_MeanThresh)),"_vmTh",gsub("-","m",as.character(inVarMean_VarmeanThresh)),".csv"),row.names = newmod,col.names = F)


## Calculating Module score on the seurat object, and visualising them

# Simple module score: log10 of 1e-4 + sum of selected gene expression divided by total gene expression in the cell 
for(i in newmod){
  #Avoid error if a module is only composed of one gene :
  if(length(eval(as.name(i)))<2){values<-log(1e-4+obj@assays$RNA@counts[eval(as.name(i)),]/colSums(obj@assays$RNA@counts))
  eval(parse(text=paste0("obj@meta.data$",i,"score<-values")))
  }else{
  values<-log(1e-4+colSums(obj@assays$RNA@counts[eval(as.name(i)),])/colSums(obj@assays$RNA@counts))
  eval(parse(text=paste0("obj@meta.data$",i,"score<-values")))
}}

obj@meta.data$UMAP_x<-obj@reductions$umap@cell.embeddings[,1]
obj@meta.data$UMAP_y<-obj@reductions$umap@cell.embeddings[,2]

# Printing the figures in one pdf:
pdf(file=paste0("ModuleFigs_",arg1,"_",nclust,"mod_mTh",gsub("-","m",as.character(inVarMean_MeanThresh)),"_vmTh",gsub("-","m",as.character(inVarMean_VarmeanThresh)),".pdf"), width=10, height=10)
pal=rev(rainbow(12))
for( j in newmod){
  p<-ggplot(obj@meta.data, aes_string(x="UMAP_x", y="UMAP_y", color=obj@meta.data[,paste0(j,"score")])) + geom_point(size=1.2)+
    ggtitle(paste("projection of module_",j," scores in Stromal cells",sep=""))+
    scale_color_gradientn(colours = pal)+theme_dark()
  print(p)
}

dev.off()


