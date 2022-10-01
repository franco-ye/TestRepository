#########installl package#########
pkg <- c("gplots","ggplot2","RColorBrewer","DESeq2","splitstackshape","factoextra","gplots","ggfortify","cluster","ggplot2","ggrepel","reshape2","easyGgplot2")
for (i in 1:length(pkg)) { if (!requireNamespace(pkg[i])) install.packages(pkg[i]) }
for (i in 1:length(pkg)) { library(pkg[i],character.only=T) }

#########workdir#########\
rm(list = ls())
options <- commandArgs(trailingOnly=T)
matrix <- options[1]
group <- options[2]
count_matrix<-read.csv('genes.count_table.matrix',header = T,row.names = 1,sep = "\t")
count_matrix = round(count_matrix, 0)
countData <- as.matrix(count_matrix)
countData <- countData[,c(7,8,9,1,2,3,4,5,6)]
colData <- read.csv('group.sh',header = F,sep = "\t",row.names = 1)
colnames(colData) <- c("Repeat","condition")
dds<- DESeqDataSetFromMatrix(countData = countData,colData = colData,design = ~condition) 
#Normalization: estimate the effective library size.
dds <- estimateSizeFactors(dds)
dds <- dds[rowSums(counts(dds)) > 1, ]
dds_sub = dds
dds_sub$condition <- droplevels(dds_sub$condition)  #?÷?droplevels(f)estimateSizeFactors(dds_sub) 
#This function obtains dispersion estimates for Negative Binomial distributed data.
cdsBlind = estimateDispersions(dds_sub)  #It doesn't changeations are similar, the rlog might perform a bit better when the size factors vary widely, and the varianceStabilizingTransformation is much faster when there are many samples.
#vsd= rlog(cdsBlind)
vsd= varianceStabilizingTransformation(cdsBlind, blind = F)   #transforms the count data to the log2 scale
#rowVars: Variance ilizingTransformation(dds_sub) 
graph <- plotPCA(vsd, intgroup=c("condition", "Repeat"), returnData=TRUE)
percentVar <- round(100 * attr(graph, "percentVar"))
cols = brewer.pal(8,"Set2")
#plot based on bacteria
#,shape= as.character(Repeat)
p <- ggplot(graph,aes(PC1,PC2,color=condition)) + 
  geom_point(size=3)+
  scale_fill_manual(values = cols) +
  #guides(fill=FALSE, shape= FALSE, colour=FALSE)+
  #stat_ellipse() + 
  #scale_fill_maaes(color =factor(graph$condition))nual(values = cols) +
  scale_color_manual(values = cols) +
  #scale_colour_brewer(palette = "Palette_Name")+
  theme_test()+
  theme(legend.text=element_text(size=12), 
        axis.text=element_text(size=12), 
        axis.title=element_text(size=14,face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  #stat_ellipse()
  ggforce::geom_mark_ellipse(aes(fill=condition,color=condition),level = 0.95 ,type = "confidence",alpha = .3, linetype = 'dashed', show.legend = FALSE)
ggsave(p,filename = "p.pdf",width = 12,height = 10)
