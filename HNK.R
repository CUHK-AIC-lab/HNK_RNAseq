library(DESeq2)         # DEG identification
library(pcaExplorer)    # PCA plot
library(ggplot2)        # plotting
library(RColorBrewer)   # color palette 
library(ggthemes)       # theme for plots
library(dplyr)          # data manipulation
library(pheatmap)       # heatmap
library(egg)            # save figures
library(fgsea)          # GSEA
library(msigdbr)        # gene sets


##raw counts file
data<-read.csv("gene_count_stringtie.csv",header=T)
data<-data[!duplicated(data$Symbol),]
rownames(data)<-data$Symbol
data<-data[,-1]
sample<-data.frame(sample=colnames(data),Label=c(rep("Sham",3),rep("SNI",3),rep("HNK",3)))
#################################################################################################################
######create the object for Deseq2 and raw counts normalization
ddsFullCountTable <- DESeqDataSetFromMatrix(countData = round(data), colData = sample, design= ~ Label)
dds <- DESeq(ddsFullCountTable)
normalized_counts<-counts(dds,normalized=TRUE)
normalized_counts_mad<-apply(normalized_counts, 1, mad)
normalized_counts<-normalized_counts[order(normalized_counts_mad,decreasing=TRUE),]
write.csv(normalized_counts,file="normalized_counts_HNK.csv")
######PCA plot (Figure 3A)
rld<-rlogTransformation(ddsFullCountTable)
p<-pcaplot(rld,intgroup="Label",ntop=1000,ellipse=FALSE)
set_panel_size(p,file="PCA.pdf",margin=unit(1,"cm"),width=unit(12,"cm"),height=unit(12,"cm"))
#################################################################################################################
######Identification of deferentially expressed genes
###SNI vs Sham
sampleA="SNI"
sampleB="Sham"
contrastV<-c("Label",sampleA,sampleB)  # sampleA / sampleB
res_SNI<-results(dds,contrast=contrastV)
head(res_SNI)
baseA<-counts(dds,normalized=TRUE)[, colData(dds)$Label == sampleA]
if (is.vector(baseA)){
  baseMeanA <-as.data.frame(baseA)
} else {
  baseMeanA <-as.data.frame(rowMeans(baseA))
}

colnames(baseMeanA) <-sampleA
baseB<-counts(dds,normalized=TRUE)[,colData(dds)$Label == sampleB]
if(is.vector(baseB)){
  baseMeanB<-as.data.frame(baseB)
} else {
  baseMeanB <-as.data.frame(rowMeans(baseB))
}

colnames(baseMeanB) <-sampleB

res_SNI<-cbind(baseMeanA, baseMeanB, as.data.frame(res_SNI))
res_SNI<-cbind(ID=rownames(res_SNI),as.data.frame(res_SNI))
res_SNI$baseMean<-rowMeans(cbind(baseA,baseB))
res_SNI$padj[is.na(res_SNI$padj)] <- 1
res_SNI<-res_SNI[order(res_SNI$padj),]
write.csv(res_SNI,"SNI_vs_Sham_DEGs.csv")  #save DEG

## SNI vs Sham vocanol plot (Figure 3B)

res_SNI$logP<--log10(res_SNI$padj)  
res_SNI$Group<-"not-significant"
res_SNI$Group[which((res_SNI$padj<0.05) & (res_SNI$log2FoldChange>=log2(1.5)))]<-"up-regulated"
res_SNI$Group[which((res_SNI$padj<0.05) & (res_SNI$log2FoldChange<= -log2(1.5)))]<-"down-regulated"
table(res_SNI$Group)
res_SNI<-res_SNI[order(res_SNI$padj),]
pdf("Vocanol_SNI.pdf",width=5,height=3)
ggscatter(res_SNI,x="log2FoldChange",y="logP",
          color="Group",
          palette=c("#2f5688","#BBBBBB","#CC0000"),
          size=0.5,
          xlab="log2FoldChange",
          ylab="-log10(Adjust P-value") + theme_base()+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)),linetype="dashed")+
  theme(plot.margin=unit(rep(1,4),'lines'))
ggpar(p,legend="right",legend.title="Group",xlim=c(-10,10))+
  font("legend.title",color="black",face="bold",size=10)+
  font("legend.text",color="black",size=9)+
  font("xy.text",size=10,face="bold")+
  font("xlab",size=10,face="bold")+
  font("ylab",size=10,face="bold")
dev.off()

###HNK vs Sham
sampleC="HNK"
sampleB="Sham"
contrastV<-c("Label",sampleC,sampleB)  # sampleC / sampleB， 注意顺序
res_HNK<-results(dds,contrast=contrastV)
head(res_HNK)
baseA<-counts(dds,normalized=TRUE)[, colData(dds)$Label == sampleA]
if (is.vector(baseA)){
  baseMeanA <-as.data.frame(baseA)
} else {
  baseMeanA <-as.data.frame(rowMeans(baseA))
}

colnames(baseMeanA) <-sampleA
baseB<-counts(dds,normalized=TRUE)[,colData(dds)$Label == sampleB]
if(is.vector(baseB)){
  baseMeanB<-as.data.frame(baseB)
} else {
  baseMeanB <-as.data.frame(rowMeans(baseB))
}

colnames(baseMeanB) <-sampleB

res_HNK<-cbind(baseMeanA, baseMeanB, as.data.frame(res_HNK))
res_HNK<-cbind(ID=rownames(res_HNK),as.data.frame(res_HNK))
res_HNK$baseMean<-rowMeans(cbind(baseA,baseB))
res_HNK$padj[is.na(res_HNK$padj)] <- 1
res_HNK<-res_HNK[order(res_HNK$padj),]
write.csv(res_HNK,"HNK_vs_Sham_DEGs.csv")  #save DEG

## HNK vs Sham vocanol plot (Figure 3B)

res_HNK$logP<--log10(res_HNK$padj)  
res_HNK$Group<-"not-significant"
res_HNK$Group[which((res_HNK$padj<0.05) & (res_HNK$log2FoldChange>=log2(1.5)))]<-"up-regulated"
res_HNK$Group[which((res_HNK$padj<0.05) & (res_HNK$log2FoldChange<= -log2(1.5)))]<-"down-regulated"
table(res_HNK$Group)
res_HNK$Label=""  
res_HNK<-res_HNK[order(res_HNK$padj),]
pdf("Vocanol_HNK.pdf",width=5,height=3)
p<-ggscatter(res_HNK,x="log2FoldChange",y="logP",
             color="Group",
             palette=c("#2f5688","#BBBBBB","#CC0000"),#"#2f5688","#BBBBBB","#CC0000"
             size=0.5,
             xlab="log2FoldChange",
             ylab="-log10(Adjust P-value") + theme_base()+
  geom_hline(yintercept = -log10(0.05),linetype="dashed")+
  geom_vline(xintercept = c(-log2(1.5),log2(1.5)),linetype="dashed")+
  theme(plot.margin=unit(rep(1,4),'lines'))
ggpar(p,legend="right",legend.title="Group",xlim=c(-10,10))+
  font("legend.title",color="black",face="bold",size=10)+
  font("legend.text",color="black",size=9)+
  font("xy.text",size=10,face="bold")+
  font("xlab",size=10,face="bold")+
  font("ylab",size=10,face="bold")
dev.off()

##scatter plot for common DEGs between SNI vs Sham and HNK vs Sham (Figure 3C)
df1<-res_SNI[,c(1,5,9)]
df2<-res_HNK[,c(1,5,9)]
colnames(df1)<-c("ID","log2FoldChange_SNI","padj_SNI")
colnames(df2)<-c("ID","log2FoldChange_HNK","padj_HNK")
df<-left_join(df1,df2,by="ID")
attach(df)
df$label[padj_SNI<0.05 & padj_HNK<0.05 & abs(log2FoldChange_SNI)>log2(1.5) & abs(log2FoldChange_HNK)>log2(1.5)]<-"DEG_Common"
df$label[padj_SNI<0.05 & padj_HNK>0.05 & abs(log2FoldChange_SNI)>log2(1.5)]<-"DEG_SNI"
df$label[padj_SNI>0.05 & padj_HNK<0.05 & abs(log2FoldChange_HNK)>log2(1.5)]<-"DEG_HNK"
df$label[padj_SNI>0.05 & padj_HNK>0.05]<-"NA"
df$label[padj_SNI<0.05 & abs(log2FoldChange_SNI)<log2(1.5)]<-"NA"
df$label[padj_HNK<0.05 & abs(log2FoldChange_HNK)<log2(1.5)]<-"NA"
df<-df[which(df$label!="NA"),]
vertical.lines<-c(-log2(1.5),log2(1.5))
p<-ggplot(df,aes(x=log2FoldChange_SNI,y=log2FoldChange_HNK, color=label))+
  geom_point(size=0.8, alpha=0.5)+
  geom_point(data=df[which(df$label=="DEG_Common"),], alpha=1,size=0.8)+
  scale_color_manual(values = c('red','black','blue'))+
  geom_text_repel(size=2.5,data=subset(df,label=="DEG_Common"),aes(log2FoldChange_SNI,log2FoldChange_HNK,label=ID))+
  geom_vline(xintercept=vertical.lines, linetype="dashed", alpha=0.2)+
  geom_hline(yintercept=vertical.lines, linetype="dashed",alpha=0.2)+
  xlim(-5,5)+ylim(-5,5)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
set_panel_size(p,file="DEG_two_treatments_scattet_plot.pdf",margin=unit(1,"cm"),width=unit(8,"cm"),height=unit(8,"cm"))

#####################################################################################################################
######fgsea
###SNI vs Sham
#retrive gobp from msigdbr package
GO_gene_sets <- msigdbr::msigdbr(
  species = "Mus musculus", 
  category = "C5" # use only GO gene sets
)
GO_gene_sets$gs_name<-gsub("GOBP_","",GO_gene_sets$gs_name) #remove prefix
GO_gene_sets2<-subset(GO_gene_sets, gs_subcat=="GO:BP") #retive only GOBP terms

#prepare GO list for GSEA
GO_list <- split(
  GO_gene_sets2$gene_symbol, 
  GO_gene_sets2$gs_name 
)

#do GSEA
set.seed(1234)
weights_SNI<-as.vector(sign(res_SNI$log2FoldChange) * -log10(res_SNI$pvalue))
names(weights_SNI)<-res_SNI$ID
fgseaRes_SNI<-fgsea(pathways=GO_list,
                stats=weights_SNI,
                eps=0.0,
                minSize=10,
                maxSize=800)
fgseaRes_SNI2<-fgseaRes_SNI
fgseaRes_SNI2$leadingEdge <- vapply(fgseaRes_SNI2$leadingEdge, paste, collapse = ", ", character(1L))
write.csv(fgseaRes_SNI2,"SNI_allGOBP_with_gene_names.csv")

### HNK vs Sham
#do GSEA
set.seed(1234)
weights_HNK<-as.vector(sign(res_HNK$log2FoldChange) * -log10(res_HNK$pvalue))
names(weights_HNK)<-res_HNK$ID
fgseaRes_HNK<-fgsea(pathways=GO_list,
                stats=weights_HNK,
                eps=0.0,
                minSize=10,
                maxSize=800)

fgseaRes_HNK2<-fgseaRes_HNK
fgseaRes_HNK2$leadingEdge <- vapply(fgseaRes_HNK2$leadingEdge, paste, collapse = ", ", character(1L))
write.csv(fgseaRes_HNK2,"HNK_allGOBP_with_gene_names.csv")


######fgseaselect pathway to plot ridgeplot
go_terms<-melt(GO_list)   
head(go_terms)
colnames(go_terms)<-c("Gene","GOBP")
selected_pathway<-c("RESPONSE_TO_PAIN","VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE",
                    "POSITIVE_REGULATION_OF_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC",
                    "PRESYNAPTIC_ENDOCYTOSIS",
                    "SYNAPTIC_SIGNALING")
#ridgeplot for SNI vs Sham (Figure 3E)
s_go_terms_SNI<-go_terms[which(go_terms$GOBP %in% selected_pathway),]
colnames(res_SNI)[1]<-"Gene"
s_go_terms_SNI<-left_join(s_go_terms_SNI,res_SNI,by="Gene")
s_fgseaRes_SNI<-fgseaRes_SNI[,c(1,2,3,5,6,7)]
colnames(s_fgseaRes_SNI)[1]<-"GOBP"
s_go_terms_SNI<-left_join(s_go_terms_SNI,s_fgseaRes_SNI,by="GOBP")

p1<-ggplot(na.omit(s_go_terms_SNI),aes(x=log2FoldChange,y=GOBP,fill=-padj.y))+
  geom_density_ridges_gradient(rel_min_height=0.01, scale=2)+
  scale_fill_viridis_c(name="padj",option="A",alpha=0.8,
                       begin=0.3,end=0.8)+
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8))+
  labs(x="log2FoldChange",title="SNI/Saline")+
  xlim(-5,5)
set_panel_size(p1,margin=unit(1,"cm"),width=unit(16,"cm"),height=unit(10,"cm"),
               file="SNI_GOBP_ridgeplot.pdf")

#ridgeplot for HNK vs Sham (Figure 3E)
s_go_terms_HNK<-go_terms[which(go_terms$GOBP %in% selected_pathway),]
colnames(res_HNK)[1]<-"Gene"
s_go_terms_HNK<-left_join(s_go_terms_HNK,res_SNI,by="Gene")
s_fgseaRes_HNK<-fgseaRes_HNK[,c(1,2,3,5,6,7)]
colnames(s_fgseaRes_HNK)[1]<-"GOBP"
s_go_terms_HNK<-left_join(s_go_terms_HNK,s_fgseaRes_HNK,by="GOBP")

p2<-ggplot(na.omit(s_go_terms_HNK),aes(x=log2FoldChange,y=GOBP,fill=-padj.y))+
  geom_density_ridges_gradient(rel_min_height=0.01, scale=2)+
  scale_fill_viridis_c(name="padj",option="A",alpha=0.8,
                       begin=0.3,end=0.8)+
  theme(axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=10, face="bold"),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8))+
  labs(x="log2FoldChange",title="HNK/Saline")+
  xlim(-5,5)
set_panel_size(p2,margin=unit(1,"cm"),width=unit(16,"cm"),height=unit(10,"cm"),
               file="HNK_GOBP_ridgeplot.pdf")

##retrieve pathway of interest and plot enrichment graph
pathway_of_interest<-"VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE"
#SNI vs Sham (Figure 3I)
p3<-plotEnrichment(GO_list[[pathway_of_interest]],
                   weights_SNI)+
  labs(title=paste0(pathway_of_interest, "\n", "\n",
                    "(NES=",round(fgseaRes_SNI$NES[which(fgseaRes_SNI$pathway==pathway_of_interest)],2), ",  ", 
                    "padj=",round(fgseaRes_SNI$padj[which(fgseaRes_SNI$pathway==pathway_of_interest)],3), ")"),
       x="Rank",y="Enrichment score")+
  theme(title=element_text(size=5,face="bold"),
        axis.text.x=element_text(size=3),
        axis.text.y=element_text(size=3),
        axis.title.x=element_text(size=6),
        axis.title.y=element_text(size=6))
set_panel_size(p3,margin=unit(1,"cm"),width=unit(6,"cm"),height=unit(4,"cm"), 
               file=paste0("SNI_GOBP_",pathway_of_interest,".pdf"))

#HNK vs Sham (Figure 3I)
p4<-plotEnrichment(GO_list[[pathway_of_interest]],
                   weights_HNK)+
  labs(title=paste0(pathway_of_interest, "\n", "\n",
                    "(NES=",round(fgseaRes_HNK$NES[which(fgseaRes_HNK$pathway==pathway_of_interest)],2), ",  ", 
                    "padj=",round(fgseaRes_HNK$padj[which(fgseaRes_HNK$pathway==pathway_of_interest)],3), ")"),
       x="Rank",y="Enrichment score")+
  theme(title=element_text(size=5,face="bold"),
        axis.text.x=element_text(size=3),
        axis.text.y=element_text(size=3),
        axis.title.x=element_text(size=6),
        axis.title.y=element_text(size=6))
set_panel_size(p4,margin=unit(1,"cm"),width=unit(6,"cm"),height=unit(4,"cm"), 
               file=paste0("HNK_GOBP_",pathway_of_interest,".pdf"))

##retrieve pathway of interest and plot heatmap (Figure 3J)
pathway_of_interest<-"VESICLE_MEDIATED_TRANSPORT_IN_SYNAPSE"
t<-subset(fgseaRes_SNI, fgseaRes_SNI$pathway == pathway_of_interest)
t2<-as.vector(unlist(t$leadingEdge))
normalized_count2<-as.data.frame(normalized_count2)[t2,]
annotation_col_selected_pathway<-data.frame(factor(sample$Label))
rownames(annotation_col_selected_pathway)<-colnames(normalized_count2)
colnames(annotation_col_selected_pathway)<-"Group"
annotation_colors<-list(Group=c(Sham='#98b600',SNI='#f39567',HNK="grey"))
pdf(paste0("heatmap_",pathway_of_interest,".pdf"),8,8)
pheatmap(normalized_count2,
         scale = "row",
         color=colorRampPalette(c("navy", "white", "red"))(1000),
         cluster_cols=FALSE,
         cluster_rows=FALSE,
         annotation_col=annotation_col_selected_pathway,
         annotation_colors=annotation_colors,
         cellwidth=10,cellheight=6, gaps_col = c(3,6),
         show_rownames=TRUE,
         show_colnames=TRUE,
         fontsize_row = 7,
         border_color = NA,
         main= pathway_of_interest
)
dev.off()

#end