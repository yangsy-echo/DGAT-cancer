required_packages <- c(
  "ggplot2", "latex2exp",  "plyr", "ggpubr","UpSetR"
)
lapply(required_packages, function(package) {
    library(package, character.only = TRUE)  
})
load("./golden_list.RData")
load("./cancer_list.RData")
load("./golden_ppi_gene.RData")
cancer <- c("BLCA","BRCA","COAD","GBM","HNSC","LUAD","STAD")
############### Figure 2A ####################
protein <- read.delim("data/protein_coding_genes.txt",header = T)
risk_list <- readRDS(file="data/risk_list.RDS")

constrained_df <- read.csv(file="/home/yanglab_data3/user/yangsy/DGAT-cancer/evaluation/enrichment/constrained_enrichment_union_p_value.csv")
CGC_OncoKB_IntOGen_df <- read.csv(file="/home/yanglab_data3/user/yangsy/review_DGAT-cancer/plot/github/CGC_OncoKB_IntOGen_df.csv")
shRNA_SCISPR_df <- read.csv(file="/home/yanglab_data3/user/yangsy/review_DGAT-cancer/plot/github/shRNA_SCISPR.csv")
plot_df1 <- rbind(CGC_OncoKB_IntOGen_df,constrained_df,shRNA_SCISPR_df)
plot_df1$geneset <- factor(plot_df1$geneset,levels = c("CGC","OncoKB","IntOGen","constraint","shRNA"))
ggplot(plot_df1,aes(y = cancer,x=-log10(p.adj),color=geneset))+
  facet_wrap(~geneset,ncol = 3)+
  geom_linerange(aes(xmin=0,xmax=-log10(p.adj)))+
  geom_point(aes(size=OR))+
  labs(y="",x=TeX("-log10(\\textit{p}{corrected})"),color="Gene set",shape="Gene set")+
  geom_vline(xintercept = -log10(0.05),linetype="dashed",color="grey20")+
  scale_color_manual(values = c("#66AF71","#EB7E95","#5EA9C2","#987284","#2a5598","#EC7E60"))+
  theme_bw()+guides(color="none")+
  theme(axis.text = element_text(size = 11.5,colour = "grey15"),
        axis.title = element_text(size = 12,colour = "black"),
        strip.text = element_text(size = 11.5,face = "bold"),
        strip.background = element_blank())
############### Figure 2B ####################
xx <- sapply(risk_list[1:7],function(x) x$gene)
xx <- Reduce(c,xx)
yy <- table(xx)
top20 <- names(yy)[order(yy,decreasing = T)[1:20]]
xx1 <- xx2 <-  matrix(nrow = 20,ncol = 9)
for(i in 1:20){
  for(j in 1:7){
    if(is.element(top20[i],risk_list[[j]]$gene)){
      xx1[i,j] <- 1
      xx2[i,j] <- risk_list[[j]]$rank[which(risk_list[[j]]$gene==top20[i],arr.ind = T)]
    }
  }
}
rownames(xx1) <- rownames(xx2) <- top20
colnames(xx1) <- colnames(xx2) <- cancer[1:9]
xx1 <- reshape2::melt(xx1)
xx2 <- reshape2::melt(xx2)
xx1$rank <- xx2$value
colnames(xx1)[1:2] <- c("gene","cancer")
xx1$gene <- factor(xx1$gene,levels = rev(top20))
# emf(file = "/home/yanglab_data3/user/yangsy/review_DGAT-cancer/plot/Figure-2B.emf",width = 7,height = 7)
ggplot(xx1[!is.na(xx1$value),],aes(x=cancer,y=gene,color=rank))+
  geom_point(shape=17,size=3)+
  scale_x_discrete(position = "top")+
  labs(x="",y="Gene",color="Rank")+
  scale_color_gradient(low = "#d62828",high = "#eae2b7")+
  theme_bw()+theme(axis.text = element_text(size = 11,colour = "black"),
                   axis.title = element_text(size=11.5,colour = "black"),
                   legend.position = "bottom")
dev.off()
############### Figure 2C ####################
approved_genes <- readRDS("/home/yanglab_data3/user/yangsy/review_DGAT-cancer/plot/github/drug_overlap.RDS")
show_drug <- data.frame(cancer=cancer[1:7],p.value=1,OR=0,overlap=0,geneset="drug")
drug.func <- function(x){
  length(intersect(x$gene,approved_genes))/nrow(x)
}
set.seed(123)
N <- length(intersect(approved_genes,protein$gene_name))
for(i in 1:7){
  aa <- length(intersect(risk_list[[i]]$gene,approved_genes))
  ff <- fisher.test(matrix(c(aa,nrow(risk_list[[i]])-aa,N,num_base-N),ncol = 2),alternative = "greater")
  show_drug$overlap[i] <- aa
  show_drug$p.value[i] <- ff$p.value
  show_drug$OR[i] <- ff$estimate
  cat(cancer[i],":",ff$p.value,"\n")
}
show_drug$p.adj <- p.adjust(show_drug$p.value,method = "bonferroni")
# emf(file = "Figure-2C.emf",width =6,height = 3.5)
ggplot(show_drug, aes(y = cancer, x = -log10(p.adj), color = geneset)) +
  geom_linerange(aes(xmin = 0, xmax = -log10(p.adj))) +
  geom_point(aes(size = OR), color = "#ec7e60") +
  labs(y = "", x = TeX("-log_{10}(\\textit{p}_{corrected})"), color = "Gene set", shape = "Gene set") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey20") +
  scale_color_manual(values = c("#ec7e60")) +
  theme_bw() +
  guides(color = "none") +
  theme(
    axis.text = element_text(size = 11.5, colour = "grey15"),
    axis.title = element_text(size = 12, colour = "black"),
    strip.text = element_text(size = 11.5, face = "bold"),
    strip.background = element_blank()
  ) +
  coord_flip()  
dev.off()
############### Figure 3A ####################
Fig3A_xx3 <- readRDS("data/Fig3A_xx3.rds")
text_df <- Fig3A_xx3[[1]]
text_df$Method <- factor(text_df$Method,levels = c("DGAT-cancer","OncodriveCLUSTL","OncodriveFML","MutSigCV","DiffMut"))
text_df$y <- (7-as.numeric(text_df$Method))*0.07+0.6
text_df$label <- paste(text_df$Method,sprintf("%0.3f", text_df$value),sep = "—")
Fig3A_xx3[[2]]$cancer <- factor(Fig3A_xx3[[2]]$cancer,levels = cancer)
Fig3A_xx3[[2]]$method <- factor(Fig3A_xx3[[2]]$method,levels = levels(text_df$Method))
# emf(file = "Figure-3A.emf",width = 12,height = 7)#14,6;12,8
ggplot(Fig3A_xx3[[2]],aes(x = recall,y=precision))+
  geom_text(data = text_df[!is.na(text_df$value),],aes(x = 0.96,y=y,label=label,color=Method),hjust=1)+
  guides(color="none")+
  geom_line(aes(color=method))+
  #scale_color_manual(values = c("#EC524B","#8AABF2","#A1676C","#629602"))+
  scale_color_manual(values = c("#EC524B","#8AABF2","#A1676C","#629602","#AF2168","#ED9660"))+
  facet_wrap(~cancer,ncol=4)+
  labs(x = "Recall",y="Precision",color="Method")+
  theme_bw()+
  theme(legend.position = "right",#c(0.9,0.2),
        legend.title = element_text(hjust = .5),
        strip.background = element_rect(fill = "white",colour = "grey"),
        strip.text = element_text(size = 11,face = "bold"),
        axis.text = element_text(size = 11,colour = "black"),
        axis.title = element_text(size = 11.5,colour = "black"))
dev.off()
############### Figure 3B ####################
orther <- Fig3A_xx3[[3]]
cancer <- names(orther)
all <- data.frame()
for (i in cancer){
    combined_df <- do.call(rbind, orther[[i]])
    combined_df$method <- rownames(combined_df)
    combined_df$cancer <- i
    all <- rbind(all,combined_df)
}
data_cleaned <- all %>%
  group_by(cancer) %>%
  distinct(threshold, precision, recall, mcc, f1_scores, specificity, .keep_all = TRUE) %>%
  ungroup()
data_cleaned <- as.data.frame(data_cleaned)
data_cleaned$method <- sub("\\..*", "", data_cleaned$method)
data_cleaned$method <- gsub("oncodriveCLUSTL2","oncodriveCLUSTL",data_cleaned$method )
data_cleaned$method <- gsub("oncodriveFML2","oncodriveFML",data_cleaned$method )
data_cleaned$method <- gsub("mutsigCV2","mutsigCV",data_cleaned$method )
data_cleaned$method <- gsub("DGAT_cancer","DGAT-cancer",data_cleaned$method )
mycol1<-c("#EC524B","#8AABF2","#A1676C","#629602","#ED9660")
data_cleaned$method <- factor(data_cleaned$method, levels = c("DGAT-cancer", "oncodriveCLUSTL", "oncodriveFML", "mutsigCV","DiffMut"))
data_cleaned$Balanced_Accuracy <- (data_cleaned$recall + data_cleaned$specificity)/2
options(repr.plot.width = 13, repr.plot.height = 6)
# emf(file = "Figure-3B_mcc.emf",width = 13,height = 6)#14,6;12,8
ggplot(data = data_cleaned, aes(x = cancer, y = mcc, fill = method)) +     #mcc	f1_scores	Balanced_Accuracy
  geom_bar(position = "dodge", stat = "identity", color = NA, width = 0.8, alpha = 0.7) + 
  scale_fill_manual(values = mycol1) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, max(data_cleaned$mcc) * 1.1)) +   #mcc	f1_scores	Balanced_Accuracy
  theme_classic() +
  labs(title = "", x = "Cancer Type", y = "mcc") +   #mcc	f1_scores	Balanced_Accuracy
  theme(
    plot.margin = unit(c(0.5, 1, 1, 1), "cm"),
    text = element_text(size = 20),  
    legend.position = "top", 
    legend.title = element_blank(), 
    legend.text = element_text(size = 16), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16), 
    axis.text.y = element_text(size = 16) 
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") 
dev.off()
############### Figure 3C ####################
novel_gene <- readRDS(file="/data/novel_gene.RDS")
upsetData=fromList(novel_gene)
outPic=paste0("novel_gene_upset.pdf") 
pdf(file=outPic,onefile = FALSE,width=14,height=7)
upset(
  upsetData,
  nsets = length(novel_gene),             
  nintersects = 20,                        
  order.by = "freq",                       
  show.numbers = "yes",                   
  number.angles = 0,                       
  point.size = 3,                         
  sets.bar.color = "#f28681",               
  main.bar.color = "#f28681",                
  matrix.color = "#bd9598",                    
  line.size = 0.8,                       
  mainbar.y.label = "Gene Intersection count", 
  sets.x.label = "Gene count",
  text.scale = 2                         
)
dev.off()
############### Figure 4A ####################
show_com <- readRDS("data/Figure_4A.RDS")
names(MutSigCV_g) <- names(oncoFML_g) <- names(oncoCLUSTL_g) <- names(DiffMut_g) <- names(DiffMut_g) <- cancer[1:7]
mutsigCV_df <- ldply(MutSigCV_g,data.frame)
mutsigCV_df$group <- "MutSigCV"
oncoCLUSTL_df <- ldply(oncoCLUSTL_g,data.frame)
oncoCLUSTL_df$group <- "OncodriveCLUSTL"
oncoFML_df <- ldply(oncoFML_g,data.frame)
oncoFML_df$group <- "OncodriveFML"
DiffMut_df <- ldply(DiffMut_g,data.frame)
DiffMut_df$group <- "DiffMut"

all_df <- rbind(mutsigCV_df[,c(".id","BYS","group","risk")],oncoFML_df[,c(".id","BYS","group","risk")],
                oncoCLUSTL_df[,c(".id","BYS","group","risk")],
                DiffMut_df[,c(".id","BYS","group","risk")])
colnames(all_df)[1] <- "cancer"
ggboxplot(all_df[!is.na(all_df$risk),],x = "cancer",y="BYS",color = "risk",fill="risk",alpha=.5)+
  facet_wrap(group~.,nrow =5,strip.position = "right")+
  ylim(0,0.011)+
  labs(x = "", y = "PP score in DGAT-cancer")+
  scale_color_manual(values = c("#6F69AC","#FD6F96"),labels=c("q>0.05","q<0.05"))+
  scale_fill_manual(values = c("#6F69AC","#FD6F96"),labels=c("q>0.05","q<0.05"))+
  theme(panel.background = element_rect(fill = "white",colour = "grey"),
        strip.text = element_text(size = 11,face = "bold"),
        axis.text = element_text(size = 11,colour = "black"),legend.position = "bottom",
        axis.title = element_text(size = 11.5),legend.key = element_rect(fill = "white",colour = "white"),
        legend.title = element_blank(), legend.box.spacing = unit(-.5,"cm"))+
  stat_compare_means(aes(group=risk),method = "wilcox.test",paired = F,label = "p.format",label.y = 0.01,
                     size=3.3)
dev.off()
############### Figure 4B ####################
