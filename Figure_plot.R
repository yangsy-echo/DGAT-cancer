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
# MutSigCV_g <- list()
# oncoFML_g <- list()
# oncoCLUSTL_g <- list()
# DiffMut_g <- list()
# TDAmut_g <- list()

# show_com <- data.frame(cancer=cancer, oncodriveFML_p=1, oncodriveCLUSTL_p=1, mutsigCV_p=1, DiffMut_g=1, DTAmut_g=1)

# for(i in 1:7){
#   prop <- read.delim(paste("/home/yanglab_data3/user/yangsy/DGAT-cancer/prop_result/pre/", cancer[i], "_prop_union_p_value.txt", sep=""), header=TRUE)
#   oncodriveFML <- read.delim(paste("/home/yanglab_data3/user/yangsy/DGAT-cancer/evaluation/oncodriveFML/", cancer[i], "_oncodriveFML.tsv", sep=""), header=TRUE)
#   oncodriveFML <- oncodriveFML[!is.na(oncodriveFML$SYMBOL), ]
#   oncodriveFML$BYS <- prop$pp[match(oncodriveFML$SYMBOL, prop$gene)]
#   oncodriveFML$risk <- oncodriveFML$Q_VALUE < 0.05
#   oncoFML_g[[i]] <- oncodriveFML
#   cat("FML:", length(unique(oncodriveFML$GENE_ID[which(oncodriveFML$Q_VALUE < 0.05)])), "\t")
  
#   oncodriveCLUSTL <- read.delim(paste("/home/yanglab_data3/user/yangsy/DGAT-cancer/evaluation/oncodriveCLUSTL/", cancer[i], "_oncodriveCLUSTL.txt", sep=""), header=TRUE)
#   oncodriveCLUSTL <- oncodriveCLUSTL[!is.na(oncodriveCLUSTL$SYMBOL), ]
#   oncodriveCLUSTL$BYS <- prop$pp[match(oncodriveCLUSTL$SYMBOL, prop$gene)]
#   oncodriveCLUSTL$risk <- oncodriveCLUSTL$Q_ANALYTICAL < 0.05
#   oncoCLUSTL_g[[i]] <- oncodriveCLUSTL
#   cat("CLUSTL:", length(unique(oncodriveCLUSTL$SYMBOL[which(oncodriveCLUSTL$Q_ANALYTICAL < 0.05)])), "\t")
  
#   mutsigCV <- read.delim(paste("/home/yanglab_data3/user/yangsy/DGAT-cancer/evaluation/mutsigCV/", cancer[i], "_mutsigCV.txt", sep=""), header=TRUE)
#   mutsigCV <- mutsigCV[!is.na(mutsigCV$gene), ]
#   mutsigCV$BYS <- prop$pp[match(mutsigCV$gene, prop$gene)]
#   mutsigCV$risk <- mutsigCV$q < 0.05
#   MutSigCV_g[[i]] <- mutsigCV
#   cat("MutSig:", length(unique(mutsigCV$gene[which(mutsigCV$q < 0.05)])), "\n")
  
#   TDAmut <- read.delim(paste("/home/yanglab_data3/user/yangsy/review_DGAT-cancer/other_method/TDAmut/", cancer[i], "_2.csv", sep=""), header=TRUE, sep=",")
#   TDAmut <- TDAmut[!is.na(TDAmut$X), ]
#   TDAmut$BYS <- prop$pp[match(TDAmut$X, prop$gene)]
#   TDAmut$risk <- TDAmut$q0 < 0.05
#   TDAmut_g[[i]] <- TDAmut
#   cat("TDAmut:", length(unique(TDAmut$X[which(TDAmut$q0 < 0.05)])), "\n")
  
#   DiffMut <- read.delim(paste("/home/yanglab_data3/user/yangsy/review_DGAT-cancer/other_method/DiffMut/Differential-Mutation-Analysis-master/Output/", cancer[i], "_mut3-DiffMut.txt", sep=""), header=TRUE, sep=" ")
#   DiffMut <- DiffMut[!is.na(DiffMut$protNames), ]
#   DiffMut$BYS <- prop$pp[match(DiffMut$protNames, prop$gene)]
#   DiffMut$risk <- DiffMut$qVal < 0.05
#   DiffMut_g[[i]] <- DiffMut
#   cat("DiffMut:", length(unique(DiffMut$protNames[which(DiffMut$qVal < 0.05)])), "\n")
  
#   t1 <- wilcox.test(oncodriveFML$BYS[which(oncodriveFML$Q_VALUE < 0.05, arr.ind=TRUE)], oncodriveFML$BYS[which(oncodriveFML$Q_VALUE >= 0.05, arr.ind=TRUE)], alternative="greater", conf.int=TRUE, na.rm=TRUE)
#   t2 <- wilcox.test(oncodriveCLUSTL$BYS[which(oncodriveCLUSTL$Q_ANALYTICAL < 0.05, arr.ind=TRUE)], oncodriveCLUSTL$BYS[which(oncodriveCLUSTL$Q_ANALYTICAL >= 0.05, arr.ind=TRUE)], alternative="greater", conf.int=TRUE, na.rm=TRUE)
#   t3 <- wilcox.test(mutsigCV$BYS[which(mutsigCV$q < 0.05, arr.ind=TRUE)], mutsigCV$BYS[which(mutsigCV$q >= 0.05, arr.ind=TRUE)], alternative="greater", conf.int=TRUE, na.rm=TRUE)
  
#   if (all(is.na(TDAmut$BYS[which(TDAmut$q0 < 0.05, arr.ind=TRUE)]))) {
#     t4 <- NA
#     cat("TDAmut all NA \n")
#   } else {
#     t4 <- wilcox.test(TDAmut$BYS[which(TDAmut$q0 < 0.05, arr.ind=TRUE)], TDAmut$BYS[which(TDAmut$q0 >= 0.05, arr.ind=TRUE)], alternative="greater", conf.int=TRUE, na.rm=TRUE)
#   }
#   t5 <- wilcox.test(DiffMut$BYS[which(DiffMut$qVal < 0.05, arr.ind=TRUE)], DiffMut$BYS[which(DiffMut$qVal >= 0.05, arr.ind=TRUE)], alternative="greater", conf.int=TRUE, na.rm=TRUE)
  
#   show_com$oncodriveFML_p[i] <- t1$p.value
#   show_com$oncodriveCLUSTL_p[i] <- t2$p.value
#   show_com$mutsigCV_p[i] <- t3$p.value
#   show_com$TDAmut_p[i] <- ifelse(is.na(t4), NA, t4$p.value)
#   show_com$DiffMut_p[i] <- t5$p.value
# }

# names(MutSigCV_g) <- names(oncoFML_g) <- names(oncoCLUSTL_g) <- names(DiffMut_g) <- names(TDAmut_g) <- names(DiffMut_g) <- cancer[1:7]
# mutsigCV_df <- ldply(MutSigCV_g,data.frame)
# mutsigCV_df$group <- "MutSigCV"
# oncoCLUSTL_df <- ldply(oncoCLUSTL_g,data.frame)
# oncoCLUSTL_df$group <- "OncodriveCLUSTL"
# oncoFML_df <- ldply(oncoFML_g,data.frame)
# oncoFML_df$group <- "OncodriveFML"
# DiffMut_df <- ldply(DiffMut_g,data.frame)
# DiffMut_df$group <- "DiffMut"
# TDAmut_df <- ldply(TDAmut_g,data.frame)
# TDAmut_df$group <- "TDAmut"

# all_df <- rbind(mutsigCV_df[,c(".id","BYS","group","risk")],oncoFML_df[,c(".id","BYS","group","risk")],
#                 oncoCLUSTL_df[,c(".id","BYS","group","risk")],TDAmut_df[,c(".id","BYS","group","risk")],
#                 DiffMut_df[,c(".id","BYS","group","risk")])
# colnames(all_df)[1] <- "cancer"
# all_df <- all_df[all_df$group != "TDAmut",]
all_df <- readRDS("/home/yanglab_data3/user/yangsy/review_DGAT-cancer/plot/github/Figure_4B.RDS")

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
############### Figure 4C ####################
# select_fea_num <- read.delim("/home/yanglab_data3/user/yangsy/DGAT-cancer/evaluation/select_fea_num", header=FALSE, row.names=1)
# compare_list <- list()
# for(i in 1:7){
#   oncodriveFML <- read.delim(paste("/home/yanglab_data3/user/yangsy/DGAT-cancer/evaluation/oncodriveFML/",cancer[i],"_oncodriveFML.tsv",sep = ""),header = T)
#   oncodriveFML <- oncodriveFML[!is.na(oncodriveFML$SYMBOL),]
#   oncodriveCLUSTL <- read.delim(paste("/home/yanglab_data3/user/yangsy/DGAT-cancer/evaluation/oncodriveCLUSTL/",cancer[i],"_oncodriveCLUSTL.txt",sep = ""),header = T)
#   oncodriveCLUSTL <- oncodriveCLUSTL[!is.na(oncodriveCLUSTL$SYMBOL),]
#   mutsigCV <- read.delim(paste("/home/yanglab_data3/user/yangsy/DGAT-cancer/evaluation/mutsigCV/",cancer[i],"_mutsigCV.txt",sep = ""),header = T)
#   mutsigCV <- mutsigCV[!is.na(mutsigCV$gene),]
#   TDAmut <- read.delim(paste("/home/yanglab_data3/user/yangsy/review_DGAT-cancer/other_method/TDAmut/",cancer[i],"_2.csv",sep = ""),header = T,sep=",")
#   TDAmut <-TDAmut[is.na(TDAmut$X),]
#   DiffMut <- read.delim(paste("/home/yanglab_data3/user/yangsy/review_DGAT-cancer/other_method/DiffMut/Differential-Mutation-Analysis-master/Output/",cancer[i],"_mut3-DiffMut.txt",sep = ""),header = T,sep=" ")
#   DiffMut <- DiffMut[!is.na(DiffMut$protNames),]
#   compare_list[[i]] <- cbind(data.frame(gene=rownames(cancer.list[[i]])),cancer.list[[i]][,as.numeric(select_fea_num[i,])])
#   compare_list[[i]]$oncodriveFML <- oncodriveFML$Q_VALUE[match(compare_list[[i]]$gene,oncodriveFML$SYMBOL)]
#   compare_list[[i]]$oncodriveCLUSTL <- oncodriveCLUSTL$Q_ANALYTICAL[match(compare_list[[i]]$gene,oncodriveCLUSTL$SYMBOL)]
#   compare_list[[i]]$mutsigCV <- mutsigCV$q[match(compare_list[[i]]$gene,mutsigCV$gene)]
#   compare_list[[i]]$TDAmut <- TDAmut$q[match(compare_list[[i]]$gene,TDAmut$X)]
#   compare_list[[i]]$DiffMut <- DiffMut$q[match(compare_list[[i]]$gene,DiffMut$protNames)]
# }
# select_fea <- Reduce(union,sapply(cancer.list,colnames))
# select_fea2 <- select_fea[-c(23,25,27)]
# select_fea2 <- select_fea2[c(1:20,22,21,23,24)]
# plot_df <- matrix(nrow = 9,ncol = 24)
# rownames(plot_df) <- cancer[1:7]
# colnames(plot_df) <- select_fea2

# for(k in 1:7){
#   xx <- compare_list[[k]]
#   cols_to_fill <- c("oncodriveFML", "oncodriveCLUSTL", "mutsigCV", "TDAmut", "DiffMut")
#   xx[cols_to_fill] <- lapply(xx[cols_to_fill], function(x) {
#     x[is.na(x)] <- 1
#     return(x)
#   })
#   xx$group <- (xx$oncodriveCLUSTL<0.05 | xx$oncodriveFML<0.05 | xx$mutsigCV<0.05 |xx$TDAmut<0.05|xx$DiffMut<0.05) & !(xx$gene%in%risk_list[[k]]$gene)
#   xx$group2 <- xx$oncodriveCLUSTL>0.05 & xx$oncodriveFML>0.05 & xx$mutsigCV>0.05 & xx$TDAmut>0.05 & xx$TDAmut>0.05 & xx$gene%in%risk_list[[k]]$gene
#   cat(cancer[k],"only predicted by other methods",sum(xx$group,na.rm = T),"\t","only predicted by DGAT-cancer",sum(xx$group2,na.rm = T),"\n")
#   xx2 <- reshape2::melt(xx[!is.na(xx$group),c(-1,-22,-23,-24)],id.vars=c("group","group2"))
#   xx2$type <- NA
#   xx2$type[xx2$group] <- "only predicted by other methods"
#   xx2$type[xx2$group2] <- "only predicted by DGAT-cancer"
#   xx2 <- xx2[!is.na(xx2$type),]
#   for(j in select_fea2){
#     if(j%in%xx2$variable & (k!=7)){
#       tt1 <- wilcox.test(xx2$value[xx2$variable==j & xx2$type=="only predicted by other methods"],xx2$value[xx2$variable==j & xx2$type=="only predicted by DGAT-cancer"],paired = F,alternative = "greater")
#       tt2 <- wilcox.test(xx2$value[xx2$variable==j & xx2$type=="only predicted by other methods"],xx2$value[xx2$variable==j & xx2$type=="only predicted by DGAT-cancer"],paired = F,alternative = "less")
#       plot_df[cancer[k],j] <- min(tt1$p.value,tt2$p.value)
#     }else if(j%in%xx2$variable & (k==7)){
#       test <- try(t.test(xx2$value[xx2$variable==j & xx2$type=="only predicted by other methods"],xx2$value[xx2$variable==j & xx2$type=="only predicted by DGAT-cancer"],paired = F,alternative = "greater"),TRUE)
#       if("try-error"%in%class(test)){
#         next
#       }else{
#         plot_df[cancer[k],j] <- test$p.value
#       }
#     }
#   }
#   plot_df[k,] <- p.adjust(plot_df[k,],method = "bonferroni")
# }
# plot_df2 <- reshape2::melt(plot_df)
# plot_df2$label <- ifelse(plot_df2$value<0.0001,"***",ifelse(plot_df2$value<0.001,"**",ifelse(plot_df2$value<0.05,"*",NA)))

# label <- sapply(select_fea2,function(x) strsplit(x,split = "_")[[1]][1])
# label[c(9,11,13,14,20:24)] <- c("M-CAP","fathmm-MKL","integrated_fitCons","GERP++","JSD score","C score","uEMD-Ex","tumor_med","normal_med")

# plot_df3 <- NULL
# for(k in 1:7){
#   xx <- compare_list[[k]]
#   xx$group <- (xx$oncodriveCLUSTL<0.05 | xx$oncodriveFML<0.05 | xx$mutsigCV<0.05|xx$DiffMut<0.05|xx$TDAmut<0.05) & !(xx$gene%in%risk_list[[k]]$gene)
#   xx$group2 <- xx$oncodriveCLUSTL>0.05 & xx$oncodriveFML>0.05 & xx$mutsigCV>0.05 & xx$DiffMut>0.05 & (xx$TDAmut>0.05 |is.na(xx$TDAmut))& xx$gene%in%risk_list[[k]]$gene
#   cat(cancer[k],"only predicted by other methods",sum(xx$group,na.rm = T),"\t","only predicted by DGAT-cancer",sum(xx$group2,na.rm = T),"\n")
#   xx2 <- reshape2::melt(xx[!is.na(xx$group),c(-1,-22,-23,-24)],id.vars=c("group","group2"))
#   xx2$type <- NA
#   xx2$type[xx2$group] <- "only predicted by other methods"
#   xx2$type[xx2$group2] <- "only predicted by DGAT-cancer"
#   xx2 <- xx2[!is.na(xx2$type) & (xx2$variable%in%select_fea2),]
#   xx2$cancer <- cancer[k]
#   plot_df3 <- rbind(plot_df3,xx2)
# }
# plot_df3$variable2 <- label[match(plot_df3$variable,names(label))]
# plot_df3$variable2 <- factor(plot_df3$variable2,levels = label)
# aggre_df3 <- plot_df2
# aggre_df3$variable2 <- label[match(aggre_df3$Var2,names(label))]
# aggre_df3$variable2 <- factor(aggre_df3$variable2,levels = label)
# colnames(aggre_df3)[1] <- "cancer"
# tt <- aggregate(value~variable2,data = plot_df3,function(x) max(x)*0.95)
# aggre_df3$y <- tt$value[match(aggre_df3$variable2,tt$variable2)]

Figure_4C <- list(aggre_df3,plot_df3) 
readRDS(file="data/Figure_4C.RDS")
ggplot(Figure_4C[[2]][Figure_4C[[2]]$variable%in%select_fea2[20:24],],aes(x = type,y=value))+
  geom_boxplot(aes(fill=type),width=0.6,outlier.size = .7)+
  geom_jitter(aes(fill=type),position = position_jitterdodge(.4),shape = 21,size=.7,alpha=.8)+
  geom_text(data = Figure_4C[[1]][Figure_4C[[1]]$variable2%in%label[20:24],],aes(x = 1.5,y=y,label=label),color="red",size=5.5)+
  facet_grid(variable2~cancer,scales = "free_y")+
  labs(x = "",y="score",fill = "")+
  scale_fill_manual(values = c("#FF8981","#00C678"),labels=c("genes missed by other methods","genes missed by DAGT-cancer"))+
  theme_bw()+
  theme(axis.text.x = element_blank(),legend.position = "bottom",
        strip.background = element_rect(fill = "white",colour = "grey"),strip.text = element_text(size = 11,face = "bold"),
        axis.text.y = element_text(size = 11,colour = "black"),legend.box.spacing = unit(-5,"mm"),
        axis.title = element_text(size = 11.5,colour = "black"))
dev.off()  
