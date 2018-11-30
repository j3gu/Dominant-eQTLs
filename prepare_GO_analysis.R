#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
f<- read.table(args[1], header = F)
colnames(f)<-c("topGene", "topSNP.beta1", "beta1","obs.p1", "adj.p1", "k1", "n1", "topSNP.beta2", "beta2","obs.p2", "adj.p2", "k2", "n2", "SNP.nums")
ref<-read.table(args[2])
#ref<-read.table("mart_gene_summary.txt")
get.geneinfo<-function(id){
        index<-which(ref$V2 == id)
        if(length(index)>0){
                return(as.character(ref[index[1],3]))
        }else{
                return(NA)
        }
}

get.genename<-function(id){
        index<-which(ref$V1 == id)
        if(length(index)>0){
                return(as.character(ref[index[1],2]))
        }else{
                return(NA)
        }
}

#divide genes into either significant or non-significant 
p2.adj<-p.adjust(f$adj.p2, method = "BH")
sig.p2<-f[which(p2.adj <= 0.05),]
p1.adj<-p.adjust(f$adj.p1, method = "BH")
#sig.p1 means genes that only show additive effect but no dominant effect
sig.p1<-f[which(p1.adj <= 0.05 & p2.adj > 0.05),]
nonsig.p2<-f[which(p2.adj > 0.05),]
p2sig.genes<-unlist(lapply(as.character(sig.p2$topGene), function(i){strsplit(i, "[.]")[[1]][1]}))
p1sig.genes<-unlist(lapply(as.character(sig.p1$topGene), function(i){strsplit(i, "[.]")[[1]][1]}))

#get the gene name and descripitons 
fg.geneinfo<-unlist(lapply(p2sig.genes, get.geneinfo))
# colnames(fg.geneinfo)<-c("Gene.name", "Gene.description")
fg.genenames<-unlist(lapply(p2sig.genes, get.genename))
bg.genenames<-unlist(lapply(p1sig.genes, get.genename))

#get the coordinates for significant genes
# expr<-read.table(args[3], header = T)
# p2.genes<-as.character(sig.p2$topGene)
# gene.coords<-t(sapply(p2.genes, function(i){
#         index<-which(expr$gene == i)
#         if(length(index)>0){
#                 output<-c(as.character(expr[index, "start"]), as.character(expr[index, "end"]), as.character(expr[index, "chr_num"]))
#         }else{
#                 return(rep(NA, 3))
#         }
# }))
# colnames(gene.coords)<-c("start","end", "chr.num")
# pre.summary<-data.frame(fg.genenames, gene.coords, sig.p2[,1:3], sig.p2$adj.p1, sig.p2[,8:9], sig.p2$adj.p2, fg.geneinfo, stringsAsFactors = F)
# #categorize genes according to its descriptions (MHC, IgG, or T-cell)
# categorize.genes<-function(i){
#         str.list<-strsplit(as.character(i[4]), "[_]")[[1]]
#         start<-as.numeric(i[1])
#         end<-as.numeric(i[2])
#         chrom<-as.character(i[3])
#         if(start >=29000000 & end <= 34000000 & chrom == "chr6"){
#                 return("MHC")
#         }else if("IgG" %in% str.list| "immunoglobulin" %in% str.list){
#                 return("IgG")
#         }else if("T-cell" %in% str.list){
#                 return("T-cell_receptor")
#         }else{
#                 return("Others")
#         }
# }
# f.subset<-cbind(as.matrix(pre.summary[2:4]), pre.summary$fg.geneinfo)
# fg.category<-as.character(unlist(apply(f.subset, 1, categorize.genes)))
# summary<-data.frame(pre.summary, fg.category)

#write.table(data.frame(fg.genenames), args[3], quote = F, row.names = F, col.names = F)
write.table(data.frame(bg.genenames), args[3], quote = F, row.names = F, col.names = F)
#write.table(summary, args[4], quote = F, row.names = F)
