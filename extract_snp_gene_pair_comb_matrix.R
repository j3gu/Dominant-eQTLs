#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
f.output<-read.table(args[1], header = T, stringsAsFactors = F)
f.chr22<-read.table(args[2], header = T)
f.bygene<-split(f.output, f.output$chr.num)

# test.meandiff <- function(inmatrix){
#         expr = as.numeric(inmatrix[291:570])
#         genotype = as.numeric(inmatrix[6:285])
#         comb.matrix.pre <- cbind(genotype, expr)
#         comb.matrix <- split(comb.matrix.pre[,2], comb.matrix.pre[,1])
#         median.matrix <- as.numeric(unlist(lapply(comb.matrix, function(i){quantile(i,0.5)})))
#         SNP.ID <- as.character(inmatrix[5])
#         chr.term <- as.character(inmatrix[287])
#         if ((abs(mean.matrix[2]-median.matrix[1])) <= abs((median.matrix[2]-mean.matrix[3]))){
#                 a<- comb.matrix[[1]]
#                 b<- comb.matrix[[2]]
#         }else{
#                 a<-comb.matrix[[2]]
#                 b<-comb.matrix[[3]]
#         }
#         results <- t.test(a,b, var.equal = FALSE, conf.level = 0.99)
#         p.value <- results$p.value
#         return(p.value)
# }

snp.gene<-lapply(f.bygene, function(i){
        chrom<-i[1,4]
        print(chrom)
        f<-read.table(sprintf("/iblm/netapp/home/j3gu/projects/GTEx/%s/snps_counts_comb/bychrom/%s_snps_counts_comb.txt", args[3],chrom), stringsAsFactors = F)
        row.num<-as.numeric(unlist(apply(i, 1,function(j){
                gene = j[5]
                snp = j[9]
                return(which(f[,1] == gene & f[,5] == snp))
        })))
        return(f[row.num,])
})
output<-c()
for(i in 1:length(snp.gene)){
        output<-rbind(output, snp.gene[[i]])
}
final.output<-data.frame(output)
colnames(final.output)<-names(f.chr22)
write.table(data.frame(final.output), args[4], quote = F, row.names = F)