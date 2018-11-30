library(peer)
#!/iblm/netapp/home/j3gu/anaconda2/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
#Apply PEER factor to remove observed and hidden covariates from gene expression matrix
f.expr <- read.table(args[1], header = TRUE)
expr.matrix<-t(f.expr[,-5:-1])
genes.info<-f.expr[,1:5]
#Covariates Matrix:
f.genotype<-read.table(args[2], header = TRUE) #pca output
genotype.cov<-f.genotype[,1:3]
f.env<-read.table(args[3], header = T)
#f.env<-read.table("/iblm/netapp/home/j3gu/data1/external/PROTECTED/GTEx/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v6.p1.c1.GRU/processed/GTEx_Subject_other_phenotypes.txt", header = T)
#Take subsets of observed covariates matrix 
subjects<-unlist(lapply(names(f.expr)[-5:-1], function(i){
        x  = strsplit(i, "[.]")[[1]][2]
        return(x)
        }))
id.name<-unlist(lapply(as.character(f.env$SUBJID), function(i){
        id = strsplit(i, "[-]")[[1]][2]
        return(id)
}))
index<-unlist(lapply(subjects, function(i){
        return(which(id.name == i))
}))
env.cov<-as.matrix(f.env[index, -2:-1])
all.covariates<-cbind(genotype.cov, env.cov)

get_simple_PEER_object<-function(K, Nmax, expr_file){
        model = PEER()
        #set data and parameters
        PEER_setNk(model, K)
        PEER_setPhenoMean(model, as.matrix(expr_file))
        PEER_setNmax_iterations(model, Nmax_iterations = Nmax)
        model
}
model = get_simple_PEER_object(K = 15, Nmax = 500, expr_file = expr.matrix) 
#set covariates
PEER_setCovariates(model, as.matrix(all.covariates))
# perform inference
PEER_update(model) 
#investigate results
# #plot the posterior variance
# pdf("../../beta_approx_WB/PEER/posterior_variance_plot.pdf")
# PEER_plotModel(model)
# dev.off()

# #factors:
# X = PEER_getX(model)
# #weights:
# W = PEER_getW(model)
# #ARD parameters
# Alpha = PEER_getAlpha(model)
#get corrected dataset:
Yc = PEER_getResiduals(model)
#get covariates:
#Cv = PEER_getCovariates(model)

#output
new.expr<-t(Yc)
colnames(new.expr)<-names(f.expr)[-5:-1]
final.df<-data.frame(genes.info, new.expr)
#df.byChr<-split(final.df, final.df$chr_num)
#head(final.df)
write.table(final.df, args[4], row.names = F, quote = F)
#write.table(final.df, "./Skin_Sun_Exposed_Lower_leg/PEER/all_chr_PEER_expr.txt", row.names = F, quote = F)

# for(i in 1:22){
#         chrom <-names(df.byChr)
#         print(chrom)
#         write.table(df.byChr[[i]], args[i+3], row.names = F, quote = F)
# }

#write.table(X, "../../beta_approx_WB/PEER/list_of_PEER_factors.txt", row.names = F, quote = F)

