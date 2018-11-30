#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
#apply linear regression model to genotype and expression data
f<- read.table(args[1], header = TRUE)
f.byGene<-split(f, f$GENE)
sample.list<-names(f)
geno.start = which(sample.list == "SNP_ID") + 1
geno.end = which(sample.list == "gene") - 1
exp.start = which(sample.list == "length") + 1
exp.end = length(sample.list)
sample_num = length(sample.list[exp.start:exp.end])
#function for linear regression model
linear.model<-function(genotype.ad, genotype.dom, expression){
        model = lm(expression ~ genotype.ad + genotype.dom)
        s = summary(model)
        if(dim(s$coefficients)[1] >=2){
                row.ad = as.vector(s$coefficients["genotype.ad",])#summary for beta1
                row.dom = as.vector(s$coefficients["genotype.dom",])#summary for beta2
                p1.val = row.ad[4]
                p2.val = row.dom[4]
                return(c(p1.val, p2.val))
        }
        else{
                return(rep(NA*2))
        }
}

#function to compute log likelihood for parameters of beta distribution 
get.lklh<-function(shape, x){
        alpha = shape[1]
        beta = shape[2]
        lklh = sum(log(dbeta(x, alpha, beta)))
        return(-1*lklh)
}

compute.mlkh<-function(min.p, input){
        maxLklh<-optim(par = c(1, 5), get.lklh, x = input)
        k<-maxLklh$par[1]
        n<-maxLklh$par[2]
        adj.pval<-pbeta(min.p, k, n)
        return(c(min.p, adj.pval, k, n))
}

center.colmeans<-function(x){
        xcenter = colMeans(x, na.rm = T)
        x-rep(xcenter, rep.int(nrow(x), ncol(x)))
}
#orthogonalize v1 with respect to v2
orthogonalize.vectors<-function(v1, v2){
        v1.proj<-v2*sum(v1*v2, na.rm = T)/sum(v2*v2, na.rm = T)
        v1.ortho<-v1-v1.proj
        return(v1.ortho)
}
#standardize vectors
standardize.v<-function(v){
        v<-as.vector(v)
        std.v<-(v-mean(v, na.rm = T))/(sd(v, na.rm = T)*sqrt(length(v[!is.na(v)])-1))
        return(std.v)
}

turn.zero<-function(cols){
        r<-as.vector(cols)
        r[is.na(r)]<-0
        return(r)
}

compute.pval<-function(r){
        r.val<-as.vector(r)
        n<-sample_num
        t = sqrt(n-2)*r.val/sqrt(1-r.val**2)
        p<-unlist(lapply(as.list(t), function(i){
                if(is.na(i) == T){
                        return(NA)
                }else{
                        pval = 2*pt(abs(i), df = n-2, lower.tail = F)
                        return(pval)  
                }
        }))
        return(p)
}

#test on f.byGene[[3]]
adjust.pvals<-function(df){
        expr<-as.vector(as.matrix(df[1, exp.start:exp.end]))
        if(sd(expr, na.rm = T) != 0){
                genotype.ad<-df[, geno.start : geno.end]
                #assign 0 to genotypes of non-reference homozygotes
                genotype.dom<-apply(genotype.ad, 1, function(i){
                        g.dom<-as.vector(as.matrix(i))
                        g.dom[g.dom == 2] = 0
                        return(g.dom)
                })
                #Permute phenotypes
                permutated.exp<-expr
                for(i in 1:1000){
                        permutated.exp<-rbind(permutated.exp, sample(expr, replace = F))
                }
                #centeralize variables
                ga.centered<-center.colmeans(t(genotype.ad))
                gd.centered<-center.colmeans(genotype.dom)
                expr.centered<-t(center.colmeans(t(permutated.exp)))
                
                merged<-rbind(ga.centered, gd.centered)
                #orthogonalize Gd with respect to Ga for every marker:
                gd.ortho<-apply(merged, 2, function(i){
                        v<-as.vector(i)
                        v1.ortho<-orthogonalize.vectors(v[382:762], v[1:381])
                })
                # #prove to be orthogonal:
                # merged2<-rbind(gd.ortho, ga.centered)
                # ortho.test<-apply(merged2, 2, function(i){
                #         v<-as.vector(i)
                #         test<-sum(v[1:381]*v[382:762])
                #         return(test)
                # })
                #Standardize gd.ortho and ga.centered
                standard.gd<-apply(gd.ortho, 2, standardize.v)
                standard.ga<-apply(ga.centered, 2, standardize.v)
                standard.expr<-t(apply(expr.centered, 1, standardize.v))
                # #test to see if being standardized:
                # sum(apply(standard.ga, 2, function(i){sum(i, na.rm = T)}), na.rm = T)
                # sum(apply(standard.gd, 2, function(i){sum(i, na.rm = T)}), na.rm = T)
                # ga.test<-apply(standard.ga, 2, function(i){sum(i**2, na.rm = T)})
                # gd.test<-apply(standard.gd, 2, function(i){sum(i**2, na.rm = T)})
                # expr.test<-apply(standard.expr, 1, function(i){sum(i**2, na.rm = T)})
                #  #Compute for correlation coefficient:
                #turn all NAs into 0 by columns
                final.ga<-apply(standard.ga, 2, turn.zero)
                final.gd<-apply(standard.gd, 2, turn.zero)
                final.expr<-t(apply(standard.expr, 1, turn.zero))
                rsquared.beta1<-final.expr %*% final.ga
                rsquared.beta2<-final.expr %*% final.gd
                #Compute for T-statistics and p-values
                p.beta1<-apply(rsquared.beta1, 2, compute.pval)
                p.beta2<-apply(rsquared.beta2, 2, compute.pval)
                #Beta estimation 
                min.p1<-apply(p.beta1, 1, min)
                min.p2<-apply(p.beta2, 1, min)
                p1.output<-compute.mlkh(min.p1[1], min.p1[-1])
                p2.output<-compute.mlkh(min.p2[1], min.p2[-1])
                #topSNP
                topSNP.p1<-as.character(df[which(p.beta1[1,] == min.p1[1])[1],5])
                topSNP.p2<-as.character(df[which(p.beta2[1,] == min.p2[1])[1],5])
                output<-c(as.character(df[1,1]), topSNP.p1, p1.output)
                output<-append(output, c(topSNP.p2, p2.output, dim(df)[1]))
                return(output)
        }else{
                return(rep(NA, 12))
        }
}
output<-t(sapply(f.byGene, adjust.pvals))
write.table(output, args[2], quote = F, row.names = F, col.names = F)


