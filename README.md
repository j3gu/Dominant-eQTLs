# dominant-eQTL
This analysis pipeline can be used to detect loci having non-additive genetic assoications with gene expression levels. 


## Theory

We applied a multiple linear regression model to identify dominant associaions between SNPs and gene expression levels. 

![image](https://user-images.githubusercontent.com/17937089/49409717-adadf480-f716-11e8-8c8b-f01a83d8e930.png)

Ga stands for the number of non-reference alleles. Ga = 0 if genotype is reference homozygous; Ga = 1 of genotype is heterozygous; Ga = 2 if genotype is non-reference homozygous. This is the variable that captures the additive effects, which is commonly used in the analysis of gene expression quantitative loci(eQTLs).

Gd captures dominant effects, where Gd = 1 if genotype is heterozygous and Gd = 0 if genotype is either reference or non-reference homozygous. 

E denotes the observed gene expression levels. Our model assumes that the noise across samples is normally distributed with 0 mean and some variance value.

Our null hypothesis is that there is no evidence for a dominant association, which is equivalent to β2 not significantly different from zero. An ANOVA model can be used instead of a multiple linear regression to detect both additive and dominant effects of the genotype; however the F-test used in the ANOVA model aims to test if at least one beta is significantly different from zero. This does not fulfill our objective, as we specifically wish to test whether there is evidence showing that β2 is not equal to zero. Alternatively, a t-test can be used to test the effect sizes of Ga and Gd separately.

Here is an illustration of what dominant eQTL looks like.
![image](https://user-images.githubusercontent.com/17937089/49410303-0088ab80-f719-11e8-8859-77401fc4cd22.png)


## Workflow

