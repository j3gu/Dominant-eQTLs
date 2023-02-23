# dominant-eQTL
This analysis pipeline can be used to detect loci having non-additive genetic assoications with gene expression levels. 


## Theory

We applied a multiple linear regression model to identify dominant associaions between SNPs and gene expression levels. 

![image](https://user-images.githubusercontent.com/17937089/49409717-adadf480-f716-11e8-8c8b-f01a83d8e930.png)

Ga stands for the number of non-reference alleles. $G_A$ = 0 if genotype is reference homozygous; Ga = 1 of genotype is heterozygous; Ga = 2 if genotype is non-reference homozygous. This is the variable that captures the additive effects, which is commonly used in the analysis of gene expression quantitative loci(eQTLs).

Gd captures dominant effects, where Gd = 1 if genotype is heterozygous and Gd = 0 if genotype is either reference or non-reference homozygous. 

E denotes the observed gene expression levels. Our model assumes that the noise across samples is normally distributed with 0 mean and some variance value.

Our null hypothesis is that there is no evidence for a dominant association, which is equivalent to β2 not significantly different from zero. An ANOVA model can be used instead of a multiple linear regression to detect both additive and dominant effects of the genotype; however the F-test used in the ANOVA model aims to test if at least one beta is significantly different from zero. This does not fulfill our objective, as we specifically wish to test whether there is evidence showing that β2 is not equal to zero. Alternatively, a t-test can be used to test the effect sizes of Ga and Gd separately.

Here is an illustration of what dominant eQTL looks like.
![image](https://user-images.githubusercontent.com/17937089/49410303-0088ab80-f719-11e8-8859-77401fc4cd22.png)


## Workflow

![image](https://user-images.githubusercontent.com/17937089/49410437-72f98b80-f719-11e8-90cc-99cbdea859e6.png)

![image](https://user-images.githubusercontent.com/17937089/49410444-83116b00-f719-11e8-9ceb-1a8ead1c4e77.png)

## Getting Started
### Dependencies
* snakemake 3.13.3
* python2.7
* numpy
* R

### Snakemake Workflow
Input files: Genotype file (genotype.txt) and RNA-Seq counts matrix file (gene_expr.txt) downloaded from GTEx project
1. Run PCA analysis on the genotype matrix
```
Rscript --vanilla pca_genotype.R genotype genotype.txt genotype_pca.txt 

```

2. Regress out sex, age, race and hidden covariates from gene expression matrix using PEER package
```
Rscript --vanilla peer_factor.R gene_expr.txt genotype_pca.txt phenotype_info.txt gene_expr_PEER.txt

```

3. Split jobs for parallele computation:1000 genes per job
 ```
 python split_into_files.py gene_expr_PEER.txt {files}_peer_expr.txt

 ```
 
4. Merge SNPs that are within +/- 100kb of gene body with gene expression matrixes by column into a large matrix
```
python snps_counts_comb_by_chr_no_filter.py snp.h5 index.h5 genotype_pca.txt {files}_peer_expr.txt {files}_snps_counts_comb.txt

```
5. Adjust p-values using beta approximation and further speed up using matrix multiplication
```
Rscript --vanilla beta_adjust_P_matrix_multip_corrected.R {files}_snps_counts_comb.txt {files}_output.txt

```
6. Merge output files
```
cat {files}_output.txt > all_chr_output.txt

```
7. extract out snp-gene pair that shows dominant effects
```
Rscript --vanilla extract_snp_gene_pair_comb_matrix.R {files}_output.txt {files}_snps_counts_comb.txt dominant_snp_gene_pair.txt

```


## Authors
Jing Gu

## Acknowledgements
* Dr. Graham McVicker Salk Institute for Biological Studies
* Patrick Fiaux UC San Diego
* Arko Sen Salk Institute for Biological Studies
* Hsiuyi Chen Salk Institute for Biological Studies
* Selene Tyndale Salk Institute for Biological Studies
