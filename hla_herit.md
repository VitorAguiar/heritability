Exploratory data analysis
================

## Sample:

  - 445 individuals from The 1000 Genomes Project (358 from 4 European
    populations + 87 Yorubans) with RNA-seq available from the Geuvadis
    consortium.

## 1\. VCF processing

### Subset and filtering

Here we select the GEUVADIS individuals, and exclude variants with
missing data or which are monomorphic in this subset of individuals.

Next we call an R script that converts VCF to GDS.

``` bash
#!/bin/bash

#PBS -l nodes=1:ppn=1
#PBS -l mem=8gb
#PBS -l walltime=24:00:00
#PBS -q short
#PBS -t 1-22
#PBS -N processVCF
#PBS -j oe
#PBS -o log/$PBS_JOBNAME

cd $PBS_O_WORKDIR

chr=$PBS_ARRAYID
samples=/raid/genevol/heritability/samples.txt
vcfin=/raid/genevol/heritability/genotypes_1000g/ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz
vcfout=./data/kgp/chr${chr}_subset.vcf

bcftools view --samples-file $samples --force-samples $vcfin |\
    bcftools view --genotype ^miss - |\
    bcftools view --min-ac=1:minor -o $vcfout -

Rscript vcf2gds.R $chr

rm $vcfout
```

The content of the vcf2gds.R script is the following:

``` r
library(SeqArray)

chr <- commandArgs(TRUE)[1]

vcf_file <- sprintf("./data/kgp/chr%s_subset.vcf", chr)
gds_file <- sprintf("./data/kgp/chr%s_tmp.gds", chr)

#convert
seqVCF2GDS(vcf_file, gds_file, fmt.import = "GT", verbose = FALSE)
```

## Analyses in R

### Packages:

``` r
library(tidyverse)
library(SeqArray) 
library(SeqVarTools) 
library(SNPRelate)
library(GENESIS)
library(Biobase)
library(cowplot)
library(ggbeeswarm)
```

### Concatenate GDS files for all chromosomes

``` r
gds_file <- "/home/vitor/heritability/data/kgp/allchrs.gds"

gds_list <- sprintf("/home/vitor/heritability/data/kgp/chr%d.gds", 1:22)

seqMerge(gds_list, gds_file)
```

### Pruning

``` r
# this step takes hours
gds <- seqOpen(gds_file)

set.seed(100)
pruned <- snpgdsLDpruning(gds, 
                          method = "corr", 
                          ld.threshold = sqrt(0.1))

prunedsnps <- unlist(pruned, use.names = FALSE)

seqSetFilter(gds, variant.id = prunedsnps)
seqExport(gds, "/home/vitor/heritability/data/kgp/allchrs_pruned.gds")
seqClose(gds)
```

### PCA on genotype data

``` r
# Close GDS and open the pruned one
gdsfmt::showfile.gds(closeall = TRUE)

gds <- seqOpen("/home/vitor/heritability/data/kgp/allchrs_pruned.gds")

pca <- snpgdsPCA(gds, num.thread = 16L)
```

    Principal Component Analysis (PCA) on genotypes:
    Calculating allele counts/frequencies ...
    # of selected variants: 2,195,733
        # of samples: 445
        # of SNVs: 2,195,733
        using 16 threads
        # of principal components: 32
    CPU capabilities: Double-Precision SSE2
    Thu Mar 11 19:40:03 2021    (internal increment: 8752)
    [..................................................]  0%, ETC: ---        [==================================================] 100%, completed, 1.3m
    Thu Mar 11 19:41:22 2021    Begin (eigenvalues and eigenvectors)
    Thu Mar 11 19:41:23 2021    Done.

``` r
pops <- read_tsv("/raid/genevol/heritability/hla_expression.tsv") %>%
    distinct(subject_id, pop) %>%
    rename(sampleid = subject_id)

pcadf <- as.data.frame(pca$eigenvect) %>%
    as_tibble() %>%
    add_column(sampleid = pca$sample.id, .before = 1) %>%
    inner_join(pops, ., by = "sampleid")
```

``` r
p1 <- ggplot(pcadf, aes(V1, V2, color = pop)) +
    geom_point() +
    scale_color_viridis_d() +
    theme_bw()

p2 <- ggplot(pcadf, aes(V2, V3, color = pop)) +
    geom_point() +
    scale_color_viridis_d() +
    theme_bw()

p3 <- ggplot(pcadf, aes(V3, V4, color = pop)) +
    geom_point() +
    scale_color_viridis_d() +
    theme_bw()

p_legend <- get_legend(p1)

plot_grid(p1 + theme(legend.position = "none"), 
          p_legend,
          p2 + theme(legend.position = "none"), 
          p3 + theme(legend.position = "none"),
          nrow = 2)
```

![](hla_herit_files/figure-gfm/pcaplot-1.png)<!-- -->

### Phenotypes (HLA-A expression in Geuvadis)

``` r
hla_expression <- read_tsv("/raid/genevol/heritability/hla_expression.tsv") %>%
    filter(gene_name == "HLA-A", subject_id %in% seqGetData(gds, "sample.id")) 
```

### Expression across populations

``` r
ggplot(hla_expression, aes(reorder(pop, tpm), tpm)) +
    geom_quasirandom(aes(color = pop), method = "smiley") +
    geom_boxplot(fill = NA, outlier.color = NA) +
    scale_color_viridis_d() +
    theme_bw() +
    labs(x = NULL, y = "TPM")
```

![](hla_herit_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### Expression across laboratories

``` r
ggplot(hla_expression, aes(reorder(lab, tpm), tpm)) +
    geom_quasirandom(aes(color = lab), method = "smiley") +
    geom_boxplot(fill = NA, outlier.color = NA) +
    theme_bw() +
    labs(x = NULL, y = "TPM")
```

![](hla_herit_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### Expression between sexes

``` r
ggplot(hla_expression, aes(reorder(sex, tpm), tpm)) +
    geom_quasirandom(aes(color = sex), method = "smiley") +
    geom_boxplot(fill = NA, outlier.color = NA) +
    theme_bw() +
    labs(x = NULL, y = "TPM")
```

![](hla_herit_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Here we got to evaluate whether such differences across groupings are
adequate for an analysis with the combined dataset, or whether they can
be accounted for by using covariates.

Before we continue, we need to transform our phenotype data.frame into
an “Annotated data frame”, which is a data structure provided by
`Biobase`. This means that we have to add metadata to our
`hla_expression` data.frame, and include a `sample.id` column as a
requirement for `GENESIS`. Since we do not have multiple samples for the
same individuals, we will just create the `sample.id` column as a copy
of `subject_id`.

``` r
metadata <- 
    data.frame(labelDescription = c("sample identifier",
                                    "subject identifier",
                                    "laboratory of RNA sequencing",
                                    "subject's sex",
                                    "subject's population",
                                    "PC 1",
                                    "PC 2",
                                    "PC 3",
                                    "expression levels in TPM"))

annotphen <- hla_expression %>%
    mutate(sample.id = subject_id) %>%
    left_join(select(pcadf, sampleid, V1:V3), 
              by = c("subject_id" = "sampleid")) %>%
    select(sample.id, subject_id, lab, sex, pop, V1, V2, V3, tpm) %>%
    as.data.frame() %>%
    AnnotatedDataFrame(metadata)

# access the metadata with the varMetadata() function
varMetadata(annotphen)
```

``` 
                       labelDescription
sample.id             sample identifier
subject_id           subject identifier
lab        laboratory of RNA sequencing
sex                       subject's sex
pop                subject's population
V1                                 PC 1
V2                                 PC 2
V3                                 PC 3
tpm            expression levels in TPM
```

``` r
# access the data with the pData() function
head(pData(annotphen))
```

``` 
  sample.id subject_id      lab    sex pop          V1            V2
1   HG00096    HG00096    UNIGE   male GBR -0.02409590  1.178854e-03
2   HG00097    HG00097     LUMC female GBR -0.02433438 -9.808199e-04
3   HG00099    HG00099     HMGU female GBR -0.02452840 -2.650025e-05
4   HG00100    HG00100 CNAG_CRG female GBR -0.02342908  1.718510e-02
5   HG00101    HG00101    UNIGE   male GBR -0.02409977  1.364934e-03
6   HG00102    HG00102    MPIMG female GBR -0.02474157 -2.891527e-03
            V3     tpm
1 -0.005208341 1604.34
2 -0.007818318 1342.20
3 -0.007551133  921.60
4  0.007413764 1392.66
5 -0.003502374 1082.02
6 -0.006582548 1053.15
```

## Fitting the NULL model

The first step in finding genetic variants which are associated with a
phenotype is preparing null model.

We fit a null model which adjusts gene expression according to
covariates such as population of origin, laboratory of sequencing, and
sex.

``` r
mod_0 <- fitNullModel(annotphen, 
                      outcome = "tpm", 
                      covars = c("pop", "sex", "lab", "V1", "V2", "V3"))
```

O output de fitNullModel tem vários elementos:

``` r
names(mod_0)
```

``` 
 [1] "family"         "hetResid"       "varComp"        "varCompCov"    
 [5] "fixef"          "betaCov"        "fitted.values"  "resid.marginal"
 [9] "logLik"         "AIC"            "workingY"       "outcome"       
[13] "model.matrix"   "group.idx"      "cholSigmaInv"   "converged"     
[17] "zeroFLAG"       "RSS"            "Ytilde"         "resid"         
[21] "CX"             "CXCXI"          "RSS0"           "sample.id"     
```

Confira se o modelo convergiu:

``` r
mod_0$converged
```

    [1] TRUE

Efeitos fixos:

``` r
mod_0$fixef
```

``` 
                Estimate  Std. Error     t value     Pr(>|t|)
(Intercept)  2718.849819   745.32292  3.64788165 2.968609e-04
popFIN       -263.170068   170.66040 -1.54206872 1.237925e-01
popGBR          4.099168    59.08472  0.06937779 9.447212e-01
popTSI         -2.599750   156.54411 -0.01660714 9.867577e-01
popYRI      -4105.038700  3719.85006 -1.10354951 2.704058e-01
sexmale       -54.378337    36.75004 -1.47968086 1.396907e-01
labHMGU      -196.782851    71.60707 -2.74809242 6.246527e-03
labICMB      -263.469927    62.29352 -4.22949201 2.863698e-05
labLUMC      -369.828292    72.52208 -5.09952683 5.115038e-07
labMPIMG     -129.932562    62.86306 -2.06691449 3.934019e-02
labUNIGE     -466.562714    54.14672 -8.61663917 1.326779e-16
labUU          51.871717    72.17204  0.71872314 4.727019e-01
V1          35175.801228 30951.90934  1.13646628 2.563944e-01
V2          -1966.705220  1996.11767 -0.98526517 3.250479e-01
V3           -119.568020   402.02886 -0.29741154 7.662960e-01
```

## Heritability

### GCTA GRM

We will use the `SNPRelate` package to compute a GRM. We will begin with
the GCTA method.

``` r
# Close GDS and open the pruned one
gdsfmt::showfile.gds(closeall = TRUE)
```

``` 
                                                            FileName ReadOnly
1 /raid/genevol/users/vitor/heritability/data/kgp/allchrs_pruned.gds     TRUE
   State
1 closed
```

``` r
pruned <- seqOpen("/home/vitor/heritability/data/kgp/allchrs_pruned.gds")

# Computar a GRM
grm_obj <- snpgdsGRM(pruned, method = "GCTA", num.thread = 16L)
```

    Genetic Relationship Matrix (GRM, GCTA):
    Calculating allele counts/frequencies ...
    # of selected variants: 2,195,733
        # of samples: 445
        # of SNVs: 2,195,733
        using 16 threads
    CPU capabilities: Double-Precision SSE2
    Thu Mar 11 19:41:32 2021    (internal increment: 8752)
    [..................................................]  0%, ETC: ---        [==================================================] 100%, completed, 1.3m
    Thu Mar 11 19:42:52 2021    Done.

We extract and rename the matrix

``` r
sample_ids <- seqGetData(pruned, "sample.id")

grm <- grm_obj$grm
rownames(grm) <- sample_ids
colnames(grm) <- sample_ids

# first 5 individuals:
grm[1:5, 1:5]
```

``` 
            HG00096      HG00097      HG00099     HG00100     HG00101
HG00096 0.956687356 0.0022034937 0.0023660411 0.001431433 0.003642691
HG00097 0.002203494 0.9754513066 0.0007980224 0.003025820 0.003069123
HG00099 0.002366041 0.0007980224 0.9609264716 0.009601814 0.005132164
HG00100 0.001431433 0.0030258202 0.0096018139 0.998632104 0.010114417
HG00101 0.003642691 0.0030691232 0.0051321635 0.010114417 0.938943766
```

Then, we build a new NULL model, now including the GRM.

``` r
mod_null <- fitNullModel(annotphen, 
                         outcome = "tpm", 
                         covars = c("lab", "sex", "pop", "V1", "V2", "V3"),
                         cov.mat = grm)
```

Now that we have a Null model adjusting expression levels for
population, sex, laboratory, population genetic structure, and the
relatedness, we can test for the association of the genetic variants
with expression levels.

The first step is to create a `SeqVarData` object including both the GDS
(genotypes) and the Annotated data.frame (phenotypes).

Then we will use the `assocTestSingle` function to assess the effect of
each variant.

``` r
# order individuals according to GDS
pData(annotphen) <- pData(annotphen) %>%
    mutate(subject_id = factor(subject_id, levels = seqGetData(pruned, "sample.id"))) %>%
    arrange(subject_id)

# create SeqVarData object and iterator
seqData <- SeqVarData(pruned, sampleData = annotphen)
iterator <- SeqVarBlockIterator(seqData, verbose = FALSE)

# test
assoc <- assocTestSingle(iterator, mod_null)

head(assoc)
```

`GENESIS` includes the `varCompCI` to compute the proportion of variance
explained (heritability) by each random effect, with a 95% CI:

``` r
varCompCI(nullmod, prop = TRUE)
```
