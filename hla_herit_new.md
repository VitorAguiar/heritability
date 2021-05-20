Heritability
================

``` r
# data processing
library(tidyverse)
library(tidytext)

# genetic data
library(SeqArray) 
library(SeqVarTools) 
library(SNPRelate)
library(GENESIS)
library(Biobase)

# plotting
library(cowplot)
library(ggbeeswarm)
library(ggsci)
library(GGally)
```

## Sample:

  - 445 individuals from The 1000 Genomes Project (358 from 4 European
    populations + 87 Yorubans) with RNA-seq available from the Geuvadis
    consortium.

## Data pre-processing

### Genotypes

We use VCF files provided by 1000 Genomes Phase III (low coverage, in
GRCh38 coordinates).

We select the GEUVADIS individuals, and exclude variants with missing
data or which are monomorphic in this subset of individuals.

Then we convert VCF to GDS, which is a data structure that optimizes
working with genotype data in R.

To perform all this processing, we run the script
`./launch_processVCF.sh`.

Then we concatenate the GDS files for each chromosome in a single GDS
with the code below:

``` r
gds_list <- sprintf("/raid/genevol/heritability/data/kgp/chr%d.gds", 1:22)

gds_file <- "/raid/genevol/heritability/data/kgp/allchrs.gds"

seqMerge(gds_list, gds_file)
```

We reduce the amount of variants by LD pruning:

``` r
# this step takes hours
gds <- seqOpen(gds_file)

set.seed(100)
pruned <- snpgdsLDpruning(gds, 
                          method = "corr", 
                          ld.threshold = sqrt(0.1))

prunedsnps <- unlist(pruned, use.names = FALSE)

seqSetFilter(gds, variant.id = prunedsnps)
seqExport(gds, "/raid/genevol/heritability/data/kgp/allchrs_pruned.gds")
seqClose(gds)
```

### Phenotypes

We use the FASTQ files provided by the Geuvadis Consortium, and estimate
expression with Salmon. We begin with the expression of HLA-A.

Transcript Per Million (TPM) values are quantile normal transformed by
`QTLtools correct`. Here is a comparison between raw and standardize
TPMs:

``` r
# raw
hla_expression <- 
    "/raid/genevol/heritability/data/geuvadis_expression.bed.gz" %>%
    read_tsv() %>%
    filter(id == "HLA-A") %>%
    pivot_longer(-(1:6), names_to = "sampleid") %>%
    select(sampleid, value)

# standardized by QTLtools
hla_expression_std <- 
    "/raid/genevol/heritability/data/geuvadis_expression_std.bed" %>%
    read_tsv() %>%
    filter(id == "HLA-A") %>%
    pivot_longer(-(1:6), names_to = "sampleid") %>%
    select(sampleid, value)
```

![](hla_herit_new_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

## Results

### PCA on genotype data

``` r
pops <- read_tsv("/raid/genevol/heritability/data/geuvadis_metadata.tsv") %>%
    distinct(sampleid, population)

# Close GDS and open the pruned one
gdsfmt::showfile.gds(closeall = TRUE)

gds <- seqOpen("/raid/genevol/heritability/data/kgp/allchrs_pruned.gds")

sample_ids <- seqGetData(gds, "sample.id")

pca <- snpgdsPCA(gds, sample.id = pops$sampleid, num.thread = 16L)
```

    Principal Component Analysis (PCA) on genotypes:
    Calculating allele counts/frequencies ...
    # of selected variants: 2,195,733
        # of samples: 445
        # of SNVs: 2,195,733
        using 16 threads
        # of principal components: 32
    CPU capabilities: Double-Precision SSE2
    Thu May 20 09:54:04 2021    (internal increment: 8752)
    [..................................................]  0%, ETC: ---        [==================================================] 100%, completed, 59s
    Thu May 20 09:55:03 2021    Begin (eigenvalues and eigenvectors)
    Thu May 20 09:55:04 2021    Done.

``` r
pcadf <- as.data.frame(pca$eigenvect) %>%
    as_tibble() %>%
    add_column(sampleid = pca$sample.id, .before = 1) %>%
    inner_join(pops, ., by = "sampleid")
```

![](hla_herit_new_files/figure-gfm/pcaplot-1.png)<!-- -->

### HLA-A expression across populations

``` r
sample_info <- "/raid/genevol/heritability/data/geuvadis_metadata.tsv" %>%
    read_tsv() %>%
    select(-ena_id)

expression_df <- left_join(hla_expression_std, sample_info, by = "sampleid") %>%
    select(sampleid, population, sex, lab, value)
```

![](hla_herit_new_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

### HLA-A expression across laboratories

![](hla_herit_new_files/figure-gfm/unnamed-chunk-7-1.png)<!-- -->

### HLA-A expression between males and females

![](hla_herit_new_files/figure-gfm/unnamed-chunk-8-1.png)<!-- -->

Here we got to evaluate whether such differences across groupings are
adequate for an analysis with the combined dataset, or whether they can
be accounted for by using covariates.

Before we continue, we need to transform our phenotype data.frame into
an “Annotated data frame”, which is a data structure provided by
`Biobase`. This means that we have to add metadata to our
`expression_df` data.frame, and include a `sample.id` column as a
requirement for `GENESIS`. Since we do not have multiple samples for the
same individuals, we will just create the `sample.id` column as a copy
of `sampleid`.

``` r
metadata <- 
    c("sample identifier",
      "subject identifier",
      "laboratory of RNA sequencing",
      "subject's sex",
      "subject's population",
      "PC 1",
      "PC 2",
      "PC 3",
      "expression levels in TPM") %>%
    data.frame(labelDescription = .) 

annotphen <- expression_df %>%
    mutate(sample.id = sampleid) %>%
    left_join(select(pcadf, sampleid, V1:V3), by = "sampleid") %>%
    select(sample.id, sampleid, lab, sex, population, V1, V2, V3, value) %>%
    as.data.frame() %>%
    AnnotatedDataFrame(metadata)

# access the metadata with the varMetadata() function
varMetadata(annotphen)
```

``` 
                       labelDescription
sample.id             sample identifier
sampleid             subject identifier
lab        laboratory of RNA sequencing
sex                       subject's sex
population         subject's population
V1                                 PC 1
V2                                 PC 2
V3                                 PC 3
value          expression levels in TPM
```

``` r
# access the data with the pData() function
head(pData(annotphen))
```

``` 
  sample.id sampleid      lab    sex population          V1            V2           V3      value
1   HG00096  HG00096    UNIGE   male        GBR -0.02409590  1.178854e-03 -0.005208341  0.0620013
2   HG00097  HG00097     LUMC female        GBR -0.02433438 -9.808199e-04 -0.007818318 -0.4860010
3   HG00099  HG00099     HMGU female        GBR -0.02452840 -2.650025e-05 -0.007551133 -1.9096400
4   HG00100  HG00100 CNAG_CRG female        GBR -0.02342908  1.718510e-02  0.007413764 -0.4235280
5   HG00101  HG00101    UNIGE   male        GBR -0.02409977  1.364934e-03 -0.003502374 -1.4375500
6   HG00102  HG00102    MPIMG female        GBR -0.02474157 -2.891527e-03 -0.006582548 -1.3484100
```

## GRM

### GCTA vs. Weir & Goudet

We will use the `SNPRelate` package to compute the GRMs.

In order to compare GRMs on different subsets of data, we will use the
same VCF which was processed to include all variable sites when taking
into account all 1KGP individuals (script
`./submit_processVCF_1KGP_all.pbs`).

``` r
# Close GDS and open the pruned one
gdsfmt::showfile.gds(closeall = TRUE)

pruned <- seqOpen("/raid/genevol/heritability/data/kgp/all_inds/allchrs_pruned.gds")

# Compute GCTA GRM with all individuals
grm_gcta_obj <- snpgdsGRM(pruned, method = "GCTA", num.thread = 16L)

# Compute Weir&Goudet GRM with all individuals
grm_wg_obj <- snpgdsGRM(pruned, method = "IndivBeta", num.thread = 16L)


# Compute GCTA GRM with Geuvadis subset
grm_gcta_geuv_obj <- snpgdsGRM(pruned, method = "GCTA", num.thread = 16L,
                               sample.id = pops$sampleid)

# Compute Weir&Goudet GRM with Geuvadis subset
grm_wg_geuv_obj <- snpgdsGRM(pruned, method = "IndivBeta", num.thread = 16L,
                             sample.id = pops$sampleid)


#Extract and rename the matrices
grm_gcta <- grm_gcta_obj$grm
rownames(grm_gcta) <- grm_gcta_obj$sample.id
colnames(grm_gcta) <- grm_gcta_obj$sample.id

grm_wg <- grm_wg_obj$grm
rownames(grm_wg) <- grm_wg_obj$sample.id
colnames(grm_wg) <- grm_wg_obj$sample.id

grm_gcta_geuv <- grm_gcta_geuv_obj$grm
rownames(grm_gcta_geuv) <- grm_gcta_geuv_obj$sample.id
colnames(grm_gcta_geuv) <- grm_gcta_geuv_obj$sample.id

grm_wg_geuv <- grm_wg_geuv_obj$grm
rownames(grm_wg_geuv) <- grm_wg_geuv_obj$sample.id
colnames(grm_wg_geuv) <- grm_wg_geuv_obj$sample.id
```

### Distribution of the GRM diagonal values

![](hla_herit_new_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

### GRM off-diagonal values

\* diagonal is set to NA so we can better see the range for off-diagonal
(much smaller) values

![](hla_herit_new_files/figure-gfm/unnamed-chunk-12-1.png)<!-- -->

### Pairs of individuals

![](hla_herit_new_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

### Fitting the Null model

The first step in finding genetic variants which are associated with a
phenotype is preparing null model.

We fit a null model which adjusts gene expression according to
covariates such as population of origin, laboratory of sequencing, and
sex. We also account for the relatedness between individuals (GRM).

``` r
mod_null_gcta <- fitNullModel(annotphen,
                        outcome = "value",
                        covars = c("lab", "sex"),
                        cov.mat = grm_gcta_geuv)
```

    [1]    0.4996713    0.4996713 -575.7362882    0.7661166
    [1]    0.65304972    0.04862065 -575.51386567    1.11813116
    [1]    0.5850034    0.1798392 -575.2525802    1.0130355
    [1]    0.5723877    0.2010636 -575.2480977    1.0002067
    [1]    0.5729969    0.2006428 -575.2480951    1.0000001

### Association testing

Now that we have a Null model adjusting expression levels for sex,
laboratory, population genetic structure, and relatedness, we can test
for the association of the genetic variants with expression levels.

The first step is to create a `SeqVarData` object including both the GDS
(genotypes) and the Annotated data.frame (phenotypes).

Then we will use the `assocTestSingle` function to assess the effect of
each variant.

``` r
# order individuals according to GDS
gdsfmt::showfile.gds(closeall = TRUE)

# it doesn't work if I use the GDS with all 1kgp individuals and filter for 
# geuvadis
pruned_geuv <- seqOpen("/raid/genevol/heritability/data/kgp/allchrs_pruned.gds")

pData(annotphen) <- pData(annotphen) %>%
    mutate(sampleid = factor(sampleid, levels = seqGetData(pruned_geuv, "sample.id"))) %>%
    arrange(sampleid)

# create SeqVarData object and iterator
seqData <- SeqVarData(pruned_geuv, sampleData = annotphen)
iterator <- SeqVarBlockIterator(seqData, variantBlock = 20000, verbose = FALSE)

# test
assoc <- assocTestSingle(iterator, mod_null_gcta)

head(assoc, 10)
```

``` 
   variant.id chr   pos allele.index n.obs        freq MAC        Score  Score.SE Score.Stat Score.pval         Est     Est.SE          PVE
1           1   1 10177            1   445 0.428089888 381 -19.89240959 12.866113 -1.5461087  0.1220783 -0.12016906 0.07772355 5.470142e-03
2           2   1 10352            1   445 0.439325843 391   2.98778594 10.226585  0.2921587  0.7701653  0.02856855 0.09778436 1.953243e-04
3           6   1 11012            1   445 0.092134831  82 -11.18542458  9.928933 -1.1265485  0.2599334 -0.11346119 0.10071576 2.904145e-03
4           8   1 13110            1   445 0.043820225  39 -10.97875301  6.993839 -1.5697749  0.1164675 -0.22445111 0.14298299 5.638886e-03
5          10   1 13118            1   445 0.143820225 128  -8.99989590 12.015513 -0.7490230  0.4538433 -0.06233800 0.08322574 1.283834e-03
6          11   1 13273            1   445 0.128089888 114  12.02161130 11.347490  1.0594071  0.2894144  0.09336048 0.08812522 2.568291e-03
7          14   1 13445            1   445 0.001123596   1   0.01991844  1.129253  0.0176386  0.9859272  0.01561971 0.88554118 7.119456e-07
8          16   1 13494            1   445 0.001123596   1   0.43793195  1.173812  0.3730851  0.7090851  0.31784047 0.85192484 3.185183e-04
9          20   1 14604            1   445 0.150561798 134 -12.86231647 11.469248 -1.1214611  0.2620916 -0.09777983 0.08718967 2.877975e-03
10         24   1 14933            1   445 0.041573034  37   4.53096304  6.926960  0.6541055  0.5130438  0.09442894 0.14436347 9.790710e-04
```

### Heritability

`GENESIS` includes the `varCompCI` to compute point estimates and
confidence intervals for each of the random effects variance component
estimates.

``` r
varCompCI(mod_null_gcta, prop = TRUE)
```

``` 
            Proportion   Lower 95  Upper 95
V_A          0.7406509  0.2295703 1.2517315
V_resid.var  0.2593491 -0.2517315 0.7704297
```

## Session Information

    Error in get(genname, envir = envir) : object 'testthat_print' not found

    ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     setting  value                       
     version  R version 4.0.3 (2020-10-10)
     os       Ubuntu 16.04.7 LTS          
     system   x86_64, linux-gnu           
     ui       X11                         
     language (EN)                        
     collate  en_US.UTF-8                 
     ctype    en_US.UTF-8                 
     tz       America/Sao_Paulo           
     date     2021-05-20                  
    
    ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     package          * version  date       lib source        
     assertthat         0.2.1    2019-03-21 [2] CRAN (R 4.0.2)
     backports          1.2.1    2020-12-09 [1] CRAN (R 4.0.2)
     beeswarm           0.2.3    2016-04-25 [1] CRAN (R 4.0.2)
     Biobase          * 2.50.0   2020-10-27 [1] Bioconductor  
     BiocGenerics     * 0.36.0   2020-10-27 [1] Bioconductor  
     Biostrings         2.58.0   2020-10-27 [1] Bioconductor  
     bit                4.0.4    2020-08-04 [1] CRAN (R 4.0.2)
     bit64              4.0.2    2020-07-30 [2] CRAN (R 4.0.2)
     bitops             1.0-6    2013-08-17 [2] CRAN (R 4.0.2)
     blob               1.2.1    2020-01-20 [2] CRAN (R 4.0.2)
     broom              0.7.5    2021-02-19 [1] CRAN (R 4.0.2)
     callr              3.5.1    2020-10-13 [1] CRAN (R 4.0.2)
     cellranger         1.1.0    2016-07-27 [2] CRAN (R 4.0.2)
     cli                2.3.0    2021-01-31 [1] CRAN (R 4.0.2)
     codetools          0.2-16   2018-12-24 [4] CRAN (R 4.0.0)
     colorspace         2.0-0    2020-11-11 [1] CRAN (R 4.0.2)
     conquer            1.0.1    2020-05-06 [2] CRAN (R 4.0.2)
     cowplot          * 1.1.1    2020-12-30 [1] CRAN (R 4.0.2)
     crayon             1.4.1    2021-02-08 [1] CRAN (R 4.0.2)
     data.table         1.13.6   2020-12-30 [1] CRAN (R 4.0.2)
     DBI                1.1.0    2019-12-15 [2] CRAN (R 4.0.2)
     dbplyr             1.4.4    2020-05-27 [2] CRAN (R 4.0.2)
     desc               1.2.0    2018-05-01 [2] CRAN (R 4.0.2)
     devtools           2.3.2    2020-09-18 [1] CRAN (R 4.0.2)
     digest             0.6.27   2020-10-24 [1] CRAN (R 4.0.2)
     DNAcopy            1.64.0   2020-10-27 [1] Bioconductor  
     dplyr            * 1.0.4    2021-02-02 [1] CRAN (R 4.0.2)
     ellipsis           0.3.1    2020-05-15 [1] CRAN (R 4.0.2)
     evaluate           0.14     2019-05-28 [2] CRAN (R 4.0.2)
     farver             2.0.3    2020-01-16 [2] CRAN (R 4.0.2)
     forcats          * 0.5.0    2020-03-01 [2] CRAN (R 4.0.2)
     foreach            1.5.0    2020-03-30 [2] CRAN (R 4.0.2)
     formula.tools      1.7.1    2018-03-01 [1] CRAN (R 4.0.2)
     fs                 1.5.0    2020-07-31 [1] CRAN (R 4.0.2)
     gdsfmt           * 1.26.1   2020-12-22 [1] Bioconductor  
     generics           0.1.0    2020-10-31 [1] CRAN (R 4.0.2)
     GENESIS          * 2.20.1   2021-01-28 [1] Bioconductor  
     GenomeInfoDb       1.26.2   2020-12-08 [1] Bioconductor  
     GenomeInfoDbData   1.2.4    2021-02-20 [1] Bioconductor  
     GenomicRanges      1.42.0   2020-10-27 [1] Bioconductor  
     GGally           * 2.1.0    2021-01-06 [1] CRAN (R 4.0.2)
     ggbeeswarm       * 0.6.0    2017-08-07 [1] CRAN (R 4.0.2)
     ggplot2          * 3.3.2    2020-06-19 [2] CRAN (R 4.0.2)
     ggsci            * 2.9      2018-05-14 [1] CRAN (R 4.0.2)
     glue               1.4.2    2020-08-27 [1] CRAN (R 4.0.2)
     gtable             0.3.0    2019-03-25 [2] CRAN (R 4.0.2)
     GWASExactHW        1.01     2013-01-05 [1] CRAN (R 4.0.2)
     GWASTools          1.36.0   2020-10-27 [1] Bioconductor  
     haven              2.3.1    2020-06-01 [2] CRAN (R 4.0.2)
     highr              0.8      2019-03-20 [2] CRAN (R 4.0.2)
     hms                0.5.3    2020-01-08 [2] CRAN (R 4.0.2)
     htmltools          0.5.1.1  2021-01-22 [1] CRAN (R 4.0.2)
     httr               1.4.2    2020-07-20 [2] CRAN (R 4.0.2)
     IRanges            2.24.1   2020-12-12 [1] Bioconductor  
     iterators          1.0.12   2019-07-26 [2] CRAN (R 4.0.2)
     janeaustenr        0.1.5    2017-06-10 [1] CRAN (R 4.0.2)
     jsonlite           1.7.2    2020-12-09 [1] CRAN (R 4.0.2)
     knitr              1.31     2021-01-27 [1] CRAN (R 4.0.2)
     labeling           0.4.2    2020-10-20 [1] CRAN (R 4.0.2)
     lattice            0.20-41  2020-04-02 [4] CRAN (R 4.0.0)
     lifecycle          1.0.0    2021-02-15 [1] CRAN (R 4.0.2)
     lmtest             0.9-38   2020-09-09 [1] CRAN (R 4.0.2)
     logistf            1.24     2020-09-16 [1] CRAN (R 4.0.2)
     lubridate          1.7.9.2  2020-11-13 [1] CRAN (R 4.0.2)
     magrittr           2.0.1    2020-11-17 [1] CRAN (R 4.0.2)
     Matrix             1.2-18   2019-11-27 [2] CRAN (R 4.0.2)
     MatrixModels       0.4-1    2015-08-22 [2] CRAN (R 4.0.2)
     matrixStats        0.58.0   2021-01-29 [1] CRAN (R 4.0.2)
     memoise            1.1.0    2017-04-21 [2] CRAN (R 4.0.2)
     mgcv               1.8-33   2020-08-27 [4] CRAN (R 4.0.2)
     mice               3.13.0   2021-01-27 [1] CRAN (R 4.0.2)
     modelr             0.1.8    2020-05-19 [1] CRAN (R 4.0.2)
     munsell            0.5.0    2018-06-12 [2] CRAN (R 4.0.2)
     nlme               3.1-149  2020-08-23 [4] CRAN (R 4.0.2)
     operator.tools     1.6.3    2017-02-28 [1] CRAN (R 4.0.2)
     pillar             1.4.7    2020-11-20 [1] CRAN (R 4.0.2)
     pkgbuild           1.1.0    2020-07-13 [2] CRAN (R 4.0.2)
     pkgconfig          2.0.3    2019-09-22 [2] CRAN (R 4.0.2)
     pkgload            1.1.0    2020-05-29 [1] CRAN (R 4.0.2)
     plyr               1.8.6    2020-03-03 [2] CRAN (R 4.0.2)
     prettyunits        1.1.1    2020-01-24 [2] CRAN (R 4.0.2)
     processx           3.4.5    2020-11-30 [1] CRAN (R 4.0.2)
     ps                 1.5.0    2020-12-05 [1] CRAN (R 4.0.2)
     purrr            * 0.3.4    2020-04-17 [2] CRAN (R 4.0.2)
     quantreg           5.61     2020-07-09 [2] CRAN (R 4.0.2)
     quantsmooth        1.56.0   2020-10-27 [1] Bioconductor  
     R6                 2.5.0    2020-10-28 [1] CRAN (R 4.0.2)
     RColorBrewer       1.1-2    2014-12-07 [2] CRAN (R 4.0.2)
     Rcpp               1.0.6    2021-01-15 [1] CRAN (R 4.0.2)
     RCurl              1.98-1.2 2020-04-18 [2] CRAN (R 4.0.2)
     readr            * 1.4.0    2020-10-05 [1] CRAN (R 4.0.2)
     readxl             1.3.1    2019-03-13 [1] CRAN (R 4.0.2)
     remotes            2.2.0    2020-07-21 [2] CRAN (R 4.0.2)
     reprex             1.0.0    2021-01-27 [1] CRAN (R 4.0.2)
     reshape            0.8.8    2018-10-23 [1] CRAN (R 4.0.2)
     rlang              0.4.10   2020-12-30 [1] CRAN (R 4.0.2)
     rmarkdown          2.7      2021-02-19 [1] CRAN (R 4.0.2)
     rprojroot          1.3-2    2018-01-03 [2] CRAN (R 4.0.2)
     RSQLite            2.2.0    2020-01-07 [2] CRAN (R 4.0.2)
     rstudioapi         0.13     2020-11-12 [1] CRAN (R 4.0.2)
     rvest              0.3.6    2020-07-25 [2] CRAN (R 4.0.2)
     S4Vectors          0.28.1   2020-12-09 [1] Bioconductor  
     sandwich           3.0-0    2020-10-02 [1] CRAN (R 4.0.2)
     scales             1.1.1    2020-05-11 [2] CRAN (R 4.0.2)
     SeqArray         * 1.30.0   2020-10-27 [1] Bioconductor  
     SeqVarTools      * 1.28.1   2020-11-20 [1] Bioconductor  
     sessioninfo        1.1.1    2018-11-05 [2] CRAN (R 4.0.2)
     SnowballC          0.7.0    2020-04-01 [1] CRAN (R 4.0.2)
     SNPRelate        * 1.24.0   2020-10-27 [1] Bioconductor  
     SparseM            1.78     2019-12-13 [2] CRAN (R 4.0.2)
     stringi            1.5.3    2020-09-09 [1] CRAN (R 4.0.2)
     stringr          * 1.4.0    2019-02-10 [2] CRAN (R 4.0.2)
     survival           3.2-7    2020-09-28 [4] CRAN (R 4.0.2)
     testthat           2.3.2    2020-03-02 [2] CRAN (R 4.0.2)
     tibble           * 3.0.6    2021-01-29 [1] CRAN (R 4.0.2)
     tidyr            * 1.1.2    2020-08-27 [1] CRAN (R 4.0.2)
     tidyselect         1.1.0    2020-05-11 [2] CRAN (R 4.0.2)
     tidytext         * 0.3.0    2021-01-06 [1] CRAN (R 4.0.2)
     tidyverse        * 1.3.0    2019-11-21 [1] CRAN (R 4.0.2)
     tokenizers         0.2.1    2018-03-29 [1] CRAN (R 4.0.2)
     usethis            2.0.1    2021-02-10 [1] CRAN (R 4.0.2)
     vctrs              0.3.6    2020-12-17 [1] CRAN (R 4.0.2)
     vipor              0.4.5    2017-03-22 [1] CRAN (R 4.0.2)
     withr              2.4.1    2021-01-26 [1] CRAN (R 4.0.2)
     xfun               0.21     2021-02-10 [1] CRAN (R 4.0.2)
     xml2               1.3.2    2020-04-23 [2] CRAN (R 4.0.2)
     XVector            0.30.0   2020-10-27 [1] Bioconductor  
     yaml               2.2.1    2020-02-01 [2] CRAN (R 4.0.2)
     zlibbioc           1.36.0   2020-10-27 [1] Bioconductor  
     zoo                1.8-8    2020-05-02 [2] CRAN (R 4.0.2)
    
    [1] /raid/genevol/users/vitor/R/x86_64-pc-linux-gnu-library/4.0
    [2] /usr/local/lib/R/site-library
    [3] /usr/lib/R/site-library
    [4] /usr/lib/R/library
