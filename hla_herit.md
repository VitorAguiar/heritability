Heritability
================

``` r
library(tidyverse)
library(SeqArray) 
library(SeqVarTools) 
library(SNPRelate)
library(GENESIS)
library(Biobase)
library(cowplot)
library(ggbeeswarm)
library(ggthemes)
```

## Sample:

-   445 individuals from The 1000 Genomes Project (358 from 4 European populations + 87 Yorubans) with RNA-seq available from the Geuvadis consortium.

## Data pre-processing

### Genotypes

We use VCF files provided by 1000 Genomes Phase III (low coverage, in GRCh38 coordinates).

We select the GEUVADIS individuals, and exclude variants with missing data or which are monomorphic in this subset of individuals.

Then we convert VCF to GDS, which is a data structure that optimizes working with genotype data in R.

To perform all this processing, we run the script `./launch_processVCF.sh`.

Then we concatenate the GDS files for each chromosome in a single GDS with the code below:

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

We use the FASTQ files provided by the Geuvadis Consortium, and estimate expression with Salmon. We begin with the expression of HLA-A.

Transcript Per Million (TPM) values are quantile normal transformed by `QTLtools correct`. Here is a comparison between raw and standardize TPMs:

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

![](hla_herit_files/figure-markdown_github/unnamed-chunk-3-1.png)

## Results

### PCA on genotype data

``` r
# Close GDS and open the pruned one
gdsfmt::showfile.gds(closeall = TRUE)

gds <- seqOpen("/raid/genevol/heritability/data/kgp/allchrs_pruned.gds")

sample_ids <- seqGetData(gds, "sample.id")

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
    Wed Mar 24 10:32:14 2021    (internal increment: 8752)

    [..................................................]  0%, ETC: ---        
    [==================================================] 100%, completed, 43s
    Wed Mar 24 10:32:57 2021    Begin (eigenvalues and eigenvectors)
    Wed Mar 24 10:32:57 2021    Done.

``` r
pops <- read_tsv("/raid/genevol/heritability/data/geuvadis_metadata.tsv") %>%
    distinct(sampleid, population)

pcadf <- as.data.frame(pca$eigenvect) %>%
    as_tibble() %>%
    add_column(sampleid = pca$sample.id, .before = 1) %>%
    inner_join(pops, ., by = "sampleid")
```

![](hla_herit_files/figure-markdown_github/pcaplot-1.png)

### HLA-A expression across populations

``` r
sample_info <- "/raid/genevol/heritability/data/geuvadis_metadata.tsv" %>%
    read_tsv() %>%
    select(-ena_id)

expression_df <- left_join(hla_expression_std, sample_info, by = "sampleid") %>%
    select(sampleid, population, sex, lab, value)

ggplot(expression_df, aes(reorder(population, value), value)) +
    geom_quasirandom(aes(color = population), method = "smiley") +
    geom_boxplot(fill = NA, outlier.color = NA) +
    scale_color_colorblind() +
    theme_bw() +
    labs(x = NULL, y = "Std TPM")
```

![](hla_herit_files/figure-markdown_github/unnamed-chunk-5-1.png)

### HLA-A expression across laboratories

``` r
ggplot(expression_df, aes(reorder(lab, value), value)) +
    geom_quasirandom(aes(color = lab), method = "smiley") +
    geom_boxplot(fill = NA, outlier.color = NA) +
    theme_bw() +
    labs(x = NULL, y = "Std TPM")
```

![](hla_herit_files/figure-markdown_github/unnamed-chunk-6-1.png)

### HLA-A expression between males and females

``` r
ggplot(expression_df, aes(reorder(sex, value), value)) +
    geom_quasirandom(aes(color = sex), method = "smiley") +
    geom_boxplot(fill = NA, outlier.color = NA) +
    theme_bw() +
    labs(x = NULL, y = "Std TPM")
```

![](hla_herit_files/figure-markdown_github/unnamed-chunk-7-1.png)

Here we got to evaluate whether such differences across groupings are adequate for an analysis with the combined dataset, or whether they can be accounted for by using covariates.

Before we continue, we need to transform our phenotype data.frame into an "Annotated data frame", which is a data structure provided by `Biobase`. This means that we have to add metadata to our `expression_df` data.frame, and include a `sample.id` column as a requirement for `GENESIS`. Since we do not have multiple samples for the same individuals, we will just create the `sample.id` column as a copy of `sampleid`.

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

``` r
# access the data with the pData() function
head(pData(annotphen))
```

      sample.id sampleid      lab    sex population          V1            V2           V3      value
    1   HG00096  HG00096    UNIGE   male        GBR -0.02409590  1.178854e-03 -0.005208341  0.0620013
    2   HG00097  HG00097     LUMC female        GBR -0.02433438 -9.808199e-04 -0.007818318 -0.4860010
    3   HG00099  HG00099     HMGU female        GBR -0.02452840 -2.650025e-05 -0.007551133 -1.9096400
    4   HG00100  HG00100 CNAG_CRG female        GBR -0.02342908  1.718510e-02  0.007413764 -0.4235280
    5   HG00101  HG00101    UNIGE   male        GBR -0.02409977  1.364934e-03 -0.003502374 -1.4375500
    6   HG00102  HG00102    MPIMG female        GBR -0.02474157 -2.891527e-03 -0.006582548 -1.3484100

## GRM

### GCTA

We will use the `SNPRelate` package to compute a GRM. We will begin with the GCTA method.

``` r
# Close GDS and open the pruned one
gdsfmt::showfile.gds(closeall = TRUE)
```

                                                    FileName ReadOnly  State
    1 /raid/genevol/heritability/data/kgp/allchrs_pruned.gds     TRUE closed

``` r
pruned <- seqOpen("/raid/genevol/heritability/data/kgp/allchrs_pruned.gds")

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
    Wed Mar 24 10:33:03 2021    (internal increment: 8752)

    [..................................................]  0%, ETC: ---        
    [==================================================] 100%, completed, 49s
    Wed Mar 24 10:33:52 2021    Done.

We extract and rename the matrix

``` r
grm <- grm_obj$grm
rownames(grm) <- sample_ids
colnames(grm) <- sample_ids

# first 5 individuals:
grm[1:5, 1:5]
```

                HG00096      HG00097      HG00099     HG00100     HG00101
    HG00096 0.956687356 0.0022034937 0.0023660411 0.001431433 0.003642691
    HG00097 0.002203494 0.9754513066 0.0007980224 0.003025820 0.003069123
    HG00099 0.002366041 0.0007980224 0.9609264716 0.009601814 0.005132164
    HG00100 0.001431433 0.0030258202 0.0096018139 0.998632104 0.010114417
    HG00101 0.003642691 0.0030691232 0.0051321635 0.010114417 0.938943766

### Distribution of the GRM diagonal values

![](hla_herit_files/figure-markdown_github/unnamed-chunk-11-1.png)

### GRM off-diagonal values

\* diagonal is set to NA so we can better see the range for off-diagonal (much smaller) values

![](hla_herit_files/figure-markdown_github/unnamed-chunk-12-1.png)

### Fitting the Null model

The first step in finding genetic variants which are associated with a phenotype is preparing null model.

We fit a null model which adjusts gene expression according to covariates such as population of origin, laboratory of sequencing, and sex. We also account for the relatedness between individuals (GRM).

``` r
mod_null <- fitNullModel(annotphen, 
                         outcome = "value", 
                         covars = c("lab", "sex", "V1", "V2"),
                         cov.mat = grm)
```

    [1]    0.4996713    0.4996713 -572.3989493    0.7609524
    [1]    0.7946567    0.1283191 -571.4559448    0.8274359
    [1]    0.85256910    0.04770545 -571.29695132    0.84968960
    [1]    0.87353098    0.01740093 -571.24256294    0.85920572
    [1]  8.826167e-01  4.040699e-03 -5.712196e+02  8.636316e-01
    [1]  8.847395e-01  8.939585e-04 -5.712143e+02  8.646985e-01
    [1]  8.852615e-01  1.185643e-04 -5.712130e+02  8.649628e-01

\*Note: if we include PC3, we get an error "the leading minor of order 445 is not positive definite"

### Association testing

Now that we have a Null model adjusting expression levels for sex, laboratory, population genetic structure, and relatedness, we can test for the association of the genetic variants with expression levels.

The first step is to create a `SeqVarData` object including both the GDS (genotypes) and the Annotated data.frame (phenotypes).

Then we will use the `assocTestSingle` function to assess the effect of each variant.

``` r
# order individuals according to GDS
pData(annotphen) <- pData(annotphen) %>%
    mutate(sampleid = factor(sampleid, levels = seqGetData(pruned, "sample.id"))) %>%
    arrange(sampleid)

# create SeqVarData object and iterator
seqData <- SeqVarData(pruned, sampleData = annotphen)
iterator <- SeqVarBlockIterator(seqData, variantBlock = 20000, verbose = FALSE)

# test
assoc <- assocTestSingle(iterator, mod_null)
```

    # of selected samples: 445

``` r
head(assoc, 10)
```

         variant.id chr     pos allele.index n.obs        freq MAC         Score   Score.SE    Score.Stat   Score.pval           Est     Est.SE          PVE
    1             1   1   10177            1   445 0.428089888 381 -1.927113e+01 11.1780661 -1.7240132578 0.0847054498 -0.1542318006 0.08946091 7.899407e-03
    2             2   1   10352            1   445 0.439325843 391  3.175980e+00  9.2765292  0.3423672198 0.7320745543  0.0369068229 0.10779894 3.115284e-04
    3             6   1   11012            1   445 0.092134831  82 -1.067291e+01  9.2196640 -1.1576241195 0.2470174620 -0.1255603375 0.10846382 3.561627e-03
    4             8   1   13110            1   445 0.043820225  39 -8.597626e+00  6.3394456 -1.3562110370 0.1750320313 -0.2139321188 0.15774250 4.888413e-03
    5            10   1   13118            1   445 0.143820225 128 -5.262152e+00 10.8740346 -0.4839189645 0.6284433807 -0.0445022460 0.09196219 6.223842e-04
    6            11   1   13273            1   445 0.128089888 114  1.059346e+01 10.1486070  1.0438336959 0.2965623388  0.1028548746 0.09853569 2.895849e-03
    7            14   1   13445            1   445 0.001123596   1  2.797633e-04  1.0215946  0.0002738496 0.9997814996  0.0002680610 0.97886190 1.993139e-10
    8            16   1   13494            1   445 0.001123596   1  3.714169e-01  1.0479487  0.3544227924 0.7230220559  0.3382062527 0.95424521 3.338540e-04
    9            20   1   14604            1   445 0.150561798 134 -9.632020e+00 10.5824489 -0.9101881917 0.3627232705 -0.0860092218 0.09449609 2.201789e-03
    10           24   1   14933            1   445 0.041573034  37  4.380638e+00  6.2722052  0.6984207150 0.4849141253  0.1113517003 0.15943356 1.296425e-03

### Heritability

`GENESIS` includes the `varCompCI` to compute the proportion of variance explained (heritability) by each random effect, with a 95% CI:

``` r
varCompCI(mod_null, prop = TRUE)
```

                  Proportion    Lower 95 Upper 95
    V_A         0.9998660865 -0.07450015 2.074232
    V_resid.var 0.0001339135 -1.07423232 1.074500

## Package versions

    Error in get(genname, envir = envir) : object 'testthat_print' not found

    ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
     setting  value                       
     version  R version 4.0.2 (2020-06-22)
     os       Ubuntu 16.04.7 LTS          
     system   x86_64, linux-gnu           
     ui       X11                         
     language en_US:en                    
     collate  en_US.UTF-8                 
     ctype    en_US.UTF-8                 
     tz       America/Sao_Paulo           
     date     2021-03-24                  

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
     ggbeeswarm       * 0.6.0    2017-08-07 [1] CRAN (R 4.0.2)
     ggplot2          * 3.3.2    2020-06-19 [2] CRAN (R 4.0.2)
     ggthemes         * 4.2.4    2021-01-20 [1] CRAN (R 4.0.2)
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
     prettyunits        1.1.1    2020-01-24 [2] CRAN (R 4.0.2)
     processx           3.4.5    2020-11-30 [1] CRAN (R 4.0.2)
     ps                 1.5.0    2020-12-05 [1] CRAN (R 4.0.2)
     purrr            * 0.3.4    2020-04-17 [2] CRAN (R 4.0.2)
     quantreg           5.61     2020-07-09 [2] CRAN (R 4.0.2)
     quantsmooth        1.56.0   2020-10-27 [1] Bioconductor  
     R6                 2.5.0    2020-10-28 [1] CRAN (R 4.0.2)
     Rcpp               1.0.6    2021-01-15 [1] CRAN (R 4.0.2)
     RCurl              1.98-1.2 2020-04-18 [2] CRAN (R 4.0.2)
     readr            * 1.4.0    2020-10-05 [1] CRAN (R 4.0.2)
     readxl             1.3.1    2019-03-13 [1] CRAN (R 4.0.2)
     remotes            2.2.0    2020-07-21 [2] CRAN (R 4.0.2)
     reprex             1.0.0    2021-01-27 [1] CRAN (R 4.0.2)
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
     SNPRelate        * 1.24.0   2020-10-27 [1] Bioconductor  
     SparseM            1.78     2019-12-13 [2] CRAN (R 4.0.2)
     stringi            1.5.3    2020-09-09 [1] CRAN (R 4.0.2)
     stringr          * 1.4.0    2019-02-10 [2] CRAN (R 4.0.2)
     survival           3.2-7    2020-09-28 [4] CRAN (R 4.0.2)
     testthat           2.3.2    2020-03-02 [2] CRAN (R 4.0.2)
     tibble           * 3.0.6    2021-01-29 [1] CRAN (R 4.0.2)
     tidyr            * 1.1.2    2020-08-27 [1] CRAN (R 4.0.2)
     tidyselect         1.1.0    2020-05-11 [2] CRAN (R 4.0.2)
     tidyverse        * 1.3.0    2019-11-21 [1] CRAN (R 4.0.2)
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
