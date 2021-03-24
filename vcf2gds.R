library(SeqArray)

chr <- commandArgs(TRUE)[1]

vcf_file <- sprintf("/raid/genevol/heritability/data/kgp/chr%s_subset.vcf", chr)
gds_file <- sprintf("/raid/genevol/heritability/data/kgp/chr%s_tmp.gds", chr)

#convert
seqVCF2GDS(vcf_file, gds_file, fmt.import = "GT", verbose = FALSE)
