library(tidyverse)

# Gene annotations form Gencode
annot <- "/home/vitor/gencode/gencode.v37.primary_assembly.annotation.gtf" %>% 
    read_tsv(comment = "#", col_names = FALSE, col_types = "c-cii-c-c")

gene_annot <- annot %>%
    filter(X3 == "gene") %>%
    transmute(chr = X1,
	      gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	      gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	      start = X4, end = X5, strand = X7)

transc_annot <- annot %>%
    filter(X3 == "transcript") %>%
    transmute(gene_id = str_extract(X9, "(?<=gene_id\\s\")[^\"]+"),
	      gene_name = str_extract(X9, "(?<=gene_name\\s\")[^\"]+"),
	      tx_id = str_extract(X9, "(?<=transcript_id\\s\")[^\"]+"))


# 1000 Genomes sample annotation
sample_ids <- readLines("/raid/genevol/geuvadis/salmon/samples.txt")

sample_info_1kgp <- "/raid/genevol/heritability/sample_annotation_1000G.xlsx" 
sample_info_geuvadis <- "/raid/genevol/heritability/sample_annotation_geuvadis.txt"

#"ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx" %>%
#    download.file(destfile = sample_info_1kgp)
#
#"http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt" %>%
#    download.file(destfile = sample_info_geuvadis)

info_1kgp_df <- readxl::read_excel(sample_info_1kgp, col_types = "text") %>%
    select(sampleid = Sample, population = Population, sex = Gender)

info_geuvadis_df <- read_tsv(sample_info_geuvadis) %>%
    select(sampleid = `Source Name`,
	   ena_id = `Comment[ENA_RUN]`,
	   lab = Performer) %>%
    distinct()

sample_info_df <- inner_join(info_1kgp_df, info_geuvadis_df, by = "sampleid") %>%
    filter(ena_id %in% sample_ids) %>%
    select(sampleid, ena_id, everything())


# Gene expression
library(furrr)

plan(multisession, workers = 10)

quants <- sprintf("/raid/genevol/geuvadis/salmon/quant/%s/quant.sf", sample_info_df$ena_id) %>%
    setNames(sample_info_df$sampleid) %>%
    future_map_dfr(~read_tsv(.) %>%
		   left_join(transc_annot, by = c("Name" = "tx_id")) %>%
		   group_by(gene_id) %>%
		   summarize(tpm = sum(TPM)) %>%
		   ungroup(), 
	           .id = "sampleid")

bed <- quants %>%
    pivot_wider(names_from = sampleid, values_from = tpm) %>%
    left_join(gene_annot, by = "gene_id") %>%
    mutate(chr = factor(chr, levels = unique(annot$X1))) %>%
    select(`#chr` = chr, start, end, id = gene_name, gid = gene_id, strd = strand, everything()) %>%
    arrange(`#chr`, start)

write_tsv(bed, "/raid/genevol/heritability/data/geuvadis_expression.bed")
write_tsv(sample_info_df, "/raid/genevol/heritability/data/geuvadis_metadata.tsv")
