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
