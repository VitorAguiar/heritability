#!/usr/bin/env bash

SAMPLES=/raid/genevol/heritability/samples.txt
BED=/raid/genevol/heritability/data/geuvadis_expression.bed
OUT=/raid/genevol/heritability/data/geuvadis_expression_std.bed

bgzip $BED && tabix -p bed $BED.gz

QTLtools correct --include-samples $SAMPLES --bed ${BED}.gz --normal --out $OUT

