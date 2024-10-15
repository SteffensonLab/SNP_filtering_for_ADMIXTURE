#!/bin/bash

#(1) convert genotype calls with low reads depth to missing

perl filter_convert_depth.pl 0.2 sample.info input.snp.vcf.gz convert.snp.vcf.gz 

#(2) SNP filtering of (1) maximum missing 0.2 (2) maximum 20% heterozygous calls 

perl filter_qua_miss_hete_nopoly.pl 0 0.2 0.2 convert.snp.vcf.gz filter.snp.vcf.gz

#(3) maf>=0.05

vcftools --gzvcf filter.snp.vcf.gz --maf 0.05 --recode --stdout | gzip -c > maf.snp.vcf.gz

#(4) to plink bed format

plink2 --vcf maf.snp.vcf.gz --threads 5 --make-bed --out all

#(5) LD prune

plink2 --bfile all --indep-pairwise 50 10 0.1
plink2 --bfile all --extract plink2.prune.in --make-bed --out LD

echo finish
