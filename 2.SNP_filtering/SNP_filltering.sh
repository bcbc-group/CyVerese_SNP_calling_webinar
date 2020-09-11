#initial pass throwing out all SNPs that are missing 60%
vcftools --gzvcf Ugibba_initial_SNP_calls.vcf.gz --max-missing 0.5 --min-alleles 2 --max-alleles 2 --min-meanDP 3 --max-meanDP 30 --maf 0.05 --recode --recode-INFO-all --out Ugibba_first_pass

#gives missing proportion of loci for each individual
vcftools --vcf Ugibba_first_pass.recode.vcf --missing-indv

#average depth for each individual
vcftools --vcf Ugibba_first_pass.recode.vcf --depth 

#observed and expected heterozygosity
vcftools --vcf Ugibba_first_pass.recode.vcf --het

#create a list of individuals with at least 50% missing data
awk '$5 > 0.5' out.imiss | cut -f1 > lowDP50.indv

#redo filtering removing any accessions with more than 50% missing data
vcftools --vcf Ugibba_first_pass.recode.vcf --max-missing 0.5  --remove lowDP50.indv --recode --recode-INFO-all --out Ready_for_plink_Ugibba

#LD prune samples, to do this will need to use plink
plink --vcf Ready_for_plink_Ugibba.recode.vcf --make-bed --allow-extra-chr --double-id --out pop_sorted 

#tests for pairwise LD in a sliding window (10kb window with a step size of 5 and r^2 of 0.5), step size and length will depend on LD decay in your system, this is an ok start though
plink --bfile pop_sorted --indep-pairwise 10 5 0.5 --allow-extra-chr --double-id --threads 8 

#This will put the plinkn output in a format that VCFtools can easily understand
awk 'FNR==NR{a[$1];next}($2 in a){print}' plink.prune.out pop_sorted.bim > excluded.txt
awk '{print $1,$4}' excluded.txt > To_prune.txt

#vcftools to prune out selected SNPs from plink
vcftools --vcf Ready_for_plink_Ugibba.recode.vcf --exclude-positions To_prune.txt --max-missing 0.5 --min-meanDP 6 --max-meanDP 30 --maf 0.05 --recode --recode-INFO-all --out Ugibba_LD_pruned_SNPs

#gives missing proportion of loci for each individual
vcftools --vcf Ugibba_LD_pruned_SNPs.recode.vcf --missing-indv

#average depth for each individual
vcftools --vcf Ugibba_LD_pruned_SNPs.recode.vcf --depth 

#observed and expected heterozygosity
vcftools --vcf Ugibba_LD_pruned_SNPs.recode.vcf --het

