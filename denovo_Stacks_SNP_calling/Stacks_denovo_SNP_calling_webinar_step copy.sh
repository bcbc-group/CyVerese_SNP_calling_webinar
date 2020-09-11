#This part of the tutorial starts with the assumption that you have already downloaded
#the sequencing files and have cleaned them

##########################################################################################
################ Read clean up, this is not as important for mapping to a reference #####
##########################################################################################

#move back up to the parent directory, the cd command may change depending on where you were last
cd SNP_calling_test
mkdir denovo_calling
mkdir de_novo_wrapper
cd denovo_calling

#run through the perl wrapper with default options
denovo_map.pl --samples ../cleaned_reads/ --popmap population_map.txt -o de_novo_wrapper/ -T 2

#export the raw SNP calls with no filtering within Stacks, subsequent filtering will occur in VCFtools
populations --batch_size 100 -P de_novo_wrapper/ -M population_map.txt -t 2 --vcf

#Do the filtering following SNP_filltering.sh
##############################################
#create a new directory for the final outputs
mkdir final_data_set

#read back in the filtered SNP file to export the necessary files for downstream analyses
populations  --batch_size 100 -V Ugibba_LD_pruned_SNPs.recode.vcf -O final_data_set/ -M population_map.txt -t 2 --fstats --vcf --treemix --plink --structure --radpainter --genepop
