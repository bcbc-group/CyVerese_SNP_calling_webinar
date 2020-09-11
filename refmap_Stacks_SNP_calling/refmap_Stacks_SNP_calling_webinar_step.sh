#Download the raw data that is needed for SNP calling. To do this we will use wget and 
#the FTP link from the European Nucleotide Archive. This an easy way to download fastq
#files via the command line. After we download each file, we will essentially rename it
#to something that has a biologically relevant meaning for us.

#first we want to create a new folder, for now on the Desktop, we will create a folder structure
#with SNP_calling_test being the main folder, and subfolder for now being Raw_files

#open up a terminal to navigate to your Desktop
cd Desktop
mkdir SNP_calling_test
cd SNP_calling_test
mkdir ../Reference_genome
mkdir Raw_files


############################
### Reference genome #######
############################

#first we are going to download the reference assembly and change the name
cd ../Reference_genome
wget ftp://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/ne/NEEC01.fasta.gz
gunzip NEEC01.fasta.gz
mv NEEC01.fasta Ugibba_reference.fasta

#while we are hear, we are going to prep this file for future use
#first create a dictionary file from the reference
gatk CreateSequenceDictionary -R Ugibba_reference.fasta -O Ugibba_reference.dict

#create an index file for mapping
samtools faidx Ugibba_reference.fasta

###############################
### Individual samples  #######
###############################

#move back up one directory to move to the raw data file directory
cd ..
cd Raw_files

#bladder rep 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/046/SRR10676746/SRR10676746.fastq.gz
mv SRR10676746.fastq.gz Ugibba_bladderR1.fastq.gz

#rhizoid rep 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/043/SRR10676743/SRR10676743.fastq.gz
mv SRR10676743.fastq.gz Ugibba_rhizoidR1.fastq.gz

#bladder rep 2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/042/SRR10676742/SRR10676742.fastq.gz
mv SRR10676742.fastq.gz Ugibba_bladderR2.fastq.gz

#stem rep 2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/040/SRR10676740/SRR10676740.fastq.gz
mv SRR10676740.fastq.gz Ugibba_stemR2.fastq.gz

#leaf rep 3
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/049/SRR10676749/SRR10676749.fastq.gz
mv SRR10676749.fastq.gz Ugibba_leafR3.fastq.gz

#bladder rep 3
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/050/SRR10676750/SRR10676750.fastq.gz
mv SRR10676750.fastq.gz Ugibba_bladderR3.fastq.gz

#stem rep 3
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/048/SRR10676748/SRR10676748.fastq.gz
mv SRR10676748.fastq.gz Ugibba_stemR3.fastq.gz

#rhizoid rep 3
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/047/SRR10676747/SRR10676747.fastq.gz
mv SRR10676747.fastq.gz Ugibba_rhizoidR3.fastq.gz

#leaf rep 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/045/SRR10676745/SRR10676745.fastq.gz
mv SRR10676745.fastq.gz Ugibba_leafR1.fastq.gz

#stem rep 1
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/044/SRR10676744/SRR10676744.fastq.gz
mv SRR10676744.fastq.gz Ugibba_stemR1.fastq.gz

#leaf rep 2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/041/SRR10676741/SRR10676741.fastq.gz
mv SRR10676741.fastq.gz Ugibba_leafR2.fastq.gz

#rhizoid rep 2
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/039/SRR10676739/SRR10676739.fastq.gz
mv SRR10676739.fastq.gz Ugibba_rhizoidR2.fastq.gz

#shoots and traps
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR504/008/SRR5046448/SRR5046448_1.fastq.gz
mv SRR5046448_1.fastq.gz Ugibba_shootstraps.fastq.gz

#vegetative pooled
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR106/051/SRR10676751/SRR10676751_1.fastq.gz
mv SRR10676751_1.fastq.gz Ugibba_vegetative.fastq.gz

##########################################################################################
################ Read clean up, this is not as important for mapping to a reference #####
##########################################################################################

#move back up to the parent directory
cd ..
mkdir ../cleaned_reads

#run fastp on each file individually to trim and remove adapters
for file in Raw_files/*.fastq.gz
do
	name=`basename $file .fastq.gz`
	echo "Running fastp on $name"
	forward=$name".fastq.gz"
	fastp -i Raw_files/$forward -o ../cleaned_reads/$forward -z 4 --adapter_sequence=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -q 20 --length_required 100 --thread 2
done


##########################################################################################
################ Read mapping to a reference genome #####
##########################################################################################
cd ../Reference_genome

#need to index the genome for read mapping
bwa index Ugibba_reference.fasta

#move back up to the parent directory and create a new folder for the mapping
cd ..
mkdir BWA
cd BWA
mkdir SAM


#Read group information starts with "@RG
#ID: is unique identifier of the samples, for now doing the sample name and the barcode info
#SM: is the sample name
#PL: is the sequencing equipment, in almost all cases this will be Illumina
#PU: is the run identifier, the lane, followed by the specific barcode of the sample
#LB: is the library count

#read mapping for each file
bwa mem -t 2 -R "@RG\tID:bladderR1\tSM:Rep1\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_bladderR1.fastq.gz > SAM/Ugibba_bladderR1.sam
bwa mem -t 2 -R "@RG\tID:rhizoidR1\tSM:Rep1\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_rhizoidR1.fastq.gz > SAM/Ugibba_rhizoidR1.sam
bwa mem -t 2 -R "@RG\tID:bladderR2\tSM:Rep2\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_bladderR2.fastq.gz > SAM/Ugibba_bladderR2.sam
bwa mem -t 2 -R "@RG\tID:stemR2\tSM:Rep2\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_stemR2.fastq.gz > SAM/Ugibba_stemR2.sam
bwa mem -t 2 -R "@RG\tID:leafR3\tSM:Rep3\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_leafR3.fastq.gz > SAM/Ugibba_leafR3.sam
bwa mem -t 2 -R "@RG\tID:bladderR3\tSM:Rep3\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_bladderR3.fastq.gz > SAM/Ugibba_bladderR3.sam
bwa mem -t 2 -R "@RG\tID:stemR3\tSM:Rep3\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_stemR3.fastq.gz > SAM/Ugibba_stemR3.sam
bwa mem -t 2 -R "@RG\tID:rhizoidR3\tSM:Rep3\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_rhizoidR3.fastq.gz > SAM/Ugibba_rhizoidR3.sam
bwa mem -t 2 -R "@RG\tID:leafR1\tSM:Rep1\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_leafR1.fastq.gz > SAM/Ugibba_leafR1.sam
bwa mem -t 2 -R "@RG\tID:stemR1\tSM:Rep1\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_stemR1.fastq.gz > SAM/Ugibba_stemR1.sam
bwa mem -t 2 -R "@RG\tID:leafR2\tSM:Rep2\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_leafR2.fastq.gz > SAM/Ugibba_leafR2.sam
bwa mem -t 2 -R "@RG\tID:rhizoidR2\tSM:Rep2\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_rhizoidR2.fastq.gz > SAM/Ugibba_rhizoidR2.sam
bwa mem -t 2 -R "@RG\tID:shootstraps\tSM:Rep\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_shootstraps.fastq.gz > SAM/Ugibba_shootstraps.sam
bwa mem -t 2 -R "@RG\tID:veg\tSM:Rep\tPL:HiSeq\tPU:HTNMKDSXX\tLB:RNA-Seq" ../Reference_genome/Ugibba_reference.fasta ../cleaned_reads/Ugibba_vegetative.fastq.gz > SAM/Ugibba_vegetative.sam

##########################################################################################
################ Convert resulting SAM files into BAM files to save space #####
##########################################################################################
mkdir BAM
mkdir sorted_bam

#convert SAM to BAM for sorting
for file in SAM/*.sam
do
	echo "Convert $file to to BAM"
	name=`basename $file .sam`
	samtools view -S -b $file > bam/$name.bam
	rm $file
done

#Sort BAM for SNP calling
for file in bam/*.bam
do
	echo "Sort $file"
	name=`basename $file .bam`
	readid=$name
	samtools sort -o sorted_bam/$readid.bam $file
	rm $file
done


##########################################################################################
################ Mark Duplicates in the resulting BAM files #####
##########################################################################################
cd ..
mkdir refmap_Stacks
cd remap_Stacks
mkdir ref_wrapper

#run the refmap perl script with default settings
ref_map.pl --samples ../BWA/sorted_bam/ --popmap population_map.txt -o ref_wrapper/ -T 2

#export the raw SNP calls with no filtering within Stacks, subsequent filtering will occur in VCFtools
populations --batch_size 1 -P de_novo_wrapper/ -M population_map.txt -t 2 --ordered-export --vcf

#Do the filtering following SNP_filltering.sh
##############################################
#create a new directory for the final outputs
mkdir final_data_set

#read back in the filtered SNP file to export the necessary files for downstream analyses
populations  --batch_size 1 -V Ugibba_LD_pruned_SNPs.recode.vcf -O final_data_set/ -M population_map.txt -t 2 --ordered-export --fstats --vcf --treemix --plink --structure --radpainter --genepop
