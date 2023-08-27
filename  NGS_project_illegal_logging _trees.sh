#!/bin/bash


#$-cwd
# This option tells gridware to change to the current directory before executing the job
# (default is the home of the user).

#$-pe serial 1
# Specify this option only for multithreaded jobs that use more than one cpu core.
# IMPORTANT! Don't use more than 4 cores to keep free space for other students!

#$-l vf=1g
# This option declares the amount of memory the job will use. Specify this value as accurate as possible.
# IMPORTANT! memory request is multiplied by the number of requested CPU cores specified with the -pe.
# Thus, you should divide the overall memory consumption of your job by the number of parallel threads.

#$-q praktikum
#$-N illegal_logging_trees



##command for the actual job

# copying data into our folder

mkdir NGS_project
cd NGS_project/
cp -r /data/proj/teaching/NGS_course/Data/projects/illegal_logging_trees .

#Quality check before trimming
cd illegal_logging_trees/fastqc_raw
mkdir fastqced_raw
fastqc -o fastqced_raw --nogroup -t 2 *.fastq.gz


#Adoptor Trimming
mkdir noadapter
cd noadapter
input_dir="/data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/fastqc_raw/"
output_dir="./"  

for number in {1..5}  
do
    r1_file="${input_dir}wood_sample_${number}_R1.fastq.gz"
    r2_file="${input_dir}wood_sample_${number}_R2.fastq.gz"
    output_r1="${output_dir}pro_wood_sample_${number}_R1.fastq.gz"
    output_r2="${output_dir}pro_wood_sample_${number}_R2.fastq.gz"
    stats_file="${output_dir}noadapters${number}.stats.txt"
    
    bbduk.sh threads=4 in="$r1_file" in2="$r2_file" out="$output_r1" out2="$output_r2" ref=/data/proj/teaching/NGS_course/Softwares/bbmap/resources/adapters.fa ktrim=r k=21 mink=11 hdist=2 tpe tbo stats="$stats_file"
done


cd -

#Removing PhiX and sequencing artifacts (nocontaminants)
mkdir nocontaminants
cd nocontaminants
input_dir="/data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/fastqc_raw/noadapter/"
output_dir="./"  

for number in {1..5}  
do
    r1_file="${input_dir}pro_wood_sample_${number}_R1.fastq.gz"
    r2_file="${input_dir}pro_wood_sample_${number}_R2.fastq.gz"
    output_r1="${output_dir}nocont_wood_sample_${number}_R1.fastq.gz"
    output_r2="${output_dir}nocont_wood_sample_${number}_R2.fastq.gz"
    stats_file="${output_dir}nocontaminants_${number}.stats.txt"
        
    bbduk.sh threads=8 in="$r1_file" in2="$r2_file" out="$output_r1" out2="$output_r2" ref=/data/proj/teaching/NGS_course/Softwares/bbmap/resources/sequencing_artifacts.fa.gz,/data/proj/teaching/NGS_course/Softwares/bbmap/resources/phix174_ill.ref.fa.gz k=31 hdist=1 stats="$stats_file"
done

cd -

#Quality trimming/filtering(sample_clean)
mkdir clean
cd clean
input_dir="/data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/fastqc_raw/nocontaminants/"
output_dir="./"  

for number in {1..5}  
do
    r1_file="${input_dir}nocont_wood_sample_${number}_R1.fastq.gz"
    r2_file="${input_dir}nocont_wood_sample_${number}_R2.fastq.gz"
    output_r1="${output_dir}clean_wood_sample_${number}_R1.fastq.gz"
    output_r2="${output_dir}clean_wood_sample_${number}_R2.fastq.gz"
        
    bbduk.sh threads=8 in="$r1_file" in2="$r2_file" out="$output_r1" out2="$output_r2" qtrim=lr trimq=21 minlength=21 maq=10 maxns=5
done


cd -

#Quality check : after trimming
mkdir fastqc_clean
fastqc -o fastqc_clean --nogroup -t 2 ~/NGS_project/illegal_logging_trees/fastqc_raw/clean/*.fastq.gz

cd ~/NGS_project/illegal_logging_trees/
#De novo assembly
mkdir spades
cd spades
mv ~/NGS_project/illegal_logging_trees/fastqc_raw/clean/*.fastq.gz .

samples=("1" "2" "3" "4" "5")
output_dir_base="spades"

for sample in "${samples[@]}"
do
    output_dir="${output_dir_base}_${sample}"
    input_r1="clean_wood_sample_${sample}_R1.fastq.gz"
    input_r2="clean_wood_sample_${sample}_R2.fastq.gz"
    
    spades.py -o "$output_dir" -1 "$input_r1" -2 "$input_r2" --only-assembler -k 77,99,127 -t 4
done

cd -
mkdir Cedrela_odorata Diospyros_rhombifolia Quercus_mongolica Quercus_robur Swietenia_mahagoni



#Download and Upload FASTA file and GFF3 file to the cluster
#Step 1 : NCBI Download FASTA and GFF3 file from NCBI (https://www.ncbi.nlm.nih.gov/nuccore/NC_043858.1, https://www.ncbi.nlm.nih.gov/nuccore/NC_040009.1, https://www.ncbi.nlm.nih.gov/nuccore/NC_046388.1, https://www.ncbi.nlm.nih.gov/nuccore/NC_039556.1, https://www.ncbi.nlm.nih.gov/nuccore/NC_037251.1)
#Step 2 After downloading them, the file name was changed "{genome_name}_Chl.fasta/gff3". ex)Quercus_mongolica_Chl.fasta
#Step 3 : Transfer the "FASTA/GFF3" files into each species' directory.


#sample 1 Global alignment = Quercus_mongolica
cd Quercus_mongolica
cp ~/NGS_project/illegal_logging_trees/spades/clean_wood_sample_1* .

bbmap.sh -Xmx8g threads=4 sam=1.3 trd ref=Quercus_mongolica_Chl.fasta nodisk in=clean_wood_sample_1_R#.fastq.gz out=sample_1_bbmap_global.sam
samtools view -@ 4 -bS sample_1_bbmap_global.sam  > sample_1_bbmap_global.bam
samtools sort  -@ 4 sample_1_bbmap_global.bam -o sample_1_bbmap_global.sorted.bam
samtools index sample_1_bbmap_global.sorted.bam
qualimap bamqc -nt 4 -c -outformat PDF:HTML -bam sample_1_bbmap_global.sorted.bam
picard MarkDuplicates INPUT=sample_1_bbmap_global.sorted.bam OUTPUT=sample_1_bbmap_global_dedup.bam METRICS_FILE=sample_1_bbmap_global_metrics.txt TAGGING_POLICY=All CREATE_INDEX=true
picard AddOrReplaceReadGroups I=sample_1_bbmap_global_dedup.bam O=sample_1_bbmap_global_RG.bam CREATE_INDEX=true RGID=sample_1 RGLB=lib1 RGPL=illumina RGPU=lane1 RGSM=sample_1


cd -

#sample 2 Global alignment = Swietenia_mahagoni
cd Swietenia_mahagoni
cp ~/NGS_project/illegal_logging_trees/spades/clean_wood_sample_2* .
bbmap.sh -Xmx8g threads=4 sam=1.3 trd ref=Swietenia_mahagoni_Chl.fasta nodisk in=clean_wood_sample_2_R#.fastq.gz out=sample_2_bbmap_global.sam
samtools view -@ 4 -bS sample_2_bbmap_global.sam  > sample_2_bbmap_global.bam
samtools sort  -@ 4 sample_2_bbmap_global.bam -o sample_2_bbmap_global.sorted.bam
samtools index sample_2_bbmap_global.sorted.bam
qualimap bamqc -nt 4 -c -outformat PDF:HTML -bam sample_2_bbmap_global.sorted.bam
picard MarkDuplicates INPUT=sample_2_bbmap_global.sorted.bam OUTPUT=sample_2_bbmap_global_dedup.bam METRICS_FILE=sample_2_bbmap_global_metrics.txt TAGGING_POLICY=All CREATE_INDEX=true
picard AddOrReplaceReadGroups I=sample_2_bbmap_global_dedup.bam O=sample_2_bbmap_global_RG.bam CREATE_INDEX=true RGID=sample_2 RGLB=lib1 RGPL=illumina RGPU=lane1 RGSM=sample_2


cd -

#sample 3 Global alignment = Quercus_robur
cd Quercus_robur
cp ~/NGS_project/illegal_logging_trees/spades/clean_wood_sample_3* .
bbmap.sh -Xmx8g threads=4 sam=1.3 trd ref=Quercus_robur_Chl.fasta nodisk in=clean_wood_sample_3_R#.fastq.gz out=sample_3_bbmap_global.sam
samtools view -@ 4 -bS sample_3_bbmap_global.sam  > sample_3_bbmap_global.bam
samtools sort  -@ 4 sample_3_bbmap_global.bam -o sample_3_bbmap_global.sorted.bam
samtools index sample_3_bbmap_global.sorted.bam
qualimap bamqc -nt 4 -c -outformat PDF:HTML -bam sample_3_bbmap_global.sorted.bam
picard MarkDuplicates INPUT=sample_3_bbmap_global.sorted.bam OUTPUT=sample_3_bbmap_global_dedup.bam METRICS_FILE=sample_3_bbmap_global_metrics.txt TAGGING_POLICY=All CREATE_INDEX=true
picard AddOrReplaceReadGroups I=sample_3_bbmap_global_dedup.bam O=sample_3_bbmap_global_RG.bam CREATE_INDEX=true RGID=sample_3 RGLB=lib1 RGPL=illumina RGPU=lane1 RGSM=sample_3


cd -

#sample 4 Global alignment = Diospyros_rhombifolia
cd Diospyros_rhombifolia
cp ~/NGS_project/illegal_logging_trees/spades/clean_wood_sample_4* .
bbmap.sh -Xmx8g threads=4 sam=1.3 trd ref=Diospyros_rhombifolia_Chl.fasta nodisk in=clean_wood_sample_4_R#.fastq.gz out=sample_4_bbmap_global.sam
samtools view -@ 4 -bS sample_4_bbmap_global.sam  > sample_4_bbmap_global.bam
samtools sort  -@ 4 sample_4_bbmap_global.bam -o sample_4_bbmap_global.sorted.bam
samtools index sample_4_bbmap_global.sorted.bam
qualimap bamqc -nt 4 -c -outformat PDF:HTML -bam sample_4_bbmap_global.sorted.bam
picard MarkDuplicates INPUT=sample_4_bbmap_global.sorted.bam OUTPUT=sample_4_bbmap_global_dedup.bam METRICS_FILE=sample_4_bbmap_global_metrics.txt TAGGING_POLICY=All CREATE_INDEX=true
picard AddOrReplaceReadGroups I=sample_4_bbmap_global_dedup.bam O=sample_4_bbmap_global_RG.bam CREATE_INDEX=true RGID=sample_4 RGLB=lib1 RGPL=illumina RGPU=lane1 RGSM=sample_4


cd -

#sample 5 Global alignment = Cedrela_odorata
cd Cedrela_odorata
cp ~/NGS_project/illegal_logging_trees/spades/clean_wood_sample_5* .
bbmap.sh -Xmx8g threads=4 sam=1.3 trd ref=Cedrela_odorata_Chl.fasta nodisk in=clean_wood_sample_5_R#.fastq.gz out=sample_5_bbmap_global.sam
samtools view -@ 4 -bS sample_5_bbmap_global.sam  > sample_5_bbmap_global.bam
samtools sort  -@ 4 sample_5_bbmap_global.bam -o sample_5_bbmap_global.sorted.bam
samtools index sample_5_bbmap_global.sorted.bam
qualimap bamqc -nt 4 -c -outformat PDF:HTML -bam sample_5_bbmap_global.sorted.bam
picard MarkDuplicates INPUT=sample_5_bbmap_global.sorted.bam OUTPUT=sample_5_bbmap_global_dedup.bam METRICS_FILE=sample_5_bbmap_global_metrics.txt TAGGING_POLICY=All CREATE_INDEX=true
picard AddOrReplaceReadGroups I=sample_5_bbmap_global_dedup.bam O=sample_5_bbmap_global_RG.bam CREATE_INDEX=true RGID=sample_5 RGLB=lib1 RGPL=illumina RGPU=lane1 RGSM=sample_5


cd -
mkdir VariantCalling
cd VariantCalling

#Call variants
#sample 1 = Quercus_mongolica
freebayes -f /data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/Quercus_mongolica/Quercus_mongolica_Chl.fasta /data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/Quercus_mongolica/sample_1_bbmap_global_RG.bam > Quercus_mongolica.raw_variant.vcf
vcfallelicprimitives -kg Quercus_mongolica.raw_variant.vcf > Quercus_mongolica.raw_variant_split.vcf
vcffilter -f 'QUAL / AO > 10' Quercus_mongolica.raw_variant_split.vcf > Quercus_mongolica.DP_Qualifilt.vcf
vcffilter -f 'QUAL / AO > 10' Quercus_mongolica.raw_variant_split.vcf > Quercus_mongolica.filt.vcf
vcftools --vcf Quercus_mongolica.filt.vcf --out Quercus_mongolica.pass --FILTER-summary

#sample 2 = Swietenia_mahagoni
freebayes -f /data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/Swietenia_mahagoni/Swietenia_mahagoni_Chl.fasta /data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/Swietenia_mahagoni/sample_2_bbmap_global_RG.bam > Swietenia_mahagoni.raw_variant.vcf
vcfallelicprimitives -kg Swietenia_mahagoni.raw_variant.vcf > Swietenia_mahagoni.raw_variant_split.vcf
vcffilter -f 'QUAL / AO > 10' Swietenia_mahagoni.raw_variant_split.vcf > Swietenia_mahagoni.DP_Qualifilt.vcf
vcffilter -f 'QUAL / AO > 10' Swietenia_mahagoni.raw_variant_split.vcf > Swietenia_mahagoni.filt.vcf
vcftools --vcf Swietenia_mahagoni.filt.vcf --out Swietenia_mahagoni.pass --FILTER-summary

#sample 3 = Quercus_robur
freebayes -f /data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/Quercus_robur/Quercus_robur_Chl.fasta /data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/Quercus_robur/sample_3_bbmap_global_RG.bam > Quercus_robur.raw_variant.vcf
vcfallelicprimitives -kg Quercus_robur.raw_variant.vcf > Quercus_robur.raw_variant_split.vcf
vcffilter -f 'QUAL / AO > 10' Quercus_robur.raw_variant_split.vcf > Quercus_robur.DP_Qualifilt.vcf
vcffilter -f 'QUAL / AO > 10' Quercus_robur.raw_variant_split.vcf > Quercus_robur.filt.vcf
vcftools --vcf Quercus_robur.filt.vcf --out Quercus_robur.pass --FILTER-summary


#sample 4 = Diospyros_rhombifolia
freebayes -f /data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/Diospyros_rhombifolia/Diospyros_rhombifolia_Chl.fasta /data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/Diospyros_rhombifolia/sample_4_bbmap_global_RG.bam > Diospyros_rhombifolia.raw_variant.vcf
vcfallelicprimitives -kg Diospyros_rhombifolia.raw_variant.vcf > Diospyros_rhombifolia.raw_variant_split.vcf
vcffilter -f 'QUAL / AO > 10' Diospyros_rhombifolia.raw_variant_split.vcf > Diospyros_rhombifolia.DP_Qualifilt.vcf
vcffilter -f 'QUAL / AO > 10' Diospyros_rhombifolia.raw_variant_split.vcf > Diospyros_rhombifolia.filt.vcf
vcftools --vcf Diospyros_rhombifolia.filt.vcf --out Diospyros_rhombifolia.pass --FILTER-summary


#sample 5 = Cedrela_odorata
freebayes -f /data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/Cedrela_odorata/Cedrela_odorata_Chl.fasta /data/proj2/home/students/prk26/NGS_project/illegal_logging_trees/Cedrela_odorata/sample_5_bbmap_global_RG.bam > Cedrela_odorata.raw_variant.vcf
vcfallelicprimitives -kg Cedrela_odorata.raw_variant.vcf > Cedrela_odorata.raw_variant_split.vcf
vcffilter -f 'QUAL / AO > 10' Cedrela_odorata.raw_variant_split.vcf > Cedrela_odorata.DP_Qualifilt.vcf
vcffilter -f 'QUAL / AO > 10' Cedrela_odorata.raw_variant_split.vcf > Cedrela_odorata.filt.vcf
vcftools --vcf Cedrela_odorata.filt.vcf --out Cedrela_odorata.pass --FILTER-summary



#Annotation
#Get SnpEff by cloning
git clone https://github.com/pcingola/SnpEff.git
cd SnpEff/config
mkdir data
cd data
mkdir Cedrela_odorata Diospyros_rhombifolia Quercus_mongolica Quercus_robur Swietenia_mahagoni

#sample 1 = Quercus_mongolica
cd Quercus_mongolica
cp ~/NGS_project/illegal_logging_trees/Quercus_mongolica/Quercus_mongolica_Chl.* .
cat Quercus_mongolica_Chl.fasta > sequences.fa
cat Quercus_mongolica_Chl.gff3 > genes.gff

cd -

#sample 2 = Swietenia_mahagoni
cd Swietenia_mahagoni
cp ~/NGS_project/illegal_logging_trees/Swietenia_mahagoni/Swietenia_mahagoni_Chl.* .
cat Swietenia_mahagoni_Chl.fasta > sequences.fa
cat Swietenia_mahagoni_Chl.gff3 > genes.gff

cd -

#sample 3 = Quercus_robur
cd Quercus_robur
cp ~/NGS_project/illegal_logging_trees/Quercus_robur/Quercus_robur_Chl.* .
cat Quercus_robur_Chl.fasta > sequences.fa
cat Quercus_robur_Chl.gff3 > genes.gff

cd -

#sample 4 = Diospyros_rhombifolia
cd Diospyros_rhombifolia
cp ~/NGS_project/illegal_logging_trees/Diospyros_rhombifolia/Diospyros_rhombifolia_Chl.* .
cat Diospyros_rhombifolia_Chl.fasta > sequences.fa
cat Diospyros_rhombifolia_Chl.gff3 > genes.gff

cd -

#sample 5 = Cedrela_odorata
cd Cedrela_odorata
cp ~/NGS_project/illegal_logging_trees/Cedrela_odorata/Cedrela_odorata_Chl.* .
cat Cedrela_odorata_Chl.fasta > sequences.fa
cat Cedrela_odorata_Chl.gff3 > genes.gff

cd ~/NGS_project/illegal_logging_trees/VariantCalling/SnpEff/

echo "Quercus_mongolica.genome : Quercus_mongolica" >> snpEff.config
snpEff build -gff3 -v Quercus_mongolica
snpEff ann Quercus_mongolica ~/NGS_project/illegal_logging_trees/VariantCalling/Quercus_mongolica.DP_Qualifilt.vcf > Quercus_mongolica_flt_ann.vcf

echo "Swietenia_mahagoni.genome : Swietenia_mahagoni" >> snpEff.config
snpEff build -gff3 -v Swietenia_mahagoni
snpEff ann Swietenia_mahagoni ~/NGS_project/illegal_logging_trees/VariantCalling/Swietenia_mahagoni.DP_Qualifilt.vcf > Swietenia_mahagoni_flt_ann.vcf

echo "Quercus_robur.genome : Quercus_robur" >> snpEff.config
snpEff build -gff3 -v Quercus_robur
snpEff ann Quercus_robur ~/NGS_project/illegal_logging_trees/VariantCalling/Quercus_robur.DP_Qualifilt.vcf > Quercus_robur_flt_ann.vcf

echo "Diospyros_rhombifolia.genome : Diospyros_rhombifolia" >> snpEff.config
snpEff build -gff3 -v Diospyros_rhombifolia
snpEff ann Diospyros_rhombifolia ~/NGS_project/illegal_logging_trees/VariantCalling/Diospyros_rhombifolia.DP_Qualifilt.vcf > Diospyros_rhombifolia_flt_ann.vcf


echo "Cedrela_odorata.genome : Cedrela_odorata" >> snpEff.config
snpEff build -gff3 -v Cedrela_odorata
snpEff ann Cedrela_odorata ~/NGS_project/illegal_logging_trees/VariantCalling/Cedrela_odorata.DP_Qualifilt.vcf > Cedrela_odorata_flt_ann.vcf














