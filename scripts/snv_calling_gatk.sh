source activate py37
cd /mnt/upenn/seq_data/scRNA_scDNA_Integration/yun/scDNA/
name=scDNA_P6593_liver_met/

mkdir "./bam/"$name
mkdir "./vcf/"$name

#rm.dup
java -Xmx4g -jar /mnt/upenn/seq_data/scRNA_scDNA_Integration/yun/software/picard.jar MarkDuplicates \
I= "/mnt/upenn/seq_data/scRNA_scDNA_Integration/scDNA/"$name"/outs/possorted_bam.bam" \
O= "./bam/"$name"/possorted_rmdup.bam" \
M= "./bam/"$name"/possorted_rmdup_marked_dup_metics.txt" \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=LENIENT

#sort
/mnt/upenn/seq_data/scRNA_scDNA_Integration/yun/software/samtools-1.10/samtools sort -o "./bam/"$name"/possorted_rmdup_sort.bam" "./bam/"$name"/possorted_rmdup.bam"

# build index
java -Xmx4g -jar /mnt/upenn/seq_data/scRNA_scDNA_Integration/yun/software/picard.jar BuildBamIndex \
I= "./bam/"$name"/possorted_rmdup_sort.bam" \
VALIDATION_STRINGENCY=LENIENT


# mutation calling
#gatk Mutect2 \
#-R /mnt/upenn/seq_data/scRNA_scDNA_Integration/yun/ref/refdata-GRCh38-1.0.0/fasta/genome.fa \
#-I "./bam/"$name"/possorted_rmdup_sort.bam" \
#-O "./vcf/"$name"/possorted_rmdup_sort.m2.vcf.gz" \
#-tumor NCI-N87 \
#--max-population-af 1 \
#--dont-use-soft-clipped-bases 


gatk --java-options "-Xmx4g" HaplotypeCaller \
-R /mnt/upenn/seq_data/scRNA_scDNA_Integration/yun/ref/refdata-GRCh38-1.0.0/fasta/genome.fa \
-I "./bam/"$name"/possorted_rmdup_sort.bam" \
-O "./vcf/"$name"/possorted_rmdup_sort.hc.vcf.gz"

## decompress file
gunzip "./vcf/"$name"/possorted_rmdup_sort.hc.vcf.gz"