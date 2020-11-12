source activate py37
cd "~/"
name=sample_name

mkdir "./bam/"$name
mkdir "./vcf/"$name

#rm.dup
java -Xmx4g -jar ~/picard.jar MarkDuplicates \
I= "possorted_bam.bam" \
O= "./bam/"$name"/possorted_rmdup.bam" \
M= "./bam/"$name"/possorted_rmdup_marked_dup_metics.txt" \
REMOVE_DUPLICATES=true \
VALIDATION_STRINGENCY=LENIENT

#sort
samtools sort -o "./bam/"$name"/possorted_rmdup_sort.bam" "./bam/"$name"/possorted_rmdup.bam"

# build index
java -Xmx4g -jar ~/picard.jar BuildBamIndex \
I= "./bam/"$name"/possorted_rmdup_sort.bam" \
VALIDATION_STRINGENCY=LENIENT


# mutation calling
#gatk Mutect2 \
#-R ~/genome.fa \
#-I "./bam/"$name"/possorted_rmdup_sort.bam" \
#-O "./vcf/"$name"/possorted_rmdup_sort.m2.vcf.gz" \
#-tumor sample_name \
#--max-population-af 1 \
#--dont-use-soft-clipped-bases 


gatk --java-options "-Xmx4g" HaplotypeCaller \
-R ~/genome.fa \
-I "./bam/"$name"/possorted_rmdup_sort.bam" \
-O "./vcf/"$name"/possorted_rmdup_sort.hc.vcf.gz"

## decompress file
gunzip "./vcf/"$name"/possorted_rmdup_sort.hc.vcf.gz"