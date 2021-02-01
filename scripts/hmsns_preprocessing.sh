## hmsns_data_preprocessing
## required tools: SRA Toolkit, FastQC, TrimGalore, bwa, samtools, picard, gatk

# For each cell
# download fastq files
  fastq-dump $SRRname


# fastqc & trimming
  mkdir fastqc
	./fastqc -o ./fastqc/ $SRRname'.fastq'
	./trim_galore $SRRname'.fastq' -o ./fastqc
	/fastqc -o ./fastqc ./$SRRname'_trimmed.fq'


## bwa alighment
# bwa index ./genome.fa
  mkdir bwa
  bwa mem ./genome.fa ./fastqc/$SRRname'_trimmed.fq' -o './bwa/'$SRRname'_trimmed.sam'

# sam to bam
  samtools view -S -b './bwa/'$SRRname'_trimmed.sam' > './bwa/'$SRRname'_trimmed.bam'


# sort
samtools sort -o './bwa/'$SRRname'_trimmed_sort.bam' './bwa/'$SRRname'_trimmed.bam'

# rmdup
java -jar ./picard.jar MarkDuplicates \
I='./bwa/'$SRRname'_trimmed_sort.bam' \
O='./bwa/'$SRRname'_trimmed_rmdup.bam' \
M='./bwa/'$SRRname'_trimmed.marked_dup_metrics.txt' \
REMOVE_DUPLICATES=true 


# add read group
java -jar ./picard.jar AddOrReplaceReadGroups \
INPUT='./bwa/'$SRRname'_trimmed_rmdup.bam' \
OUTPUT='./bwa/'$SRRname'_trimmed_rmdup_rg.bam' \
RGID=$SRRname \
RGLB=genome \
RGPL=ILLUMINA \
RGPU=machine \
RGSM=$SRRname


samtools sort -o './bwa/'$SRRname'_trimmed_rmdup_rg_sort.bam' './bwa/'$SRRname'_trimmed_rmdup_rg.bam'

# build index
java -Xmx4g -jar ./picard.jar BuildBamIndex \
I= './bwa/'$SRRname'_trimmed_rmdup_rg_sort.bam' \
VALIDATION_STRINGENCY=LENIENT


# rm
rm './bwa/'$SRRname'_trimmed.sam'
rm './bwa/'$SRRname'_trimmed.bam'
rm './bwa/'$SRRname'_trimmed_sort.bam'
rm './bwa/'$SRRname'_trimmed_rmdup.bam'
rm './bwa/'$SRRname'_trimmed_rmdup_rg.bam'


## snp calling using HaplotypeCaller
mkdir vcf
./gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R ./genome.fa \
   -I './bwa/'$SRRname'_trimmed_rmdup_rg_sort.bam' \
   -O './vcf/'$SRRname'.output.g.vcf.gz' \
   -ERC GVCF

# conslidate vcf files
cd ./vcf
for ii in {1..22};
do
  ./gatk --java-options "-Xmx4g"\
       GenomicsDBImport \
       --genomicsdb-workspace-path './cvcf'$chrr \
       -L $chrr \
       --sample-name-map ./sample_map.txt \
       --tmp-dir=. \
       --reader-threads 5
done


# genotype
for ii in {1..22};
do
  ./gatk --java-options "-Xmx4g" GenotypeGVCFs \
    -R ./genome.fa \
    -V gendb://'./cvcf'$chrr \
    -O 'genotype.chr'$chrr'.vcf.gz'
done



