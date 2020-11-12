# --------------------------------------------------
# Running VarTrix
# --------------------------------------------------
# Geneate SNP by cell matrix from the bam file
# Inputs: vcf file, bam file, fasta file, barcodes.tsv 
# VarTrix only allows parsing one chromosome in each run. 
# This (python) function can help to generate separate vcf files for each chromosome: 
# https://github.com/seasoncloud/Basic_CNV_SNV/blob/main/scripts/vcf_sep.py
# Install VarTrix (https://github.com/10XGenomics/vartrix)

cd ~/
mkdir chr

for chr in {1..22};
do
	mkdir chr/chr"$chr"_matrix
	./vartrix_linux -v "~/vcf/chr"$chr".vcf" -b ~/sample.bam -f ~/ref.fa -c ~/barcodes.tsv -s "coverage"	
	mv out_matrix.mtx chr/chr"$chr"_matrix/out_matrix.mtx
	mv ref_matrix.mtx chr/chr"$chr"_matrix/ref_matrix.mtx
done