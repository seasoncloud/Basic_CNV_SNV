# --------------------------------------------------
# Running VarTrix
# --------------------------------------------------
# Geneate SNP by cell matrix from the bam file
# Inputs: vcf file, bam file, fasta file, barcodes.tsv 
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