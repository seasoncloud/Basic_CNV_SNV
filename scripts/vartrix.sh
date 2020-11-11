# Install VarTrix

mkdir chr

for chr in {1..22};
do
	mkdir chr/chr"$chr"_matrix
	/mnt/upenn/seq_data/scRNA_scDNA_Integration/yun/software/vartrix_linux -v "~/chr"$chr".vcf" -b ~/sample.bam -f ~/ref.fa -c ~/barcodes.tsv -s "coverage"	
	mv out_matrix.mtx chr/chr"$chr"_matrix/out_matrix.mtx
	mv ref_matrix.mtx chr/chr"$chr"_matrix/ref_matrix.mtx
done