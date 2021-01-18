#' Used for generating combined alt_all.rds, ref_all.rds and var_all.rds
#'
#' @param dir_path Directory where the 'chr' directory is located and the output rds files
#' @param vcf_path Path for the 'vcf_sub' directory
#' @param barcodes_path Path to the barcodes.tsv file with one column, each row is a cell barcode.
#' @param samplename Sample name of the data
#' @param plot_stat Whether to generate the statistic plot for each cell and each SNP in the directory
#' 
#' export
# example:Cbn_matrix(dir_path = "./", vcf_path="./vcf_sub/")
Cbn_matrix=function(dir_path='./',vcf_path='./vcf_sub/',barcodes=NULL, samplename='Sample', plot_stat=TRUE){
## setting path and name
vcf_names=list.files(vcf_path)
vcf_cus=sapply(strsplit(vcf_names,"chr"),'[',1)[1]
assay='scDNAseq'
barcode=read.table(paste0(dir_path, "barcodes.tsv"), sep='\t', stringsAsFactors = F)

##combine matrices for each chromosome
alt_all=NULL
ref_all=NULL
var_all=NULL

for(ii in 1:22){
  alt=readMM(paste0(dir_path,"/chr/chr",as.character(ii),"_matrix/out_matrix.mtx"))
  ref=readMM(paste0(dir_path,"/chr/chr",as.character(ii),"_matrix/ref_matrix.mtx"))
  var=read.table(paste0(vcf_path,"/", vcf_cus
                        ,"chr",as.character(ii),".vcf"), header=F, stringsAsFactors = F)
  r0=rowSums(alt)
  ind0=r0>0


  ###

  var_all=rbind(var_all, var[ind0,])
  alt_all=rbind(alt_all, alt[ind0,])
  ref_all=rbind(ref_all, ref[ind0,])
  cat(paste0('chr',ii," "))
}

cat("\n")
rm(alt)
rm(ref)
rm(var)

message(paste0("Total SNP number:",as.character(nrow(alt_all))))
message(paste0("Total cell number:",as.character(ncol(alt_all))))


colnames(alt_all)=barcode[,1]
colnames(ref_all)=barcode[,1]

message(paste0("Save REF matrix, ALT matrix, and variant info to the path:", dir_path))
saveRDS(alt_all,paste0(dir_path,"alt_all.rds"))
saveRDS(ref_all,paste0(dir_path,"ref_all.rds"))
saveRDS(var_all,paste0(dir_path,"var_all.rds"))


total_all=alt_all+ref_all

## calculating VAF for each SNP by pooling all cells
af=rowSums(alt_all)/rowSums(total_all)
af[is.na(af)]=0


## plot the statistics
if(plot_stat==TRUE){
pdf(paste0(dir_path,"/statistics_","_raw.pdf" ))
par(mfrow=c(2,1))
hist(colSums(total_all), main=paste0(samplename," ",assay," coverage (",as.character(dim(total_all)[2])," cells across ",  as.character(dim(total_all)[1]), " SNPs)"), xlab="coverage of individul cells", ylab='frequency', breaks = 300)
hist(rowSums(total_all), main=paste0(samplename," ",assay," coverage (",as.character(dim(total_all)[1])," SNPs across ",  as.character(dim(total_all)[2]), " cells)"), xlab="coverage of individul SNPs", xlim=c(0,100), ylab='frequency', breaks = 1000000)
hist(af, 100, main=paste0(samplename, " ", assay, " VAF for each SNP"))
dev.off()
}
}




