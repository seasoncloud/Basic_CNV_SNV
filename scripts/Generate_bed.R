#' Help generate customized bin_bed file of the BED forma as the input for the Gen_bin_cell_atac.R function
#'
#' @param chr_size A 22x2 matrix. The first column is 'chr' and the second column is the size (bps), Each row is a chromosome (from 1 to 22).
#' The matrix for hg19 and GRCh38 can be direcly downloaded from ./data-raw/ directory.
#' @param bin_res Numeric. The fixed bin size for generating the bed file across the chromosomes.
#' @param out_path The path to save the output bed file.
#' 
#' @return A bin_bed file across the chromosomes. 
#'
#' @export
## generate bed files for the bins
Generate_bed=function(chr_size=NULL, bin_res=200000, out_path='./'){
chr_bin=NULL
for(ii in 1:22){
  sub=chr_size[which(chr_size[,1]==paste0('chr',ii)),, drop=F]
  sub=as.numeric(c(0,sub[,2]) )
  minn=min(sub)
  maxx=max(sub)
  rr=maxx-minn
  nn=ceiling(rr/bin_res)
  mm=cbind(rep(ii, nn), (0:(nn-1))*bin_res, c((1:(nn-1))*bin_res, maxx))
  chr_bin=rbind(chr_bin, mm)
  print(ii)
}

return(chr_bin)
write.table(chr_bin,paste0(out_path, "./chrbin", bin_res, ".bed"), sep='\t', col.names=F, row.names=F, quote=F)

}



# For example, if the reference genome is hg19, to generate a bin_bed file with bin size=200k,
# chr_size=read.table("./data-raw/sizes.cellranger-atac-hg19-1.2.0.txt", stringsAsFactors = F, sep='\t')
# chr_size=chr_size[1:22,] 
### read the function
bin_bed=Generate_bed(chr_size = chr_size, bin_res = 200000)
# 