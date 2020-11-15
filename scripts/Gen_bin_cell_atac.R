#' Generate bin by cell matrix for scATAC-seq data
#'
#' @param bin_bed A matrix of the BED format. The first three columns are 'chr', "start site", "end site". Each row is a bin region.
#' @param barcodes A matrix/ data.frame with barcodes for each cell in the first column.
#' @param path_to_fragments The path to the "fragments.tsv.gz" file with the same format as that from the Cell Ranger software.
#' @param out_path The path to save the output matrix.
#' 
#' @import rtracklayer
#' @return A Alleloscope object including the necessary information.
#'
#' @export
## generate bed files for the bins
Gen_bin_cell_atac=function(bin_bed=NULL, barcodes=NULL, path_to_fragments="./fragments.tsv.gz", out_path="./" ){
  # generate bin by cell matrix from fragment file
  barcodes=barcodes$V1
  
  cat("Read fragment file...\n")
  fragments <- import.bed(path_to_fragments, extraCols = c( "type"="character", "type"="integer"))
  colnames(mcols(fragments)) <- c("barcode", "dup_counts")
  
  cat(paste0("Total ", length(fragments)," fragments.\n"))
  #length(unique(fragments$barcode))
  fragments_incell <- fragments[fragments$barcode %in% barcodes]
  
  cat(paste0("Total ", length(fragments_incell)," fragments in cell.\n")) # 236637826
  
  ## check bin bed chr column
  if(grepl("chr",bin_bed[1,1])){
    bin_bed[,1]=gsub("chr","",bin_bed[,1])
  }else{
    bin_bed=bin_bed
  }
  
  chr200k=bin_bed
  chr200k=chr200k[order(as.numeric(chr200k[,1]), as.numeric(chr200k[,2])),]
  bins=paste0('chr',chr200k[,1],':',chr200k[,2],"_", chr200k[,3])
  query=GRanges(paste0('chr',chr200k[,1]), IRanges(chr200k[,2]+1,chr200k[,3]))
  
  ov=findOverlaps(query, fragments_incell )
  ov=as.matrix(ov)
  tmp=fragments_incell$barcode[ov[,2]]
  ov=cbind(ov,match(tmp, barcodes))
  
  cat("Generate bin-by-cell matrix...")
  mm=table(ov[,1],ov[,3])
  colnames(mm)=barcodes[as.numeric(colnames(mm))]
  rownames(mm)=bins[as.numeric(rownames(mm))]
  
  message("The bin-by-cell matrix has beed successfully generated!")
  saveRDS(mm,paste0(out_path,"/bin_cell_atac_fragments.rds"))
  return(mm)
}