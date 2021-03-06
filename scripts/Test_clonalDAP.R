#' Identify DAPs between cells in two clones with CNA adjustment for scATAC-seq data
#'
#' @param Y1k Integer. Number of reads observed in peak k for clone 1
#' @param Y2k Integer. Number of reads observed ub peak k for clone 2
#' @param N1k Integer. Number of reads observed across all peaks for all cells in clone 1
#' @param N2k Integer. Number of reads observed across all peaks for all cells in clone 2
#' @param f1k Required if cna_adj=TRUE. Proportion of DNA elements in the peak k region across the genome for clone 1.
#' @param f1k Required if cna_adj=TRUE. Proportion of DNA elements in the peak k region across the genome for clone 2
#' @param cna_adj Logical(TRUE/FALSE). Whether or not to adjust for the underlying copy number when performing the test. 
#'
#' @return A list including "llrstat": test statistic of the test; "pval": p-value of the test; "ratio":(Y1k/N1k)/(Y2k/N2k)
#' 
#' @export
Test_clonalDAP=function(Y1k=NULL, Y2k=NULL, N1k=NULL, N2k=NULL, f1k=NULL, f2k=NULL, cna_adj=FALSE){
  if(cna_adj==TRUE){
    if(is.null(f1k)|is.null(f2k)){
      stop("f1k and f2k are required to adjust for copy number.")
    }
  }
  if(cna_adj==TRUE){
    b4ac=(f1k*N1k+f2k*N2k+f1k*Y2k+f2k*Y1k)^2-4*(f1k*f2k*(N1k+N2k)*(Y1k+Y2k))
    pp=((f1k*N1k+f2k*N2k+f1k*Y2k+f2k*Y1k)-sqrt(b4ac))/(2*(N1k+N2k)*f1k*f2k)
    llrstat=2*(Y1k*log((Y1k+0.000000001)/N1k)+(N1k-Y1k)*log(1-Y1k/N1k)+Y2k*log((Y2k+0.000000001)/N2k)+(N2k-Y2k)*log(1-Y2k/N2k)
               -Y1k*log(f1k*pp)-(N1k-Y1k)*log(1-f1k*pp)-Y2k*log(f2k*pp)-(N2k-Y2k)*log(1-f2k*pp))
  }else{
    pp=(Y1k+Y2k)/(N1k+N2k)
    llrstat=2*(Y1k*log((Y1k+0.000000001)/N1k)+(N1k-Y1k)*log(1-Y1k/N1k)+Y2k*log((Y2k+0.000000001)/N2k)+(N2k-Y2k)*log(1-Y2k/N2k)
               -Y1k*log(pp)-(N1k-Y1k)*log(1-pp)-Y2k*log(pp)-(N2k-Y2k)*log(1-pp))
  }
  pval=pchisq(q=llrstat, df=1, lower.tail = F)
  ratio=(Y1k/N1k)/((Y2k+0.000000001)/N2k)
  return(list(llrstat=llrstat, pval=pval, ratio=ratio, cna_adj=cna_adj))
}
