#' Identify DAPs between cells in two clones with CNA adjustment for scATAC-seq data
#'
#' @param Y1k Number of reads observed for peak k clone 1
#' @param Y2k Number of reads observed for peak k in clone 2
#' @param N1k Number of reads observed across all peaks for all cells in clone 1
#' @param N2k Number of reads observed across all peaks for all cells in clone 2
#' @param f1k Proportion of DNA elements in the peak k region across the genome for clone 1
#' @param f1k Proportion of DNA elements in the peak k region across the genome for clone 2
#'
#' @return A list including "llrstat": test statistic of the test; "pval": p-value of the test; "ratio":(Y1k/N1k)/(Y2k/N2k)
#' 
#' @export
Test_clonalDAP_adj=function(Y1k=NULL, Y2k=NULL, N1k=NULL, N2k=NULL, f1k=NULL, f2k=NULL){
    b4ac=(f1k*N1k+f2k*N2k+f1k*Y2k+f2k*Y1k)^2-4*(f1k*f2k*(N1k+N2k)*(Y1k+Y2k))
    pp=((f1k*N1k+f2k*N2k+f1k*Y2k+f2k*Y1k)-sqrt(b4ac))/(2*(N1k+N2k)*f1k*f2k)
    llrstat=2*(Y1k*log((Y1k+0.000000001)/N1k)+(N1k-Y1k)*log(1-Y1k/N1k)+Y2k*log((Y2k+0.000000001)/N2k)+(N2k-Y2k)*log(1-Y2k/N2k)
               -Y1k*log(f1k*pp)-(N1k-Y1k)*log(1-f1k*pp)-Y2k*log(f2k*pp)-(N2k-Y2k)*log(1-f2k*pp))
    pval=pchisq(q=llrstat, df=1, lower.tail = F)
    ratio=(Y1k/N1k)/((Y2k+0.000000001)/N2k)
    return(list(llrstat=llrstat, pval=pval, ratio=ratio))
}
