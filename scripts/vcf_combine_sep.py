# --------------------------------------------------
# The python script combine two vcf files and splits a single vcf file into separate vcf files for chromosome 1 to 22.
# --------------------------------------------------
# Geneate SNP by cell matrix from the bam file
# Param<vcf_path1>: The path to the first vcf file.
# Param<vcf_path2>: The path to the second vcf file.
# Param<out_dir>: The directory path for the output vcf files for chr1-22.
# Param<chr_in>: If the input vcf file with the "chr" labeling.
# Param<chr_out>: If the output vcf files with the "chr" labeling.

import argparse

__author__ = 'Chi-Yun Wu'

# str2bool function from "https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse"
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def process_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-vp1','--vcf_path1',dest='vcf_path1', help='The path to the first input vcf file',)
    parser.add_argument('-vp2','--vcf_path2',dest='vcf_path2', help='The path to the second input vcf file',)
    parser.add_argument('-od','--out_dir',dest='out_dir', help='The directory path for the output vcf files for chr1-22.',)
    parser.add_argument('-in','--chr_in',dest='chr_in', type=str2bool,nargs='?',const=True,default=True, help='True/False: if the input vcf file with the "chr" labeling',)
    parser.add_argument('-out','--chr_out',dest='chr_out',type=str2bool,nargs='?',const=True,default=True, help='True/False: if the output vcf file with the "chr" labeling',)
    return parser

def vcf_cbn_split(vcf_path1= None,vcf_path2= None, out_dir= None, chr_in= None, chr_out=None ):

# import packages
    import pandas as pd
    import io
    import os

# Input vcf file1
    vcf_head1=[]
    vcf_content1=[]


    with open(vcf_path1, mode='r') as vcf1:
        for ll in vcf1:
            if(ll[0]=='#' and ll[1]=='#'):
                if(chr_in==True and chr_out==True):
                    vcf_head1.append(ll)
                else:
                    lls=ll.replace("chr","")##
                    vcf_head1.append(lls) ##
            else:
                lls=ll.strip('\n')
                lls=lls.split('\t')
                vcf_content1.append(lls)
                
# Input vcf file2
    vcf_head2=[]
    vcf_content2=[]


    with open(vcf_path2, mode='r') as vcf2:
        for ll in vcf2:
            if(ll[0]=='#' and ll[1]=='#'):
                if(chr_in==True and chr_out==True):
                    vcf_head2.append(ll)
                else:
                    lls=ll.replace("chr","")##
                    vcf_head2.append(lls) ##
            else:
                lls=ll.strip('\n')
                lls=lls.split('\t')
                vcf_content2.append(lls)

# select heterozygous sites
    df1=pd.DataFrame(vcf_content1)
    df1.columns = df1.iloc[0]
    df1=df1[1:]

    df2=pd.DataFrame(vcf_content2)
    df2.columns = df2.iloc[0]
    df2=df2[1:]
    
## select SNVs with >=1 reads
    ref1=[]
    alt1=[]
    info=df1.iloc[:,9].tolist()
    for ii in range(len(info)):
        rd=info[ii].split(':')[1].split(',')
        ref1.append(rd[0])
        alt1.append(rd[1])

    df1['ref']=ref1
    df1['alt']=alt1
    
    
    ref2=[]
    alt2=[]
    info=df2.iloc[:,9].tolist()
    for ii in range(len(info)):
        rd=info[ii].split(':')[1].split(',')
        ref2.append(rd[0])
        alt2.append(rd[1])

    df2['ref']=ref2
    df2['alt']=alt2
    
# nochr    
    nochr1=[]
    CHROM=df1.iloc[:,0].tolist()
    for cc in CHROM:
        if(chr_in==True):
            ccn=cc.replace('chr',"")
        else:
            ccn=cc
        nochr1.append(ccn)

    df1['#CHROM']=nochr1


    df1_sub=df1.loc[ (df1['ref']!='0') & (df1['alt']!='0')]
    
    
    nochr2=[]
    CHROM=df2.iloc[:,0].tolist()
    for cc in CHROM:
        if(chr_in==True):
            ccn=cc.replace('chr',"")
        else:
            ccn=cc
        nochr2.append(ccn)

    df2['#CHROM']=nochr2


    df2_sub=df2.loc[ (df2['ref']!='0') & (df2['alt']!='0')]

# subset the vcf files for chr1-22
    chr=[]
    for ii in range(22):
        chr.append(str(ii+1))

    df1_sub2=df1_sub.loc[ df1_sub['#CHROM'].isin(chr)]
#    if(chr_out==True):
#        df1_sub2['#CHROM']=['chr' + s for s in df1_sub2['#CHROM']]


    df2_sub2=df2_sub.loc[ df2_sub['#CHROM'].isin(chr)]
#    if(chr_out==True):
#        df2_sub2['#CHROM']=['chr' + s for s in df2_sub2['#CHROM']]


# combind the two vcf files
    df1_sub2=df1_sub2.rename(columns={df1_sub2.columns[9]:'Sample'})
    df2_sub2=df2_sub2.rename(columns={df2_sub2.columns[9]:'Sample'})
    frames = [df1_sub2, df2_sub2]
    df_all=pd.concat(frames)
    df_all['#CHROM'] = df_all['#CHROM'].astype(int)
    df_all['POS'] = df_all['POS'].astype(int)
    df_all=df_all.sort_values(['#CHROM','POS'], ascending=[True, True])

    df_all["#CHROM"]=df_all["#CHROM"].astype(str)
    df_all["POS"]=df_all["POS"].astype(str)

    tmp = df_all[['#CHROM','POS']]
    tmp2=tmp.duplicated(keep='first')
    #tmp3=tmp2.value_counts()
    df_all=df_all[~tmp2]


    if(chr_out==True):
        df_all['#CHROM']=['chr' + s for s in df_all['#CHROM']]

# save the filtered vcf file
    with open(out_dir+"/filtered.vcf",'w') as ww:
        for ll in vcf_head1:
            ww.write(ll)

    df_all.to_csv(out_dir+"/filtered.vcf", index=None, sep='\t', mode='a')


#VCF separate

    df=df_all

    for ii in range(1,23):
        exportfile=out_dir+"/chr"+str(ii)+'.vcf'
        if(chr_out==True):
            vcf_chr=df.loc[df['#CHROM']=='chr'+str(ii)]
        else:
            vcf_chr=df.loc[df['#CHROM']==str(ii)]
        with open(exportfile,'w') as ww:
            for ll in vcf_head1:
                ww.write(ll)    
        vcf_chr.to_csv(exportfile, sep='\t', index=None, mode='a')


def main(args=None):
    parser = process_parser()
    args = parser.parse_args()
    vcf_cbn_split(vcf_path1= args.vcf_path1, vcf_path2= args.vcf_path2, out_dir= args.out_dir, chr_in= args.chr_in, chr_out= args.chr_out )
    

if __name__ == '__main__':
    main()
    
