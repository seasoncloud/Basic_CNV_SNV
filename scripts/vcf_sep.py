# --------------------------------------------------
# The python script splits a single vcf file into separate vcf files for chromosome 1 to 22.
# --------------------------------------------------
# Geneate SNP by cell matrix from the bam file
# Param<vcf_path>: The path to the vcf file
# Param<out_dir>: The directory path for the output vcf files for chr1-22.
# Param<chr_in>: If the input vcf file with the "chr" labeling
# Param<chr_out>: If the output vcf files with the "chr" labeling


def vcf_split(vcf_path= None, out_dir= None, chr_in= True, chr_out=True ):

# import packages
    import pandas as pd
    import io
    import os

# Input vcf file
    vcf_head=[]
    vcf_content=[]


    with open(vcf_path, mode='r', chr) as vcf:
        for ll in vcf:
            if(ll[0]=='#' and ll[1]=='#'):
                if(chr_in==True and chr_out==True):
                    vcf_head.append(ll)
                else:
                    lls=ll.replace("chr","")##
                    vcf_head.append(lls) ##
            else:
                lls=ll.strip('\n')
                lls=lls.split('\t')
                vcf_content.append(lls)

# select heterozygous sites
    df=pd.DataFrame(vcf_content)
    df.columns = df.iloc[0]
    df=df[1:]

    ref=[]
    alt=[]
    info=df.iloc[:,9].tolist()
    for ii in range(len(info)):
        rd=info[ii].split(':')[1].split(',')
        ref.append(rd[0])
        alt.append(rd[1])

    df['ref']=ref
    df['alt']=alt

    nochr=[]
    CHROM=df.iloc[:,0].tolist()
    for cc in CHROM:
        if(chr_in==True):
            ccn=cc.split('hr')[1]
        else:
            ccn=cc
        nochr.append(ccn)

    df['#CHROM']=nochr


    df_sub=df.loc[ (df['ref']!='0') & (df['alt']!='0')]

# subset the vcf files for chr1-22
    chr=[]
    for ii in range(22):
        chr.append(str(ii+1))

    df_sub2=df_sub.loc[ df_sub['#CHROM'].isin(chr)]

    if(chr_out==True):
        df_sub2['#CHROM']=['chr' + s for s in nochr]
    


# save the filtered vcf file
    with open(out_dir+"/filtered.vcf",'w') as ww:
        for ll in vcf_head:
            ww.write(ll)

    df_sub2.to_csv(out_dir+"/filtered.vcf", index=None, sep='\t', mode='a')


#VCF separate

    df=df_sub2


    for ii in range(1,23):
        exportfile=out_dir+"/chr"+str(ii)+'.vcf'
        if(chr_out==True):
            vcf_chr=df.loc[df['#CHROM']=='chr'+str(ii)]
        else:
            vcf_chr=df.loc[df['#CHROM']==str(ii)]
        with open(exportfile,'w') as ww:
            for ll in vcf_head:
                ww.write(ll)    
        vcf_chr.to_csv(exportfile, sep='\t', index=None, mode='a')


