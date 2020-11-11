import pandas as pd
import io
import os

tmp_dir="~/"


vcf_name=tmp_dir+"sample.vcf"

os.makedirs(tmp_dir+"/vcf_sub/")

vcf_head=[]
vcf_content=[]

with open(vcf_name, mode='r') as vcf:
    for ll in vcf:
        if(ll[0]=='#' and ll[1]=='#'):
            #lls=ll.replace("chr","")##
            #vcf_head.append(lls) ##
            vcf_head.append(ll)
        else:
            lls=ll.strip('\n')
            lls=lls.split('\t')
            vcf_content.append(lls)


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
   ccn=cc.split('hr')[1]
   nochr.append(ccn)

df['#CHROM']=nochr


## select heterozygous
df_sub=df.loc[ (df['ref']!='0') & (df['alt']!='0')]

chr=[]
for ii in range(22):
    chr.append("chr"+str(ii+1))



df_sub2=df_sub.loc[ df_sub['#CHROM'].isin(chr)]


# save the filtered vcf file
with open(tmp_dir+"/vcf_sub/sample_filtered.vcf",'w') as ww:
    for ll in vcf_head:
        ww.write(ll)

df_sub2.to_csv(tmp_dir+"/vcf_sub/sample_filtered.vcf", index=None, sep='\t', mode='a')


#VCF separate

df=df_sub2
df.columns = df.iloc[0]
df=df[1:]

for ii in range(1,23):
    exportfile=tmp_dir+"/vcf_sub/sample_filtered.chr"+str(ii)+'.vcf'
    vcf_chr=df.loc[df['#CHROM']=='chr'+str(ii)]
    with open(exportfile,'w') as ww:
        for ll in vcf_head:
            ww.write(ll)    
    vcf_chr.to_csv(exportfile, sep='\t', index=None, mode='a')


