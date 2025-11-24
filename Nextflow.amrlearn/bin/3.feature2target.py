#!/usr/bin/env python3

##*************************************************************************##
##   Step3. pre-processing the data to create a feature-to-target table    ##          
##*************************************************************************##

import pandas as pd
import sys

#if the input arguments not 4, showing the usage.
if len(sys.argv)!=4: 
    print("Usage:python3 feature2target.py <vcf_snp_count.txt> <antibiotics table> <outfile> \n\n e.g., python3 feature2target.py vcf_snp_count.txt Antibiotics.txt feature2target.txt ")
    sys.exit()
df = pd.read_csv(sys.argv[1], sep = '\t',dtype=str, header = 0)

# hard code
df1 = df.drop(df.columns[1:6], axis=1).set_index('locus_tag').transpose() #reset the index
df1.reset_index(level=0, inplace = True)
df1 = df1.rename(columns={"index": "locus_tag"})
df2 = pd.read_csv(sys.argv[2], sep = '\t',dtype=str, header = 0)
first_col = df2.columns[0]
# if the antibiotics.txt already have the R, I, S classification, then user can uncomment this code.
#df2 = df2.fillna(0).replace('0','NA').replace(['R','I','S'],[1,0,-1])

df3 = pd.merge(df2, df1, left_on=first_col, right_on='locus_tag', how='right').drop(first_col, axis=1)
df3.to_csv(sys.argv[3],index=False,sep='\t')