import sys,os
import pandas as pd
for files in os.listdir('./'):
	if files.startswith('HTSeqCount'):
		df = pd.read_csv(files,header=None,index_col=0,sep='\t')
		try:
			res[files.split('HTSeqCount_')[1].split('.out')[0]] = df.iloc[:,0]
		except:
			res = df
			res.columns = [files.split('HTSeqCount_')[1].split('.out')[0]]

res.iloc[0:-5,:].to_csv("Read_counts.txt", index=True, header=True,sep="\t")

