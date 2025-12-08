import pandas as pd
import sys

df1 = pd.read_csv(sys.argv[1], sep='\t')
df2 = pd.read_csv(sys.argv[2], sep='\t')

merged = pd.merge(df1, df2, on='SampleID', how='outer')

merged.to_csv(sys.argv[3], sep='\t', index=False)