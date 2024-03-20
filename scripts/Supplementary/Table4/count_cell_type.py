import pandas as pd
df = pd.read_csv("/storage/chentemp/zz4/adult_dev_compare/results/Annotation/ALL.csv")

df.celltype.value_counts().to_csv("/storage/chentemp/zz4/adult_dev_compare/Supplementary/Table4/celltype_meta.csv")