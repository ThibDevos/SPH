import numpy as np
import pandas as pd
import sys
pd.set_option('display.float_format', lambda x: '%.8f' % x)
file1 = sys.argv[1]

warmup = 500

f1 = pd.read_csv(file1, sep=" ", header=None, names=["t", "construction", "calcul", "total"])
f1 = f1.iloc[warmup:]


print(file1)
df_describe_1 = pd.DataFrame(f1)
df_describe_1.describe()
print(df_describe_1.describe()["construction"]["mean"])
print(df_describe_1.describe()["calcul"]["mean"])
print(df_describe_1.describe()["total"]["mean"])


