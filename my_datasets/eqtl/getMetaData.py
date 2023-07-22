import json

import pandas as pd

filename = "/Users/bfpbr5/Documents/GitHub/jaspar_no_seq.json"
df = pd.read_json(filename,encoding="utf-8", orient="records")
df['class'] = df['class'].apply(lambda row: row[0] if row else '')
df['family'] = df['family'].apply(lambda row: row[0] if row else '')
classes = df['family'].unique()
class_sample = []
for item in list(classes):
    sample = list(df['id'].loc[df['class'] == item])[0]
    class_sample.append(sample)

# filename = "/Users/bfpbr5/Documents/GitHub/errors.json"
# df = pd.read_json(filename,encoding="utf-8", orient="records")
# df['class'] = df['class'].apply(lambda row: row[0] if row else '')
# df['family'] = df['family'].apply(lambda row: row[0] if row else '')
# classes = df['class'].unique()
# class_sample = []
# for item in list(classes):
#     sample = list(df['id'].loc[df['class'] == item])[0]
#     class_sample.append(sample)