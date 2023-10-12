import pandas as pd
import os
import sys
from Bio import SeqIO

unitig_colors_mtx_path = sys.argv[1]
split_path = sys.argv[2]
frequency_filter = float(sys.argv[3])
outputpath = sys.argv[4]
assemblies = sys.argv[5:]

filebase_to_id = {}

for assembly in assemblies:
    assembly_basename = os.path.basename(assembly)
    basename_woext = '.'.join(assembly_basename.rsplit('.', 1)[:-1])
    for record in SeqIO.parse(assembly, 'fasta'):
        filebase_to_id[record.id] = basename_woext

unitig_color_mtx = pd.read_csv(unitig_colors_mtx_path, delim_whitespace=True)
unitig_color_mtx.columns = unitig_color_mtx.columns.map(lambda x: filebase_to_id[x.split(':')[2]])

with open(split_path, 'r') as file:
    split_values = file.read().splitlines()

def get_binary_vector(X, Y, filter=0.9):
    # Filter the DataFrame X to keep only columns in Y
    filtered_X = X[Y]

    # Calculate the row-wise sum and compare with 90% of the length of Y
    result = (filtered_X.sum(axis=1) / len(Y)) >= filter

    # Convert boolean series to integer
    result = result.astype(int)

    return result

binary_vector = get_binary_vector(unitig_color_mtx, split_values, filter=frequency_filter)

indices = binary_vector[binary_vector == 1].index.tolist()

with open(outputpath, 'w') as f:
    f.writelines(f"{index}\n" for index in indices)
