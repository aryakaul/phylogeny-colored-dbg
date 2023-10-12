import pandas as pd
import sys

segfile = sys.argv[1]
seqfile = sys.argv[2]
output = sys.argv[3]

# Step 1: Read cf_seg and get unique unitig identifiers
with open(segfile, "r") as file:
    identifiers = [line.split()[0] for line in file]

# Dictionary to store our data
data_dict = {}

# Step 2: Read cf_seq and populate the data_dict
with open(seqfile, "r") as file:
    for line in file:
        parts = line.strip().split()
        ref_seq = parts[0]
        associated_ids = [x[:-1] for x in parts[1:]]  # Removing '+' or '-'
        if ref_seq not in data_dict:
            data_dict[ref_seq] = {id_: 0 for id_ in identifiers}

        # Update the dictionary for the associated IDs
        for id_ in associated_ids:
            if id_ in data_dict[ref_seq]:  # To ensure the ID exists in our dictionary's inner dict
                data_dict[ref_seq][id_] = 1

# Convert the dictionary to a DataFrame
df = pd.DataFrame.from_dict(data_dict, orient='index')

df.transpose().to_csv(output, sep='\t')
