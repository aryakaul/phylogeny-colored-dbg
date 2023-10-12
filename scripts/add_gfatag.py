import sys
import gzip
import seaborn as sns

gfa_in = sys.argv[1]
colors_in = sys.argv[2]
gfa_out = sys.argv[3]

dict_of_colors = {}
dict_unitig_to_color = {}

with open(colors_in, 'r') as f:
    header=next(f)
    for lines in f:
        fields = lines.rstrip().split()
        unitig, pres_abs_vector = fields[0], ' '.join(fields[1:])

        if pres_abs_vector not in dict_of_colors:
            dict_of_colors[pres_abs_vector] = len(dict_of_colors)

        dict_unitig_to_color[fields[0]] = dict_of_colors[pres_abs_vector]

max_colors = len(dict_of_colors)
palette = sns.color_palette("deep", max_colors).as_hex()

for unitig, index in dict_unitig_to_color.items():
    dict_unitig_to_color[unitig] = palette[index]

with gzip.open(gfa_in, 'rt') as f, gzip.open(gfa_out, 'wt') as o:
    o.write("#\t" + header)
    o.write("#\t" + "-".join(dict_of_colors.keys()) + "\n")
    o.write("#\t" + "-".join(palette) + "\n")

    for i in f:
        if not i.startswith('S'):
            o.write(i)
        else:
            segment_line = i.rstrip().split()
            color = dict_unitig_to_color[segment_line[1]]
            segment_line.append(f"CB:Z:{color}")
            o.write("\t".join(segment_line) + "\n")
