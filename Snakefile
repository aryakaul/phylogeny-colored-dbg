import os
import sys
import glob
import argparse

# ╭───────────────────────────────────────────────────────────────────────────╮
#   SETUP
# ╰───────────────────────────────────────────────────────────────────────────╯

datasets = set(
    # ["chlamydia", "mtbc-highquality", "rickettsia", "mycoplasma", "testdata"]
    ["mtbc-highquality"]
)

# Assigning obtained values
DATASET_CHOSEN = config["DATA"]
if DATASET_CHOSEN not in datasets:
    sys.error(f"Dataset chosen is not available from the options: {datasets}")
    sys.exit(2)
MAIN_DIR = os.path.join(os.getcwd(), "data", DATASET_CHOSEN)
KMERSIZE = config["KMER"]
FREQUENCYTHRESH = config["FREQUENCY"]


def main_dir(child):
    return os.path.join(MAIN_DIR, child)


sample_dirs = glob.glob(f"{MAIN_DIR}/assemblies/*")
samples = [os.path.basename(acc) for acc in sample_dirs]

# -----------------------------------------------------------------------------

# ╭───────────────────────────────────────────────────────────────────────────╮
#   RULES
# ╰───────────────────────────────────────────────────────────────────────────╯
#
# RULE ALL
# -----------------------------------------------------------------------------


def enumerate_allsplits(samples, outputdir):
    number_of_possible_splits = len(samples) * 2 - 3
    parts = ["1", "2"]
    possiblesplits = list(range(0, number_of_possible_splits))
    return expand(
        os.path.join(outputdir, "unitigs_in_split{SPLITNO}_part{PART}.txt"),
        SPLITNO=possiblesplits,
        PART=parts,
    )


rule all:
    input:
        #allsplits=enumerate_allsplits(samples, main_dir("trees/phylogenetic_splits")),
        cdbg=main_dir(f"colored-debruijngraph/colored_dbg_k{KMERSIZE}_bifrost.gfa.gz"),
        binarymatrix = main_dir("unitig_sample_matrix/unitigs.colors.tsv"),
        unitig_presence_absencesplits=enumerate_allsplits(samples, main_dir(f"present_unitigs_at_frequency_{FREQUENCYTHRESH}"))


# RULE: Make colored, compacted de bruijn graph using Bifrost
# -----------------------------------------------------------------------------


rule run_bifrost:
    output:
        cdbg=main_dir(f"colored-debruijngraph/colored_dbg_k{KMERSIZE}_bifrost.gfa.gz"),
        colors=main_dir(f"colored-debruijngraph/colored_dbg_k{KMERSIZE}_bifrost.color.bfg")
    params:
        kmersize=KMERSIZE,
    conda:
        "./envs/bifrost.yml"
    threads: 32
    shell:
        """
        OUTPUTDIR=$(dirname {output.cdbg})
        mkdir -p $OUTPUTDIR
        find $OUTPUTDIR/../assemblies -type f > $OUTPUTDIR/list_of_files.txt
        Bifrost build \
                -t {threads} \
                -k {params.kmersize} \
                -c \
                -r $OUTPUTDIR/list_of_files.txt \
                -o $OUTPUTDIR/$(basename {output.cdbg} .gfa.gz)
        """


# RULE: Generate all phylogenetic splits on the tree
# -----------------------------------------------------------------------------


rule split_tree:
    input:
        tree=main_dir("trees/sample_phylogeny.nwk"),
    output:
        enumerate_allsplits(samples, main_dir("trees/phylogenetic_splits")),
    params:
        outputdir=main_dir("trees/phylogenetic_splits"),
    conda:
        "./envs/ete3.yml"
    shell:
        """
        mkdir -p {params.outputdir}
        python scripts/split_tree.py \
                {input.tree} \
                {params.outputdir}
        """


# RULE: Generate unitig by sample matrix
# -----------------------------------------------------------------------------

rule unitig_sample_matrix:
    input:
        cdbg=main_dir(f"colored-debruijngraph/colored_dbg_k{KMERSIZE}_bifrost.gfa.gz"),
        colors=main_dir(f"colored-debruijngraph/colored_dbg_k{KMERSIZE}_bifrost.color.bfg")
    output:
        binarymatrix = main_dir("unitig_sample_matrix/unitigs.colors.tsv")
    params:
        outputdir=main_dir("unitig_sample_matrix"),
    conda:
        "./envs/bifrost.yml"
    shell:
        """
        mkdir -p {params.outputdir}
        zcat {input.cdbg} | \
                awk '{{if ($1=="S") {{print ">" $2 "\\n" $3}}}}' \
                > {params.outputdir}/allunitigs.fasta

        Bifrost query \
                -v \
                -t {threads} \
                -e 1.0 \
                -g {input.cdbg} \
                -C {input.colors} \
                -q {params.outputdir}/allunitigs.fasta \
                -o {params.outputdir}/unitigs.colors
        """


# RULE: Generate binary column vectors for each part of each phylogenetic 
## split and the provided frequency threshold
# -----------------------------------------------------------------------------

rule unitig_splits_presenceabsence:
    input:
        binarymatrix = main_dir("unitig_sample_matrix/unitigs.colors.tsv"),
        split_considered = main_dir("trees/phylogenetic_splits/split{SPLITNO}_part{PART}.txt")
    output:
        unitigs_present = main_dir(f"present_unitigs_at_frequency_{FREQUENCYTHRESH}/unitigs_in_split{{SPLITNO}}_part{{PART}}.txt")
    params:
        freq = FREQUENCYTHRESH
    conda:
        "./envs/pandas.yml"
    shell:
        """
        mkdir -p $(dirname {output.unitigs_present})
        python scripts/generate_unitig_presenceabsence.py \
            {input.binarymatrix} \
            {input.split_considered} \
            {params.freq} \
            {output.unitigs_present}
        """
