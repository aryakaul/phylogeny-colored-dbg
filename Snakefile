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

def enumerate_all_splits(samples, outputdir, kmersize):
    number_of_possible_splits = len(samples) * 2 - 3
    possiblesplits = list(range(0, number_of_possible_splits))
    parts=["1", "2"]
    return expand(
        os.path.join(
            outputdir, f"unitigs_in_split{{SPLITNO}}_part{{PART}}.txt"
        ),
        SPLITNO=possiblesplits,
        PART=parts
    )

def enumerate_all_phylogenycDBGs(samples, outputdir, kmersize):
    number_of_possible_splits = len(samples) * 2 - 3
    possiblesplits = list(range(0, number_of_possible_splits))
    return expand(
        os.path.join(
            outputdir, f"split{{SPLITNO}}_k{kmersize}_cuttlefish.gfa1"
        ),
        SPLITNO=possiblesplits,
    )


rule all:
    input:
        seqfile=main_dir(f"colored-debruijngraph/colored_dbg_k{KMERSIZE}_cuttlefish.cf_seq"),
        binarymatrix=main_dir("unitig_sample_matrix/unitigs.colors.tsv"),
        phylogenyCdbgs=enumerate_all_phylogenycDBGs(
                samples,
                main_dir(f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/gfas"),
                KMERSIZE
        )


# RULE: Make colored, compacted de bruijn graph using Cuttlefish
# -----------------------------------------------------------------------------


rule run_cuttlefish:
    output:
        seqfile=main_dir(f"colored-debruijngraph/colored_dbg_k{KMERSIZE}_cuttlefish.cf_seq"),
        segfile=main_dir(f"colored-debruijngraph/colored_dbg_k{KMERSIZE}_cuttlefish.cf_seg"),
    params:
        kmersize=KMERSIZE,
        tempdir="/n/scratch3/users/a/ak586/cuttlefish/tmp"
    conda:
        "./envs/cuttlefish.yml"
    threads: 32
    resources:
        mem_gb=50
    shell:
        """
        ulimit -n 2048
        OUTPUTDIR=$(dirname {output.seqfile})
        mkdir -p $OUTPUTDIR
        find $OUTPUTDIR/../assemblies -type f > $OUTPUTDIR/list_of_files.txt
        cuttlefish build \
                -l $OUTPUTDIR/list_of_files.txt \
                -t {threads} \
                -k {params.kmersize} \
                -f 3 \
                -m 45 \
                -w {params.tempdir} \
                -o $OUTPUTDIR/colored_dbg_k{params.kmersize}_cuttlefish
        """


# RULE: Generate all phylogenetic splits on the tree
# -----------------------------------------------------------------------------


def enumerate_allsplitsandparts(samples, outputdir):
    number_of_possible_splits = len(samples) * 2 - 3
    parts = ["1", "2"]
    possiblesplits = list(range(0, number_of_possible_splits))
    return expand(
        os.path.join(outputdir, "unitigs_in_split{SPLITNO}_part{PART}.txt"),
        SPLITNO=possiblesplits,
        PART=parts,
    )


rule split_tree:
    input:
        tree=main_dir("trees/sample_phylogeny.nwk"),
    output:
        enumerate_allsplitsandparts(samples, main_dir("trees/phylogenetic_splits")),
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
        seqfile=main_dir(f"colored-debruijngraph/colored_dbg_k{KMERSIZE}_cuttlefish.cf_seq"),
        segfile=main_dir(f"colored-debruijngraph/colored_dbg_k{KMERSIZE}_cuttlefish.cf_seg"),
    output:
        binarymatrix=main_dir("unitig_sample_matrix/unitigs.colors.tsv"),
    conda:
        "./envs/pandas.yml"
    shell:
        """
        python ./scripts/generate_unitig_colormatrix_cuttlefish.py \
                {input.segfile} \
                {input.seqfile} \
                {output}
        """


# RULE: Generate binary column vectors for each part of each phylogenetic
## split and the provided frequency threshold
# -----------------------------------------------------------------------------


rule unitig_splits_presenceabsence:
    input:
        binarymatrix=main_dir("unitig_sample_matrix/unitigs.colors.tsv"),
        split_considered=main_dir(
            "trees/phylogenetic_splits/split{SPLITNO}_part{PART}.txt"
        ),
    output:
        unitigs_present=main_dir(
            f"present_unitigs_at_frequency_{FREQUENCYTHRESH}/unitigs_in_split{{SPLITNO}}_part{{PART}}.txt"
        ),
    params:
        freq=FREQUENCYTHRESH,
    conda:
        "./envs/pandas.yml"
    shell:
        """
        OUTPUTDIR=$(dirname {output.unitigs_present})
        mkdir -p $OUTPUTDIR

        python scripts/generate_unitig_presenceabsence_cuttlefish.py \
            {input.binarymatrix} \
            {input.split_considered} \
            {params.freq} \
            {output.unitigs_present} \
            $(find $OUTPUTDIR/../assemblies -type f | tr "\\n" " " )
        """


# RULE: For each part of each split make a fasta file containing only the
## unitigs present there. ### Concern --> this is going to increase memory
## requirements substantially?
# -----------------------------------------------------------------------------


rule generate_fasta_splitpart:
    input:
        unitigs_present=main_dir(
            f"present_unitigs_at_frequency_{FREQUENCYTHRESH}/unitigs_in_split{{SPLITNO}}_part{{PART}}.txt"
        ),
        segfile=main_dir(f"colored-debruijngraph/colored_dbg_k{KMERSIZE}_cuttlefish.cf_seg"),
    output:
        fasta_splitpart=main_dir(
            f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/fastas/split{{SPLITNO}}_part{{PART}}.fasta"
        ),
    params:
        freq=FREQUENCYTHRESH,
    conda:
        "./envs/seqkit.yml"
    threads: 12
    shell:
        """
        mkdir -p $(dirname {output.fasta_splitpart})
        cat {input.unitigs_present} | \
                parallel -j {threads} --no-notice grep -w {{}} {input.segfile} | \
                awk '{{printf ">%s\\n%s\\n", $1, $2}}'  \
                > {output}
        """


# RULE: Rerun Cuttlefish for these fasta files!
# -----------------------------------------------------------------------------


rule create_gfas_withcuttlefish:
    input:
        fasta_splitpart1=main_dir(
            f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/fastas/split{{SPLITNO}}_part1.fasta"
        ),
        fasta_splitpart2=main_dir(
            f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/fastas/split{{SPLITNO}}_part2.fasta"
        ),
    output:
        gfa=main_dir(f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/gfas/split{{SPLITNO}}_k{KMERSIZE}_cuttlefish.gfa1"),
        json=main_dir(f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/gfas/split{{SPLITNO}}_k{KMERSIZE}_cuttlefish.json"),
    params:
        kmersize=KMERSIZE,
        tempdir = "/n/scratch3/users/a/ak586/tmp_cuttlefish"
    conda:
        "./envs/cuttlefish.yml"
    threads: 32
    resources:
        mem_gb=30
    shell:
        """
        ulimit -n 2048
        OUTPUTDIR=$(dirname {output.gfa})
        mkdir -p $OUTPUTDIR
        mkdir -p {params.tempdir}
        OUTPUT=$OUTPUTDIR/$(basename {output.gfa} .gfa1)
        cuttlefish build \
                -s {input.fasta_splitpart1} \
                -s {input.fasta_splitpart2} \
                -t {threads} \
                -k {params.kmersize} \
                -m 28 \
                -f 1 \
                -w {params.tempdir} \
                -o $OUTPUT
        """


# RULE: Generate the unitig by color matrix for each of these files
# -----------------------------------------------------------------------------


rule unitig_split_matrix:
    input:
        pcdbg=main_dir(
            f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/gfas/split{{SPLITNO}}_k{KMERSIZE}_bifrost.gfa.gz"
        ),
        colors=main_dir(
            f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/gfas/split{{SPLITNO}}_k{KMERSIZE}_bifrost.color.bfg"
        ),
    output:
        unitigs_fasta_split=main_dir(
            f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/gfas/split{{SPLITNO}}_unitigs.fasta"
        ),
        binarymatrix_split=main_dir(
            f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/gfas/split{{SPLITNO}}_unitigs.colors.tsv"
        ),
    params:
        thresh=0.8,
    threads: 12
    conda:
        "./envs/bifrost.yml"
    shell:
        """
        zcat {input.pcdbg} | \
                awk '{{if ($1=="S") {{print ">" $2 "\\n" $3}}}}' \
                > {output.unitigs_fasta_split}

        Bifrost query \
                -v \
                -t {threads} \
                -e {params.thresh} \
                -g {input.pcdbg} \
                -C {input.colors} \
                -q {output.unitigs_fasta_split} \
                -o $(dirname {output.binarymatrix_split})/$(basename {output.binarymatrix_split} .tsv)
        """


rule add_color_gfatag:
    input:
        pcdbg=main_dir(
            f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/gfas/split{{SPLITNO}}_k{KMERSIZE}_bifrost.gfa.gz"
        ),
        colors=main_dir(
            f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/gfas/split{{SPLITNO}}_unitigs.colors.tsv"
        ),
    output:
        updated_pcdbg=main_dir(
            f"phylogenyColoredDeBruijnGraphs/frequency_{FREQUENCYTHRESH}/gfas/split{{SPLITNO}}_k{KMERSIZE}_bifrost_colorsadded.gfa.gz"
        ),
    threads: 12
    conda:
        "./envs/pandas.yml"
    shell:
        """
        python ./scripts/add_gfatag.py \
                {input.pcdbg} \
                {input.colors} \
                {output}
        """
