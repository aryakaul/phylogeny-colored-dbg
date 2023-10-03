usage() {
    echo "
    Snakemake runner.
        -d, --data              Name of the dataset to use. Default is 'testdata' containing 2 M. bovis genomes and 2 M. Tb genomes
        -k, --kmersize          Length of k-mer to use. Default is 31
        -t, --tree              Unrooted phylogenetic tree corresponding to the samples. Each leaf should have the same basename
                                as the genome file basename
        -f, --frequency         Frequency threshold to classify a given unitig as appearing in a split. Default is 0.9.
    "
}

run_snakemake() {
    logdir="logs/$(basename $DATA)"
    sed "s|~~LOGSDIR~~|$logdir|g" ./o2_profile/config_template.yaml > ./o2_profile/config.yaml
    mkdir -p $logdir
    snakemake \
        -p \
        --conda-frontend mamba \
        --profile o2_profile \
        --cluster-status ./o2_profile/status-sacct.sh \
        --config DATA\=$DATA KMER\=$KMERSIZE FREQUENCY\=$FREQUENCY \
        -k \
        --reason \
        --verbose \

}

main() {
    DATA="testdata"
    KMERSIZE=31
    FREQUENCY=0.9
    for i in "$@"; do
    case $i in
        -d|--data)
        DATA="$2"
        shift
        shift
        ;;
        -t|--tree)
        TREE="$2"
        shift
        shift
        ;;
        -k|--kmersize)
        KMERSIZE="$2"
        shift
        shift
        ;;
        -f|--frequency)
        FREQUENCY="$2"
        shift
        shift
        ;;
        -h|--help)
        usage
        exit
        ;;
        *)
        ;;
    esac
    done
    run_snakemake
}
main "$@"
