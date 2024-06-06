#!/bin/bash

genome=""
num_cores=1
NT=""

display_help() {
    echo "Usage: $0 -g <genome.fasta> -n <NT_path> -p <num_cores>"
    echo "Options:"
    echo "  -g <genome.fasta>: Path to the genome FASTA file"
    echo "  -n <NT_path>: Path to the NCBI nucleotide NT database (with TXDB files inside also)"
    echo "  -p <num_cores>: Number of CPU cores to use for parallelization (the default is 1)"
    echo "  -h: Display help message"
}

while getopts ":g:n:p:h" opt; do
    case $opt in
        g)
            genome="$OPTARG"
            ;;
        n)
            NT="$OPTARG"
            ;;
        p)
            num_cores="$OPTARG"
            ;;
        h)
            display_help
            exit 0
            ;;
        :)
            echo "Error: Option -$OPTARG requires an argument."
            exit 1
            ;;
        \?)
            echo "Error: Invalid option -$OPTARG"
            exit 1
            ;;
    esac
done

if [[ -z $genome || -z $NT ]]; then
    echo "Error: Some options are missing."
    display_help
    exit 1
fi

# Export BLASTDB variable
export BLASTDB="$NT"

# Split the genome FASTA into individual sequences
num_sequences=$(grep -c ">" "$genome")
seqkit split -i -p "$num_sequences" "$genome" || exit 1

# Run BLAST in parallel for each sequence subset
ls "${genome}.split"/*.fasta | parallel -j "$num_cores" \
    'blastn \
    -task megablast \
    -query {} \
    -db nt \
    -outfmt "6 qseqid sacc staxids sskingdoms sscinames scomnames evalue bitscore" \
    -max_target_seqs 1 \
    -max_hsps 1 \
    -num_threads 32 \
    -evalue 1e-25' \
    > $(basename "$genome").vs.nt.mts1.hsp1.1e25.megablast.out

# Clean up temp files
rm -rf "${genome}.split"