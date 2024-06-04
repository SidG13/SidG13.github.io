---
title: "Increasing the efficiency of BLAST"
sub_title: "BLASTing multi-fasta queries against multi-fasta databases, with taxonomy information"
categories:
  - Bioinformatics
elements:
  - genomes
  - blast
  - NCBI
  - parallelization
  - taxonomy
---
Below, I provide methods and a script to show you how to set up a BLAST search to identify contaminant sequences (of non-eukaryotic origin) in a eukaryotic genome assembly. You can easily modify this procedure for any use case however (such as different search parameters, databases etc.) by changing the script provided. The main point is to make big BLAST tasks faster by using GNU `parallel`.

**Make sure you have plenty of available space on your machine (several Tb) before uncompressing any files!**

<br>

## Large sequence databases
A problem that arises every once in a while is the need to blast very large query sets, such as a pre-assembly genome or transcriptome, against NCBI databases like the `nt` database, which as of the time of this post, when uncompressed is just shy of 2 Tb! 

Blasting a large genome/transcriptome against this database can take days or weeks. We can greatly speed up this process with GNU `parallel`. 

<br>

## Dependancies

[blast executables](https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html "blast") for the set of main blast algorithms. It can be easily installed on most Linux machines with `conda` or `mamba` (preferred nowadays). For example:
```shell
mamba install -c bioconda blast
```
<br>

[seqkit](https://bioinf.shenwei.me/seqkit/ "seqkit") is used to split a fasta into parts. It too can be installed with `conda` or `mamba`:

```shell
mamba install -c bioconda seqkit
```
<br>

[GNU parallel](https://www.gnu.org/software/parallel/) lets you distribute jobs across different threads. Install the appropriate version for your OS.

<br>

## Downloading and setting up appropriate files

We'll not only need to download the `nt` database for this, but also a few taxonomy files from NCBI. This will allow us to access and parse taxonomy-related information from our BLAST search. All these accessory files get dumped into the `nt/` directory.


**nt.gz**
```shell
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz # Remember, this is a huge file. This may take a while.
gunzip nt.gz # This will take up a lot of space!
```

**Taxonomy files**
```shell
## First, download and move the accession2taxid files (connects accessions from NT to an NCBI taxonomy ID)

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
gunzip nucl_gb.accession2taxid.gz
sed '1d' nucl_gb.accession2taxid | awk '{print $2" "$3}' > taxidmapfile # this needs to be done to convert the taxidmapfile to the correct format
mv taxidmapfile nt/

## Next, download and uncompress other taxonomy files. These encode lots of taxonomy information that can be accessed later through blast search parameters.

wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
mv taxdump.tar.gz nt/
tar -xvzf nt/taxdump.tar.gz
```

**Making the blast database**
```shell
makeblastdb -in nt/ -dbtype nucl -parse_seqids -taxid_map taxidmapfile
```

Now we're ready to run BLAST!

<br>

## Running parallel blast, including taxonomy outputs

By default, the BLAST algorithm will take a single query sequence and perform the local alignment against the entire database. If we have a file with thousands of queries, like a pre-assembly genome or transcriptome, this will take a long time. The basis of this script is to split up a BLAST problem into several individual parts using `seqkit` and `parallel`. If you have the computational power, different threads will partition the BLAST task to run in parallel, speeding it up greatly. 

Because we have also set up taxonomy information in the `nt/` database, we will be able to extract taxonomy-relevant information from our search. Below are some of the options I use when querying, a full list is found [here](https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/ "Blast options").

| Option | Definition | 
|:--------|:-------:|
| sacc  | subject accession   | 
| staxids   | unique subject taxonomy ID(s)   |
| sskingdoms   | unique subject super kingdom(s)   | 
| sscinames   | unique subject scientific name(s)   | 
| scomnames   | unique subject common name(s)   | 
|====================

<br>
#### Parallel BLAST script

The actual BLAST command is near the end of the script, and its options can be modified in any way you choose.

[parallel_genome_blasting.sh][1]

[1]:{{ site.url }}{{ site.baseurl }}/assets/blog_files/parallel_genome_blasting.sh

```bash
#!/bin/bash

# Default values
genome=""
num_cores=1
NT=""

# Function to display help message
display_help() {
    echo "Usage: $0 -g <genome.fasta> -n <NT_path> -p <num_cores>"
    echo "Options:"
    echo "  -g <genome.fasta>: Path to the genome FASTA file"
    echo "  -n <NT_path>: Path to the NCBI nucleotide NT database (with TXDB files inside also)"
    echo "  -p <num_cores>: Number of CPU cores to use for parallelization (default: 1)"
    echo "  -h: Display this help message"
}

# Parse command-line options
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

# Check if required options are provided
if [[ -z $genome || -z $NT ]]; then
    echo "Error: Missing required options."
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

# Clean up temporary files
rm -rf "${genome}.split"
```

## Output

An example output found in the `megablast.out` is shown below. The columns correspond to options set with the `-outfmt` argument.

```bash
Scaffold_31_HRSCAF_176	MT070612	8732	Eukaryota	Crotalus durissus terrificus	tropical rattlesnake	0.0	6948
Scaffold_3_HRSCAF_30	OX365972	103942	Eukaryota	Vipera ursinii	Vipera ursinii	0.0	4106
Scaffold_24_HRSCAF_168	MT019621	8730	Eukaryota	Crotalus atrox	western diamondback rattlesnake	0.0	8527
Scaffold_19_HRSCAF_159	OX365970	103942	Eukaryota	Vipera ursinii	Vipera ursinii	0.0	8347
Scaffold_22_HRSCAF_162	OX365981	103942	Eukaryota	Vipera ursinii	Vipera ursinii	0.0	3099
```

<br>
Quickly <kbd>Control + F</kbd> or <kbd>⌘ + F</kbd> this file to see if there are any non-eukaryotic sequences in your assembly. If so, any non-eukaryotic contaminant sequences that were assembled can easily be identified and removed:
```bash
awk -F'\t' '$4 == "Eukaryota" {print $1}' megablast.out > eukaryotic_sequences.txt
seqkit grep -f eukaryotic_sequences.txt genome.fasta > no_contaminants_genome.fasta
```