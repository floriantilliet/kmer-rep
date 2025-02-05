------------------------------------------------------------------------

# K-mergency

------------------------------------------------------------------------

## Description

This tool is designed to estimate the compression rate of repetitive elements within a genome assembly. It offers two key functionalities:

1.  It generates graphs that visualize the compression of repeated sequences, focusing on the compression of highly frequent kmers. It also annotates different k-mers based on their sequence, providing a clear overview of their distribution.
2.  It calculates a metric for assessing the compression of repeated kmers, evaluates the overall genome compression, and offers a probabilistic analysis of repetition types based on their frequency. \
    Additionally, if the user has previously performed a Merqury analysis, the tool can refine the genome completeness estimate based on Merqury results, providing a more accurate assessment.

------------------------------------------------------------------------

## Prerequisites

-   Bash version 4.0 or higher
-   Java Runtime Environment (JRE)
-   R with ggplot2 (recommended R 4.0.3+)
-   Jellyfish 2.3.0 or higher
-   GenomeScope2.0

------------------------------------------------------------------------

## Installation

No installation required. Simply clone the repository and run the script.

``` bash
git clone https://github.com/floriantilliet/kmer-rep
cd kmer-rep
```

------------------------------------------------------------------------

## Usage

To run the script, use the following command:

``` bash
./kmergency.sh -i "path_to_reads_fasta path_to_assembly_fasta path_to_annotation_fasta" [-o output] [-f true|false] [-v true|false] [-m path_to_merqury_results_dir] [-h]
```

-   `-h`: Display the help message.
-   `-f`: Force deletion of the output directory if it already exists (default: false).
-   `-m`: Use Merqury results to calculate the completeness value (requires the path to the completeness folder).
-   `-v`: Enable verbose mode for additional log output (default: true).
-   `-o`: Specify a name for the output directory (default: K_MERGENCY).
-   `-i`: Provide the 3 input FASTA files (reads, assembly, and annotation files in this exact order, separated by spaces, between quotes). These files can be in a compressed format.

------------------------------------------------------------------------

## Steps performed by the script:

**K-mer Counting**: The script first uses Jellyfish to count the k-mers in the input reads and assembly files, storing the results in .jf files.\
**K-mer Histogram**: It generates histograms of the k-mer frequencies using the jellyfish histo command.\
**GenomeScope Analysis**: The tool runs GenomeScope on the histograms to assess genome completeness and calculate the average k-mer coverage.\
**K-mer Dumping**: The script dumps k-mers from the reads and the assembly into .txt files. **Merge of Reads and Assembly Dumps**: It merges the k-mer dumps from both the reads and assembly, annotating repeated k-mers from the annotations file.\
**Repetition Statistics**: The tool computes statistics on the repetition rates of k-mers and estimates the genome compression based on the repetition levels.\
**Merqury Completion (Optional)**: If Merqury results are provided, the script calculates a refined genome completeness estimate using Merqury data. **Report Generation**: An HTML report is generated using RMarkdown, containing plots and statistics from the analysis.

------------------------------------------------------------------------

## Outputs

`1_jellyfish_reads_output/`: Contains the Jellyfish count and histogram files for the reads. `2_genomescope_output/`: Contains the GenomeScope output used for genome completeness analysis.\
`3_jellyfish_assembly_output/`: Contains the Jellyfish count file for the assembly.\
`4_dumps/`: Contains the dumps of k-mers from the reads and the assembly. `5_annotated_repeated_kmers/`: Annotates repeated k-mers found in the annotations file. `6_repetitions_stats/`: Contains statistics on k-mer repetitions, including compression and genome assembly estimation.\
`report_<output>.html`: An HTML report with visualizations of the k-mer compression and repetition analysis.

------------------------------------------------------------------------

## Example Workflow
