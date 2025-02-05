---
---
---

# K-mergency

------------------------------------------------------------------------

## Description

This tool is designed to estimate the compression rate of repetitive elements within a genome assembly.
It offers two key functionalities:

1.  It generates graphs that visualize the compression of repeated sequences, focusing on the compression of highly frequent kmers. It also annotates different k-mers based on their sequence, providing a clear overview of their distribution.
2.  It calculates a metric for assessing the compression of repeated kmers, evaluates the overall genome compression, and offers a probabilistic analysis of repetition types based on their frequency.

Additionally, if the user has previously performed a Merqury analysis, the tool can refine the genome completeness estimate based on Merqury results, providing a more accurate assessment.

## Features

-   **K-mer analysis**: Counts, histograms, and dumps of k-mers from input reads and genome assemblies using Jellyfish.
-   **Compression estimation**: Computes compression of repeated k-mers in the assembly compared to the raw reads.
-   **Genome compression metrics**: Generates statistics on the repetition rates of k-mers in both the genome and the assembly.
-   **Integration with Merqury**: Optionally refines the genome completeness estimate using Merqury results.
-   **Visualization**: Generates an HTML report with graphs to visualize the compression and distribution of repeated sequences.
-   **Repeat annotation**: Annotates repeated k-mers based on their sequence and categorizes them into different repetition types.

------------------------------------------------------------------------

## Prerequisites

-   Bash version 4.0 or higher
-   Java Runtime Environment (JRE)
-   R with ggplot2 (recommended R 4.0.3+)

------------------------------------------------------------------------

## Installation

No installation required.
Simply clone the repository and run the script.

``` bash
git clone https://github.com/floriantilliet/kmer-rep
cd kmer-rep
```

------------------------------------------------------------------------

## Usage

To run the script, use the following command:

``` bash
sbatch run_kmergency.sh -i "reads_fasta assembly_fasta annotation_fasta" -o output_directory -f true|false -v true|false -m merqury_results_dir
```

------------------------------------------------------------------------

## Options

-   `-h`: Display the help message.
-   `-f`: Force deletion of the output directory if it already exists (default: false).
-   `-m`: Use Merqury results to calculate the completeness value (requires the path to the completeness folder).
-   `-v`: Enable verbose mode for additional log output (default: true).
-   `-o`: Specify the output directory (default: K_MERGENCY).
-   `-i`: Provide the input FASTA files (reads, assembly, and annotation files in this exact order, separated by spaces). These files can be compressed.

------------------------------------------------------------------------

## Steps performed by the script:

**K-mer Counting**: The script first uses Jellyfish to count the k-mers in the input reads and assembly files, storing the results in .jf files.
**K-mer Histogram**: It generates histograms of the k-mer frequencies using the jellyfish histo command.
**GenomeScope Analysis**: The tool runs GenomeScope on the histograms to assess genome completeness and calculate the average k-mer coverage.
**K-mer Dumping**: The script dumps k-mers from the reads and the assembly into .txt files.
**Merge of Reads and Assembly Dumps**: It merges the k-mer dumps from both the reads and assembly, annotating repeated k-mers from the annotations file.
**Repetition Statistics**: The tool computes statistics on the repetition rates of k-mers and estimates the genome compression based on the repetition levels.
**Merqury Completion (Optional)**: If Merqury results are provided, the script calculates a refined genome completeness estimate using Merqury data.
**Report Generation**: An HTML report is generated using RMarkdown, containing plots and statistics from the analysis.

------------------------------------------------------------------------

## R Script in the Pipeline

**Purpose** The pipeline includes an R script that is responsible for generating an HTML report based on the k-mer analysis.
This report contains visualizations and statistics about the compression of repeated k-mers in the genome assembly.

**Steps in the R Script** After completing the k-mer and genome analyses, the following R script is executed:

``` bash
Rscript -e "rmarkdown::render('$DIR/../report.rmd', params = list(dir='$DIR', output='$OUTPUT', depth='$DEPTH'), output_dir = '$DIR/', output_file = 'report_$OUTPUT.html')"
```

This command generates an HTML report from an R Markdown file `(report.rmd)` located in the directory specified by the `$DIR` variable.
The parameters passed to the R script include: `dir`: Path to the working directory `output`: Name of the output directory `depth`: Genome depth value (calculated earlier in the script)

The generated HTML report will be saved in the working directory `($DIR/)` as `report_$OUTPUT.html`.

------------------------------------------------------------------------

## Example Workflow

1.  Prepare your input files:
    -   reads_fasta: The raw reads in FASTA format (compressed or uncompressed).
    -   assembly_fasta: The genome assembly in FASTA format (compressed or uncompressed).
    -   annotation_fasta: The annotation file in FASTA format (compressed or uncompressed).
2.  Run the script with SLURM on a cluster (adjust the parameters as needed):

``` bash
sbatch run_kmergency.sh -i "reads.fasta assembly.fasta annotations.fasta" -o my_output -v true
```

3.  Once the script has completed, check the output directory for the results:
    -   A report in HTML format (report_my_output.html).
    -   K-mer statistics and repetition metrics.
    -   Optional Merqury completion values (if Merqury results were provided).

------------------------------------------------------------------------

## Outputs

`1_jellyfish_reads_output/`: Contains the Jellyfish count and histogram files for the reads.
`2_genomescope_output/`: Contains the GenomeScope output used for genome completeness analysis.
`4_dumps/`: Contains the dumps of k-mers from the reads and the assembly.
`5_annotated_repeated_kmers/`: Annotates repeated k-mers found in the annotations file.
`6_repetitions_stats/`: Contains statistics on k-mer repetitions, including compression and genome assembly estimation.
`report_<output>.html`: An HTML report with visualizations of the k-mer compression and repetition analysis.

------------------------------------------------------------------------

## Citing K-mergency

Là je ne sais pas si on doit demander à ce que Klopp et nous soyons cités lors de l'utilisation ?
