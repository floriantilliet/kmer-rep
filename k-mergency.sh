#!/bin/bash

#SBATCH --job-name=K-mergency
#SBATCH --output=K-mergency.out
#SBATCH --error=K-mergency.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=256G


module purge
module load bioinfo/Jellyfish/2.3.0
module load statistics/R/4.2.2
module load bioinfo/GenomeScope2.0/eca7b88

OUTPUT=K_MERGENCY
VERBOSE=true
FORCE=false
MERQURY=false

while getopts "hfm:v:o:i:" option
do
        case $option in
                h)
                        echo "Usage: $0 [-h] [-f true|false] [-m merqury_results_dir] [-v true|false] [-o output_dir] -i \"reads_fasta assembly_fasta annotation_fasta\""
                        echo ""
                        echo "Options:"
                        echo "  -h  Display this help message."
                        echo "  -f  Force the deletion of the output directory if it already exists. Default: false."
                        echo "  -m  Use Merqury results to calculate the completion value. Requires the path to the completeness folder."
                        echo "  -v  Verbose mode. Default: true."
                        echo "  -o  Output directory. Default: K_MERGENCY."
                        echo "  -i  Input fasta files. Requires exactly 3 arguments separated by spaces: reads_fasta assembly_fasta annotation_fasta."
                        echo "     The paths to the 3 input files have to be written between quotes, seperated by white spaces, and in this exact order."
                        echo "     The input fasta files can be compressed."
                        exit 0
                        ;;
                f)
                        if [[ "$OPTARG" != "true" && "$OPTARG" != "false" ]];
                        then
                        echo "Error: -f option requires an argument of either 'true' or 'false'."
                        exit 1
                        fi
                        FORCE="$OPTARG"
                        ;;
                m)
                        MERQURY=true
                        COMPLETENESS="$OPTARG"
                        ;;
                v)
                        echo "$OPTARG"
                        if [[ "$OPTARG" != "true" && "$OPTARG" != "false" ]];
                        then
                        echo "Error: -v option requires an argument of either 'true' or 'false'."
                        exit 1
                        fi
                        VERBOSE="$OPTARG"
                        ;;
                o)
                        OUTPUT="$OPTARG"
                        ;;
                i)
                        IFS=' ' read -r READS ASSEMBLY ANNOTATIONS <<< "$OPTARG"
                        if [[ -z "$ANNOTATIONS" || ! -z "$(echo $OPTARG | awk '{print $4}')" ]];
                        then
                        echo "Error: -i option requires exactly 3 substrings separated by spaces."
                        exit 1
                        fi
                        ;;
        esac
done

if (file $READS | grep -q "gzip compressed data");
then
    READS="<(zcat $READS)"
fi
if (file $ASSEMBLY | grep -q "gzip compressed data");
then
    ASSEMBLY="<(zcat $ASSEMBLY)"
fi
if (file $ANNOTATIONS | grep -q "gzip compressed data");
then
    ANNOTATIONS="<(zcat $ANNOTATIONS)"
fi

if $FORCE = true;
then
    rm -rf $OUTPUT
fi

DIR=$(pwd)/$OUTPUT

mkdir -p $DIR/1_jellyfish_reads_output/

if $VERBOSE = true;
then
    echo "Starting jellyfish count of k-mers in the reads file: $READS"
    echo ""
fi

if [[ ! -s "$DIR/1_jellyfish_reads_output/reads.jf" ]];
then
jellyfish count -m 21 -s 1G -t 16 -C $READS -o $DIR/1_jellyfish_reads_output/reads.jf
fi

if $VERBOSE = true;
then
    echo "Jellyfish count of k-mers in the reads is done, output file: $DIR/1_jellyfish_reads_output/reads.jf"
    echo ""
    echo "Starting jellyfish histo of k-mers in the reads"
    echo ""
fi

if [[ ! -s "$DIR/1_jellyfish_reads_output/reads.histo" ]];
then
jellyfish histo -t 16 $DIR/1_jellyfish_reads_output/reads.jf > $DIR/1_jellyfish_reads_output/reads.histo
fi

if $VERBOSE = true;
then
    echo "Jellyfish histo of k-mers in the reads is done, output file: $DIR/1_jellyfish_reads_output/reads.histo"
    echo ""
    echo "Starting GenomeScope analysis"
    echo ""
fi

if [[ ! -s "$DIR/2_genomescope_output" ]];
then
Rscript /usr/local/bioinfo/src/GenomeScope2.0/genomescope2.0-eca7b88/genomescope.R -i $DIR/1_jellyfish_reads_output/reads.histo -o $DIR/2_genomescope_output -k 21
fi

COV=$(grep "kmercov" $DIR/2_genomescope_output/model.txt | awk '{print $2}' | grep -Eo '[0-9.]+e[+-][0-9]+' | awk '{printf "%.0f\n", $1}')

DEPTH=$(echo $COV | awk '{a=$0*2; printf(a"\n")}')

THRESHOLD=$(echo $COV | awk '{a=$0*20; printf(a"\n")}')

if $VERBOSE = true;
then
    echo "GenomeScope analysis is done, output file: $DIR/2_genomescope_output"
    echo ""
    echo "Starting jellyfish dump of k-mers in the reads"
    echo ""
fi

if [[ ! -s "$DIR/4_dumps/reads_dump_$OUTPUT.txt" ]];
then
mkdir -p $DIR/4_dumps/

jellyfish dump -c -t --lower-count=$THRESHOLD $DIR/1_jellyfish_reads_output/reads.jf | sort -k1,1 > $DIR/4_dumps/reads_dump_$OUTPUT.txt
fi

if $VERBOSE = true;
then
    echo "Jellyfish dump of k-mers in the reads is done, output file: $DIR/4_dumps/reads_dump_$OUTPUT.txt"
    echo ""
    echo "Starting jellyfish count of k-mers in the assembly file: $ASSEMBLY"
    echo ""
fi

if [[ ! -s "$DIR/3_jellyfish_assembly_output/assembly_$OUTPUT.jf" ]];
then
mkdir -p $DIR/3_jellyfish_assembly_output/

jellyfish count -m 21 -s 100M -t 16 -C $ASSEMBLY -o $DIR/3_jellyfish_assembly_output/assembly_$OUTPUT.jf
fi

if $VERBOSE = true;
then
    echo "Jellyfish count of k-mers in the assembly is done, output file: $DIR/3_jellyfish_assembly_output/assembly_$OUTPUT.jf"
    echo ""
    echo "Starting jellyfish dump of k-mers in the assembly"
    echo ""
fi

if [[ ! -s "$DIR/4_dumps/assembly_dump_$OUTPUT.txt" ]];
then
jellyfish dump -t -c --lower-count=2 $DIR/3_jellyfish_assembly_output/assembly_$OUTPUT.jf | sort -k1,1 > $DIR/4_dumps/assembly_dump_$OUTPUT.txt
fi

if $VERBOSE = true;
then
    echo "Jellyfish dump of k-mers in the assembly is done, output file: $DIR/4_dumps/assembly_dump_$OUTPUT.txt"
    echo ""
    echo "Starting dump merge of k-mers in the reads and assembly"
    echo ""
fi

if [[ ! -s "$DIR/4_dumps/merged_dump_$OUTPUT.txt" ]];
then
join -j 1 -t $'\t' -a 1 $DIR/4_dumps/reads_dump_$OUTPUT.txt $DIR/4_dumps/assembly_dump_$OUTPUT.txt > $DIR/4_dumps/temp_merged_dump_$OUTPUT.txt

awk -F'\t' '{ if (NF < 3) print $0, 1; else print $0 }' OFS='\t' $DIR/4_dumps/temp_merged_dump_$OUTPUT.txt > $DIR/4_dumps/merged_dump_$OUTPUT.txt
fi

if $VERBOSE = true;
then
    echo "Dump merge of k-mers in the reads and assembly is done, output file: $DIR/4_dumps/merged_dump_$OUTPUT.txt"
    echo ""
    echo "Starting annotation of repeated k-mers"
    echo ""
fi

if [[ ! -s "$DIR/5_annotated_repeated_kmers/annotation_kmers_$OUTPUT.txt" ]];
then

    mkdir -p $DIR/5_annotated_repeated_kmers

    awk '
    /^>/ {
        if ($0 !~ /#Unknown/ && $0 !~ /\?/) {
            split($0, tmp1, "#")
            full_repeat_type = tmp1[2]
            split(full_repeat_type, tmp2, "/")
            repeat_type = tmp2[1]
        } else {
            repeat_type = ""
        }
    }
    !/^>/ && repeat_type != "" {
        seq = $0
        for (i = 1; i <= length(seq) - 20; i++) {
            kmer = substr(seq, i, 21)
            if (kmer in dict) {
                if (dict[kmer] != repeat_type) {
                    delete dict[kmer]
                }
            } else {
                dict[kmer] = repeat_type
            }
        }
    }
    END {
        for (k in dict) {
            if (dict[k] != "") {
                print k "\t" dict[k]
            }
        }

    }' "$ANNOTATIONS" | sort -k1,1 > $DIR/5_annotated_repeated_kmers/annotation_kmers_$OUTPUT.txt

fi

if $VERBOSE = true;
then
    echo "Annotation of repeated k-mers is done, output file: $DIR/5_annotated_repeated_kmers/annotation_kmers_$OUTPUT.txt"
    echo ""
    echo "Starting merge of dump and annotation of repeated k-mers"
    echo ""
fi

if [[ ! -s "$DIR/5_annotated_repeated_kmers/merged_dump_annotated_$OUTPUT.txt" ]];
then
join -1 1 -2 1 -t $'\t' $DIR/4_dumps/merged_dump_$OUTPUT.txt $DIR/5_annotated_repeated_kmers/annotation_kmers_$OUTPUT.txt > $DIR/5_annotated_repeated_kmers/merged_dump_annotated_$OUTPUT.txt
fi

if $VERBOSE = true;
then
    echo "Merge of dump and annotation of repeated k-mers is done, ouput file: $DIR/5_annotated_repeated_kmers/merged_dump_annotated_$OUTPUT.txt"
    echo ""
    echo "Starting statistics on repetitions"
    echo ""
fi

if [[ ! -s "$DIR/6_repetitions_stats/repetition_stats_$OUTPUT.txt" ]];
then

mkdir -p $DIR/6_repetitions_stats

awk -F'\t' -v depth=$DEPTH '
BEGIN {
    count_10000 = 0; count_200 = 0; count_10 = 0;
    compression_10 = 0; compression_200 = 0; compression_10000 = 0;
    observed_repeated_kmer = 0;
}

{
# normalize the reads kmer count with the depth to allow its comparison with the assembly kmer count
expected_count = $2 / depth;

# compression of a given kmer in the assembly compared to the reads
compression = $3 / expected_count;

total_expected_count += $2;
total_observed_count += $3;

if (expected_count > 10) {
    expected_repeated_kmer += $2;

    count_10 += 1*expected_count;
    if (compression >= 1) {
        compression_10 += 1*expected_count;
    }
    else{
        compression_10 += compression*expected_count;
    }
}

if (expected_count > 200) {
    count_200 += 1*expected_count;
    if (compression >= 1) {
        compression_200 += 1*expected_count;
    }
    else{
        compression_200 += compression*expected_count;
    }
}

if (expected_count > 10000) {
        count_10000 += 1*expected_count;
        if (compression >= 1) {
            compression_10000 += 1*expected_count;
        }
        else{
            compression_10000 += compression*expected_count;
        }
}
}

END {

    sum_compression_10 = 100 - 100 * compression_10/count_10;

    repeated_kmers_in_genome = 100 * expected_repeated_kmer/84000000000;
    
    repeated_kmers_in_assembly = 100 * observed_repeated_kmer/total_observed_count;

    print "Compression of k-mers repeated at least 10 times in the genome:",sum_compression_10, "%";
    print "This represents all the k-mers we considered sufficiently present to be repeated. \n";

    print "Compression of k-mers repeated at least 200 times in the genome:",100 - 100 * compression_200/count_200, "%";
    print "This represents highly repeated k-mers, which are therefore unlikely to be related to gene duplications or highly conserved elements in the genome (which are most often repeated fewer than 100 times). \n";

    print "Compression of k-mers repeated at least 10000 times in the genome:", 100 - 100 * compression_10000/count_10000, "%";
    print "This represents highly repeated k-mers generally found in transposons or satellite DNA (the only ones often repeated tens of thousands of times). \n";
 
    print "Percentage of repeated kmers in the genome assembly:", repeated_kmers_in_assembly, "%";

    print "Estimation of assembly genome compression", 100 - ((50 * sum_compression_10) / 100), "%" ;
}' $DIR/4_dumps/merged_dump_$OUTPUT.txt > $DIR/6_repetitions_stats/repetition_stats_$OUTPUT.txt

fi

if $VERBOSE = true;
then
    echo "Statistics on repetitions is done, ouput file: $DIR/6_repetitions_stats/repetition_stats_$OUTPUT.txt"
    echo ""
    echo "Starting html report generation"
    echo ""
fi

Rscript -e "rmarkdown::render('$DIR/../report.rmd', params = list(dir='$DIR', output='$OUTPUT', depth='$DEPTH'), output_dir = '$DIR/', output_file = 'report_$OUTPUT.html')"

if $VERBOSE = true;
then
    echo "Html report generation is done, output file: $DIR/report_$OUTPUT.html"
    echo ""
fi

if $MERQURY = true;
then
    if $VERBOSE = true;
    then
        echo "Starting Merqury completion value calculation"
        echo ""
    fi
    completeness_files=($COMPLETENESS/*.completeness.stats)
    if [ ${#completeness_files[@]} -eq 0 ];
    then
        echo "No completeness stats files found!"
        exit 1
    fi

    merqury_value=$(awk -F'\t' '{print $5}' "${completeness_files[0]}")
    sum_above_10DEPTH=0
    total_sum=0

    while read -r col1 col2; do
        total_sum=$(echo "$total_sum + $col2" | bc -l)
        
        if (( $(echo "$col2 > 10 * $DEPTH" | bc -l) )); then
            sum_above_10DEPTH=$(echo "$sum_above_10DEPTH + $col2" | bc -l)
        fi
    done < "$DIR/1_jellyfish_reads_output/reads.histo"


    if (( $(echo "$total_sum != 0" | bc -l) )); then
        percent_repeted=$(echo "$sum_above_10DEPTH / $total_sum" | bc -l)
        echo "Percent of reptition is : $percent_repeted"
    else
        echo "The total sum is equal to zero, impossible to calculate."
    fi
    
    sum_compression_20=$(awk -F': ' '/Compression of k-mers repeated at least 20 times in the genome/ {gsub(/%/, "", $2); print $2}' $DIR/6_repetitions_stats/repetition_stats_$OUTPUT.txt)
    completion=$(echo "$merqury_value - ($sum_compression_20 * $percent_repeted) / 100" | bc -l)
    completion=$(printf "%.3f" $completion)
  {
    echo "Percent of repeated kmers: $percent_repeted %"
    echo "Merqury completion value: $merqury_value %"
    echo "Completion value corrected: $completion %"
  } >> $DIR/6_repetitions_stats/repetition_stats.txt
    if $VERBOSE = true;
    then
        echo "Merqury completion value calculation is done, output file: $DIR/6_repetitions_stats/repetition_stats.txt"
        echo ""
    fi
fi

echo "The pipeline was successful"
