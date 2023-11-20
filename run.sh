#!/bin/bash

# Default values
signal_bam_file=""
control_bam_file=""
bam_file=""
output_path=""
chrom_size=""
subfam_size=""
rmsk_path=""
cage_window=""

# Function to display usage information
usage() {
    echo "Usage: $0 [--signal_bam_file <signal_bam_file>] [--control_bam_file <control_bam_file>] [--bam_file <bam_file>] --output_path <output_path> --chrom_size <chrom_size> --subfam_size <subfam_size> --rmsk_path <rmsk_path> [--cage_window <cage_window>]"
    exit 1
}

# Check if any arguments are provided
if [ "$#" -eq 0 ]; then
    usage
fi

# Parse command-line arguments
while [ "$#" -gt 0 ]; do
    case "$1" in
        --signal_bam_file)
            shift
            signal_bam_file="$1"
            ;;
        --control_bam_file)
            shift
            control_bam_file="$1"
            ;;
        --bam_file)
            shift
            bam_file="$1"
            ;;
        --output_path)
            shift
            output_path="$1"
            ;;
        --chrom_size)
            shift
            chrom_size="$1"
            ;;
        --subfam_size)
            shift
            subfam_size="$1"
            ;;
        --rmsk_path)
            shift
            rmsk_path="$1"
            ;;
        --cage_window)
            shift
            cage_window="$1"
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
    shift
done

# Check for mutually exclusive conditions
if [[ -n "$bam_file" && (-n "$signal_bam_file" || -n "$control_bam_file") ]]; then
    echo "Error: Provide either --bam_file or --signal_bam_file and --control_bam_file, not both."
    usage
fi

# Check if required arguments are provided
if [ -z "$output_path" ] || [ -z "$chrom_size" ] || [ -z "$subfam_size" ] || [ -z "$rmsk_path" ]; then
    echo "Missing required arguments."
    usage
fi

# Display the entered values
if [ -n "$bam_file" ]; then
    echo "BAM File: $bam_file"
else
    echo "Signal BAM File: $signal_bam_file"
    echo "Control BAM File: $control_bam_file"
fi

echo "Output Path: $output_path"
echo "Chrom Size: $chrom_size"
echo "Subfam Size: $subfam_size"
echo "RepeatMasker Path: $rmsk_path"

# Check if cage_window is provided
if [ -n "$cage_window" ]; then
    echo "CAGE Window: $cage_window"
fi
echo ""----------------------------------------""

# Remove trailing slashes from folder path
output_path=$(echo "$output_path" | sed 's:/*$::')

# Display the entered values
cd ./iteres
echo "Iteres processing item: $bam_file"
bam=$(basename -- "$bam_file")
bam="${bam%.*}"
mkdir ./iteres_output

# Cage-Seq
if [ -n "$cage_window" ]; then
    echo "./iteres stat -w -o ./iteres_output/$bam/$bam -W $cage_window $chrom_size $subfam_size $rmsk_path $bam_file"
    echo "./iteres filter -o ./iteres_output/$bam/$bam -W $cage_window $chrom_size $subfam_size $rmsk_path $bam_file"

    python zarrScript.py --local --CAGE -op $output_path/"$item"/"$item".zarr -lp ./iteres_output/"$item"/"$item"_ALL.iteres.loci -ssp ./iteres_output/"$item"/"$item".iteres.subfamily.stat --cageWigPlus ./iteres_output/"$item"/"$item"_+.iteres.wig  --cageWigMinus ./iteres_output/"$item"/"$item"_-.iteres.wig --cageWigUniPlus ./iteres_output/"$item"/"$item"_+.iteres.unique.wig --cageWigUniMinus ./iteres_output/"$item"/"$item"_-.iteres.unique.wig --subfamily_size_filepath "$subfam_size"
elif [ -n "$bam_file" ]; then
    # Atac-Seq and Dnase-Seq
      echo "./iteres stat -w -o ./iteres_output/$bam/$bam $chrom_size $subfam_size $rmsk_path $bam_file"
      echo "./iteres filter -o ./iteres_output/$bam/$bam $chrom_size $subfam_size $rmsk_path $bam_file"

      python zarrScript.py --local -op $output_path/"$item"/"$item".zarr -lp ./iteres_output/"$item"/"$item"_ALL.iteres.loci -ssp ./iteres_output/"$item"/"$item".iteres.subfamily.stat -awp ./iteres_output/"$item"/"$item".iteres.wig -uwp ./iteres_output/"$item"/"$item".iteres.unique.wig --subfamily_size_filepath "$subfam_size"
else
    # ChIP-Seq
      echo "Signal BAM File: $signal_bam_file"
      signal_bam=$(basename -- "$signal_bam_file")
      signal_bam="${signal_bam%.*}"
      echo "./iteres stat -w -o ./iteres_output/$signal_bam/$signal_bam $chrom_size $subfam_size $rmsk_path $signal_bam_file"
      echo "./iteres filter -o ./iteres_output/$signal_bam/$signal_bam $chrom_size $subfam_size $rmsk_path $signal_bam_file"

      echo "Control BAM File: $control_bam_file"
      control_bam=$(basename -- "$control_bam_file")
      control_bam="${control_bam%.*}"
      echo "./iteres stat -w -o ./iteres_output/$control_bam/$control_bam $chrom_size $subfam_size $rmsk_path $control_bam_file"
      echo "./iteres filter -o ./iteres_output/$control_bam/$control_bam $chrom_size $subfam_size $rmsk_path $control_bam_file"

      echo "Generating Zarr"
      mkdir ./zarr_output
      python zarrScript.py --EXPERIMENT --local -op "$output_path/${signal_bam_file}.zarr"\
            --signal_subfam_stat "./iteres_output/$sig/$sig.iteres.subfamily.stat" --control_subfam_stat "./iteres_output/$con/$con.iteres.subfamily.stat" \
            --signal_all_wig_path "./iteres_output/$sig/$sig.iteres.wig" --control_all_wig_path "./iteres_output/$con/$con.iteres.wig" \
            --signal_uni_wig_path "./iteres_output/$sig/$sig.iteres.unique.wig" --control_uni_wig_path "./iteres_output/$con/$con.iteres.unique.wig" \
            --signal_loci_path "./iteres_output/$sig/${sig}_ALL.iteres.loci" --control_loci_path "./iteres_output/$con/${con}_ALL.iteres.loci" \
            --subfamily_size_filepath "$subfam_size"
fi

echo "Done Iteres Processing item: $bam_file"
echo "-- "
