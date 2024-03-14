# Repeat Browser Data Processing

---

In this repo, we provide a pipeline to process the transposable elements(TE),
which is used for getting alignment statistics of transposable elements and 
save the processed data into a proper file format that suitable for our repeat browser. 

This pipeline includes a modified [iteres](https://epigenome.wustl.edu/iteres/) pipeline and 
a python script to convert the analysis into [zarr](https://zarr.dev/) format. 
The output zarr files can be uploaded into our [repeat browser](https://repeatbrowser.org/) for the visualization.


## Prerequisites

---

### 1. Compile the iteres and downloaded related files
In the top directory, run 
```bash
git clone clone Jiawei-Shen/Repeat-Browser_data_processing
cd Repeat-Browser_data_processing
make
```
And here's some related files you may need:

Repeat size file: [here](https://epigenome.wustl.edu/iteres/download/hg19/subfam.size) (length of consensus sequence of repeat subfamily)

hg38 Repeat annotation: [download from UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz)

hg38 Chromosome size file: [Full](https://epigenome.wustl.edu/iteres/download/hg38/hg38_full.size)  [Lite](https://epigenome.wustl.edu/iteres/download/hg38/hg38_lite.size) (without supercontigs)

### 2. Install the required packages of python
- [Python](https://www.python.org/) (version 3.6.13 recommended)
- [pip](https://pip.pypa.io/en/stable/installation/) (Python package installer)

```bash
# Install dependencies from requirements.txt
pip install -r requirements.txt
```

## Implementation

---

### Step 0. Align the reads

If you already have a bam file, this step can be omitted. If not, here's some tutorials to align the reads in different scenarios.

#### (1). ChIP-Seq data
We recommend users to use [BWA](https://github.com/lh3/bwa) to align the ChIP-Seq data.
```bash
# change the path to the bwa folder
./bwa index ref.fa read-se.fq.gz | gzip -3 > aln-se.sam.gz
```

#### (2). CAGE-Seq data
We recommend users to use [STAR](https://github.com/alexdobin/STAR/tree/master) to align the CAGE-Seq data.
Since we are focused on the multireads, it will have some differences from the default settings of STAR. 
```bash
STAR --chimSegmentMin 100  
    --outFilterMultimapNmax 100 
    --winAnchorMultimapNmax 100 
    --alignEndsType EndToEnd 
    --alignEndsProtrude 100 DiscordantPair 
    --outFilterScoreMinOverLread 0.4 
    --outFilterMatchNminOverLread 0.4 
    --outSAMtype BAM Unsorted 
    --outSAMattributes All 
    --outSAMstrandField intronMotif 
    --outSAMattrIHstart 0 
    --readFilesCommand zcat 
    --chimOutType WithinBAM SoftClip
```
We consult the [SQuIRE](https://github.com/wyang17/SQuIRE) repository for guidance on the STAR parameters related to handling CAGE-Seq multireads.

### Step 1. Run the pipeline
You can choose to run the bash file or run the whole pipeline by the bash file **run.sh**. 

Here are some sample files for human you may need:
```bash
--chrom_size
--subfam_size 
--rmsk_path 
```

**Repeat size file**: [here](https://epigenome.wustl.edu/iteres/download/hg19/subfam.size) (length of consensus sequence of repeat subfamily)

**hg19**: Repeat annotation: [download from UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/rmsk.txt.gz)

**hg19**: Chromosome size file: [Full](https://epigenome.wustl.edu/iteres/download/hg19/hg19_full.size)    [Lite](https://epigenome.wustl.edu/iteres/download/hg19/hg19_lite.size) (without supercontigs)

**hg38**: Repeat annotation: [download from UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/rmsk.txt.gz)

**hg38**: Chromosome size file: [Full](https://epigenome.wustl.edu/iteres/download/hg38/hg38_full.size)    [Lite](https://epigenome.wustl.edu/iteres/download/hg38/hg38_lite.size) (without supercontigs)

For the bash file:
#### (1). The default scenario
In this scenario, you only have one bam file to process. And this bam file is not from CAGE-Seq
```bash
bash run.sh --bam_file /path/to/your/bam_file 
            --output_path /path/to/output 
            --chrom_size /path/to/chrom_size_file 
            --subfam_size /path/to/subfam_size_file 
            --rmsk_path /path/to/rmsk_file 
```

#### (2). The data is from ChIP-Seq
In this case, it will have two bam files. One is the signal bam file, another is IgG control bam file.
```bash
bash run.sh --signal_bam_file /path/to/signal.bam 
            --control_bam_file /path/to/control.bam
            --output_path /path/to/output 
            --chrom_size /path/to/chrom_size_file 
            --subfam_size /path/to/subfam_size_file 
            --rmsk_path /path/to/rmsk_file 
```

#### (3). The data is from CAGE-Seq
In this case, you will have to set the length of cage_window, 
which is the length of basepairs segments around 5' end during our process. 

The default value of cage_window is 20, which meaning the segment we select is from 20 bp in front of 5' end to 20 bp behind it.

```bash
bash run.sh --bam_file /path/to/your/bam_file 
            --output_path /path/to/output 
            --chrom_size /path/to/chrom_size_file 
            --subfam_size /path/to/subfam_size_file 
            --rmsk_path /path/to/rmsk_file 
            --cage_window 20
```

We provide some sample files in the Prerequisites section. 