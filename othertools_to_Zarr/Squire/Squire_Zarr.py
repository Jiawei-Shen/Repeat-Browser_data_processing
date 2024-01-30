import zarr
import s3fs
import time
import argparse
import json
import numpy as np
import pandas as pd


# The length of wig data depends on the max length of the subfamily data. The rest part of the shorter wig data is padding -1.
def read_TEcounts_bysubfamily(filename):
    TE_set = {}
    with open(filename) as file:
        for i, line in enumerate(file):
            if i > 0:
                TE_name = line.split('\t')[11].split(':')[0]
                TE_start = int(line.split('\t')[9])
                TE_stop = int(line.split('\t')[10])
                TE_chr = line.split('\t')[8]
                FPKM = float(line.split('\t')[4])
                if TE_name in TE_set.keys():
                    TE_set[TE_name].append({"TE_start": TE_start, "TE_stop": TE_stop, "TE_chr": TE_chr, "FPKM": FPKM})
                else:
                    TE_set[TE_name] = [{"TE_start": TE_start, "TE_stop": TE_stop, "TE_chr": TE_chr, "FPKM": FPKM}]

    return TE_set


def rmsk_line_read(line):
    chr = line.split('\t')[5]
    genoStart = line.split('\t')[6]
    genoEnd = line.split('\t')[7]
    repName = line.split('\t')[10]
    repClass = line.split('\t')[11]
    strand = line.split('\t')[9]
    repEnd = line.split('\t')[14]
    if strand == '+':
        repStart = line.split('\t')[13]
        repLeft = line.split('\t')[15]
    else:
        repStart = line.split('\t')[15]
        repLeft = line.split('\t')[13]
    return chr, int(genoStart), int(genoEnd), repName, repClass, int(repStart), int(repEnd), abs(int(repLeft))


def file_lines_count(file_name):
    from itertools import (takewhile, repeat)
    buffer = 1024 * 1024
    with open(file_name) as f:
        buf_gen = takewhile(lambda x: x, (f.read(buffer) for _ in repeat(None)))
        return sum(buf.count('\n') for buf in buf_gen)


def calculate_consensus_dict(df, type):
    df_len = len(df)
    consenus_dict = {}

    for index, row in df.iterrows():
        if type == "all":
            add_number = round(row['tot_counts'])
        elif type == "uni":
            add_number = round(row['uniq_counts'])

        if row['repName'] not in consenus_dict:
            consensus_length = row['repEnd'] - row['repLeft'] if row['repLeft'] <= 0 else row['repEnd'] - row[
                'repStart']
            consenus_dict[row['repName']] = np.zeros(consensus_length)
            if row['repLeft'] <= 0:
                consenus_dict[row['repName']][row['repStart']: row['repEnd']] += add_number
            else:
                consenus_dict[row['repName']][row['repLeft']: row['repEnd']] += add_number
        else:
            if row['repLeft'] <= 0:
                consenus_dict[row['repName']][row['repStart']: row['repEnd']] += add_number
            else:
                consenus_dict[row['repName']][row['repLeft']: row['repEnd']] += add_number
        print(f"Processing {type} data... {round(index/df_len, 3)*100}%", end='\r', flush=True)
    return consenus_dict


def cal_TEcounts_consensus(filename, rmsk):
    start_time = time.time()
    chr_list = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX",
                "chrY"]
    consensus_map_all = {}
    consensus_map_unique = {}
    df = pd.read_table(filename, '\t')
    # Filter rows where 'tx_chr' is not in the allowed_chr list
    df = df[df['tx_chr'].isin(chr_list)]

    # Define column names based on the provided information
    # rmsk_column_names = [
    #     'bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns',
    #     'chr', 'genoStart', 'genoEnd', 'genoLeft', 'strand',
    #     'repName', 'repClass', 'repFamily', 'repStart', 'repEnd',
    #     'repLeft', 'id'
    # ]
    rmsk_column_names = [
        'bin', 'swScore', 'milliDiv', 'milliDel', 'milliIns',
        'TE_chr', 'TE_start', 'TE_stop', 'genoLeft', 'strand',
        'repName', 'repClass', 'repFamily', 'repStart', 'repEnd',
        'repLeft', 'id'
    ]
    # Read the rmsk file into a DataFrame
    df_rmsk = pd.read_table(rmsk, sep='\t', header=None, names=rmsk_column_names)
    df_rmsk_demo = df_rmsk.head(10000)

    # Merge DataFrames based on the specified conditions
    # merged_table = pd.concat([df, df_rmsk_demo[(df['TE_chr'] == df_rmsk_demo['chr']) &
    #                                            (df['TE_start'] >= df_rmsk_demo['genoStart']) &
    #                                            (df['TE_stop'] <= df_rmsk_demo['genoEnd'])]])
    merged_df = pd.merge(df, df_rmsk, on=['TE_chr', 'TE_start', 'TE_stop'])

    consenus_dict_uni = calculate_consensus_dict(merged_df, 'uni')
    consenus_dict_all = calculate_consensus_dict(merged_df, 'all')

    return consenus_dict_all, consenus_dict_uni


def read_rmsk_bysubfamily(filename):
    rep_set = {}
    with open(filename) as file:
        for i, line in enumerate(file):
            rep_name = line.split('\t')[10]
            rep_chr = line.split('\t')[5]
            geno_start = int(line.split('\t')[6])
            geno_end = int(line.split('\t')[7])
            if line.split('\t')[9] == '-':
                rep_start = int(line.split('\t')[15])
            else:
                rep_start = int(line.split('\t')[13])
            rep_end = int(line.split('\t')[14])
            rep_left = int(line.split('\t')[15])

            if rep_name in rep_set.keys():
                if rep_chr in rep_set[rep_name].keys():
                    rep_set[rep_name][rep_chr].append({"rep_start": rep_start, "rep_end": rep_end, "rep_left": rep_left,
                                                       "geno_start": geno_start, "geno_end": geno_end})
                else:
                    rep_set[rep_name][rep_chr] = [{"rep_start": rep_start, "rep_end": rep_end, "rep_left": rep_left,
                                                   "geno_start": geno_start, "geno_end": geno_end}]
            else:
                rep_set[rep_name] = {rep_chr: [{"rep_start": rep_start, "rep_end": rep_end, "rep_left": rep_left,
                                                "geno_start": geno_start, "geno_end": geno_end}]}

    return rep_set


def read_consensus_length(subfamily_size_filepath):
    consensus_len = {}
    with open(subfamily_size_filepath) as file:
        for i, line in enumerate(file):
            if i > 0:
                subfam_name = line.split('\t')[0]
                consensus_len[subfam_name] = int(line.split('\t')[1])
    return consensus_len


def heatmap_data(subFcounts_filepath, consensus_length_data, TE_Set):
    all_tmp_set = {}
    uni_tmp_set = {}
    with open(subFcounts_filepath) as file:
        uni_read_counts = 0
        all_read_counts = 0
        for i, line in enumerate(file):
            if i > 0:
                subfamily_name = line.split('\t')[2].split(':')[0]
                if subfamily_name in consensus_length_data.keys():
                    TE_segment = TE_Set[subfamily_name]
                    subfamily_len = sum(list(map(lambda x: x['TE_stop'] - x['TE_start'], TE_segment)))
                    multi_counts = float(line.split('\t')[7])
                    unique_counts = int(line.split('\t')[5])
                    print('{} consensus len is {}.'.format(subfamily_name, subfamily_len))
                    all_tmp_set[subfamily_name] = (multi_counts + unique_counts) / subfamily_len * 1e6
                    uni_tmp_set[subfamily_name] = unique_counts / subfamily_len * 1e6
                    all_read_counts += multi_counts + unique_counts
                    uni_read_counts += unique_counts
                    # read_counts = int(line.split('\t')[7])
                    # all_FPKM_set = {key: value / read_counts for key, value in all_tmp_set.items()}
                    # uni_FPKM_set = {key: value / read_counts for key, value in uni_tmp_set.items()}
                else:
                    print("{} Not Detected!!!".format(subfamily_name), end='\r', flush=True)
                    # FPKM_set[subfamily_name] = float(line.split('\t')[4])
        all_FPKM_set = {key: value / (all_read_counts / 1e3) for key, value in all_tmp_set.items()}
        uni_FPKM_set = {key: value / (uni_read_counts / 1e3) for key, value in uni_tmp_set.items()}

    # with open('1101test_36.txt', 'w') as f:
    #     subfam_list = ['LTR12', 'LTR12E', 'LTR12D', 'LTR12C', 'THE1C', 'THE1B', 'THE1D', 'LTR12B', 'LTR12F']
    #     for key, value in all_FPKM_set.items():
    #         if key in subfam_list:
    #             f.writelines('key: {}, value: {}.\n'.format(key, value))

    return all_FPKM_set, uni_FPKM_set


def genomview_loci_data(TE):
    data_dict = {}
    for subfamily in TE.keys():
        for TE_info in TE[subfamily]:
            # data format is [chrome, Start, End, RPKM]
            data = [TE_info["TE_chr"], TE_info["TE_start"], TE_info["TE_stop"], TE_info["FPKM"]]
            if subfamily in data_dict:
                data_dict[subfamily].append(data)
            else:
                data_dict[subfamily] = [data]
    return data_dict


def create_RNA_zarr(
        TEcounts_filepath='./data/jiawei_1101_zarr_input/SRR3498333_squire_count/SRR3498333.fastq_TEcounts.txt',
        rmsk_filepath='./data/RNA-seq pipeline/gene_annotation/hg38_rmsk.txt',
        subfamily_size_filepath='./data/RNA-seq pipeline/gene_annotation/hg38TEsubF.size.txt',
        subFcounts_filepath='./data/jiawei_1101_zarr_input/SRR3498333_squire_count/SRR3498333.fastq_subFcounts.txt',
        # consensus_json_all_filepath='consensusView_all.json', consensus_json_uni_filepath='consensusView_unique.json',
        parameters={}, local=True, local_path='./data/example0126.zarr',
        s3_path='s3://repeatbrowsers/![](../../Downloads/3218908-200.png)Zarr_data/example_RNA.zarr'):

    consensus_json_all_set, consensus_json_uni_set = cal_TEcounts_consensus(TEcounts_filepath, rmsk_filepath)
    TE_set = read_TEcounts_bysubfamily(TEcounts_filepath)
    rmsk = read_rmsk_bysubfamily(rmsk_filepath)
    consensus_length_set = read_consensus_length(subfamily_size_filepath)
    TE_consensus_set = {}
    not_in_subfam = []
    for TE_name in TE_set.keys():
        if TE_name not in consensus_length_set.keys():
            # print("{} Not in!!!".format(TE_name))
            not_in_subfam.append(TE_name)
        else:
            TE_consensus_set[TE_name] = consensus_length_set[TE_name]

    if not_in_subfam:
        [TE_set.pop(item) for item in not_in_subfam]

    all_FPKM_data, uni_FPKM_data = heatmap_data(subFcounts_filepath, TE_consensus_set, TE_set)
    # consensus_set = consensus_map(rmsk, TE_set, TE_consensus_set)

    # with open(consensus_json_all_filepath, 'r') as f:
    #     consensus_json_all_set = json.load(f)
    # with open(consensus_json_uni_filepath, 'r') as f:
    #     consensus_json_uni_set = json.load(f)

    loci_dict = genomview_loci_data(TE_set)
    # subfam_range = set()
    # with open(subfamily_size_filepath) as file:
    #     for i, line in enumerate(file):
    #         if i > 0:
    #             sf = line.split('\t')[0]
    #             subfam_range.add(sf)
    # loci_dict = {k: v for k, v in org_loci_dict.items() if k in subfam_range}
    loci_chunk_sequence = ['repname', 'chrome', 'start', 'end', 'FPKM']

    subfamily = [key for key in TE_set.keys()]
    subfamily.insert(0, '#subfamily')
    all_reads_FPKM = [value for value in all_FPKM_data.values()]
    all_reads_FPKM.insert(0, 'all_reads_FPKM')
    uni_reads_FPKM = [value for value in uni_FPKM_data.values()]
    uni_reads_FPKM.insert(0, 'uni_reads_FPKM')

    stat_subfamily = np.array([subfamily, all_reads_FPKM, uni_reads_FPKM])
    all_wig_data = []
    uni_wig_data = []
    max_len = max(list(map(len, consensus_json_all_set.values())))
    for key in consensus_json_all_set:
        d_all = np.array(consensus_json_all_set[key], dtype='i4')
        d_uni = np.array(consensus_json_uni_set[key], dtype='i4')
        a_all = np.pad(d_all, (0, max_len - len(d_all)), 'constant', constant_values=(0, -1))
        a_uni = np.pad(d_uni, (0, max_len - len(d_uni)), 'constant', constant_values=(0, -1))
        all_wig_data.append(a_all)
        uni_wig_data.append(a_uni)
    uni_wig_data = np.array(uni_wig_data)
    all_wig_data = np.array(all_wig_data)
    all_wig_keys = list(consensus_json_all_set.keys())
    uni_wig_keys = list(consensus_json_uni_set.keys())

    print('Building the Zarr format!')
    if local:
        root = zarr.open(local_path, mode='w')
    else:
        # s3_path = 's3://repeatbrowsers/Zarr_data/example.zarr'
        # Initilize the S3 file system
        s3 = s3fs.S3FileSystem()
        store = s3fs.S3Map(root=s3_path, s3=s3, check=False)
        root = zarr.group(store=store, overwrite=True)

    subfam_attrs = subfamily[1:]
    root.attrs['subfamily_stat'] = subfam_attrs
    subfam_data = np.array([all_reads_FPKM[1:], uni_reads_FPKM[1:]], dtype='float')
    subfam_stat = root.create_dataset('subfam_stat', data=subfam_data, compressor=zarr.Zlib(level=1),
                                      chunks=(1, stat_subfamily.shape[1]))

    uni_bigwig = root.create_dataset('uni_bigwig', data=uni_wig_data, compressor=zarr.Zlib(level=1),
                                     chunks=(uni_wig_data.shape[0], uni_wig_data.shape[1]))
    all_bigwig = root.create_dataset('all_bigwig', data=all_wig_data, compressor=zarr.Zlib(level=1),
                                     chunks=(all_wig_data.shape[0], all_wig_data.shape[1]))
    root.attrs['uni_wig'] = uni_wig_keys
    root.attrs['all_wig'] = all_wig_keys

    for key in loci_dict:
        data = np.array(loci_dict[key]).flatten()
        root.create_dataset('loci_{}'.format(key), data=data, compressor=zarr.Zlib(level=1), chunks=(len(data) or 1))
    root.attrs['loci'] = loci_chunk_sequence

    if parameters['id'] == 'unknown':
        parameters['id'] = TEcounts_filepath.split('/')[-1].split('.')[0].split('_')[0]
    if parameters['Assay'] == 'unknown':
        parameters['Assay'] = 'RNAseq'
    root.attrs['Parameters'] = parameters

    if local:
        zarr.consolidate_metadata(local_path)
    else:
        zarr.consolidate_metadata(store)

    return root


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process the SQuIRE output files to the .Zarr format.')

    parser.add_argument('--Biosample', '-b', default='unknown', type=str, help="The biosample type of your input file.")
    parser.add_argument('--Tissue', '-t', default='unknown', type=str, help="The tissue type of your input file.")
    parser.add_argument('--Assay', '-a', default='unknown', type=str, help="The assay type of your input file.")
    parser.add_argument('--Control', '-c', default='unknown', type=str, help="The control id of your input file.")
    parser.add_argument('--Organism', default='unknown', type=str, help="The organism of your input file.")
    parser.add_argument('--id', default='unknown', type=str, help="The file id (accession) of your input file.")

    parser.add_argument('--subFcounts',
                        default='./SRR3498327_squire_count/SRR3498327.fastq_subFcounts.txt',
                        type=str, help="The path of subFcounts file (All) output by SQUIRE.")
    parser.add_argument('--TEcounts',
                        default='./SRR3498327_squire_count/SRR3498327.fastq_TEcounts.txt',
                        type=str, help="The path of TEcounts file (All) output by SQUIRE.")

    parser.add_argument('--subFsize',
                        default='./RNA-seq pipeline/gene_annotation/hg38TEsubF.size.txt', type=str,
                        help="The path of subfamily size file (All) output by SQUIRE.")
    parser.add_argument('--rmsk', default='./RNA-seq pipeline/gene_annotation/hg38_rmsk.txt',
                        type=str, help="The file path of rmsk file.")
    parser.add_argument('--output', '-o', default='./example.zarr', type=str,
                        help="The local path of output files, which must end with *.zarr")

    args = parser.parse_args()
    parameters = {"Biosample": args.Biosample, "Assay": args.Assay, "Control": args.Control, "Organism": args.Organism,
                  "id": args.id, "Tissue": args.Tissue}

    root = create_RNA_zarr(TEcounts_filepath=args.TEcounts,
                           rmsk_filepath=args.rmsk,
                           subfamily_size_filepath=args.subFsize,
                           subFcounts_filepath=args.subFcounts,
                           # consensus_json_all_filepath=args.consensus_json_all_filepath,
                           # consensus_json_uni_filepath=args.consensus_json_uni_filepath,
                           parameters=parameters, local_path=args.output)
    print(root.info)