# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.
import zarr
import s3fs
import time
import argparse
import json
import numpy as np


def read_statfile_byline(filename):
    data = []
    subfamily = []
    all_reads_RPKM = []
    uni_reads_RPKM = []

    with open(filename) as file:
        for line in file:
            data.append(np.array(line.replace('\n', '').split('\t')))
            subfamily.append(line.replace('\n', '').split('\t')[0])
            all_reads_RPKM.append(line.replace('\n', '').split('\t')[8])  # subfam_all_reads_RPKM
            uni_reads_RPKM.append(line.replace('\n', '').split('\t')[10])  # subfam_unique_reads_RPKM

    return_array = np.array([subfamily, all_reads_RPKM, uni_reads_RPKM])
    return return_array


def read_locifile_byline(filename='./testAligned.sortedByCoord.out_ALL.iteres.loci'):
    repname = []
    RPKM = []
    start = []
    end = []
    chrome = []

    with open(filename) as file:
        for i, line in enumerate(file):
            if i > 0:
                chrome.append(line.replace('\n', '').split('\t')[0])
                start.append(int(line.replace('\n', '').split('\t')[1]))
                end.append(int(line.replace('\n', '').split('\t')[2]))
                repname.append(line.replace('\n', '').split('\t')[4])
                RPKM.append(float(line.replace('\n', '').split('\t')[8]))

    return_array = np.array([repname, chrome, start, end, RPKM])
    return return_array


def read_locifile_bysubfamily(filename='./testAligned.sortedByCoord.out_ALL.iteres.loci',
                              range='/Users/jiaweishen/PycharmProjects/Zarr_format/data/RNA-seq pipeline/gene_annotation/hg38TEsubF.size.txt'):
    data_dict = {}

    with open(filename) as file:
        for i, line in enumerate(file):
            if i > 0:
                subfam_name = line.replace('\n', '').split('\t')[4]
                if subfam_name in data_dict:
                    # data format is [chrome, Start, End, RPKM]
                    data = [line.replace('\n', '').split('\t')[0], int(line.replace('\n', '').split('\t')[1]),
                            int(line.replace('\n', '').split('\t')[2]), float(line.replace('\n', '').split('\t')[8])]
                    data_dict[subfam_name].append(data)
                else:
                    data_dict[subfam_name] = []

    return data_dict


def read_wigfile_bysubfamily(filename):
    subfam_data = []
    subfam_name = ''
    data = {}
    with open(filename) as file:
        for line in file:
            line_data = np.array(line.replace('\n', '').split(' '))
            if len(line_data) > 2:
                subfam_name = line_data[1].split('=')[1]
                data[subfam_name] = []
            else:
                data[subfam_name].append(int(line_data[0]))

    max_len = max(list(map(len, data.values())))
    for key in data:
        d = np.array(data[key], dtype='i4')
        a = np.pad(d, (0, max_len - len(d)), 'constant', constant_values=(0, -1))
        subfam_data.append(a)

    return list(data.keys()), np.array(subfam_data)


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


# rep_set: {'repeat_name':
#               {'chr': [{"rep_start": rep_start, "rep_end": rep_end, "rep_left": rep_left,
#                           "geno_start": geno_start, "geno_end": geno_end},
#                           ...]}}


def read_consensus_length(subfamily_size_filepath):
    consensus_len = {}
    with open(subfamily_size_filepath) as file:
        for i, line in enumerate(file):
            if i > 0:
                subfam_name = line.split('\t')[0]
                consensus_len[subfam_name] = int(line.split('\t')[1])
    return consensus_len


def consensus_map(rmsk, TE, consensus_set):
    TE_consensus_set = {}
    for TE_name in TE.keys():
        if TE_name not in consensus_set.keys():
            print("{} Not in!!!".format(TE_name))
        else:
            TE_consensus_set[TE_name] = np.zeros(consensus_set[TE_name])
    time_start = time.time()
    print('consensus mapping start time counting!')
    for i, subfamily in enumerate(TE.keys()):
        print('NOW: {}/{}.'.format(i, len(TE)))
        if subfamily in TE_consensus_set.keys():
            for TE_info in TE[subfamily]:
                chromosome = TE_info["TE_chr"]
                for rmsk_info in rmsk[subfamily][chromosome]:
                    if TE_info["TE_start"] == rmsk_info["geno_start"] and TE_info["TE_stop"] == rmsk_info["geno_end"]:
                        te_start = rmsk_info["rep_start"]
                        te_end = rmsk_info["rep_end"] if rmsk_info["rep_end"] < len(
                            TE_consensus_set[subfamily]) else len(
                            TE_consensus_set[subfamily])
                        TE_consensus_set[subfamily][te_start:te_end] += 1
    time_end = time.time()
    print('consensus mapping time cost: ', time_end - time_start, 's')
    return TE_consensus_set


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


# def heatmap_data(subFcounts_filepath, consensus_length_data):
#     all_tmp_set = {}
#     uni_tmp_set = {}
#     FPKM_set = {}
#     with open(subFcounts_filepath) as file:
#         read_counts = 0
#         for i, line in enumerate(file):
#             if i > 0:
#                 subfamily_name = line.split('\t')[2].split(':')[0]
#                 if subfamily_name in consensus_length_data.keys():
#                     subfamily_len = consensus_length_data[subfamily_name]
#                     total_counts = float(line.split('\t')[6])
#                     unique_counts = int(line.split('\t')[5])
#                     print('{} consensus len is {}.'.format(subfamily_name, subfamily_len))
#                     all_tmp_set[subfamily_name] = total_counts / subfamily_len * 1e6
#                     uni_tmp_set[subfamily_name] = unique_counts / subfamily_len * 1e6
#                     read_counts += int(line.split('\t')[7])
#                     # read_counts = int(line.split('\t')[7])
#                     # all_FPKM_set = {key: value / read_counts for key, value in all_tmp_set.items()}
#                     # uni_FPKM_set = {key: value / read_counts for key, value in uni_tmp_set.items()}
#                 else:
#                     print("{} Not in!!!".format(subfamily_name))
#                     # FPKM_set[subfamily_name] = float(line.split('\t')[4])
#         all_FPKM_set = {key: value / read_counts for key, value in all_tmp_set.items()}
#         uni_FPKM_set = {key: value / read_counts for key, value in uni_tmp_set.items()}
#     return all_FPKM_set, uni_FPKM_set

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
                    print("{} Not in!!!".format(subfamily_name))
                    # FPKM_set[subfamily_name] = float(line.split('\t')[4])
        all_FPKM_set = {key: value / (all_read_counts / 1e3) for key, value in all_tmp_set.items()}
        uni_FPKM_set = {key: value / (uni_read_counts / 1e3) for key, value in uni_tmp_set.items()}

    # with open('1101test_36.txt', 'w') as f:
    #     subfam_list = ['LTR12', 'LTR12E', 'LTR12D', 'LTR12C', 'THE1C', 'THE1B', 'THE1D', 'LTR12B', 'LTR12F']
    #     for key, value in all_FPKM_set.items():
    #         if key in subfam_list:
    #             f.writelines('key: {}, value: {}.\n'.format(key, value))

    return all_FPKM_set, uni_FPKM_set


def calculate_ratio(signal_stat_subfamily, control_stat_subfamily):
    marker_ALL_RPKM = dict(zip(signal_stat_subfamily[0][1:], signal_stat_subfamily[1][1:]))
    marker_UNI_RPKM = dict(zip(signal_stat_subfamily[0][1:], signal_stat_subfamily[2][1:]))

    control_ALL_RPKM = dict(zip(control_stat_subfamily[0][1:], control_stat_subfamily[1][1:]))
    control_UNI_RPKM = dict(zip(control_stat_subfamily[0][1:], control_stat_subfamily[2][1:]))

    value_ALL_RPKM = {}
    value_UNI_RPKM = {}
    for subfam_name in marker_ALL_RPKM:
        if subfam_name in control_ALL_RPKM.keys():
            value_ALL_RPKM[subfam_name] = float(marker_ALL_RPKM[subfam_name]) if float(
                control_ALL_RPKM[subfam_name]) == 0 else float(marker_ALL_RPKM[subfam_name]) / float(
                control_ALL_RPKM[subfam_name])
        else:
            value_ALL_RPKM[subfam_name] = float(marker_ALL_RPKM[subfam_name])

    for subfam_name in marker_UNI_RPKM:
        if subfam_name in control_UNI_RPKM.keys():
            value_UNI_RPKM[subfam_name] = float(marker_UNI_RPKM[subfam_name]) if float(
                control_UNI_RPKM[subfam_name]) == 0 else float(marker_UNI_RPKM[subfam_name]) / float(
                control_UNI_RPKM[subfam_name])
        else:
            value_UNI_RPKM[subfam_name] = float(marker_UNI_RPKM[subfam_name])

    value_stat_subfamily = [['#subfamily'], ['all_reads_RPKM'], ['unique_reads_RPKM']]
    value_stat_subfamily[0].extend(list(value_ALL_RPKM.keys()))
    value_stat_subfamily[1].extend(list(value_ALL_RPKM.values()))
    value_stat_subfamily[2].extend(list(value_UNI_RPKM.values()))

    return np.array(value_stat_subfamily)


# def parse_parameters(type, file_path):
#
#     if type == 'RNA':
#         if ".json" in file_path.split('/')[-1]:
#             ...
#         elif ".txt" in file_path.split('/')[-1]:
#             with open(file_path) as file:
#                 read_counts = 0
#                 paras = {}
#                 for i, line in enumerate(file):
#                         TEcounts_filepath = line.split('\t')[0],
#                         rmsk_filepath = line.split('\t')[1],
#                         subfamily_size_filepath = line.split('\t')[2],
#                         subFcounts_filepath = line.split('\t')[3],
#                         subfamily_name = line.split('\t')[2].split(':')[0]
#
#     elif type == 'EXPERIMENT':
#         ...


def create_iteres_zarr(
        subfam_stat_path='default_with_unmapped_bam_and_readsAligned.sortedByCoord.out.iteres.subfamily.stat',
        all_wig_path='out_all.wig', uni_wig_path='out_unique.wig',
        loci_path="testAligned.sortedByCoord.out_ALL.iteres.loci",
        parameters={}, s3_path='s3://repeatbrowsers/Zarr_data/example.zarr',
        range='/Users/jiaweishen/PycharmProjects/Zarr_format/data/RNA-seq pipeline/gene_annotation/hg38TEsubF.size.txt',
        local=False, local_path='./data/example.zarr'):
    # stat_family = read_statfile_byline(filename=
    # 'default_with_unmapped_bam_and_readsAligned.sortedByCoord.out.iteres.family.stat')
    stat_subfamily = read_statfile_byline(
        filename=subfam_stat_path)
    # stat_class = read_statfile_byline(
    #     filename='default_with_unmapped_bam_and_readsAligned.sortedByCoord.out.iteres.class.stat')

    all_keys, all_wig = read_wigfile_bysubfamily(all_wig_path)
    uni_keys, uni_wig = read_wigfile_bysubfamily(uni_wig_path)

    # loci_data = read_locifile_byline('testAligned.sortedByCoord.out_ALL.iteres.loci')
    org_loci_dict = read_locifile_bysubfamily(loci_path)
    subfam_range = set()
    with open(range) as file:
        for i, line in enumerate(file):
            if i > 0:
                sf = line.split('\t')[0]
                subfam_range.add(sf)

    loci_dict = {k: v for k, v in org_loci_dict.items() if k in subfam_range}
    loci_chunk_sequence = ['repname', 'chrome', 'start', 'end', 'RPKM']

    time_start = time.time()
    print('start iteres time counting!')
    if local:
        root = zarr.open(local_path, mode='w')
    else:
        # s3_path = 's3://repeatbrowsers/Zarr_data/example.zarr'
        # Initilize the S3 file system
        s3 = s3fs.S3FileSystem()
        store = s3fs.S3Map(root=s3_path, s3=s3, check=False)
        root = zarr.group(store=store, overwrite=True)

    time_end = time.time()
    print('time cost: ', time_end - time_start, 's')

    subfam_attrs = list(stat_subfamily[0][1:])
    root.attrs['subfamily_stat'] = subfam_attrs
    subfam_data = np.array([stat_subfamily[1][1:], stat_subfamily[2][1:]], dtype='float')
    subfam_stat = root.create_dataset('subfam_stat', data=subfam_data, compressor=zarr.Zlib(level=1),
                                      chunks=(1, stat_subfamily.shape[1]))

    # fam_stat = stat.create_dataset('fam_stat', data=stat_family)
    # class_stat = stat.create_dataset('class_stat', data=stat_class)

    time_end = time.time()
    print('time cost: ', time_end - time_start, 's')

    uni_bigwig = root.create_dataset('uni_bigwig', data=uni_wig, compressor=zarr.Zlib(level=1),
                                     chunks=(uni_wig.shape[0], uni_wig.shape[1]))
    all_bigwig = root.create_dataset('all_bigwig', data=all_wig, compressor=zarr.Zlib(level=1),
                                     chunks=(all_wig.shape[0], all_wig.shape[1]))
    root.attrs['uni_wig'] = uni_keys
    root.attrs['all_wig'] = all_keys

    for key in loci_dict:
        data = np.array(loci_dict[key]).flatten()
        root.create_dataset('loci_{}'.format(key), data=data, compressor=zarr.Zlib(level=1), chunks=(len(data) or 1))
    # loci = root.create_dataset('loci', data=loci_data, compressor=zarr.Zlib(level=1), chunks=(1, loci_data.shape[1]))
    root.attrs['loci'] = loci_chunk_sequence

    # uni_bigwig = bigwig.create_group('uni_bigwig')
    # for uni_key, uni_value in uni_wig.items():
    #     uni_bigwig.create_dataset('{}'.format(uni_key), data=np.array(uni_value))
    #
    # all_bigwig = bigwig.create_group('all_bigwig')
    # for all_key, all_value in all_wig.items():
    #     all_bigwig.create_dataset('{}'.format(all_key), data=np.array(all_value))
    if parameters['id'] == 'unknown':
        parameters['id'] = subfam_stat_path.split('/')[-1].split('.')[0]
    root.attrs['Parameters'] = parameters

    if local:
        zarr.consolidate_metadata(local_path)
    else:
        zarr.consolidate_metadata(store)

    time_end = time.time()
    print('time cost: ', time_end - time_start, 's')

    return root


def create_iteres_zarr_cage(
        subfam_stat_path='default_with_unmapped_bam_and_readsAligned.sortedByCoord.out.iteres.subfamily.stat',
        all_wig_path_plus='out_all.wig', all_wig_path_minus='out_all.wig',
        uni_wig_path_plus='out_unique.wig', uni_wig_path_minus='out_unique.wig',
        loci_path="testAligned.sortedByCoord.out_ALL.iteres.loci",
        parameters={}, s3_path='s3://repeatbrowsers/Zarr_data/example.zarr',
        range='/Users/jiaweishen/PycharmProjects/Zarr_format/data/RNA-seq pipeline/gene_annotation/hg38TEsubF.size.txt',
        local=False, local_path='./data/example.zarr'):
    stat_subfamily = read_statfile_byline(filename=subfam_stat_path)
    all_plus_keys, all_plus_wig = read_wigfile_bysubfamily(all_wig_path_plus)
    all_minus_keys, all_minus_wig = read_wigfile_bysubfamily(all_wig_path_minus)
    uni_plus_keys, uni_plus_wig = read_wigfile_bysubfamily(uni_wig_path_plus)
    uni_minus_keys, uni_minus_wig = read_wigfile_bysubfamily(uni_wig_path_minus)

    # loci_data = read_locifile_byline('testAligned.sortedByCoord.out_ALL.iteres.loci')
    org_loci_dict = read_locifile_bysubfamily(loci_path)
    subfam_range = set()
    with open(range) as file:
        for i, line in enumerate(file):
            if i > 0:
                sf = line.split('\t')[0]
                subfam_range.add(sf)

    loci_dict = {k: v for k, v in org_loci_dict.items() if k in subfam_range}
    loci_chunk_sequence = ['repname', 'chrome', 'start', 'end', 'RPKM']

    time_start = time.time()
    print('start iteres time counting!')
    if local:
        root = zarr.open(local_path, mode='w')
    else:
        # s3_path = 's3://repeatbrowsers/Zarr_data/example.zarr'
        # Initilize the S3 file system
        s3 = s3fs.S3FileSystem()
        store = s3fs.S3Map(root=s3_path, s3=s3, check=False)
        root = zarr.group(store=store, overwrite=True)

    time_end = time.time()
    print('time cost: ', time_end - time_start, 's')

    subfam_attrs = list(stat_subfamily[0][1:])
    root.attrs['subfamily_stat'] = subfam_attrs
    subfam_data = np.array([stat_subfamily[1][1:], stat_subfamily[2][1:]], dtype='float')
    subfam_stat = root.create_dataset('subfam_stat', data=subfam_data, compressor=zarr.Zlib(level=1),
                                      chunks=(1, stat_subfamily.shape[1]))

    # fam_stat = stat.create_dataset('fam_stat', data=stat_family)
    # class_stat = stat.create_dataset('class_stat', data=stat_class)

    time_end = time.time()
    print('time cost: ', time_end - time_start, 's')

    uni_bigwig_plus = root.create_dataset('uni+_bigwig', data=uni_plus_wig, compressor=zarr.Zlib(level=1),
                                          chunks=(uni_plus_wig.shape[0], uni_plus_wig.shape[1]))
    uni_bigwig_minus = root.create_dataset('uni-_bigwig', data=uni_minus_wig, compressor=zarr.Zlib(level=1),
                                           chunks=(uni_minus_wig.shape[0], uni_minus_wig.shape[1]))
    all_bigwig_plus = root.create_dataset('all+_bigwig', data=all_plus_wig, compressor=zarr.Zlib(level=1),
                                          chunks=(all_plus_wig.shape[0], all_plus_wig.shape[1]))
    all_bigwig_minus = root.create_dataset('all-_bigwig', data=all_minus_wig, compressor=zarr.Zlib(level=1),
                                           chunks=(all_minus_wig.shape[0], all_minus_wig.shape[1]))
    root.attrs['uni+_wig'] = uni_plus_keys
    root.attrs['uni-_wig'] = uni_minus_keys
    root.attrs['all+_wig'] = all_plus_keys
    root.attrs['all-_wig'] = all_minus_keys

    for key in loci_dict:
        data = np.array(loci_dict[key]).flatten()
        root.create_dataset('loci_{}'.format(key), data=data, compressor=zarr.Zlib(level=1), chunks=(len(data) or 1))
    # loci = root.create_dataset('loci', data=loci_data, compressor=zarr.Zlib(level=1), chunks=(1, loci_data.shape[1]))
    root.attrs['loci'] = loci_chunk_sequence

    # uni_bigwig = bigwig.create_group('uni_bigwig')
    # for uni_key, uni_value in uni_wig.items():
    #     uni_bigwig.create_dataset('{}'.format(uni_key), data=np.array(uni_value))
    #
    # all_bigwig = bigwig.create_group('all_bigwig')
    # for all_key, all_value in all_wig.items():
    #     all_bigwig.create_dataset('{}'.format(all_key), data=np.array(all_value))
    if parameters['id'] == 'unknown':
        parameters['id'] = subfam_stat_path.split('/')[-1].split('.')[0]
    if parameters['Assay'] == 'unknown':
        parameters['Assay'] = 'RNA-seq'
    root.attrs['Parameters'] = parameters

    if local:
        zarr.consolidate_metadata(local_path)
    else:
        zarr.consolidate_metadata(store)

    time_end = time.time()
    print('time cost: ', time_end - time_start, 's')

    return root


# def create_RNA_zarr(TEcounts_filepath='./data/jiawei_1101_zarr_input/SRR3498333_squire_count/SRR3498333.fastq_TEcounts.txt',
#                     rmsk_filepath='./data/RNA-seq pipeline/gene_annotation/hg38_rmsk.txt',
#                     subfamily_size_filepath='./data/RNA-seq pipeline/gene_annotation/hg38TEsubF.size.txt',
#                     subFcounts_filepath='./data/jiawei_1101_zarr_input/SRR3498333_squire_count/SRR3498333.fastq_subFcounts.txt',
#                     consensus_json_all_filepath='consensusView_all.json', consensus_json_uni_filepath='consensusView_all.json',
#                     parameters={}, local=True, local_path='./data/example.zarr',
#                     s3_path='s3://repeatbrowsers/![](../../Downloads/3218908-200.png)Zarr_data/example_RNA.zarr'):
#     TE_set = read_TEcounts_bysubfamily(TEcounts_filepath)
#     rmsk = read_rmsk_bysubfamily(rmsk_filepath)
#     consensus_length_set = read_consensus_length(subfamily_size_filepath)
#     TE_consensus_set = {}
#     not_in_subfam = []
#     for TE_name in TE_set.keys():
#         if TE_name not in consensus_length_set.keys():
#             print("{} Not in!!!".format(TE_name))
#             not_in_subfam.append(TE_name)
#         else:
#             TE_consensus_set[TE_name] = consensus_length_set[TE_name]
#
#     if not_in_subfam:
#         [TE_set.pop(item) for item in not_in_subfam]
#
#     all_FPKM_data, uni_FPKM_data = heatmap_data(subFcounts_filepath, TE_consensus_set, TE_set)
#     consensus_set = consensus_map(rmsk, TE_set, TE_consensus_set)
#
#     with open(consensus_json_all_filepath, 'r') as f:
#         consensus_json_all_set = json.load(f)
#     with open(consensus_json_uni_filepath, 'r') as f:
#         consensus_json_uni_set = json.load(f)
#
#     loci_dict = genomview_loci_data(TE_set)
#     # subfam_range = set()
#     # with open(subfamily_size_filepath) as file:
#     #     for i, line in enumerate(file):
#     #         if i > 0:
#     #             sf = line.split('\t')[0]
#     #             subfam_range.add(sf)
#     # loci_dict = {k: v for k, v in org_loci_dict.items() if k in subfam_range}
#     loci_chunk_sequence = ['repname', 'chrome', 'start', 'end', 'FPKM']
#
#     subfamily = [key for key in TE_set.keys()]
#     subfamily.insert(0, '#subfamily')
#     all_reads_FPKM = [value for value in all_FPKM_data.values()]
#     all_reads_FPKM.insert(0, 'all_reads_FPKM')
#     uni_reads_FPKM = [value for value in uni_FPKM_data.values()]
#     uni_reads_FPKM.insert(0, 'uni_reads_FPKM')
#
#     stat_subfamily = np.array([subfamily, all_reads_FPKM, uni_reads_FPKM])
#     all_wig_data = []
#     uni_wig_data = []
#     max_len = max(list(map(len, consensus_set.values())))
#     for key in consensus_set:
#         d = np.array(consensus_set[key], dtype='i4')
#         a = np.pad(d, (0, max_len - len(d)), 'constant', constant_values=(0, -1))
#         all_wig_data.append(a)
#         uni_wig_data.append(a * 0)
#     uni_wig_data = np.array(uni_wig_data)
#     all_wig_data = np.array(all_wig_data)
#     all_wig_keys = list(consensus_set.keys())
#     uni_wig_keys = list(consensus_set.keys())
#
#     time_start = time.time()
#     print('start RNA pipeline time counting!')
#     if local:
#         root = zarr.open(local_path, mode='w')
#     else:
#         # s3_path = 's3://repeatbrowsers/Zarr_data/example.zarr'
#         # Initilize the S3 file system
#         s3 = s3fs.S3FileSystem()
#         store = s3fs.S3Map(root=s3_path, s3=s3, check=False)
#         root = zarr.group(store=store, overwrite=True)
#
#     time_end = time.time()
#     print('time cost: ', time_end - time_start, 's')
#
#     subfam_attrs = subfamily[1:]
#     root.attrs['subfamily_stat'] = subfam_attrs
#     subfam_data = np.array([all_reads_FPKM[1:], uni_reads_FPKM[1:]], dtype='float')
#     subfam_stat = root.create_dataset('subfam_stat', data=subfam_data, compressor=zarr.Zlib(level=1),
#                                       chunks=(1, stat_subfamily.shape[1]))
#
#     time_end = time.time()
#     print('time cost: ', time_end - time_start, 's')
#
#     uni_bigwig = root.create_dataset('uni_bigwig', data=uni_wig_data, compressor=zarr.Zlib(level=1),
#                                      chunks=(uni_wig_data.shape[0], uni_wig_data.shape[1]))
#     all_bigwig = root.create_dataset('all_bigwig', data=all_wig_data, compressor=zarr.Zlib(level=1),
#                                      chunks=(all_wig_data.shape[0], all_wig_data.shape[1]))
#     root.attrs['uni_wig'] = uni_wig_keys
#     root.attrs['all_wig'] = all_wig_keys
#
#     for key in loci_dict:
#         data = np.array(loci_dict[key]).flatten()
#         root.create_dataset('loci_{}'.format(key), data=data, compressor=zarr.Zlib(level=1), chunks=(len(data) or 1))
#     root.attrs['loci'] = loci_chunk_sequence
#
#     if parameters['id'] == 'unknown':
#         parameters['id'] = TEcounts_filepath.split('/')[-1].split('.')[0].split('_')[0]
#     root.attrs['Parameters'] = parameters
#
#     if local:
#         zarr.consolidate_metadata(local_path)
#     else:
#         zarr.consolidate_metadata(store)
#
#     time_end = time.time()
#     print('time cost: ', time_end - time_start, 's')
#
#     return root

def create_RNA_zarr(
        TEcounts_filepath='./data/jiawei_1101_zarr_input/SRR3498333_squire_count/SRR3498333.fastq_TEcounts.txt',
        rmsk_filepath='./data/RNA-seq pipeline/gene_annotation/hg38_rmsk.txt',
        subfamily_size_filepath='./data/RNA-seq pipeline/gene_annotation/hg38TEsubF.size.txt',
        subFcounts_filepath='./data/jiawei_1101_zarr_input/SRR3498333_squire_count/SRR3498333.fastq_subFcounts.txt',
        consensus_json_all_filepath='consensusView_all.json', consensus_json_uni_filepath='consensusView_unique.json',
        parameters={}, local=True, local_path='./data/example.zarr',
        s3_path='s3://repeatbrowsers/![](../../Downloads/3218908-200.png)Zarr_data/example_RNA.zarr'):
    TE_set = read_TEcounts_bysubfamily(TEcounts_filepath)
    rmsk = read_rmsk_bysubfamily(rmsk_filepath)
    consensus_length_set = read_consensus_length(subfamily_size_filepath)
    TE_consensus_set = {}
    not_in_subfam = []
    for TE_name in TE_set.keys():
        if TE_name not in consensus_length_set.keys():
            print("{} Not in!!!".format(TE_name))
            not_in_subfam.append(TE_name)
        else:
            TE_consensus_set[TE_name] = consensus_length_set[TE_name]

    if not_in_subfam:
        [TE_set.pop(item) for item in not_in_subfam]

    all_FPKM_data, uni_FPKM_data = heatmap_data(subFcounts_filepath, TE_consensus_set, TE_set)
    # consensus_set = consensus_map(rmsk, TE_set, TE_consensus_set)

    with open(consensus_json_all_filepath, 'r') as f:
        consensus_json_all_set = json.load(f)
    with open(consensus_json_uni_filepath, 'r') as f:
        consensus_json_uni_set = json.load(f)

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

    time_start = time.time()
    print('start RNA pipeline time counting!')
    if local:
        root = zarr.open(local_path, mode='w')
    else:
        # s3_path = 's3://repeatbrowsers/Zarr_data/example.zarr'
        # Initilize the S3 file system
        s3 = s3fs.S3FileSystem()
        store = s3fs.S3Map(root=s3_path, s3=s3, check=False)
        root = zarr.group(store=store, overwrite=True)

    time_end = time.time()
    print('time cost: ', time_end - time_start, 's')

    subfam_attrs = subfamily[1:]
    root.attrs['subfamily_stat'] = subfam_attrs
    subfam_data = np.array([all_reads_FPKM[1:], uni_reads_FPKM[1:]], dtype='float')
    subfam_stat = root.create_dataset('subfam_stat', data=subfam_data, compressor=zarr.Zlib(level=1),
                                      chunks=(1, stat_subfamily.shape[1]))

    time_end = time.time()
    print('time cost: ', time_end - time_start, 's')

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

    time_end = time.time()
    print('time cost: ', time_end - time_start, 's')

    return root


def create_Experiment_zarr(
        signal_subfam_stat='./ENCFF076NFT_signal_iteres/stat/ENCFF076NFT.sorted.iteres.subfamily.stat',
        control_subfam_stat='./ENCFF610EFV_control_iteres/stat/ENCFF610EFV.sorted.iteres.subfamily.stat',
        signal_all_wig_path='./ENCFF076NFT_signal_iteres/stat/ENCFF076NFT.sorted.iteres.wig',
        signal_uni_wig_path='./ENCFF076NFT_signal_iteres/stat/ENCFF076NFT.sorted.iteres.unique.wig',
        control_all_wig_path='./ENCFF610EFV_control_iteres/stat/ENCFF610EFV.sorted.iteres.wig',
        control_uni_wig_path='./ENCFF610EFV_control_iteres/stat/ENCFF610EFV.sorted.iteres.unique.wig',
        signal_loci_path='./ENCFF076NFT_signal_iteres/filter/ENCFF076NFT.sorted_ALL.iteres.loci',
        control_loci_path='./ENCFF610EFV_control_iteres/filter/ENCFF610EFV.sorted_ALL.iteres.loci',
        range='/Users/jiaweishen/PycharmProjects/Zarr_format/data/RNA-seq pipeline/gene_annotation/hg38TEsubF.size.txt',
        parameters={}, s3_path='s3://repeatbrowsers/Zarr_data/example.zarr',
        local=False, local_path='./data/Experiment.zarr'):
    signal_stat_subfamily = read_statfile_byline(filename=signal_subfam_stat)
    control_stat_subfamily = read_statfile_byline(filename=control_subfam_stat)
    value_stat_subfamily = calculate_ratio(signal_stat_subfamily, control_stat_subfamily)

    signal_all_keys, signal_all_wig = read_wigfile_bysubfamily(signal_all_wig_path)
    signal_uni_keys, signal_uni_wig = read_wigfile_bysubfamily(signal_uni_wig_path)

    control_all_keys, control_all_wig = read_wigfile_bysubfamily(control_all_wig_path)
    control_uni_keys, control_uni_wig = read_wigfile_bysubfamily(control_uni_wig_path)

    signal_loci_dict = read_locifile_bysubfamily(signal_loci_path)
    control_loci_dict = read_locifile_bysubfamily(control_loci_path)

    subfam_range = set()
    with open(range) as file:
        for i, line in enumerate(file):
            if i > 0:
                sf = line.split('\t')[0]
                subfam_range.add(sf)
    signal_loci_dict = {k: v for k, v in signal_loci_dict.items() if k in subfam_range}
    control_loci_dict = {k: v for k, v in control_loci_dict.items() if k in subfam_range}

    loci_chunk_sequence = ['repname', 'chrome', 'start', 'end', 'RPKM']

    time_start = time.time()
    print('start Experiment Zarr Processing.')
    if local:
        root = zarr.open(local_path, mode='w')
    else:
        # s3_path = 's3://repeatbrowsers/Zarr_data/example.zarr'
        # Initilize the S3 file system
        s3 = s3fs.S3FileSystem()
        store = s3fs.S3Map(root=s3_path, s3=s3, check=False)
        root = zarr.group(store=store, overwrite=True)

    time_end = time.time()
    print('Time cost so far: ', time_end - time_start, 's')

    subfam_attrs = list(value_stat_subfamily[0][1:])
    root.attrs['subfamily_stat'] = subfam_attrs
    subfam_data = np.array([value_stat_subfamily[1][1:], value_stat_subfamily[2][1:]], dtype='float')  # 1:ALL and 2:UNI
    subfam_stat = root.create_dataset('subfam_stat', data=subfam_data, compressor=zarr.Zlib(level=1),
                                      chunks=(1, value_stat_subfamily.shape[1]))

    signal_uni_bigwig = root.create_dataset('signal_uni_bigwig', data=signal_uni_wig, compressor=zarr.Zlib(level=1),
                                            chunks=(signal_uni_wig.shape[0], signal_uni_wig.shape[1]))
    signal_all_bigwig = root.create_dataset('signal_all_bigwig', data=signal_all_wig, compressor=zarr.Zlib(level=1),
                                            chunks=(signal_all_wig.shape[0], signal_all_wig.shape[1]))
    root.attrs['signal_uni_wig'] = signal_uni_keys
    root.attrs['signal_all_wig'] = signal_all_keys

    control_uni_bigwig = root.create_dataset('control_uni_bigwig', data=control_uni_wig, compressor=zarr.Zlib(level=1),
                                             chunks=(control_uni_wig.shape[0], control_uni_wig.shape[1]))
    control_all_bigwig = root.create_dataset('control_all_bigwig', data=control_all_wig, compressor=zarr.Zlib(level=1),
                                             chunks=(control_all_wig.shape[0], control_all_wig.shape[1]))
    root.attrs['control_uni_wig'] = control_uni_keys
    root.attrs['control_all_wig'] = control_all_keys

    time_end = time.time()
    print('Time cost so far: ', time_end - time_start, 's')

    for key in signal_loci_dict:
        data = np.array(signal_loci_dict[key]).flatten()
        root.create_dataset('signal_loci_{}'.format(key), data=data, compressor=zarr.Zlib(level=1),
                            chunks=(len(data) or 1))
    for key in control_loci_dict:
        data = np.array(control_loci_dict[key]).flatten()
        root.create_dataset('control_loci_{}'.format(key), data=data, compressor=zarr.Zlib(level=1),
                            chunks=(len(data) or 1))
    root.attrs['loci'] = loci_chunk_sequence

    if parameters['id'] == 'unknown':
        parameters['id'] = signal_subfam_stat.split('/')[-1].split('.')[0]
    root.attrs['Parameters'] = parameters

    if local:
        zarr.consolidate_metadata(local_path)
    else:
        zarr.consolidate_metadata(store)

    time_end = time.time()
    print('Time cost so far: ', time_end - time_start, 's')
    return root


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process the input files to the .Zarr format.')

    parser.add_argument('--s3_path', '-s', default='s3://repeatbrowsers/Zarr_data/example.zarr', type=str,
                        help="s3 path for AWS S3 bucket")
    parser.add_argument('--Biosample', '-b', default='unknown', type=str, help="The biosample type of your input file.")
    parser.add_argument('--Tissue', '-t', default='unknown', type=str, help="The tissue type of your input file.")
    parser.add_argument('--Assay', '-a', default='unknown', type=str, help="The assay type of your input file.")
    parser.add_argument('--Control', '-c', default='unknown', type=str, help="The control id of your input file.")
    parser.add_argument('--Organism', '-o', default='unknown', type=str, help="The organism of your input file.")
    parser.add_argument('--id', default='unknown', type=str, help="The file id (accession) of your input file.")

    parser.add_argument('--loci_path', '-lp', default='testAligned.sortedByCoord.out_ALL.iteres.loci', type=str,
                        help="The path of loci file output by iteres.")
    parser.add_argument('--subfamily_stat_path', '-ssp',
                        default='default_with_unmapped_bam_and_readsAligned.sortedByCoord.out.iteres.subfamily.stat',
                        type=str, help="The path of subfamily stat file output by iteres.")
    parser.add_argument('--all_wig_path', '-awp', default='out_all.wig', type=str,
                        help="The path of .wig file (All) output by iteres.")
    parser.add_argument('--uni_wig_path', '-uwp', default='out_unique.wig', type=str,
                        help="The path of .wig file (uni) output by iteres.")
    parser.add_argument("--local", action="store_true", help="Create .Zarr file locally.")
    parser.add_argument('--output_path', '-op', default='./example.zarr', type=str,
                        help="The local path of output files, which must end with *.zarr")

    parser.add_argument("--RNA", action="store_true", help="The input files come from RNA pipeline.")
    parser.add_argument("--CAGE", action="store_true", help="The input files come from CAGE pipeline.")
    parser.add_argument('--rmsk_filepath', '-rmsk', default='./data/RNA-seq pipeline/gene_annotation/hg38_rmsk.txt',
                        type=str, help="The path of rmsk file.")
    parser.add_argument('--TEcounts_filepath',
                        default='./data/jiawei_1101_zarr_input/SRR3498336_squire_count/SRR3498336.fastq_TEcounts.txt',
                        type=str,
                        help="The path of TEcounts file (All) output by SQUIRE.")
    parser.add_argument('--subFcounts_filepath',
                        default='./data/jiawei_1101_zarr_input/SRR3498336_squire_count/SRR3498336.fastq_subFcounts.txt',
                        type=str,
                        help="The path of subFcounts file (All) output by SQUIRE.")
    parser.add_argument('--subfamily_size_filepath',
                        default='./data/RNA-seq pipeline/gene_annotation/hg38TEsubF.size.txt', type=str,
                        help="The path of subfamily size file (All) output by SQUIRE.")

    parser.add_argument('--signal_subfam_stat',
                        default='./ENCFF076NFT_signal_iteres/stat/ENCFF076NFT.sorted.iteres.subfamily.stat', type=str,
                        help="The path of signal's subfamily stat file output by iteres.")
    parser.add_argument('--control_subfam_stat',
                        default='./ENCFF610EFV_control_iteres/stat/ENCFF610EFV.sorted.iteres.subfamily.stat', type=str,
                        help="The path of control's subfamily stat file output by iteres.")
    parser.add_argument('--signal_all_wig_path',
                        default='./ENCFF076NFT_signal_iteres/stat/ENCFF076NFT.sorted.iteres.wig', type=str,
                        help="The path of signal's .wig file (All) output by iteres.")
    parser.add_argument('--signal_uni_wig_path',
                        default='./ENCFF076NFT_signal_iteres/stat/ENCFF076NFT.sorted.iteres.unique.wig', type=str,
                        help="The path of signal's .wig file (uni) output by iteres.")
    parser.add_argument('--control_all_wig_path',
                        default='./ENCFF610EFV_control_iteres/stat/ENCFF610EFV.sorted.iteres.wig', type=str,
                        help="The path of control's .wig file (All) output by iteres.")
    parser.add_argument('--control_uni_wig_path',
                        default='./ENCFF610EFV_control_iteres/stat/ENCFF610EFV.sorted.iteres.unique.wig', type=str,
                        help="The path of control's .wig file (uni) output by iteres.")
    parser.add_argument('--signal_loci_path',
                        default='./ENCFF076NFT_signal_iteres/filter/ENCFF076NFT.sorted_ALL.iteres.loci', type=str,
                        help="The signal's path of loci file output by iteres.")
    parser.add_argument('--control_loci_path',
                        default='./ENCFF610EFV_control_iteres/filter/ENCFF610EFV.sorted_ALL.iteres.loci', type=str,
                        help="The control's path of loci file output by iteres.")
    parser.add_argument("--EXPERIMENT", action="store_true", help="The input should contain both signal and control.")

    parser.add_argument('--consensus_json_all_filepath', default='consensusView_all.json', type=str,
                        help="The consensus json file's path of the output by the script (ALL).")
    parser.add_argument('--consensus_json_uni_filepath', default='consensusView_unique.json', type=str,
                        help="The consensus json file's path of the output by the script (UNIQUE).")

    parser.add_argument('--cageWigPlus', default='cageWigPlus.json', type=str,
                        help="The wig file output by the iteres-cage, the + strand.")
    parser.add_argument('--cageWigMinus', default='cageWigMinus.json', type=str,
                        help="The wig file output by the iteres-cage, the - strand.")
    parser.add_argument('--cageWigUniPlus', default='cageWigUniPlus.json', type=str,
                        help="The wig file output by the iteres-cage, the + strand for the unique reads.")
    parser.add_argument('--cageWigUniMinus', default='cageWigUniMinus.json', type=str,
                        help="The wig file output by the iteres-cage, the - strand for the unique reads.")

    # parser.add_argument("--file_input_RNA",
    #                     default='./FILE_RNA.txt', type=str,
    #                     help="The file_input_RNA is a .txt or .json file contains the multiple inputs' info.")

    args = parser.parse_args()
    parameters = {"Biosample": args.Biosample, "Assay": args.Assay, "Control": args.Control, "Organism": args.Organism,
                  "id": args.id, "Tissue": args.Tissue, "Mode": "Experiment" if args.EXPERIMENT else "File"}

    if args.RNA:
        root = create_RNA_zarr(TEcounts_filepath=args.TEcounts_filepath,
                               rmsk_filepath=args.rmsk_filepath,
                               subfamily_size_filepath=args.subfamily_size_filepath,
                               subFcounts_filepath=args.subFcounts_filepath,
                               consensus_json_all_filepath=args.consensus_json_all_filepath,
                               consensus_json_uni_filepath=args.consensus_json_uni_filepath,
                               parameters=parameters, local=args.local, local_path=args.output_path,
                               s3_path=args.s3_path)
    elif args.EXPERIMENT:
        root = create_Experiment_zarr(signal_subfam_stat=args.signal_subfam_stat,
                                      control_subfam_stat=args.control_subfam_stat,
                                      signal_all_wig_path=args.signal_all_wig_path,
                                      signal_uni_wig_path=args.signal_uni_wig_path,
                                      control_all_wig_path=args.control_all_wig_path,
                                      control_uni_wig_path=args.control_uni_wig_path,
                                      signal_loci_path=args.signal_loci_path,
                                      control_loci_path=args.control_loci_path,
                                      range=args.subfamily_size_filepath,
                                      parameters=parameters, s3_path=args.s3_path,
                                      local=args.local, local_path=args.output_path)
    elif args.CAGE:
        root = create_iteres_zarr_cage(
            subfam_stat_path=args.subfamily_stat_path,
            all_wig_path_plus=args.cageWigPlus, all_wig_path_minus=args.cageWigMinus,
            uni_wig_path_plus=args.cageWigUniPlus, uni_wig_path_minus=args.cageWigUniMinus,
            loci_path=args.loci_path,
            parameters=parameters, s3_path=args.s3_path,
            range=args.subfamily_size_filepath,
            local=args.local, local_path=args.output_path)
    else:
        root = create_iteres_zarr(subfam_stat_path=args.subfamily_stat_path, loci_path=args.loci_path,
                                  range=args.subfamily_size_filepath,
                                  all_wig_path=args.all_wig_path, uni_wig_path=args.uni_wig_path, parameters=parameters,
                                  s3_path=args.s3_path, local=args.local, local_path=args.output_path)
        # zarr.save_group('data/example.zarr', root)
    print(root.info)

    # # s3_path = 's3://repeatbrowser/stat/example.zarr'
    # # s3_path = 's3://repeatbrowsers/Zarr_data/example.zarr'
    # s3_path = args.s3_path
    # # Initilize the S3 file system
    # s3 = s3fs.S3FileSystem()
    # store = s3fs.S3Map(root=s3_path, s3=s3, check=False)
    # # Read Zarr file
    # z2 = zarr.group(store=store)
    # # z2 = zarr.open_group('data/example.zarr', 'r')
    # print(z2.info)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
