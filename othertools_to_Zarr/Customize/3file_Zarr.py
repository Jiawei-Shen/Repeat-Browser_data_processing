import pandas as pd
import zarr
import json
import numpy as np
import argparse


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


def read_json(file_path):
    # Open the JSON file
    with open(file_path, 'r') as json_file:
        # Load the JSON data into a dictionary
        data_dict = json.load(json_file)

        # Print or process the loaded dictionary
        print(data_dict)
    return data_dict


def generate_Zarr(parameters, heatmap="./sample_files/heatmap.csv",
                  uni_wig="./sample_files/ENCFF076NFT.sorted.iteres.unique.wig",
                  all_wig="./sample_files/ENCFF076NFT.sorted.iteres.wig",
                  loci_file="./sample_files/data_list.json",
                  local_path='./sample.zarr'):
    heatmap = pd.read_csv(heatmap)
    all_keys, all_wig = read_wigfile_bysubfamily(all_wig)
    uni_keys, uni_wig = read_wigfile_bysubfamily(uni_wig)
    loci_dict = read_json(loci_file)

    root = zarr.open(local_path, mode='w')

    heatmap_names = np.concatenate((['#subfamily'], heatmap["TE"].values))
    heatmap_all_value = np.concatenate((['all_value'], heatmap['all_value'].values))
    heatmap_uni_value = np.concatenate((['uni_value'], heatmap['uni_value'].values))
    subfam_data = np.transpose(heatmap[['all_value', 'uni_value']].values)
    stat_subfamily = np.stack((heatmap_names, heatmap_all_value, heatmap_uni_value), axis=0)

    subfam_stat = root.create_dataset('subfam_stat', data=subfam_data, compressor=zarr.Zlib(level=1),
                                      chunks=(1, stat_subfamily.shape[1]))
    uni_bigwig = root.create_dataset('uni_bigwig', data=uni_wig, compressor=zarr.Zlib(level=1),
                                     chunks=(uni_wig.shape[0], uni_wig.shape[1]))
    all_bigwig = root.create_dataset('all_bigwig', data=all_wig, compressor=zarr.Zlib(level=1),
                                     chunks=(all_wig.shape[0], all_wig.shape[1]))
    root.attrs['uni_wig'] = uni_keys
    root.attrs['all_wig'] = all_keys

    for key in loci_dict:
        data = np.array(loci_dict[key]).flatten()
        root.create_dataset('loci_{}'.format(key), data=data, compressor=zarr.Zlib(level=1), chunks=(len(data) or 1))
    root.attrs['loci'] = ['repname', 'chrome', 'start', 'end', 'FPKM']

    if parameters['id'] == 'unknown':
        parameters['id'] = heatmap.split('/')[-1].split('.')[0].split('_')[0]
    root.attrs['Parameters'] = parameters

    zarr.consolidate_metadata(local_path)

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process the SQuIRE output files to the .Zarr format.')

    parser.add_argument('--Biosample', '-b', default='unknown', type=str, help="The biosample type of your input file.")
    parser.add_argument('--Tissue', '-t', default='unknown', type=str, help="The tissue type of your input file.")
    parser.add_argument('--Assay', '-a', default='unknown', type=str, help="The assay type of your input file.")
    parser.add_argument('--Control', '-c', default='unknown', type=str, help="The control id of your input file.")
    parser.add_argument('--Organism', default='unknown', type=str, help="The organism of your input file.")
    parser.add_argument('--id', default='unknown', type=str, help="The file id (accession) of your input file.")

    parser.add_argument('--heatmap',
                        default="./sample_files/heatmap.csv",
                        type=str, help="The file path of heatmap data, it should be a formated csv file.")
    parser.add_argument('--uniWig',
                        default='./sample_files/ENCFF076NFT.sorted.iteres.unique.wig',
                        type=str, help="The path of wig files (uni) of the consensus info.")
    parser.add_argument('--allWig',
                        default='./sample_files/ENCFF076NFT.sorted.iteres.wig', type=str,
                        help="The path of subfamily size file (All) output by SQUIRE.")
    parser.add_argument('--lociFile', default='./sample_files/data_list.json',
                        type=str, help="The file path of rmsk file.")
    parser.add_argument('--output', '-o', default='./example.zarr', type=str,
                        help="The local path of output files, which must end with *.zarr")


    args = parser.parse_args()
    parameters = {"Biosample": args.Biosample, "Assay": args.Assay, "Control": args.Control, "Organism": args.Organism,
                  "id": args.id, "Tissue": args.Tissue}

    generate_Zarr(parameters, heatmap=args.heatmap,
                  uni_wig=args.uniWig,
                  all_wig=args.allWig,
                  loci_file=args.lociFile,
                  local_path=args.output)
