# importing libraries
import pandas as pd
import os
import time
import threading
import requests
import json
import argparse

def parse_args():
    p = argparse.ArgumentParser(
        description="Take ENCODE tsv and download BED files with highest FRiP scores"
    )

    p.add_argument(
        "-i", "--input",
        required=True,
        help="Path to input ENCODE TSV file"
    )

    p.add_argument(
        "-m", "--min_peaks",
        required=False,
        default=0,
        help="Minimum number of peaks to be considered a valid, high quality file"
    )

    p.add_argument(
        "-t", "--type",
        required=True,
        choices=["histone", "TF"],
        help="Specifying histone or TF TSV file"
    )

    p.add_argument(
        "-p", "--parallel",
        required=False,
        default=8,
        type=int,
        help="Number of parallel threads to use for pipeline"
    )

    return p.parse_args()

def split_into_chunks(lst, n):
    avg_len = len(lst)//n
    remainder = len(lst)%n
    chunks = []
    start = 0
    
    for i in range(n):
        # Each chunk gets an additional element if there's a remainder
        end = start + avg_len + (1 if i < remainder else 0)
        chunks.append(lst[start:end])
        start = end

    return chunks

def download_metadata(accession_list, metadata_dict):
    headers = {'accept': 'application/json'}
    counter = 0
    for accession in accession_list:
        if counter % 10 == 0:
            time.sleep(1)
        url = 'https://www.encodeproject.org/files/' + accession
        response = requests.get(url, headers=headers).json()

        #save as json file
        with open( '../data/metadata/' + accession + '.json', 'w') as f:
            json.dump(response, f)

        metadata_dict[accession] = response
        counter += 1

def metadata_run(ENCODE_tsv, n):
    """
    Downloads metadata JSONS from samples specified in the input TSV file from the ENCODE api.
    """
    # Create the output directory "../data/metadata" . if already exists, delete all files in it
    os.makedirs("../data/metadata", exist_ok=True)
    # Delete all files in the directory
    for filename in os.listdir("../data/metadata"):
        file_path = os.path.join("../data/metadata", filename)
        try:
            if os.path.isfile(file_path):
                os.remove(file_path)
        except Exception as e:
            print(f"Error deleting file {file_path}: {e}")

    # Read the TSV file
    ENCODE_df = pd.read_csv(ENCODE_tsv, sep='\t', skiprows=1)
    ENCODE_df = ENCODE_df[['Dataset', 'Accession', 'Target label']]

    #split list into n groups for parallel processing
    ENCODE_df = ENCODE_df['Accession'].tolist()
    bed_chunks = split_into_chunks(ENCODE_df, n)
    metadata_dict = {}
    threads = []

    for chunk in bed_chunks:
        thread = threading.Thread(target=download_metadata, args=(chunk, metadata_dict))
        threads.append(thread)
        thread.start()
    for thread in threads:
        thread.join()

def read_metadata_TF(min_peaks):
    """
    Downloads BED file with highest frip for each target. Excludes files with less than threshold # of peaks.
    """

    # Create the output directory if it doesn't exist
    os.makedirs("../data/bed", exist_ok=True)
    bed_json_path = '../data/metadata/'

    # Get data from JSON files
    BED_frip_dict = {}
    BED_dataset_dict = {}
    BED_target_dict = {}
    BED_reproducible_peaks = {}

    for json_file in os.listdir(bed_json_path):
        with open(bed_json_path + json_file) as f:
            data = json.load(f)
            if len(data['quality_metrics']) != 0:
                if 'frip' in data['quality_metrics'][0]:
                    frip = data['quality_metrics'][0]['frip']
                elif len(data['quality_metrics']) == 2 and 'frip' in data['quality_metrics'][1]:
                    frip = data['quality_metrics'][1]['frip']
                else:
                    frip = None
            else:
                frip = None

            if len(data['quality_metrics']) != 0:
                if 'reproducible_peaks' in data['quality_metrics'][0]:
                    reproducible_peaks = data['quality_metrics'][0]['reproducible_peaks']
                elif len(data['quality_metrics']) == 2 and 'reproducible_peaks' in data['quality_metrics'][1]:
                    reproducible_peaks = data['quality_metrics'][1]['reproducible_peaks']
                else:
                    reproducible_peaks = None
        
                BED_frip_dict[json_file] = frip

            BED_reproducible_peaks[json_file] = reproducible_peaks
            dataset = data['dataset']
            BED_dataset_dict[json_file] = dataset
            target = data['target']['label']
            BED_target_dict[json_file] = target
    
    # make into dataframe
    #save to df
    df_final = pd.DataFrame.from_dict(BED_frip_dict, orient='index', columns=['frip'])

    #add reproducible_peaks column
    df_final['reproducible_peaks'] = df_final.index
    df_final['reproducible_peaks'] = df_final['reproducible_peaks'].map(BED_reproducible_peaks)

    #add dataset column
    df_final['dataset'] = df_final.index
    df_final['dataset'] = df_final['dataset'].map(BED_dataset_dict)
    #keep element 2 of dataset values using / as delimiter
    df_final['dataset'] = df_final['dataset'].str.split('/').str[2]

    #add target column
    df_final['target'] = df_final.index
    df_final['target'] = df_final['target'].map(BED_target_dict)

    #add accession column which is index
    df_final['accession'] = df_final.index
    df_final['accession'] = df_final['accession'].str[:-5]

    #reset index
    df_final.reset_index(drop=True, inplace=True)

    #order columns as accession, target, dataset, frip
    df_final = df_final[['accession', 'target', 'dataset', 'frip', 'reproducible_peaks']]

    #order by frip
    df_final = df_final.sort_values(by='frip', ascending=False)

    #keep first instance of target
    df_final = df_final.drop_duplicates(subset='target', keep='first')

    #keep only rows with reproducible_peaks > min_peaks
    df_final = df_final[df_final['reproducible_peaks'] > min_peaks]

    return df_final

def read_metadata_histone(min_peaks):
    # Create the output directory if it doesn't exist
    os.makedirs("../data/bed", exist_ok=True)
    bed_json_path = '../data/metadata/'

    BED_frip_dict = {}
    BED_dataset_dict = {}
    BED_target_dict = {}

    for json_file in os.listdir(bed_json_path):
        with open(bed_json_path + json_file) as f:
            data = json.load(f)
            if len(data['quality_metrics']) != 0:
                if 'frip' in data['quality_metrics'][0]:
                    frip = data['quality_metrics'][0]['frip']
                elif len(data['quality_metrics']) == 2 and 'frip' in data['quality_metrics'][1]:
                    frip = data['quality_metrics'][1]['frip']
                else:
                    frip = None
            else:
                frip = None

            dataset = data['dataset']
            target = data['target']['label']

            BED_frip_dict[json_file] = frip
            BED_dataset_dict[json_file] = dataset
            BED_target_dict[json_file] = target

    #save frip to df
    df_final = pd.DataFrame.from_dict(BED_frip_dict, orient='index', columns=['frip'])

    #add target column
    df_final['target'] = df_final.index
    df_final['target'] = df_final['target'].map(BED_target_dict)

    #add dataset column
    df_final['dataset'] = df_final.index
    df_final['dataset'] = df_final['dataset'].map(BED_dataset_dict)
    df_final['dataset'] = df_final['dataset'].str.split('/').str[2]

    #order by frip
    df_final = df_final.sort_values(by='frip', ascending=False)

    #keep first instance of target
    df_final = df_final.drop_duplicates(subset='target', keep='first')

    #parse index to get accession and store in column
    df_final['accession'] = df_final.index
    df_final['accession'] = df_final['accession'].str[:-5]

    #reset index
    df_final.reset_index(drop=True, inplace=True)

    #order columns as accession,target,dataset,frip
    df_final = df_final[['accession', 'target', 'dataset', 'frip']]

    #TODO: add peaks columns and check if > min_peaks
    
    return df_final

#download BED files and name them according to target format: "https://www.encodeproject.org/files/ENCFF636FWF/@@download/ENCFF636FWF.bed.gz"
def download_BED_file_helper(target, accession):
    os.system(f"wget -q -O ../data/bed/{target}.bed.gz 'https://www.encodeproject.org/files/{accession}/@@download/{accession}.bed.gz'")
    print(f"Downloaded {target}.bed")

def download_bed_file(accession_df, threads):
    """
    Make accession dictionary with target as key and accession as value. Download bed files from ENCODE.
    """
    #keep first 3 columns
    accession_df = accession_df[['accession', 'target', 'dataset']]

    #make dict (key = target, value = accession) from pandas dataframe
    accession_dict = accession_df.set_index('target')['accession'].to_dict()

    active = []  # list of currently running Thread objects

    for target, accession in accession_dict.items():
        # start a new thread for this download
        t = threading.Thread(target=download_BED_file_helper, args=(target, accession), daemon=True)
        t.start()
        active.append(t)

        # if we've hit our limit, wait for the oldest to finish
        if len(active) >= threads:
            active[0].join()
            active.pop(0)

    # wait for any remaining threads
    for t in active:
        t.join()

    print("All BED files downloaded")


def main():
    args = parse_args()
    input_file = args.input
    min_peaks = int(args.min_peaks)
    parallel = args.parallel
    type = args.type

    # Download metadata
    print("Downloading metadata...")
    metadata_run(input_file, parallel)
    print("Metadata downloaded.")

    # Read metadata and download BED files
    if type == "TF":
        accession_df = read_metadata_TF(min_peaks)
        print("TF metadata read.")
    elif type == "histone":
        accession_df = read_metadata_histone(min_peaks)
        print("Histone metadata read.")
    else:
        raise ValueError("Invalid type specified. Use 'histone' or 'TF'.")
    
    # Download BED files
    download_bed_file(accession_df, parallel)


if __name__ == "__main__":
    main()