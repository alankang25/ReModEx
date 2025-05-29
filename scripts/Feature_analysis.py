#TODO: Clean up interim files and add comments

#importing libraries
import pandas as pd
import numpy as np
import os
import argparse
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
from sklearn.linear_model import Ridge

#--------------block for parsing arguments-------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Take ENCODE tsv and download BED files with highest FRiP scores"
    )

    p.add_argument(
        "-c", "--control_ATAC",
        type=str,
        required=True,
        help="Path to the directory containing control ATAC-seq bed and bigwig file (e.g., ./control_peaks/)"
    )

    p.add_argument(
        "-t", "--treatment_ATAC",
        type=str,
        required=True,
        help="Path to the directory containing treatment ATAC-seq bed and bigwig file (e.g., ./treatment_peaks/)"
    )
    args = p.parse_args()
    if not os.path.exists(args.control_ATAC):
        raise FileNotFoundError(f"Control ATAC-seq directory {args.control_ATAC} does not exist.")
    if not os.path.exists(args.treatment_ATAC):
        raise FileNotFoundError(f"Treatment ATAC-seq directory {args.treatment_ATAC} does not exist.")
    return args

#--------------block for making feature matrix-------------------
def make_feature_matrix(bed_list, feature_matrix):
    """
    Overlaps the ATAC-seq peaks with the features in the bed_list and returns a binary feature matrix.
    Parameters:
    bed_list: list of bed files containing features
    feature_matrix: pandas DataFrame with features as columns and ATAC-seq peaks as rows with target as log2 fold change
    """
    #calculate overlap between ATAC-seq peaks and features using bedtools overlap
    #if ./data/overlap_beds exists, delete it and create a new one
    if os.path.exists("../data/output/overlap_beds"):
        os.system("rm -r ../data/output/overlap_beds")
    os.makedirs("../data/output/overlap_beds", exist_ok=True)

    for bed_file in bed_list:
        #get the name of the bed file without the path
        bed_name = os.path.basename(bed_file).split('.')[0]
        #use bedtools intersect to get the overlap between ATAC-seq peaks and features
        os.system(f"bedtools intersect -a ../data/output/final_merged_peaks.bed -b {bed_file} -wa -wb > ../data/output/overlap_beds/{bed_name}_overlap.bed")

    #read the overlap files and create a binary feature matrix
    feature_matrix = pd.DataFrame(columns=['peak_name'] + [os.path.basename(f).split('.')[0] for f in bed_list])
    feature_matrix['peak_name'] = pd.read_csv("../data/output/final_merged_peaks.bed", sep="\t", header=None)[3]
    feature_matrix.set_index('peak_name', inplace=True)
    for bed_file in bed_list:
        #get the name of the bed file without the path
        bed_name = os.path.basename(bed_file).split('.')[0]
        #read the overlap file
        overlap_file = f"../data/output/overlap_beds/{bed_name}_overlap.bed"
        if os.path.exists(overlap_file):
            overlap_df = pd.read_csv(overlap_file, sep="\t", header=None)
            #create a binary column for the feature
            feature_matrix[bed_name] = 0
            #set the value to 1 for the peaks that overlap with the feature
            feature_matrix.loc[overlap_df[3], bed_name] = 1

    #fill NaN values with 0
    feature_matrix.fillna(0, inplace=True)
    #save the feature matrix to a csv file
    feature_matrix.to_csv("../data/output/feature_matrix.csv")

    return feature_matrix

def modify_args(control_ATAC, treatment_ATAC):
    """
    Get bed and bigwig files from control and treatment ATAC-seq directories and return them.
    """

    #check if paths end with "/"
    if not control_ATAC.endswith('/'):
        control_ATAC += '/'
    if not treatment_ATAC.endswith('/'):
        treatment_ATAC += '/'

    print(f"Control ATAC-seq directory: {control_ATAC}")
    print(f"Treatment ATAC-seq directory: {treatment_ATAC}")

    #get bed and bigwig files from control ATAC-seq directory
    control_bed = [f for f in os.listdir(control_ATAC) if f.endswith('.bed') or f.endswith('.bed.gz')]
    if not control_bed:
        raise FileNotFoundError(f"No bed files found in control ATAC-seq directory {control_ATAC}.")
    control_bed = os.path.join(control_ATAC, control_bed[0])  # take the first bed file found
    control_bigwig = [f for f in os.listdir(control_ATAC) if f.endswith('.bw') or f.endswith('.bigwig')]
    if not control_bigwig:
        raise FileNotFoundError(f"No bigwig files found in control ATAC-seq directory {control_ATAC}.")
    control_bigwig = os.path.join(control_ATAC, control_bigwig[0])  # take the first bigwig file found
    print(f"Control ATAC-seq bed file: {control_bed}")
    print(f"Control ATAC-seq bigwig file: {control_bigwig}")

    #get bed and bigwig files from treatment ATAC-seq directory
    treatment_bed = [f for f in os.listdir(treatment_ATAC) if f.endswith('.bed') or f.endswith('.bed.gz')]
    if not treatment_bed:
        raise FileNotFoundError(f"No bed files found in treatment ATAC-seq directory {treatment_ATAC}.")
    treatment_bed = os.path.join(treatment_ATAC, treatment_bed[0])  # take the first bed file found
    treatment_bigwig = [f for f in os.listdir(treatment_ATAC) if f.endswith('.bw') or f.endswith('.bigwig')]
    if not treatment_bigwig:
        raise FileNotFoundError(f"No bigwig files found in treatment ATAC-seq directory {treatment_ATAC}.")
    treatment_bigwig = os.path.join(treatment_ATAC, treatment_bigwig[0])  # take the first bigwig file found
    print(f"Treatment ATAC-seq bed file: {treatment_bed}")
    print(f"Treatment ATAC-seq bigwig file: {treatment_bigwig}")

    #return bed and bigwig files
    return control_bed, treatment_bed, control_bigwig, treatment_bigwig


def make_peak_set(control_bed, treatment_bed):
    """
    Make a set of ATAC-seq peaks from control and treatment ATAC-seq bed files.
    Parameters:
    control_bed: path to the control ATAC-seq bed file
    treatment_bed: path to the treatment ATAC-seq bed file
    Returns: merged set of ATAC-seq peaks labeled with peak number
    """

    #bedtools merge two bed files
    os.system("cat " + control_bed + " " + treatment_bed + " > ../data/output/cat_peaks.bed")
    os.system("sort -k1,1 -k2,2n ../data/output/cat_peaks.bed > ../data/output/cat_peaks_sorted.bed")
    os.system("bedtools merge -i ../data/output/cat_peaks_sorted.bed > ../data/output/merged_peaks.bed")

    #read merged peaks and label them with peak number
    merged_peaks = pd.read_csv("../data/output/merged_peaks.bed", sep="\t", header=None)
    merged_peaks.columns = ['chrom', 'start', 'end']
    #peak name is Peak_{index}
    merged_peaks['peak_name'] = ['Peak_' + str(i) for i in range(len(merged_peaks))]
    #save merged peaks to bed file
    merged_peaks.to_csv("../data/output/final_merged_peaks.bed", sep="\t", header=False, index=False)

def calculate_log2_fold_change(control_bigwig, treatment_bigwig, merged_peaks="../data/output/final_merged_peaks.bed"):
    #TODO: retire this method and use differential peak calling
    #use bigwigaverage over bed to calculate log2 fold change for ATAC-seq peaks
    os.system(f"/mnt/disk/share/software/ucsc_utilities/bigWigAverageOverBed {control_bigwig} {merged_peaks} ../data/output/control_ATAC_peaks.bed")
    os.system(f"/mnt/disk/share/software/ucsc_utilities/bigWigAverageOverBed {treatment_bigwig} {merged_peaks} ../data/output/treatment_ATAC_peaks.bed")

    #read control and treatment ATAC-seq peaks
    control_peaks = pd.read_csv("../data/output/control_ATAC_peaks.bed", sep="\t", header=None)
    treatment_peaks = pd.read_csv("../data/output/treatment_ATAC_peaks.bed", sep="\t", header=None)

    #make dictionary with column 0 as key and column 5 as value
    control_dict = dict(zip(control_peaks[0], control_peaks[5]))
    treatment_dict = dict(zip(treatment_peaks[0], treatment_peaks[5]))

    #make dataframe with length of merged peaks
    merged_peaks_df = pd.read_csv(merged_peaks, sep="\t", header=None)
    merged_peaks_df.columns = ['chrom', 'start', 'end', 'peak_name']
    merged_peaks_df['control'] = merged_peaks_df['peak_name'].map(control_dict)
    merged_peaks_df['treatment'] = merged_peaks_df['peak_name'].map(treatment_dict)

    #calculate log2 fold change
    merged_peaks_df['target'] = np.log2(merged_peaks_df['treatment'] / merged_peaks_df['control'])

    return merged_peaks_df[['peak_name', 'target']]

def random_forest(feature_matrix, target_column='target'):

    return

def main():
    # Define the list of bed files containing features
    bed_path = '../data/bed/'
    bed_list = [os.path.join(bed_path, f) for f in os.listdir(bed_path) if f.endswith('.bed') or f.endswith('.bed.gz')]

    # Parse command line arguments
    args = parse_args()
    control_ATAC = args.control_ATAC
    treatment_ATAC = args.treatment_ATAC

    #save bed and bigwig file paths as variables
    control_bed, treatment_bed, control_bigwig, treatment_bigwig = modify_args(control_ATAC, treatment_ATAC)

    # Make peak set
    make_peak_set(control_bed, treatment_bed)

    # calculate log2 fold change for ATAC-seq peaks
    feature_df = calculate_log2_fold_change(control_bigwig, treatment_bigwig)

    # Make feature matrix
    feature_matrix = make_feature_matrix(bed_list, feature_df)

    #random forest classifier


    return 


if __name__ == "__main__":
    main()