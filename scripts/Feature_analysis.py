#TODO: Clean up interim files and add comments
#TODO: Add feature where you can exclude certain features (i.e. remove rows with certain features as 1)
#TODO: Add comments that specify the equations used for feature importance calculations in sklearn

#importing libraries
import pandas as pd
import numpy as np
import os
import argparse

#----------------block for parsing arguments---------------------
def parse_args():
    p = argparse.ArgumentParser(
        description="Take ENCODE tsv and download BED files with highest FRiP scores"
    )
    
    p.add_argument(
        "-x", "--exclude_features",
        type=str,
        nargs='*',
        default=[],
        help="List of features to exclude from the feature matrix (e.g., ['feature1', 'feature2'])."
        " Features should be specified as the name of the bed file without the .bed extension."
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
#--------------------------end of block--------------------------
#--------------block for making feature matrix-------------------
def make_feature_matrix(bed_list, feature_matrix_input):
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
    feature_matrix['target'] = feature_matrix_input.set_index('peak_name')['target']
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
#--------------------------end of block--------------------------

#-------------------block for machine learning-------------------
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc, precision_recall_curve
from sklearn.linear_model import Ridge
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler

def random_forest(feature_matrix, target_column='target'):
    """
    Train a Random Forest classifier on the feature matrix.
    Parameters:
    feature_matrix: pandas DataFrame with features as columns and ATAC-seq peaks as rows with target as log2 fold change
    target_column: column name for the target variable (default is 'target')
    """
    #remove nan values from the feature matrix
    feature_matrix = feature_matrix.dropna()

    # 1 if log2 fold change is less than -1, 0 if other
    feature_matrix[target_column] = (feature_matrix[target_column] < -1).astype(int)

    #balance the dataset by undersampling
    class_counts = feature_matrix[target_column].value_counts()
    min_class_count = class_counts.min()
    balanced_data = pd.concat([
        feature_matrix[feature_matrix[target_column] == cls].sample(min_class_count, random_state=3)
        for cls in class_counts.index
    ])
    # Shuffle the balanced data
    balanced_data = balanced_data.sample(frac=1, random_state=3).reset_index(drop=True)
    # Create the feature matrix
    feature_matrix = balanced_data

    # Split the data into features and target
    X = feature_matrix.drop(columns=[target_column])
    y = feature_matrix[target_column]

    # Split the data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=3)

    # Initialize the Random Forest classifier
    rf_classifier = RandomForestClassifier(n_estimators=100, random_state=3)
    # Fit the classifier to the training data
    rf_classifier.fit(X_train, y_train)

    # Make predictions on the test set
    y_pred = rf_classifier.predict(X_test)
    # Calculate accuracy
    accuracy = accuracy_score(y_test, y_pred)
    print(f"Random Forest Classifier Accuracy: {accuracy:.2f}")
    # Print classification report
    print("Classification Report:")
    print(classification_report(y_test, y_pred))
    # Print confusion matrix
    print("Confusion Matrix:")
    print(confusion_matrix(y_test, y_pred))

    #---------------block for plotting performance metrics-------------------
    # Plot ROC and Precision-Recall curves on the same figuer side by side
    fpr, tpr, _ = roc_curve(y_test, rf_classifier.predict_proba(X_test)[:, 1])
    roc_auc = auc(fpr, tpr)
    precision, recall, _ = precision_recall_curve(y_test, rf_classifier.predict_proba(X_test)[:, 1])
    pr_auc = auc(recall, precision)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    ax1.plot(fpr, tpr, color='orange', label=f'ROC curve (area = {roc_auc:.2f})')
    ax1.plot([0, 1], [0, 1], color='grey', linestyle='--', alpha=0.5)
    ax1.set_xlabel('False Positive Rate')
    ax1.set_ylabel('True Positive Rate')
    ax1.set_title('Receiver Operating Characteristic (ROC) Curve')
    ax1.legend(loc='lower right')
    ax2.plot(recall, precision, color='green', label=f'Precision-Recall curve (area = {pr_auc:.2f})')
    ax2.plot([0, 1], [0.5, 0.5], color='grey', linestyle='--', alpha=0.5)
    ax2.set_xlabel('Recall')
    ax2.set_ylabel('Precision')
    ax2.set_title('Precision-Recall Curve')
    ax2.legend(loc='lower left')
    plt.tight_layout()
    plt.savefig("../data/output/roc_pr_curves.pdf", format='pdf')
    #--------------------------end of block--------------------------------
    
    return rf_classifier

def random_forest_feature_importance(rf_classifier):
    """
    Calculate and print feature importance from the Random Forest classifier.
    Parameters:
    rf_classifier: trained Random Forest classifier
    """
    # Get feature importances
    importances = rf_classifier.feature_importances_
    # Create a DataFrame for feature importances
    feature_importance_df = pd.DataFrame({
        'feature': rf_classifier.feature_names_in_,
        'importance': importances
    })
    # Sort by importance
    feature_importance_df = feature_importance_df.sort_values(by='importance', ascending=False)
    
    print("Feature Importances:")
    print(feature_importance_df)

    #---------------block for plotting MDI feature importances-------------------
    plt.figure(figsize=(6, 8))
    #plot top 30 features
    feature_importance_plot = feature_importance_df.head(30)
    sns.barplot(x='importance', y='feature', data=feature_importance_plot, color=sns.xkcd_rgb['dark peach'])
    plt.title('MDI Feature Importances from RF Classifier')
    plt.xlabel('Mean Decrease in Impurity (MDI)')
    plt.ylabel('Feature')
    plt.tight_layout()
    plt.savefig("../data/output/feature_importances_rf.pdf", format='pdf')
    #--------------------------end of block--------------------------------

    return feature_importance_df

def ridge_regression(feature_matrix, target_column='target'):
    """
    Train a Ridge regularized linear regression model on the feature matrix.
    Parameters:
    feature_matrix: pandas DataFrame with features as columns and ATAC-seq peaks as rows with target as log2 fold change
    target_column: column name for the target variable (default is 'target')
    """
    
    #remove inf and nan values from the feature matrix
    feature_matrix = feature_matrix.replace([np.inf, -np.inf], np.nan).dropna()

    # Split the data into features and target
    X_df = feature_matrix.drop(columns=[target_column])
    y = feature_matrix[target_column]

    scaler = StandardScaler()
    X = scaler.fit_transform(X_df)

    #subsample 70% of the data and run ridge regression 10 times with random seed --> save the average coefficients with standard deviation
    indices = [feature_matrix.sample(frac=0.7, random_state=42 + i).index for i in range(10)]

    coefficients = {col: [] for col in X_df.columns}

    # Loop through each set of sampled indices
    for index_list in indices:
        # Convert index_list to positional indices for numpy array X
        pos_indices = feature_matrix.index.get_indexer(index_list)
        temp_X = X[pos_indices, :]
        temp_y = y.iloc[pos_indices]
        
        # Fit ridge regression (adjust alpha as needed)
        ridge_model = Ridge(alpha=1.0)
        ridge_model.fit(temp_X, temp_y)
        
        # Store the coefficient for each feature
        for j, col in enumerate(X_df.columns):
            coefficients[col].append(ridge_model.coef_[j])

    #calculate the mean and standard deviation of the coefficients for each feature
    mean_coefficients = {col: np.mean(coefs) for col, coefs in coefficients.items()}
    std_coefficients  = {col: np.std(coefs)  for col, coefs in coefficients.items()}

    #sort the features by the mean coefficient in descending order
    sorted_features = sorted(mean_coefficients.items(), key=lambda x: x[1], reverse=True)

    # Print coefficients
    feature_importance_df = pd.DataFrame(sorted_features, columns=['feature', 'importance'])
    feature_importance_df['std'] = feature_importance_df['feature'].map(std_coefficients)

    return feature_importance_df

def plot_ridge(feature_importance_df):
    """
    Calculate feature importance from Ridge regression coefficient df
    """
    # Create a DataFrame for feature importance
    feature_importance_df = pd.DataFrame({
        'feature': feature_importance_df['feature'],
        'importance': feature_importance_df['importance'],
        'std': feature_importance_df['std']
    })

    #---------------block for plotting Ridge feature importances-------------------
    plot_df = feature_importance_df[['feature', 'importance']]
    plt.figure(figsize=(6, 20))
    # plot all values in a heatmap with a color gradient (zero as white in coolwarm)
    sns.heatmap(
        plot_df.set_index('feature'),
        annot=True,
        cmap='coolwarm',
        center=0,
        yticklabels=plot_df['feature']
    )
    plt.title('Feature Importances from Ridge Regression')
    plt.ylabel('Feature')
    plt.xlabel('Importance')
    plt.tight_layout()
    plt.savefig("../data/output/feature_importances_ridge.pdf", format='pdf')
    #--------------------------end of block--------------------------------

    return feature_importance_df
#--------------------------end of block--------------------------

def remove_features(feature_matrix, exclude_features):
    """
    Remove rows with 1 in specified columns from the feature matrix.
    Parameters:
    feature_matrix: pandas DataFrame with features as columns and ATAC-seq peaks as rows
    exclude_features: list of features to exclude from the feature matrix
    """
    #drop rows where any of the specified features are 1
    for feature in exclude_features:
        if feature in feature_matrix.columns:
            feature_matrix = feature_matrix[feature_matrix[feature] != 1]
        else:
            print(f"Feature '{feature}' not found in the feature matrix. Skipping.")
    #reset index
    feature_matrix.reset_index(drop=True, inplace=True)
    
    return feature_matrix

def main():
    # Define the list of bed files containing features
    bed_path = '../data/bed/'
    bed_list = [os.path.join(bed_path, f) for f in os.listdir(bed_path) if f.endswith('.bed') or f.endswith('.bed.gz')]

    # Parse command line arguments
    args = parse_args()
    control_ATAC = args.control_ATAC
    treatment_ATAC = args.treatment_ATAC
    exclude_features = args.exclude_features
    print(f"Excluding features: {exclude_features}")

    #save bed and bigwig file paths as variables
    control_bed, treatment_bed, control_bigwig, treatment_bigwig = modify_args(control_ATAC, treatment_ATAC)

    # Make peak set
    make_peak_set(control_bed, treatment_bed)

    # calculate log2 fold change for ATAC-seq peaks
    feature_df = calculate_log2_fold_change(control_bigwig, treatment_bigwig)

    # Make feature matrix
    feature_matrix = make_feature_matrix(bed_list, feature_df)

    #TODO: exclude features are are specified through command line arguments --> drop rows with 1 in specified columns
    if exclude_features:
        feature_matrix = remove_features(feature_matrix, exclude_features)

    #random forest classifier
    rf_classifier = random_forest(feature_matrix)
    rf_feature_importance_df = random_forest_feature_importance(rf_classifier)

    #ridge regularized linear regression
    ridge_df = ridge_regression(feature_matrix)
    ridge_feature_importance_df = plot_ridge(ridge_df)

    #concat feature importances --> column names are 'feature', 'MDI_importance', 'Ridge_coef'
    feature_importance_df = pd.concat([
        rf_feature_importance_df.rename(columns={'importance': 'MDI_importance'}),
        ridge_feature_importance_df.rename(columns={'importance': 'Ridge_coef'})
    ], axis=1)
    feature_importance_df = feature_importance_df.loc[:, ~feature_importance_df.columns.duplicated()]
    feature_importance_df.to_csv("../data/output/feature_importance.csv", index=False)

    return 


if __name__ == "__main__":
    main()