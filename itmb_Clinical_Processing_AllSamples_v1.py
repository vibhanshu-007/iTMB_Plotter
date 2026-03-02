# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 12:43:52 2024

Main script TMB_Plotter (Satya Prakash K)
Below Script is to run for All Samples (Vyomesh J)
"""

import pandas as pd
import os
import sys

# Get the input file from the first argument
if len(sys.argv) < 2:
    print("Please provide the input file as the first argument.")
    sys.exit(1)

input_file = sys.argv[1]

# Check the file extension and read accordingly
file_extension = os.path.splitext(input_file)[-1].lower()

# Read the input file based on its extension
if file_extension == '.xlsx':
    df = pd.read_excel(input_file, engine='openpyxl',header=0)
elif file_extension == '.csv':
    df = pd.read_csv(input_file)
else:
    print("Unsupported file format. Please provide a .csv or .xlsx file.")
    sys.exit(1)
 
print(df)
# Iterate over each row of the DataFrame
for index, row in df.iterrows():
    batch = row['Batch']
    sample_name = row['Sample_Name']
    cancer_type = row['Broad_Category_Cancer_Type']
    tmb_score = row['TMB_Score']
    
    # Check if TMB_Score is not missing or 'NA'
    if pd.isna(tmb_score):
        print(f"Skipping row {index + 1} due to missing TMB_Score.")
        continue

    # Combine Batch and Sample_Name into a single string
    sample_combined = f"{sample_name}_{batch}"
                
    # Prepare the command string to call script_2.py with the arguments
    command = f"python /home/basecare/Patient_samples/01_iTMB_Graph/Common_Files/itmb_plotter2.py --file /home/basecare/Patient_samples/01_iTMB_Graph/Common_Files/itmb_final.tsv --score {tmb_score} --cancer \"{cancer_type}\" --sample {sample_combined}"
    
    # Execute the command using os.system
    exit_code = os.system(command)
    
    # Check the exit code and print appropriate messages
    if exit_code == 0:
        print(f"Successfully processed row {index + 1}: Sample {sample_combined}")
    else:
        print(f"Error processing row {index + 1}. Exit code: {exit_code}")


### TMB Percentile Mapping (Added Below) - Karthikasri - 07-August-2025 ###


# Load percentile mapping data (Sheet2 of the same Excel input)
try:
    if file_extension == '.xlsx':
        df_percentiles = pd.read_excel(input_file, sheet_name='Sheet2', engine='openpyxl')
    else:
        print("Percentile mapping only supported for .xlsx input files with Sheet2.")
        sys.exit(1)
except Exception as e:
    print(f"Failed to read Sheet2 for percentile mapping: {e}")
    sys.exit(1)

# Define percentiles columns for reference
percentiles_columns = [
    "5th", "10th", "15th", "20th", "25th", "30th", "35th", "40th", "45th", "50th",
    "55th", "60th", "65th", "70th", "75th", "80th", "85th", "90th", "95th", "100th"
]

# Function to determine TMB_Percentile
def get_percentile(cancer_type, tmb_score, percentiles_df):
    cancer_row = percentiles_df[percentiles_df["Cancer"] == cancer_type]
    if cancer_row.empty or pd.isna(tmb_score):
        return None
    for col in percentiles_columns:
        if tmb_score <= cancer_row[col].values[0]:
            return int(col.replace("th", ""))
    return 100

# Add TMB_Percentile column
df["TMB_Percentile"] = df.apply(
    lambda row: get_percentile(row["Broad_Category_Cancer_Type"], row["TMB_Score"], df_percentiles), axis=1
)

# Save the updated DataFrame with percentiles
output_file_path = "TMB_Percentile_Output_Final.xlsx"
df.to_excel(output_file_path, index=False)

print(f"\nTMB Percentile output saved to: {output_file_path}")
