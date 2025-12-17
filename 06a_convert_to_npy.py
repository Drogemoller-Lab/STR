import os
import numpy as np

# Change working directory
os.chdir("/home/projects/hearing_loss/clsaARHL_SA/str/phenotype")

# Define file paths (now they are relative to the current working directory)
file_pairs = [
    ("met_merge_filt_males.txt", "met_merge_filt_males.npy"),
    ("met_merge_filt_females.txt", "met_merge_filt_females.npy"),
    ("sen_merge_filt_males.txt", "sen_merge_filt_males.npy"),
    ("sen_merge_filt_females.txt", "sen_merge_filt_females.npy"),

]
# Process each file
for input_file, output_file in file_pairs:
    try:
        # Load the text file as strings
        data = np.loadtxt(input_file, delimiter=' ', dtype=str)

        # Replace 'NA' with NaN
        data[data == 'NA'] = np.nan

        # Convert to float
        data = data.astype(float)

        # Save as .npy
        np.save(output_file, data)
        print(f"Data saved as {output_file}")
    except Exception as e:
        print(f"Error processing {input_file}: {e}")
