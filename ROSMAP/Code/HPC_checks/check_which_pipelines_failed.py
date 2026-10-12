import os
import pandas as pd

# Define the main folder and major cell type filenames
main_folder = "test"
cell_types = [
    "Immune_memory_usage.txt",
    "Inhibitory_memory_usage.txt",
    "OPCs_memory_usage.txt",
    "Oligodendrocytes_memory_usage.txt",
    "Vasculature_memory_usage.txt",
    "Astrocytes_memory_usage.txt",
    "Excitatory_memory_usage.txt",
]

# Initialize a list to store the results
results = []

# Loop through each set folder in the main folder
for set_folder in sorted(os.listdir(main_folder)):
    set_path = os.path.join(main_folder, set_folder)

    # Check if it's a directory
    if os.path.isdir(set_path):
        # Record the set name and check for the presence of each file
        set_results = {"Set": set_folder}
        for cell_type in cell_types:
            file_path = os.path.join(set_path, cell_type)
            set_results[cell_type] = os.path.isfile(file_path)  # True if the file exists, False otherwise
        results.append(set_results)

# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Replace True/False with "Ran"/"Not Ran"
results_df.replace({True: "Ran", False: "Not Ran"}, inplace=True)

# Save the summary to a CSV file
output_file = os.path.join(main_folder, "pipeline_run_summary.csv")
results_df.to_csv(output_file, index=False)

print(f"Summary table saved to {output_file}")