import os
import pandas as pd

# config file
configfile: workflow.source_path("../../config/config.yaml")

# Define output directory
output_dir = config["output_dir"]

# Assembly name
asm_name = config["asm_name"]

# add file separator to output_dir if missing 
if not output_dir.endswith(os.path.sep):
    output_dir += os.path.sep

# load sample sheet
sample_sheet = pd.read_csv(config["sample_sheet"], dtype=str).set_index(
    ["sample"], drop=False
)

# get reads for trinity
R1 = sample_sheet["fastq_1"].dropna().tolist()
R2 = sample_sheet["fastq_2"].dropna().tolist()

# get memory for trinity
def get_mem_mb_trinity(wildcards, attempt):
    return 128000 + attempt * 32000
