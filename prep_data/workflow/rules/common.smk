import os
import pandas as pd

# config file
configfile: workflow.source_path("../../config/config.yaml")

# Define output directory
output_dir = config["output_dir"]

# add file separator to output_dir if missing 
if not output_dir.endswith(os.path.sep):
    output_dir += os.path.sep

# load sample sheet
sample_sheet = pd.read_csv(config["sample_sheet"], dtype=str).set_index(
    ["sample"], drop=False
)
# get sample names
samples = sorted(set(sample_sheet.index))

# get fastqs for fastp
def get_fastq(wildcards):
    """Get fastq files of given sample"""
    fastqs = sample_sheet.loc[wildcards.sample, ["fastq_1", "fastq_2"]].dropna()
    return [fastqs.fastq_1, fastqs.fastq_2]

# get memory for seal
def get_mem_mb_seal(wildcards, attempt):
    return 10000 + attempt * 10000