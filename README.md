# A.glyphus_Bd_pre-exposure_RNAseq

Contains code for RNAseq analysis of *Atelopus glyphus* for the project: **Prior-exposing highly susceptible harlequin toads to Batrachochytrium dendrobatidis worsens disease outcomes upon re-exposure**. 

The repo contains the following code directories:
1. **prep_data**: snakemake pipeline to trim reads and remove contaminant and rRNA sequences.
2. **tr_ass**: snakemake pipeline to assemble de novo transcriptome, annotate transcripts and run QC.
3. **quant**: snakemake pipeline to quantify transcript expression.
4. **DEG_analysis**: R code for differential expression analyses.
5. **data**: Metadata for R analyses