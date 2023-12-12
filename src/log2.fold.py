# !pip install pydeseq2
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import pandas as pd
import numpy as np

def tpm_log2foldchange_table(samples_db_csv,refSTjude,geneId_col_name):
    """
    samples_db_csv: csv file with sample names as column_names, rows with GEneId. For each geneId, the TPM values for each sample are shown in the table
    refSTjude: csv file with sample names as column_names, rows with GEneId. For each geneId, the TPM values for each sample are shown in the table.
    """
    # Load the data
    samples_db = pd.read_csv(samples_db_csv)
    ref_st_jude = pd.read_csv(refSTjude)
    # Set geneId as index
    samples_db = samples_db.set_index(geneId_col_name)
    ref_st_jude = ref_st_jude.set_index(geneId_col_name)

    # Merge the dataframes
    merged_tpm = pd.concat([samples_db, ref_st_jude], axis=1)
    # converting tpms to integers for deseq2
    merged_tpm = merged_tpm.round().astype(int)
    # Filter out genes with no expression in all samples
    merged_tpm = merged_tpm[merged_tpm.sum(axis=1) > 0]

    #Set Conditions for Deseq object. Labeled as sample samples in samples_db_csv and as control samples in refSTjude
    conditions = ['sample'] * len(samples_db.columns) + ['control'] * len(ref_st_jude.columns)
    col_data = pd.DataFrame({'condition': conditions}, index=merged_tpm.columns)
    #Transpose df to have genes as columns and samples as rows
    merged_tpm = merged_tpm.T
    #Create Deseq object
    dds = DeseqDataSet(counts=merged_tpm,metadata=col_data,design_factors="condition")
    dds.deseq2()
    stat_res = DeseqStats(dds, n_cpus=8, contrast=('condition', 'sample', 'control'))
    stat_res.summary()
    res = stat_res.results_df
    # stat_res.run_deseq()
    # res = stat_res.get_deseq_result()
    # results = stat_res.deseq_result
    # For now just return the DESq2 results that also contains log2FoldChange values

    return res
    # return merged_tpm, col_data,dds,res


# results = tpm_log2foldchange_table("/rez/samples_db.csv", "/rez/refSTjude.csv", "GeneID")