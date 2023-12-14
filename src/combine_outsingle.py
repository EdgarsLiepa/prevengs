import pandas as pd
import argparse

def combine_outsingle(file_zvalues, file_pvalues, output_file=False,gene_id_in_index=True):
    # Read the CSV files. Šī konkrēti priekš nerediģētiem outsingle failiem. Izkomentētās līnijas ir tad, ja
    # ja ir kkādi citi faili, bet tad funkcijā arī jānorāda gene_id kolonnu nosaukumi

    # if gene_id_in_index:
    df1 = pd.read_csv(file_zvalues, index_col=0, sep='\t').reset_index().rename(columns={'index': 'gene_id'})
    # print(df1)
    df2 = pd.read_csv(file_pvalues, index_col=0,sep='\t').reset_index().rename(columns={'index': 'gene_id'})
    combined_df = pd.merge(df1, df2, on='gene_id',suffixes=('_zscore', '_pvalue'))
    # else:
    #     df1 = pd.read_csv(file1,sep='\t')
    #     df2 = pd.read_csv(file2,sep='\t')
    #     combined_df = pd.merge(df1, df2, left_on=gene_id_col1, right_on=gene_id_col2,suffixes=('_zscore', '_pvalue'))

    # Write the combined data to a new CSV file
    if output_file:
        combined_df.to_csv('outsingle- z and p values.csv', index=False)
        print(f'Combined data written to outsingle- z and p values.csv')
    return combined_df

def main(args):
    combine_outsingle(args.path_to_zvalues_file,args.path_to_pvalues_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine outsingle results")
    
    parser.add_argument("path_to_zvalues_file", help="The Z-score file to process.")
    parser.add_argument("path_to_pvalues_file", help="The Z-score p value file")
    parser.add_argument("out_dir", help="The directory to save the results in.")

    args = parser.parse_args()
    
    main(args)