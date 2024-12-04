import pyranges as pr
import pandas as pd


class FilterTranscripts:
    def __init__(self, gtf_path, output_gtf_path):
        """
        Initialize the FilterTranscripts class.

        Input files:
        -----------
        gtf_path : str
            File path to input GTF file.
        output_gtf_path : str
            File path to output the filtered GTF file.
        """
        # Read the input GTF file into a pandas DataFrame
        self.gtf_pd = pr.read_gtf(gtf_path, as_df=True)
        # Store the output file path
        self.output_gtf = output_gtf_path

    def filter(self,
               gene_types=['lncRNA', 'protein_coding'],
               min_isoforms=2,
               transcript_column='transcript_id',
               gene_id_column='gene_id'):
        """
        Filter genes based on specified criteria.

        Input parameters:
        -----------
        gene_types : list, optional
            List of gene types to keep (default: lncRNA and protein_coding).
        min_isoforms : int, optional
            Minimum number of unique isoforms to keep per gene (default: 2).
        transcript_column : str, optional
            Column name for transcript identification.
        gene_id_column : str, optional
            Column name for gene identification.
        """
        # Validate required columns
        required_columns = [transcript_column, 'gene_type', gene_id_column, 'Feature']
        for col in required_columns:
            if col not in self.gtf_pd.columns:
                raise ValueError(f"Column '{col}' not found in the DataFrame")

        # Filter by gene types and feature type
        type_filtered_df = self.gtf_pd[
            (self.gtf_pd['gene_type'].isin(gene_types)) &
            (self.gtf_pd['Feature'] == 'transcript')
        ]

        # Count unique isoforms per gene
        isoform_counts = type_filtered_df.groupby(gene_id_column)[transcript_column].nunique()

        # Select genes with sufficient isoforms
        multi_isoform_genes = isoform_counts[isoform_counts >= min_isoforms].index

        # Filter the DataFrame based on the selected genes
        self.filtered_df = type_filtered_df[
            type_filtered_df[gene_id_column].isin(multi_isoform_genes)
        ]

    def save_filtered_gtf(self):
        """
        Save the filtered DataFrame as a GTF file to the specified output path.
        """
        if not hasattr(self, 'filtered_df'):
            raise ValueError("Filtered DataFrame not found. Run the 'filter' method first.")

        # Convert the filtered DataFrame to PyRanges for GTF output
        filtered_pr = pr.PyRanges(self.filtered_df)
        filtered_pr.to_gtf(self.output_gtf)

        print(f"Output filtered GTF: {self.output_gtf}")