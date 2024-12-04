import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

class TranscriptAnalyzer:
    def __init__(self, dataframe):
        """
        -----------
        dataframe : pd.DataFrame
            DataFrame containing transcript-level TPM data
            Expected columns: transcript_id, gene_id, followed by sample columns
        """
        # Ensure numeric processing
        self.original_df = dataframe.copy()
        self.processed_df = dataframe.copy()
        
        # Convert all sample columns to numeric, coercing errors to NaN
        sample_columns = self.processed_df.columns[2:]
        for col in sample_columns:
            self.processed_df[col] = pd.to_numeric(self.processed_df[col], errors='coerce')
        
        print(self.processed_df.head())
        
    def filter_low_counts(self, min_tpm_threshold=1):
        """
        Filter out transcripts with low total expression across all samples
        
        Parameters:
        -----------
        min_tpm_threshold : float, optional (default=1)
            Minimum total TPM threshold for a transcript to be retained
        
        Returns:
        --------
        self : TranscriptAnalyzer
            Returns self for method chaining
        """
        # Identify sample columns (exclude first two ID columns)
        sample_columns = self.processed_df.columns[2:]
        
        # Compute total TPM across all sample columns, handling potential NaNs
        self.processed_df['total_tpm'] = self.processed_df[sample_columns].sum(axis=1)
        
        # Filter transcripts above threshold
        self.processed_df = self.processed_df[self.processed_df['total_tpm'] > min_tpm_threshold]
        
        # Group by gene_id and filter genes with at least 2 transcripts
        gene_transcript_counts = self.processed_df.groupby('gene_id').size()
        valid_genes = gene_transcript_counts[gene_transcript_counts >= 2].index
        
        self.processed_df = self.processed_df[self.processed_df['gene_id'].isin(valid_genes)]
        
        # Drop the temporary total_tpm column
        self.processed_df = self.processed_df.drop(columns=['total_tpm'])
        
        return self
    
    def compute_transcript_fractions(self):
        """
        Convert TPM to fraction of total gene expression for each sample
        
        Returns:
        --------
        self : TranscriptAnalyzer
            Returns self for method chaining
        """
        # Identify sample columns (exclude first two ID columns)
        sample_columns = self.processed_df.columns[2:]
        
        # Robust method to compute gene-level totals
        def safe_gene_total(gene_group):
            # For each sample, sum the TPM values for the gene's transcripts
            return gene_group[sample_columns].sum()
        
        # Compute gene-level totals
        gene_totals = self.processed_df.groupby('gene_id').apply(safe_gene_total)
        
        # Create a copy of the processed dataframe to store fractions
        fraction_df = self.processed_df.copy()
        
        # Compute fractions for each transcript in each sample
        for sample in sample_columns:
            # Divide transcript TPM by total gene TPM for that sample
            # Use numpy.divide with a small epsilon to avoid division by zero
            fraction_df[sample] = fraction_df.apply(
                lambda row: np.divide(
                    row[sample], 
                    gene_totals.loc[row['gene_id'], sample], 
                    out=np.zeros_like(row[sample], dtype=float), 
                    where=gene_totals.loc[row['gene_id'], sample] != 0
                ), 
                axis=1
            )
        
        self.processed_df = fraction_df
        self.processed_df.to_csv("TPM_transformed.csv",sep="\t")
        print("Writing fractions to CSV...")
        return self
    def differential_transcript_usage(self, condition_map, method='wilcoxon'):
        """
        Perform statistical analysis to detect differential transcript usage
        
        Parameters:
        -----------
        condition_map : dict
            Mapping of sample names to their conditions (e.g., {'Normal1': 'Normal', 'Tumor1': 'Tumor'})
        method : str, optional (default='wilcoxon')
            Statistical test method. Options: 'wilcoxon' or 't-test'
        
        Returns:
        --------
        results_df : pd.DataFrame
            DataFrame with statistical results for each transcript
        """
        # Validate condition map covers all samples
        sample_columns = self.processed_df.columns[2:]
        if set(condition_map.keys()) != set(sample_columns):
            raise ValueError("Condition map must exactly match sample columns")
        
        # Separate normal and tumor samples
        normal_samples = [col for col, cond in condition_map.items() if cond == 'Normal']
        tumor_samples = [col for col, cond in condition_map.items() if cond == 'Tumor']
        
        # Compute mean fractions for each condition
        results = []
        for _, row in self.processed_df.iterrows():
            # Extract fractions, removing any NaN values
            normal_fractions = row[normal_samples]
            tumor_fractions = row[tumor_samples]
            
            # Remove NaN values before statistical test
            normal_fractions = pd.to_numeric(normal_fractions, errors='coerce')
            tumor_fractions = pd.to_numeric(tumor_fractions, errors='coerce')

            normal_fractions = normal_fractions[normal_fractions.notna()]
            tumor_fractions = tumor_fractions[tumor_fractions.notna()]
            
            # Skip if either group is empty
            if len(normal_fractions) == 0 or len(tumor_fractions) == 0:
                continue
            
            # Perform statistical test based on selected method
            try:
                if method == 'wilcoxon':
                    statistic, p_value = stats.mannwhitneyu(
                        normal_fractions, tumor_fractions, 
                        alternative='two-sided'
                    )
                    test_name = 'Wilcoxon Rank-Sum'
                elif method == 't-test':
                    statistic, p_value = stats.ttest_ind(
                        normal_fractions, tumor_fractions, 
                        equal_var=False
                    )
                    test_name = "Welch's t-test"
                else:
                    raise ValueError("Method must be 'wilcoxon' or 't-test'")
                
                results.append({
                    'transcript_id': row['transcript_id'],
                    'gene_id': row['gene_id'],
                    'mean_normal_fraction': normal_fractions.mean(),
                    'mean_tumor_fraction': tumor_fractions.mean(),
                    'median_normal_fraction': np.median(normal_fractions),
                    'median_tumor_fraction': np.median(tumor_fractions),
                    'log2_fold_change': np.log2(
                        (tumor_fractions.mean() + 1e-10) / (normal_fractions.mean() + 1e-10)
                    ),
                    'test_statistic': statistic,
                    'p_value': p_value,
                    'test_method': test_name
                })
            except Exception as e:
                print(f"Skipping transcript due to error: {e}")
                continue
        
        # Convert to DataFrame and apply FDR correction
        results_df = pd.DataFrame(results)
        
        # Only apply FDR correction if there are results
        if not results_df.empty:
            _, results_df['adjusted_p_value'] = fdrcorrection(results_df['p_value'])
            # Sort by adjusted p-value
            results_df = results_df.sort_values('adjusted_p_value')
        
        return results_df
    