import pandas as pd
import argparse

def main():
    """
    Main wrapper script for the TranscriptAnalyzer
    """
    parser = argparse.ArgumentParser(description="Process transcript TPM data.")
    parser.add_argument("input_csv", type=str, help="Path to the input CSV file.")
    parser.add_argument("output_tsv", type=str, help="Path to the output TSV file.")
    
    args = parser.parse_args()
    
    # Load the input CSV file
    try:
        dataframe = pd.read_csv(args.input_csv)
    except Exception as e:
        print(f"Error loading input CSV: {e}")
        return
    
    # Process the data
    analyzer = TranscriptAnalyzer(dataframe)
    analyzer.filter_low_counts().compute_transcript_fractions()

    res = analyzer.differential_transcript_usage()
    
    # Save the results to the output file
    try:
        res.to_csv(args.output_tsv)
        print(f"Processed data saved to {args.output_tsv}")
    except Exception as e:
        print(f"Error saving output TSV: {e}")

if __name__ == "__main__":
    main()