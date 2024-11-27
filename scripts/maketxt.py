def convert_to_fasta(input_file, output_file):
    """
    Convert a line-by-line sequence file into a proper FASTA file.
    
    Parameters:
    input_file (str): Path to the input file with sequences.
    output_file (str): Path to the output FASTA file.
    """
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            for line in infile:
                line = line.strip()
                if line.startswith(">"):  # Header line in FASTA
                    outfile.write(f"{line}\n")
                else:  # Sequence line
                    outfile.write(f"{line}\n")
        print(f"FASTA conversion completed. Output saved to {output_file}")
    except Exception as e:
        print(f"An error occurred: {e}")

# File paths
input_path = "../files/sim_reads_new.fa"  # Replace with your input file path
output_path = "../files/sim_reads_new.txt"  # Replace with your desired output file path

# Run the conversion
convert_to_fasta(input_path, output_path)