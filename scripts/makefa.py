# Input and output file paths
input_file = "../files/sim_reads.txt"  # Replace with the path to your input file
output_file = "../files/sim_reads.fa"  # Output FASTA file

# Open input and output files
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    # Read each line from the input file
    for i, line in enumerate(infile, start=1):
        # Clean up the line to remove any extra whitespace
        read = line.strip()
        if read:  # Ensure the line isn't empty
            # Write a header line and the read to the FASTA file
            outfile.write(f">read_{i}\n{read}\n")

print(f"FASTA file created: {output_file}")