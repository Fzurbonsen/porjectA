import matplotlib.pyplot as plt
import numpy as np
import re
import os

def extract_cigar_strings(filename):
    gssw_cigars = []
    gwfa_cigars = []

    with open(filename, 'r') as file:
        for line in file:
            if line.startswith('gssw:'):
                parts = line.split()
                if len(parts) > 2:
                    gssw_cigars.append(parts[2])
            elif line.startswith('gwfa:'):
                parts = line.split()
                if len(parts) > 2:
                    gwfa_cigars.append(parts[2])

    return gssw_cigars, gwfa_cigars

def parse_cigar(cigar):
    return re.findall(r'(\d+)([MID])', cigar)

def plot_cigar_matrix(cigar, ref, read, file_path):
    operations = parse_cigar(cigar)

    ref_index = 0
    read_index = 0

    alignment_matrix = []

    for length, op in operations:
        length = int(length)
        if op == 'M':  # Match or mismatch
            for i in range(length):
                alignment_matrix.append((ref_index, read_index, 'M'))
                ref_index += 1
                read_index += 1
        elif op == 'I':  # Insertion to the reference
            for i in range(length):
                alignment_matrix.append((None, read_index, 'I'))
                read_index += 1
        elif op == 'D':  # Deletion from the reference
            for i in range(length):
                alignment_matrix.append((ref_index, None, 'D'))
                ref_index += 1

    max_ref_index = max((i for i, j, op in alignment_matrix if i is not None), default=0)
    max_read_index = max((j for i, j, op in alignment_matrix if j is not None), default=0)

    matrix = np.full((max_ref_index + 1, max_read_index + 1), ' ', dtype='<U1')

    for ref_i, read_i, op in alignment_matrix:
        if ref_i is not None and read_i is not None:
            matrix[ref_i, read_i] = '\\'  # Match/Mismatch

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.imshow(matrix != '\\', cmap='binary_r', interpolation='none', aspect='auto')

    # ax.set_xticks(np.arange(max_read_index + 1))
    # ax.set_yticks(np.arange(max_ref_index + 1))

    # ax.set_xticklabels(list(read[:max_read_index + 1]))
    # ax.set_yticklabels(list(ref[:max_ref_index + 1]))

    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    plt.xticks(rotation=90)

    # Create directory if it doesn't exist
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    # Save the figure
    plt.savefig(file_path)
    plt.close()

# Example usage
# cigar1 = "1M1D5M2D1M2D10M1D3M1D2M1D3M4D6M1D1M1D1M2D4M1D1M4D2M1D1M1D7M1D14M1D2M1D3M1D4M1D3M1D3M1D1M1D3M1D1M1D3M1D1M2D4M4D2M1D12M3D5M3D2M1D1M4D5M1D1M1D11M1D1M1D2M3D5M3D1M1D8M1D2M1D1M1D10M2D7M1D6M1D4M1D3M2D20M2D4M1D7M2D3M1D7M2D9M2D3M1D3M1D2M1D1M2D3M1D2M1D5M1D3M1D3M1D1M1D5M1D4M2D6M1D1M1D5M5D6M1D1M1D12M1D3M1D1M1D7M1D5M2D1M6D5M1D1M1D3M1D1M1D5M3D4M1D11M2D2M1D8M1D5M1D14M1D3M1D1M1D4M1D6M2D1M1D10M1D6M1D3M1D5M2D4M1D2M1D9M1D2M1D1M1D5M2D5M2D3M1D4M1D2M1D18M1D1M2D5M1D2M1D3M1D10M1D2M2D2M3D1M1D1M1D10M1D6M1D8M3D5M1D6M2D2M1D2M1D3M4D6M1D1M1D1M2D5M4D11M1D1M1D5M1D7M1D4M1D5M1D2M3D3M1D2M1D4M1D4M2D6M3D8M1D2M2D3M3D1M1D3M3D2M1D3M1D6M2D10M1D2M1D3M3D3M1D1M1D3M1D1M3D13M1D6M1D3M1D2M1D2M1D7M2D1M1D3M2D8M1D3M1D5M1D7M1D5M1D6M1D12M1D7M3D1M1D2M1D2M1D7M1D2M1D3M1D1M1D1M1D3M2D1M1D3M2D1M1D6M1D6M2D4M2D5M2D2M1D1M1D12M1D3M1D1M1D1M2D1M1D1M2D2M1D1M1D7M1D3M1D3M1D2M1D8M1D5M1D5M1D3M3D1M3D1M1D14M2D4M1D12M1D1M1D11M1D2M1D1M1D2M1D1M1D3M1D3M3D4M1D1M1D3M1D2M1D1M1D6M1D10M1D1M1D1M1D3M1D1M1D10M1D10M1D5M1D3M1D4M1D2M1D5M2D1M1D2M2D2M2D3M1D2M1D4M1D7M"
# plot_cigar_matrix(cigar1, "ref1", "read1", "./file1.png")

# cigar2 = "2D1M1D1M3D1M1D4M5D1M4D1M5D1M5D1M4D1M1D1M13D1M4D1M1D2M1D1M7D2M1D1M2D5M2D1M13D1M2D2M3D1M5D1M1D2M5D1M2D2M4D1M2D2M2D1M1D1M1D3M3D2M1D1M6D1M1D1M6D1M3D2M3D46M1D1M1I47M1D135M1D1M2D1M1D2M1D1M1D2M4D2M4D1M1D1M7D1M1D1M4D1M1D4M4D2M1D1M1D1M1D3M1D1M3D2M8D1M3D1M6D1M6D2M4D1M1D3M3D2M3D1M2D1M2D1M4D1M3D1M2D1M3D1M5D1M3D1M2D3M3D1M7D1M1D1M1D2M2D2M3D1M3D1M1D1M4D1M2D5M1D1M1D1M3D1M1D1M5D1M1D2M2D1M2D2M6D1M3D1M6D2M2D1M6D1M2D9M8D1M1D2M2D1M1D3M2D1M8D2M4D1M1D1M2D3M1D2M5D1M1D1M2D3M1D2M13D1M3D1M3D1M1D1M4D2M3D3M1D1M2D1M2D1M1D1M1D2M4D1M1D1M2D1M2D1M1D1M2D3M11D1M1D4M1D1M2D1M4D2M5D1M1D1M6D1M4D1M3D1M6D3M4D1M1D1M6D2M1D1M6D2M1D1M3D1M2D1M4D1M3D2M5D1M3D2M1D2M6D26M1I2M1D13M1I53M2I5M1D2M1D124M1I1M1D73M4D1M4D1M7D1M1D1M4D1M4D1M3D1M2D1M1D1M5D1M5D1M5D2M7D1M1D4M2D1M1D1M5D2M1D3M1D1M1D1M3D1M2D1M2D2M1D5M2D1M3D6M11D1M3D1M1D1M3D1M1D1M5D1M1D1M2D2M1D1M1D1M3D1M1D2M2D1M5D1M1D1M2D1M18D1M9D2M1D2M1D2M3D2M4D1M2D1M11D1M7D1M2D1M4D2M2D1M5D1M2D1M6D1M3D1M6D1M1D1M4D1M3D2M1D1M3D1M2D1M3D1M6D1M3D1M11D1M1D2M10D3M2D2M6D1M1D2M1D3M2D1M2D3M3D5M9D1M3D2M1D1M1D2M4D2M6D1M1D1M2D8M6D1M1D1M1D1M6D1M12D2M1D1M5D1M1D1M6D19M3D1M18D1M5D1M1D1M4D1M2D1M3D1M5D1M5D4M3D2M2D1M1D2M11D2M1D2M6D1M4D1M1D2M1D1M11D2M6D1M12D1M1D7M2D1M1D2M3D2M3D1M1D1M1D3M4D1M2D2M3D2M2D2M3D1M2D2M7D1M6D2M1D1M11D2M3D1M5D1M2D2M1D1M2D3M4D3M6D1M4D1M4D1M54D"
# plot_cigar_matrix(cigar2, "ref2", "read2", "./file2.png")

gssw_cigars, gwfa_cigars = extract_cigar_strings("./cigars.txt")

for i in range(len(gssw_cigars)):
    plot_cigar_matrix(gssw_cigars[i], "", "", f"./img/gssw/gssw_cigar{i}.png")
    plot_cigar_matrix(gwfa_cigars[i], "", "", f"./img/gwfa/gwfa_cigar{i}.png")