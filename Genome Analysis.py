import os

# Function to find motifs in a DNA sequence
def find_motifs(dna_sequence, motif):
    motifs = []
    sequence_length = len(dna_sequence)
    motif_length = len(motif)

    for i in range(sequence_length - motif_length + 1):
        if dna_sequence[i:i + motif_length] == motif:
            motifs.append(i)

    return motifs

# Function to find the reverse complement of a DNA sequence
def reverse_complement(dna_sequence):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(dna_sequence))

# Function to find open reading frames (ORFs) in a DNA sequence
def find_orfs(dna_sequence, start_codon="ATG", stop_codons=["TAA", "TAG", "TGA"]):
    orfs = []
    i = 0

    while i < len(dna_sequence):
        start = dna_sequence.find(start_codon, i)
        if start == -1:
            break
        end = len(dna_sequence)

        for stop_codon in stop_codons:
            stop = dna_sequence.find(stop_codon, start + 3)
            if stop != -1 and (stop - start) % 3 == 0:
                end = min(end, stop)

        if end - start > 3:  # ORF must be at least 3 codons long
            orfs.append((start, end + 3))
        i = start + 1

    return orfs

# Function to translate an ORF into amino acids
def translate_orf(orf):
    genetic_code = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }

    amino_acids = []
    for i in range(0, len(orf), 3):
        codon = orf[i:i+3]
        amino_acid = genetic_code.get(codon, 'X')  # 'X' for unknown
        amino_acids.append(amino_acid)

    return ''.join(amino_acids)

# Main program
if __name__ == "__main__":
    print("Welcome to the Even More Advanced Genome Sequence Analysis Tool")

    # Get the script's directory
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Construct the full file path for 'dna_sequence.txt'
    file_path = os.path.join(script_dir, "dna_sequence.txt")

    try:
        # Open 'dna_sequence.txt'
        with open(file_path, "r") as file:
            dna_sequence = file.read()

        # Input motif to search for
        motif = "AGC"

        # Find motifs
        motifs = find_motifs(dna_sequence, motif)

        if len(motifs) > 0:
            print(f"The motif '{motif}' was found at positions: {motifs}")
        else:
            print(f"The motif '{motif}' was not found in the DNA sequence.")

        # Find ORFs
        orfs = find_orfs(dna_sequence)
        if len(orfs) > 0:
            print("Open Reading Frames (ORFs) found:")
            for i, (start, end) in enumerate(orfs):
                orf_sequence = dna_sequence[start:end]
                amino_acids = translate_orf(orf_sequence)
                print(f"ORF {i + 1}: Start position: {start}, End position: {end}, Amino Acids: {amino_acids}")
        else:
            print("No ORFs found in the DNA sequence.")
    except FileNotFoundError:
        print("Error: 'dna_sequence.txt' not found. Make sure the file exists in the script's directory.")
