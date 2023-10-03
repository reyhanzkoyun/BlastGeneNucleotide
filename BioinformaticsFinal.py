import matplotlib.pyplot as plt
import seaborn as sns
from Bio.Blast import NCBIWWW, NCBIXML
from tabulate import tabulate


def read_file_contents(file_path):
    """
    Read the contents of a file and return as a string.
    """
    try:
        with open(file_path, 'r') as file:
            file_contents = file.read()
        return file_contents
    except IOError:
        print(f"Error: Unable to read file {file_path}")
        return ""


def calculate_differences(seq1, seq2):
    """
    Calculate the differences in dinucleotide pairs between two sequences.
    """
    differences = 0
    for i in range(len(seq1) - 1):
        pair1 = seq1[i:i + 2]
        pair2 = seq2[i:i + 2]
        if pair1 != pair2:
            differences += 1
    return differences


def plot_dinucleotide_counts(seq, title):
    """
    Plot the counts of dinucleotide pairs in a sequence.
    """
    dinucleotides = [seq[i:i + 2] for i in range(len(seq) - 1)]
    counts = [dinucleotides.count(pair) for pair in set(dinucleotides)]

    plt.figure(figsize=(8, 6))
    sns.barplot(x=list(set(dinucleotides)), y=counts)
    plt.xlabel("Dinucleotide Pairs")
    plt.ylabel("Count")
    plt.title(title)
    plt.xticks(rotation=45)
    plt.show()


def run_blast(sequence):
    """
    Perform a BLAST search using the given sequence.
    """
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_records = NCBIXML.parse(result_handle)
    return blast_records


fasta_file1 = r"C:\Users\roman\PycharmProjects\pythonProject3\NBPF10.fasta"
fasta_file2 = r"C:\Users\roman\PycharmProjects\pythonProject3\NBPF11.fasta"

# Read the contents of the FASTA files
gene1_seq = read_file_contents(fasta_file1)
gene2_seq = read_file_contents(fasta_file2)

# Calculate the differences in dinucleotide pairs between the sequences
differences = calculate_differences(gene1_seq, gene2_seq)

# Display the number of differences
print("Number of differences without BLAST:", differences)

# Create two separate graphs for comparison

# Graph 1: NBPF10 Dinucleotide Counts
plot_dinucleotide_counts(gene1_seq, "NBPF10 Dinucleotide Counts")

# Graph 2: NBPF11 Dinucleotide Counts
plot_dinucleotide_counts(gene2_seq, "NBPF11 Dinucleotide Counts")

# Perform BLAST search for gene1_seq
blast_results = run_blast(gene1_seq)
for blast_record in blast_results:
    # Process the BLAST results as needed
    # Example: Print the alignment title and the first HSP (High-scoring Segment Pair)
    print("Alignment title:", blast_record.alignments[0].title)
    print("HSP bit score:", blast_record.alignments[0].hsps[0].bits)
