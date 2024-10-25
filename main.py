from Bio.Seq import Seq
from Bio import SeqIO


dna_sequence = Seq('ATGAGCCAGAAATACCTCAACCAGATTCGAGAGGAAACTACTCAGGCAACCGAACTGTCTGCCTTCCTGGAGAAACAGTTTGGGGATCAAGGCAAATACGTAGCCGATGTTTCAGAGCGAACATCAACGGACTACACCGGCGACTCCCAGACCAGCATTCAGAACTAA')

rna_sequence = dna_sequence.transcribe()
protein = rna_sequence.translate()

print(f'DNA sequence {dna_sequence}')
print(f'RNA sequence: {rna_sequence}')
print(f'Protein: {protein}')


amino_acid = 'T' # Threonine
amino_acid_count = protein.count(amino_acid)
amino_acid_percentage = round(amino_acid_count / ((len(protein) - 1)) * 100)

print(f'Count of amino acid {amino_acid} in protein: {amino_acid_count}')
print(f'Percentage of {amino_acid} in protein: {amino_acid_percentage}%')


glutamine_codons = ['AAC', 'AGC']
glutamine_codons_count = 0
j = 0


while j < len(rna_sequence) - 1:
    if rna_sequence[j] + rna_sequence[j + 1] + rna_sequence[j + 2] in glutamine_codons:
        glutamine_codons_count += 1

    j += 3

print(f'Count of glutamine codons in mRNA: {glutamine_codons_count}')


# File convert
SeqIO.convert('sequence.fastq',  'fastq', 'sequence.fasta', 'fasta')
