def read_fasta(fasta_file_path):
    sequences = {}
    current_chrom = None
    sequence_parts = []

    with open(fasta_file_path, 'r') as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                if current_chrom is not None:
                    sequences[current_chrom] = ''.join(sequence_parts)

                current_chrom = line[1:].split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line.strip().upper()) 

        if current_chrom is not None:
            sequences[current_chrom] = ''.join(sequence_parts)

    return sequences

# Filtriranje malih slova i pretvaranje 'N' u 'A'
def filter_sequence(sequence):
    #no_lower_case = ''.join(char for char in sequence if not char.islower())
    #modified_seq = no_lower_case.replace('N', 'A')
    modified_seq = sequence.replace('N', 'A')
    return modified_seq

def extract_sequences(bed_file_path, fasta_file_path, output_fasta_path):
    ref = read_fasta(fasta_file_path)
    sequences = []

    with open(bed_file_path, 'r') as bed_file:
        for line in bed_file:
            fields = line.strip().split()

            if fields[8] == "+" and fields[9] == "+": 
                chrom1 = fields[0]
                chrom2 = fields[3]

                # izvlacimo sekvence koje su za 30% (orginalne duzine) po strani duze
                start1 = int(fields[1]) - 1 
                end1 = int(fields[2]) 
                org_dif1 = 1.6*(end1 - start1) # da dodatno sirimo sekvencu ako je prosireni dio kasnije izbacen u modifikacijama
                dif = 0.3*(end1 - start1)
                start1 -= dif
                end1 += dif

                start2 = int(fields[4]) - 1
                end2 = int(fields[5]) 
                org_dif2 = 1.6*(end2 - start2)
                dif = 0.3*(end2 - start2)
                start2 -= dif
                end2 += dif

                if chrom1 == chrom2 and chrom1 == 'chr10':  
                    seq1 = ref[chrom1][int(start1):int(end1)]  
                    seq2 = ref[chrom2][int(start2):int(end2)]

                    modified_seq1 = filter_sequence(seq1)
                    modified_seq2 = filter_sequence(seq2)

                    """cnt = 0
                    while len(modified_seq1) < org_dif1:
                        seq1 = ref[chrom1][int(start1 - cnt):int(end1 + cnt)]
                        modified_seq1 = filter_sequence(seq1)
                        cnt += 1

                    cnt = 0
                    while len(modified_seq2) < org_dif2:
                        seq2 = ref[chrom2][int(start2 - cnt):int(end2 + cnt)]
                        modified_seq2 = filter_sequence(seq2)
                        cnt += 1"""

                    header1 = f">{chrom1}:{int(start1) + 1}-{int(end1)}"
                    sequences.append((header1, modified_seq1))

                    header2 = f">{chrom2}:{int(start2) + 1}-{int(end2)}"
                    sequences.append((header2, modified_seq2))

    with open(output_fasta_path, 'w') as output_fasta:
        for header, seq in sequences:
            output_fasta.write(f"{header}\n{seq}\n")

bed_file_path = 'final.bed'
fasta_file_path = 'hg38.fa'
output_fasta_path = 'filtered_sequences_extended.fa'

extract_sequences(bed_file_path, fasta_file_path, output_fasta_path)
