from pyfaidx import Fasta

ref = Fasta('hg38.fa')

def filter_sequence(sequence):
    modified_seq = sequence.replace('N', 'A')
    return modified_seq

# Extract kromozoma 10 iz hg38.fa
def extract_chromosome_10(output_fasta_path):
    chrom = 'chr10'  
    orig_seq = ref[chrom][:].seq  # Cijeli chr10
    
    seq = []
    index = [i for i in range(len(orig_seq))]
    for ci, c in enumerate(orig_seq):
        if c.isupper():
            index[ci] = len(seq)
            seq.append(c)
    seq = ''.join(seq) # String samo velikih slova chr10

    with open(output_fasta_path, 'w') as output_fasta:
        output_fasta.write(f">{chrom}\n")
        
        # Zapisujemo chr10 u chr10_filtered.fa
        for i in range(0, len(seq), 60):
            output_fasta.write(seq[i:i+60] + '\n')
    
    return seq, orig_seq, index

def extract_duplications(seq, orig_seq, index, bed_file, output_fasta_file, extended_duplications_file, extended = 0):
    sequences = []
    extended_sequences = []

    with open(bed_file_path, 'r') as bed_file:
        for line in bed_file:
            fields = line.strip().split()

            if fields[8] == "+" and fields[9] == "+": 
                chrom1 = fields[0]
                chrom2 = fields[3]

                start1 = int(fields[1]) - 1 
                end1 = int(fields[2]) 
                start2 = int(fields[4]) - 1
                end2 = int(fields[5]) 

                if chrom1 == chrom2 and chrom1 == 'chr10': 
                    # Da bi pristupili ispravnim vrijednostima indexa, granice moraju biti velika slova
                    while not orig_seq[start1].isupper():
                        start1 -= 1
                    while not orig_seq[end1].isupper():
                        end1 += 1
                    while not orig_seq[start2].isupper():
                        start2 -= 1
                    while not orig_seq[end2].isupper():
                        end2 += 1
                    
                    # Koordinate granica u sekvenci samo velikih slova
                    start1 = index[start1]
                    end1 = index[end1]
                    start2 = index[start2]
                    end2 = index[end2]
                    
                    seq1 = seq[int(start1):int(end1)]  
                    seq2 = seq[int(start2):int(end2)]

                    n_count1 = seq1.count('N')
                    n_count2 = seq2.count('N')

                    if n_count1 <= 50 and n_count2 <= 50:
                        modified_seq1 = filter_sequence(seq1)
                        modified_seq2 = filter_sequence(seq2)

                        header1 = f">{chrom1}:{start1 + 1}-{end1}"
                        sequences.append((header1, modified_seq1))

                        header2 = f">{chrom2}:{start2 + 1}-{end2}"
                        sequences.append((header2, modified_seq2))

                        # Preracunavamo koordinate za prosirenu sekvencu 
                        dif = extended/100 * (end1 - start1) # za 30% sirenje, dif = 0.3*len(seq1)
                        start11 = int(start1 - dif)
                        end11 = int(end1 + dif)

                        dif = extended/100 * (end2 - start2)
                        start22 = int(start2 - dif)
                        end22 = int(end2 + dif) 

                        seq11 = seq[start11:end11]  
                        seq22 = seq[start22:end22]

                        modified_seq11 = filter_sequence(seq11)
                        modified_seq22 = filter_sequence(seq22)

                        header11 = f">{chrom1}:{start11 + 1}-{end11}"
                        extended_sequences.append((header11, modified_seq11))

                        header22 = f">{chrom2}:{start22 + 1}-{end22}"
                        extended_sequences.append((header22, modified_seq22))

    with open(output_fasta_file, 'w') as output_fasta:
        for header, seq in sequences:
            output_fasta.write(f"{header}\n{seq}\n")
    
    with open(extended_duplications_file, 'w') as output_fasta:
        for header, seq in extended_sequences:
            output_fasta.write(f"{header}\n{seq}\n")

if __name__ == "__main__":
    output_fasta_path = 'chr10_filtered.fa'
    seq, orig_seq, index = extract_chromosome_10(output_fasta_path)

    bed_file_path = 'final.bed'
    duplications_file = 'filtered_sequencespom.fa'
    extended_duplications_file = 'filtered_sequences_extended.fa'

    # Izvlacimo duplikacije bez malih slova i duplikacije prosirene za 30%
    extract_duplications(seq, orig_seq, index, bed_file_path, duplications_file, extended_duplications_file, 30)

