from python import math
from itertools import groupby, combinations
from collections import Counter, defaultdict
from python import Levenshtein as lev
from python import mmh3
from python import seahash
from python import Bio.SeqIO
import re
import time

def Pi(x):
    return mmh3.hash(x)

def f(x):
    """if isinstance(x, str):
        x = bytes(x, 'utf-8')
    return seahash.hash(x)"""
    return mmh3.hash(x) 

def find_anchors(s, T, Pi):
    A = [0]
    q = 3
    r = math.floor((len(s) - q + 1 - T) / (2 * T + 2))
    
    h = [Pi(s[i:i + q]) for i in range(len(s) - q + 1)]

    for i in range(r, len(s) - q + 1 - r):
        if all(h[i] < h[j] for j in range(i - r, i + r + 1) if j != i):
            A.append(i)

    A.append(len(s))
    return A

def partition_string(s, T, Pi):
    P = []
    A = find_anchors(s, T, Pi)

    for i in range(len(A) - 1):
        P.append((A[i], A[i + 1] - A[i]))

    return P

def Minjoin(S, K, T):
    O = set()
    C = []
    A = []

    for s in S:
        P = partition_string(s, T, Pi)
        for pos, length in P:
            A.append((f(s[pos:pos + length]), pos, length, S.index(s)))

    A.sort(key=lambda x: x[0])
    grouped = defaultdict(list)
    for key, group in groupby(A, key=lambda x: x[0]):
        grouped[key].extend(group)

    for G in grouped.values():
        for (hi, posi, leni, i), (hj, posj, lenj, j) in combinations(G, 2):
            if i != j and abs(len(S[i]) - len(S[j])) <= K:
                if abs(posi - posj) + abs((len(S[i]) - posi) - (len(S[j]) - posj)) <= K:
                    C.append((i, j))

    count = Counter(C)
    threshold = T / 20
    C = [pair for pair in C if count[pair] >= threshold]

    for (i, j) in C:
        if lev.distance(S[i], S[j]) <= K:
            O.add((i, j))

    return O

def split_into_windows(sequence, window_size=1000, step_size=500):
    windows = []
    for start in range(0, len(sequence) - window_size + 1, step_size):
        end = start + window_size
        windows.append(sequence[start:end])
    return windows

def compare_windows(seq1, seq2, windows1, windows2, K, T):
    results = []
    for win1 in windows1:
        for win2 in windows2:
            comparison_results = Minjoin([win1, win2], K, T)
            if comparison_results:
                results.append(comparison_results)
    return results

def process_pair(pair, K, T):
        (chrom1, start1, end1, sequence1) = pair[0]
        (chrom2, start2, end2, sequence2) = pair[1]

        if chrom1 == chrom2:
            windows1 = split_into_windows(sequence1)
            windows2 = split_into_windows(sequence2)
            comparisons = compare_windows(sequence1, sequence2, windows1, windows2, K, T)
            return comparisons
        
        return None

def compare_adjacent_sequences(sequences, K, T):
    results = []
    #pairs = [(sequences[i], sequences[i + 1]) for i in range(0, len(sequences) - 1, 2)]
    pairs = [(sequences[i], sequences[i + 1]) for i in range(0, 20, 2)]
    
    for pair in pairs:
        result = process_pair(pair, K, T)
        if result:
            results.append(result)
    
    return results

if __name__ == "__main__":
    start_time = time.time()
    output_fasta_path = 'filtered_sequences.fa'
    
    # Pattern za header
    pattern = re.compile(r"^(chr\d+):(\d+)-(\d+)$")
    sequences = []

    for record in SeqIO.parse(output_fasta_path, "fasta"):
        match = pattern.match(record.id)
        if match:
            chrom, start, end = match.groups()
            sequence = str(record.seq)

            sequences.append((chrom, int(start), int(end), sequence))
        else:
            print(f"Header se ne poklapa s zadanim formatom: {record.id}")

    K = 150.
    T = 20 + K / 8.
    
    comparison_results = compare_adjacent_sequences(sequences, K, T)

    cnt = 0
    for result in comparison_results:
        #chrom1, start1, end1 = result['sequence1']
        #chrom2, start2, end2 = result['sequence2']
        print(len(result))
    
    end_time = time.time()
    print(end_time - start_time)
