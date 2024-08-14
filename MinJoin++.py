from itertools import groupby, combinations
from collections import Counter, defaultdict
from python import Levenshtein as lev
from python import mmh3
from python import Bio.SeqIO
import re
import time

def Pi(x):
    return mmh3.hash(x)

def sea_swap64(x):
    return (((x << 56u64) & 0xff00000000000000u64) | 
            ((x << 40u64) & 0x00ff000000000000u64) | 
            ((x << 24u64) & 0x0000ff0000000000u64) |
            ((x << 8u64)  & 0x000000ff00000000u64) |
            ((x >> 8u64)  & 0x00000000ff000000u64) |
            ((x >> 24u64) & 0x0000000000ff0000u64) |
            ((x >> 40u64) & 0x000000000000ff00u64) |
            ((x >> 56u64) & 0x00000000000000ffu64))

def diffuse(val):
    val *= 0x6eed0e9da4d94a4fu64
    a = val >> 32u64
    b = val >> 60u64
    val ^= a >> b
    val *= 0x6eed0e9da4d94a4fu64
    return val

def f(key, seed = 0):
    a = 0x16f11fe89b0d677cu64 ^ u64(seed)
    b = 0xb480a793d8e6c86cu64
    c = 0x6fe2e5aaf078ebc9u64
    d = 0x14f994a4c5259381u64
    pad = [u64(0)]

    p = Ptr[u64](key.ptr)
    i, l = 0, len(key)
    while l >= 32:
        a ^= p[i]
        b ^= p[i + 1]
        c ^= p[i + 2]
        d ^= p[i + 3]
        a = diffuse(a)
        b = diffuse(b)
        c = diffuse(c)
        d = diffuse(d)
        i += 4
        l -= 32
    if 25 <= l <= 31:
        a ^= p[i]
        b ^= p[i + 1]
        c ^= p[i + 2]
        str.memcpy(pad.arr.ptr.as_byte(), (p + 3).as_byte(), l - 24)
        d ^= pad[0]
        a = diffuse(a)
        b = diffuse(b)
        c = diffuse(c)
        d = diffuse(d)
    elif l == 24:
        a ^= p[i]
        b ^= p[i + 1]
        c ^= p[i + 2]
        a = diffuse(a)
        b = diffuse(b)
        c = diffuse(c)
    elif 17 <= l <= 23:
        a ^= p[i]
        b ^= p[i + 1]
        str.memcpy(pad.arr.ptr.as_byte(), (p + 2).as_byte(), l - 16)
        c ^= pad[0]
        a = diffuse(a)
        b = diffuse(b)
        c = diffuse(c)
    elif l == 16:
        a ^= p[i]
        b ^= p[i + 1]
        a = diffuse(a)
        b = diffuse(b)
    elif 9 <= l <= 15:
        a ^= p[i]
        str.memcpy(pad.arr.ptr.as_byte(), (p + 1).as_byte(), l - 8)
        b ^= pad[0]
        a = diffuse(a)
        b = diffuse(b)
    elif l == 8:
        a ^= p[i]
        a = diffuse(a)
    elif 1 <= l <= 7:
        str.memcpy(pad.arr.ptr.as_byte(), p.as_byte(), l)
        a ^= pad[0]
        a = diffuse(a)
        a ^= b
        c ^= d
        a ^= c
        a ^= u64(len(key))
    return sea_swap64(diffuse(a))

def floor(x):
    return int(x) if x >= 0 else int(x) - 1 if x != int(x) else int(x)

def find_anchors(s, T, Pi):
    A = [0]
    q = 3
    r = floor((len(s) - q + 1 - T) / (2 * T + 2))
    
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

def compare_windows(windows1, windows2, K, T):
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
            comparisons = compare_windows(windows1, windows2, K, T)
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
    #print(len(sequences))

    cnt = 0
    for result in comparison_results:
        #chrom1, start1, end1 = result['sequence1']
        #chrom2, start2, end2 = result['sequence2']
        print(len(result))

    end_time = time.time()
    print(end_time - start_time)

