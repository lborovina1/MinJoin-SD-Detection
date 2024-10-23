from itertools import groupby, combinations
from collections import Counter, defaultdict
import re
import time

def Pi(x, seed = 0):
    m, r = 0xc6a4a7935bd1e995, 47
    h = seed ^ (8*m)

    k = 0
    for char in x:
        k = (k * 256 + ord(char)) & 0xFFFFFFFF

    k *= m
    k ^= k >> r
    k *= m
    h ^= k
    h *= m
    h ^= h >> r
    h *= m
    h ^= h >> r
    return h

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

def expand(chunks1, chunks2, hashed_partitions, i, j, num_same_segments, K, threshold, q = 8, r = 10, percentage_error = 15):
    first_cp = chunks1[i] # kopija prvog segmenta
    second_cp = chunks2[j] # kopija drugog segmenta
    counter = 0
    expand_left = True # za granicne slucajeve
    expand_right = True
    changed = False

    # indeksi segmenata oko duplikacija na koje sirimo
    curr_before_first = i - 1
    curr_after_first = i + 1
    curr_before_second = j - 1
    curr_after_second = j + 1

    if curr_before_first == 0 and (curr_after_first + 1)*1000 < len(chunks1):
        i += 1
        curr_before_first += 1
        curr_after_first += 1

    if curr_before_second == 0 and (curr_after_second + 1)*1000 < len(chunks2):
        j += 1
        curr_before_second += 1
        curr_after_second += 1

    expanded_left = -1 # uzimamo zadnji element prethodnog segmenta
    expanded_right = 0 # uzimamo prvi element sljedeceg segmenta

    # najmanji i najveci pivoti po poziciji
    filtered_first = [elem for elem in hashed_partitions if elem[-1] == i]
    filtered_second = [elem for elem in hashed_partitions if elem[-1] == j]
    first_sorted = sorted(filtered_first, key=lambda x: x[1]) # sortiramo po pozicijama pivota
    second_sorted = sorted(filtered_second, key=lambda x: x[1])

    smallest_pivot_first = first_sorted[0] # sortiramo po poziciji u segmentima
    largest_pivot_first = first_sorted[-1]
    if len(first_sorted) > 1:
        second_largest_first = first_sorted[-2] # preuzima najveci ako on izgubi status pivota
    else:
        second_largest_first = first_sorted[-1]
    smallest_pivot_second = second_sorted[0]
    largest_pivot_second = second_sorted[-1]
    if len(first_sorted) > 1:
        second_largest_second = second_sorted[-2] # preuzima najveci ako on izgubi status pivota
    else:
        second_largest_second = second_sorted[-1]

    # racuna u odnosu na pocetak prosirene sekvence -> za globalnu poziciju mora dodati na pocetak prosirenog segmenta
    first_before_last_it = (i * 1000 + ((curr_before_first - i + 1) * 1000 + expanded_left), (i + 1) * 1000 - 1 + (curr_after_first - i - 1) * 1000 + expanded_right)
    second_before_last_it = (j * 1000 + ((curr_before_second - j + 1) * 1000 + expanded_left), (j + 1) * 1000 - 1 + (curr_after_second - j - 1) * 1000 + expanded_right)

    while num_same_segments >= int(threshold):
        # indeksi koliko je do sada prosireno za povratak iz funkcije
        first_before_last_it = (i * 1000 + ((curr_before_first - i + 1) * 1000 + expanded_left), (i + 1) * 1000 - 1 + (curr_after_first - i - 1) * 1000 + expanded_right)
        second_before_last_it = (j * 1000 + ((curr_before_second - j + 1) * 1000 + expanded_left), (j + 1) * 1000 - 1 + (curr_after_second - j - 1) * 1000 + expanded_right)

        if counter % 2 == 0 and expand_left == True: #jednom siri lijevo, drugi put desno
            #ne siri ako ne moze vise lijevo -> pocetak cijelog hromozoma
            if curr_before_first < 0 or curr_before_second < 0:
                expand_left = False
                continue

            # sirimo obje duplikacije za jedan element lijevo
            first_cp = chunks1[curr_before_first][expanded_left] + first_cp
            second_cp = chunks2[curr_before_second][expanded_left] + second_cp

            q_mer_first = f(first_cp[0:q]) #novonastali q-meri
            q_mer_second = f(second_cp[0:q])

            pos_expanded_left = expanded_left + (curr_before_first - i + 1) * 1000 # npr -3005 ako je prosireno preko tri segmenta prije toga i uslo u cetvrti
            if abs(pos_expanded_left - smallest_pivot_first[1]) <= r: # spada u okolinu najmanjeg trenutnog pivota
                pivot_pos = abs(curr_before_first - (i - 1)) * 1000 + smallest_pivot_first[1] + abs(expanded_left)
                if q_mer_first < f(first_cp[pivot_pos:(pivot_pos + q)]): # poredi s q-merom najmanjeg pivota
                    index_smallest = hashed_partitions.index(smallest_pivot_first)
                    length = abs(pos_expanded_left - smallest_pivot_first[1]) + smallest_pivot_first[2]
                    if smallest_pivot_first[1] > 0:
                        length += 1
                    smallest_pivot_first = (f(first_cp[0:length]), pos_expanded_left, length, 0, i)
                    hashed_partitions[index_smallest] = smallest_pivot_first 
                    changed = True

            else: # nije u okolini najmanjeg -> provjerava je li najmanji u svojoj okolini
                cnt = 1
                while cnt <= r:
                    hashed = f(first_cp[cnt:(cnt + q)])
                    if q_mer_first >= hashed:
                        break
                    cnt += 1

                if cnt == 21: # najmanji u svojoj okolini
                    length = abs(pos_expanded_left - smallest_pivot_first[1])
                    if smallest_pivot_first[1] > 0:
                        length += 1
                    smallest_pivot_first = (f(first_cp[0:length]), pos_expanded_left, length, 0, i) # on postaje najmanji lijevi pivot
                    hashed_partitions.append(smallest_pivot_first)
                    changed = True

            # ----- na isti nacin siri drugi segment lijevo -----
            pos_expanded_left = expanded_left + (curr_before_second - j + 1) * 1000
            if abs(pos_expanded_left - smallest_pivot_second[1]) <= r: # spada u okolinu najmanjeg trenutnog pivota
                pivot_pos = abs(curr_before_second - (j - 1)) * 1000 + smallest_pivot_second[1] + abs(expanded_left) # u odnosu na pocetak stringa kada se pocetak gleda kao 0
                if q_mer_second < f(second_cp[pivot_pos:(pivot_pos + q)]): # mijenja trenutnog najmanjeg pivota
                    index_smallest = hashed_partitions.index(smallest_pivot_second)
                    length = abs(pos_expanded_left - smallest_pivot_second[1]) + smallest_pivot_second[2]
                    if smallest_pivot_second[1] > 0:
                        length += 1
                    smallest_pivot_second = (f(second_cp[0:length]), pos_expanded_left, length, 1, j)
                    hashed_partitions[index_smallest] = smallest_pivot_second
                    changed = True

            else: # nije u okolini najmanjeg -> provjerava jel najmanji u svojoj okolini
                cnt = 1
                while cnt <= r:
                    hashed = f(second_cp[cnt:(cnt + q)])
                    if q_mer_second >= hashed:
                        break
                    cnt += 1

                if cnt == 21: # najmanji u svojoj okolini
                    length = abs(pos_expanded_left - smallest_pivot_second[1])
                    if smallest_pivot_second[1] > 0:
                        length += 1
                    smallest_pivot_second = (f(second_cp[0:length]), pos_expanded_left, length, 1, j) # on postaje najmanji lijevi pivot
                    hashed_partitions.append(smallest_pivot_second)
                    changed = True

            expanded_left -= 1
            if expanded_left == -1001:
                expanded_left = -1
                curr_before_first -= 1
                curr_before_second -= 1

        # -------------------------------------------------- sirenje desno ----------------------------------------------------
        elif counter % 2 != 0 and expand_right == True:
            if curr_after_first >= len(chunks1) or curr_after_second >= len(chunks2):
                expand_right = False
                continue

            # sirimo obje duplikacije za jedan element desno
            first_cp = first_cp + chunks1[curr_after_first][expanded_right] 
            second_cp = second_cp + chunks2[curr_after_second][expanded_right] 

            q_mer_first = f(first_cp[(len(first_cp) - q):len(first_cp)]) # novonastali q-meri
            q_mer_second = f(second_cp[(len(second_cp) - q):len(second_cp)])

            pos_expanded_right = expanded_right + (curr_after_first - i) * 1000 - q + 1 # npr 1001 ako je prosireno za 4 slova desno, a q-mer duzine 4 (0-1000 originalni string)
            if abs(pos_expanded_right - largest_pivot_first[1]) <= r: # spada u okolinu najmanjeg trenutnog pivota
                pivot_pos = abs(curr_before_first - (i - 1)) * 1000 + abs(expanded_left) + largest_pivot_first[1] # u odnosu na pocetak stringa kada se pocetak gleda kao 0
                if q_mer_first < f(first_cp[pivot_pos:(pivot_pos + q)]): # poredi s q-merom najveceg pivota
                    index_largest = hashed_partitions.index(largest_pivot_first)
                    index_second_largest = hashed_partitions.index(second_largest_first)

                    length_second_largest = second_largest_first[2] + largest_pivot_first[2] # predzadnjem pivotu se pripaja largest_pivot_first jer je on izgubio status pivota
                    pos_second_largest = abs(curr_before_first - (i - 1)) * 1000 + abs(expanded_left) + second_largest_first[1] # u odnosu na pocetak = 0

                    largest_pivot_first = (q_mer_first, pos_expanded_right, q, 0, i)
                    second_largest_first = (f(first_cp[pos_second_largest:pos_second_largest + length_second_largest]), pos_second_largest, length_second_largest, 0, i)
                    hashed_partitions[index_largest] = largest_pivot_first 
                    hashed_partitions[index_second_largest] = second_largest_first
                    changed = True

            else: # nije u okolini najveceg -> provjerava je li najmanji u svojoj okolini
                cnt = 1
                while cnt <= r:
                    hashed = f(first_cp[(len(first_cp) - q - cnt):(len(first_cp) - cnt)])
                    if q_mer_first >= hashed:
                        break
                    cnt += 1

                if cnt == 21: # najmanji u svojoj okolini
                    largest_pivot_first = (q_mer_first, pos_expanded_right, q, 0, i) # on postaje najmanji lijevi pivot
                    hashed_partitions.append(largest_pivot_first)
                    changed = True
            
            # ----- na isti nacin siri drugi segment desno -----
            pos_expanded_right = expanded_right + (curr_after_second - j) * 1000 - q + 1 # npr 1001 ako je prosireno za 4 slova desno, a q-mer duzine 4 (0-1000 originalni string)
            if abs(pos_expanded_right - largest_pivot_second[1]) <= r: # spada u okolinu najmanjeg trenutnog pivota
                pivot_pos = abs(curr_before_second - (j - 1)) * 1000 + abs(expanded_left) + largest_pivot_second[1] # u odnosu na pocetak stringa kada se pocetak gleda kao 0
                if q_mer_second < f(second_cp[pivot_pos:(pivot_pos + q)]): # poredi s q-merom najveceg pivota
                    index_largest = hashed_partitions.index(largest_pivot_second)
                    index_second_largest = hashed_partitions.index(second_largest_second)

                    length_second_largest = second_largest_second[2] + largest_pivot_second[2] # predzadnjem pivotu se pripaja largest_pivot_first jer je on izgubio status pivota
                    pos_second_largest = abs(curr_before_second - (j - 1)) * 1000 + abs(expanded_left) + second_largest_second[1] # u odnosu na pocetak = 0

                    largest_pivot_second = (q_mer_second, pos_expanded_right, q, 1, j)
                    second_largest_second = (f(second_cp[pos_second_largest:pos_second_largest + length_second_largest]), pos_second_largest, length_second_largest, 1, j)
                    hashed_partitions[index_largest] = largest_pivot_second
                    hashed_partitions[index_second_largest] = second_largest_second
                    changed = True

            else: # nije u okolini najveceg -> provjerava je li najmanji u svojoj okolini
                cnt = 1
                while cnt <= r:
                    hashed = f(second_cp[(len(second_cp) - q - cnt):(len(second_cp) - cnt)])
                    if q_mer_second >= hashed:
                        break
                    cnt += 1

                if cnt == 21: # najmanji u svojoj okolini
                    largest_pivot_second = (q_mer_second, pos_expanded_right, q, 1, j) # on postaje najmanji lijevi pivot
                    hashed_partitions.append(largest_pivot_second)
                    changed = True

            expanded_right += 1
            if expanded_right == 1000:
                expanded_right = 0
                curr_after_first += 1
                curr_after_second += 1

        # update thresholda jer se povecala duzina stringova 
        K = percentage_error / 100 * len(first_cp)
        T = 20 + K / 8.
        threshold = T / 20

        # provjeravamo poklapanja segmenata 
        if changed:
            num_same_segments = 0
            hashed_partitions.sort(key=lambda x: x[0])
            grouped = defaultdict(list)
            for key, group in groupby(hashed_partitions, key=lambda x: x[0]): # grupisemo uzastopne elemente koji imaju istu hashiranu vrijednost
                grouped[key].extend(group)

            for G in grouped.values(): # za svaku grupu
                for (hi, posi, leni, i, x), (hj, posj, lenj, j, y) in combinations(G, 2): # za svaki par unutar jedne grupe; ne gleda se raspored -> svaki par uzima samo jednom
                    if i != j:
                        if posi > 0:
                            first_arg = (curr_after_first - 1) * 1000 + expanded_right - posi
                        else:
                            first_arg = abs(posi) + 1 + (curr_after_first - 1) * 1000 + expanded_right

                        if posj > 0:
                            second_arg = (curr_after_second - 1) * 1000 + expanded_right - posj
                        else:
                            second_arg = abs(posj) + 1 + (curr_after_second - 1) * 1000 + expanded_right

                        if abs(posi - posj) + abs(first_arg - second_arg) <= K:
                            num_same_segments += 1
            changed = False

        counter += 1
        if expand_left == False or expand_right == False:
            return None # dosao do granice jedne od prosirenih sekvenci -> ne moze se u manje od 10% poklapati s dokazanim

    # zanemarujemo zadnju iteraciju petlje i vracamo do tad prosireno iz funkcije
    return (first_before_last_it, second_before_last_it)
    

def floor(x):
    return int(x) if x >= 0 else int(x) - 1 if x != int(x) else int(x)

def find_anchors(s, T, Pi):
    A = [0]
    q = 8
    #r = floor((len(s) - q + 1 - T) / (2 * T + 2))
    r = 10
    
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

def Minjoin(S, K, T, i, j):
    O = set()
    C = []
    A = []

    P = partition_string(S[0], T, Pi)
    for pos, length in P:
        A.append((f(S[0][pos:pos + length]), pos, length, 0, i))

    P = partition_string(S[1], T, Pi)
    for pos, length in P:
        A.append((f(S[1][pos:pos + length]), pos, length, 1, j))

    A.sort(key=lambda x: x[0])
    grouped = defaultdict(list)
    for key, group in groupby(A, key=lambda x: x[0]): # grupisemo uzastopne elemente koji imaju istu hashiranu vrijednost
        grouped[key].extend(group)

    for G in grouped.values(): # za svaku grupu
        for (hi, posi, leni, i, x), (hj, posj, lenj, j, y) in combinations(G, 2): # za svaki par unutar jedne grupe; ne gleda se raspored -> svaki par uzima samo jednom
            if i != j and abs(len(S[0]) - len(S[1])) <= K:
                if abs(posi - posj) + abs((len(S[0]) - posi) - (len(S[1]) - posj)) <= K:
                    C.append((i, j))

    count = Counter(C)
    threshold = T / 20

    for (i, j) in C:
        if count[(i, j)] >= int(threshold):
            O.add((i, j))

    return (O, A, count[(i, j)])

def split_into_windows(sequence, window_size=1000, step_size=1000):
    windows = []
    if len(sequence) < window_size:
        windows.append(sequence)
        return windows
    
    for start in range(0, len(sequence) - window_size, step_size):
        end = start + window_size
        windows.append(sequence[start:end])

    return windows

def compare_windows(windows1, windows2, K, T):
    results = []
    for i, win1 in enumerate(windows1):
        for j, win2 in enumerate(windows2):
            comparison_results = (Minjoin([win1, win2], K, T, i, j))[0]
            if comparison_results:
                results.append(comparison_results)
    return results

def process_pair(pair, K, T):
        # poredi dvije sekvence
        (chrom1, start1, end1, sequence1) = pair[0]
        (chrom2, start2, end2, sequence2) = pair[1]

        if chrom1 == chrom2:
            windows1 = split_into_windows(sequence1)
            windows2 = split_into_windows(sequence2)
            comparisons = compare_windows(windows1, windows2, K, T)
            return comparisons
        
        return None

def find_windows(pair, org_pair, K, T):
    windows1 = split_into_windows(pair[0][-1], 1000, 1000)
    windows2 = split_into_windows(pair[1][-1], 1000, 1000)

    if not windows1 or not windows2:
        return None

    for i, win1 in enumerate(windows1):
        for j, win2 in enumerate(windows2):
            if i * 1000 + pair[0][1] < org_pair[0][1] or j * 1000 + pair[1][1] < org_pair[1][1]: 
                continue
            comparison_results = Minjoin([win1, win2], K, T, i, j)
            if comparison_results[2] > 1: # Prepoznao barem jedan prozor koji zadovoljava MinJoin
                # i je pozicija poklapajuceg prozora u skupini prozora za extended_pairs[0] -> njegova pozicija u cijelom kromozomu extended_pairs[0][1] + i*1000
                return (windows1, windows2, i, j, comparison_results[1], comparison_results[2])
    
    return None

if __name__ == "__main__":
    start_time = time.time()
    seq_file = 'filtered_sequencespom.fa'
    seq_extended = 'filtered_sequences_extended.fa'
    sequences = []
    sequences_extended = []
    pattern = re.compile(r"^>(chr\d+):(\d+)-(\d+)$")

    # Otvaramo dokument s duplikacijama
    with open(seq_file, 'r') as f:
        for data in f:
            seq = next(f, None)
            data = data.replace('\n', '')
            seq = seq.replace('\n', '')

            match = pattern.match(data)
            if match:
                chrom, start, end = match.groups()
                sequences.append((chrom, int(start), int(end), seq))
    # Prosjecna duzina duplikacija 13256.6

    # Otvaramo dokument s prosirenim duplikacijama za poredjenje
    with open(seq_extended, 'r') as f:
        for data in f:
            seq = next(f, None)
            data = data.replace('\n', '')
            seq = seq.replace('\n', '')

            match = pattern.match(data)
            if match:
                chrom, start, end = match.groups()
                sequences_extended.append((chrom, int(start), int(end), seq))
    
    K = 150.
    T = 20 + K / 8.
    threshold = T / 20

    pairs = [(sequences[i], sequences[i + 1]) for i in range(0, len(sequences) - 1, 10)]
    pairs_extended = [(sequences_extended[i], sequences_extended[i + 1]) for i in range(0, len(sequences_extended) - 1, 10)]

    cnt_pass_minjoin = 0
    cnt_found_window = 0
    cnt_false = 0
    cnt_true = 0
    cnt_not_enough = 0
    cnt_too_much = 0
    cnt_left_first_not_enough = 0
    cnt_left_first_too_much = 0
    cnt_right_first_not_enough = 0
    cnt_right_first_too_much = 0

    cnt_left_second_not_enough = 0
    cnt_left_second_too_much = 0
    cnt_right_second_not_enough = 0
    cnt_right_second_too_much = 0

    for i, pair in enumerate(pairs):
        result = process_pair(pair, K, T)
        if len(result) != 0: # MinJoin prepoznao duplikaciju 
            cnt_pass_minjoin += 1
            # Ako je uslov da moraju granice odstupati najvise 15% ukupne duzine sekvenci -> ako moze da siri a dostigao granicu stringa -> ne zadovoljava
            win = find_windows(pairs_extended[i], pair, K, T)
            if win is not None: # Pronadjen barem jedan prozor duzine 1000 koji zadovoljava MinJoin
                cnt_found_window += 1
                windows1 = win[0]
                windows2 = win[1]
                num_i = win[2]
                j = win[3]
                A = win[4]
                len_C = win[5]
                borders = expand(windows1, windows2, A, num_i, j, len_C, K, threshold)

                if borders is None: # Pokusao siriti van granica 
                    cnt_too_much += 1
                    cnt_false += 1

                else:
                    if borders[0][0] + pairs_extended[i][0][1] - pair[0][1] > 0:
                        cnt_left_first_not_enough += 1
                        cnt_not_enough += 1
                    else:
                        cnt_left_first_too_much += 1
                        #cnt_too_much += 1
                    if pair[0][2] - borders[0][1] - pairs_extended[i][0][1] > 0:
                        cnt_right_first_not_enough += 1
                        cnt_not_enough += 1
                    else:
                        cnt_right_first_too_much += 1
                        #cnt_too_much += 1
                    if borders[1][0] + pairs_extended[i][1][1] - pair[1][1] > 0:
                        cnt_left_second_not_enough += 1
                        cnt_not_enough += 1
                    else:
                        cnt_left_second_too_much += 1
                        #cnt_too_much += 1
                    if pair[1][2] - borders[1][1] - pairs_extended[i][1][1] > 0:
                        cnt_right_second_not_enough += 1
                        cnt_not_enough += 1
                    else:
                        cnt_right_second_too_much += 1
                        #cnt_too_much += 1

                    if abs(borders[0][0] + pairs_extended[i][0][1] - pair[0][1]) <= 0.3 * len(pair[0][-1]) and abs(borders[0][1] + pairs_extended[i][0][1] - pair[0][2]) <= 0.3 * len(pair[0][-1]) and abs(borders[1][0] + pairs_extended[i][1][1] - pair[1][1]) <= 0.3 * len(pair[1][-1]) and abs(borders[1][1] + pairs_extended[i][1][1] - pair[1][2]) <= 0.3 * len(pair[1][-1]) or abs(borders[0][0] + pairs_extended[i][0][1] - pair[1][1]) <= 0.3 * len(pair[1][-1]) and abs(borders[0][1] + pairs_extended[i][0][1] - pair[1][2]) <= 0.3 * len(pair[1][-1]) and abs(borders[1][0] + pairs_extended[i][1][1] - pair[0][1]) <= 0.3 * len(pair[0][-1]) and abs(borders[1][1] + pairs_extended[i][1][1] - pair[0][2]) <= 0.3 * len(pair[0][-1]):
                        cnt_true += 1
                    else:
                        cnt_false += 1
    

    print(len(pairs), cnt_pass_minjoin, cnt_found_window, cnt_false, cnt_true)  
    print(cnt_left_first_not_enough, cnt_left_first_too_much, cnt_right_first_not_enough, cnt_right_first_too_much, cnt_left_second_not_enough, cnt_left_second_too_much, cnt_right_second_not_enough, cnt_right_second_too_much)
    print(cnt_too_much)

    end_time = time.time()
    print(end_time - start_time)


