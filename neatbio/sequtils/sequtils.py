aa3_to_1_dict = {'Ala': 'A',
                     'Cys': 'C',
                     'Asp': 'D',
                     'Glu': 'E',
                     'Phe': 'F',
                     'Gly': 'G',
                     'His': 'H',
                     'Ile': 'I',
                     'Lys': 'K',
                     'Leu': 'L',
                     'Met': 'M',
                     'Asn': 'N',
                     'Pro': 'P',
                     'Gln': 'Q',
                     'Arg': 'R',
                     'Ser': 'S',
                     'Thr': 'T',
                     'Val': 'V',
                     'Trp': 'W',
                     'Tyr': 'Y'}
    

def gc_content(seq):
    result = float(seq.count('G') + seq.count('C')) / len(seq) * 100
    return round(result, 4)

def at_content(seq):
    result = float(seq.count('A') + seq.count('T')) / len(seq) * 100
    return round(result, 4)

def __get_key(val, my_dict):
    for key, value in my_dict.items():
        if val == value:
            return key

def __get_value(val, my_dict):
    for key, value in my_dict.items():
        if val == key:
            return value
        
def convert_1to3(seq):
    term_list = [__get_key(n, aa3_to_1_dict) for n in seq]

    return ''.join(term_list)
    
def __kmers(self, k=3):
    triplet_list = [self.seq[i:i+k] for i in range(0, len(self.seq), k)]
    return triplet_list

def convert_3to1(seq):
    term_list = [__get_value(n, aa3_to_1_dict) for n in __kmers(seq)]
    return ''.join(term_list)

def hamming_distance(lhs, rhs):
    return len([(x,y) for x,y in zip(lhs,rhs) if x != y])

def occurence(main_seq, sub_seq):
    start = 0
    indices = []
    while True:
        start = main_seq.find(sub_seq, start)
        if start > 0:
            indices.append(start)
        else:
            break
        start += 1
    return indices