import sys

CodonTable = {
    # 'M' - START, '*' - STOP
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "TGT": "C",
    "TGC": "C",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "TTT": "F",
    "TTC": "F",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
    "CAT": "H",
    "CAC": "H",
    "ATA": "I",
    "ATT": "I",
    "ATC": "I",
    "AAA": "K",
    "AAG": "K",
    "TTA": "L",
    "TTG": "L",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "ATG": "M",
    "AAT": "N",
    "AAC": "N",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AGA": "R",
    "AGG": "R",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "AGT": "S",
    "AGC": "S",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "TGG": "W",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGA": "*",
}


class Sequence(object):
    """Create a Valid Sequence for DNA, RNA

    example: seq = Sequence('ATGC')

    """

    def __init__(self, seq=None):
        super(Sequence, self).__init__()
        self.seq = seq
        if not isinstance(self.__validate_seq(seq), str):
            raise TypeError(
                "The sequence data given to a Sequence() should be a string (not another Sequence object)"
                "nor a Non Valid Nucleotide [A, T, G, C, U]"
            )

    def __repr__(self):
        return f"Sequence(seq='{self.seq}')"

    def __str__(self):
        return self.seq

    def __validate_seq(self, seq):
        base_nucleotide = ["A", "T", "G", "C", "U"]
        real_seq = seq.upper()
        for base in real_seq:
            if base not in base_nucleotide:
                return False

        return real_seq

    def __len__(self):
        return len(self.seq)

    def __contains__(self, sub_char):
        return sub_char in str(self)

    def __getitem__(self, index):
        if isinstance(index, int):
            return self.seq[index]

    # Basic Fxn - count, find, index

    def count(self, subseq, start=0, end=sys.maxsize):
        """Return the count of the number of nucleotide in a sequence."""
        return str(self).count(subseq, start, end)

    def find(self, subseq, start=0, end=sys.maxsize):
        """Find the position of a nucleotide in a sequence."""
        return str(self).find(subseq, start, end)

    def rfind(self, subseq, start=0, end=sys.maxsize):
        """Find the position of a nucleotide in a sequence from right to left."""
        return str(self).rfind(subseq, start, end)

    def index(self, subseq, start=0, end=sys.maxsize):
        """Find the index of a nucleotide in a sequence."""
        return str(self).index(subseq, start, end)

    def rindex(self, subseq, start=0, end=sys.maxsize):
        """Find the position of a nucleotide in a sequence from right to left."""
        return str(self).rindex(subseq, start, end)

    ### Main Fxn
    def get_frequency(self):
        """Get the Frequency of a nucleotide in a sequence."""
        base_dict = {"A": 0, "T": 0, "G": 0, "C": 0}

        for n in self.seq:
            if self.__validate_seq(n) != False:
                base_dict[n] += 1
            else:
                return "NucleotideError: {} is not a nucleotide [A, T, G, C]".format(n)
        return base_dict

    @property
    def gc(self):
        """Get GC content of a sequence"""
        result = float(self.seq.count("G") + self.seq.count("C")) / len(self.seq) * 100
        return round(result, 4)

    @property
    def at(self):
        """Get AT content of a sequence"""
        result = float(self.seq.count("A") + self.seq.count("T")) / len(self.seq) * 100
        return round(result, 4)

    def complement(self):
        """Return complement of sequence"""
        base_pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}
        result = [base_pairs[n] for n in self.seq if n in base_pairs.keys()]

        return "".join(result)

    def reverse_complement(self):
        """Return reverse complement of sequence"""
        base_pairs = {"A": "T", "T": "A", "G": "C", "C": "G"}
        result = [base_pairs[n] for n in self.seq if n in base_pairs.keys()]

        return "".join(result)[::-1]

    def transcribe(self):
        """Transcribe sequence into mRNA"""
        mrna = self.seq.replace("T", "U")
        return mrna

    def translate(self, start_pos=0):
        """Translate DNA sequence into Amino Acid"""
        amino_acid_list = [
            CodonTable[self.seq[pos : pos + 3]]
            for pos in range(start_pos, len(self.seq) - 2, 3)
        ]

        return "".join(amino_acid_list)
