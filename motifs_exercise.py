lecture_dna = [
    "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
    "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
    "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
    "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
    "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
    "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
]
import random
from Bio import motifs
from Bio.Seq import Seq

def count_matrix(motifs):
    A_list=[]
    C_list=[]
    G_list=[]
    T_list=[]
    for i in range(len(motifs[0])):
        A_count = 0
        C_count = 0
        G_count = 0
        T_count = 0
        for motif in motifs:
            if motif[i] == "A":
                A_count=A_count+1
            elif motif[i] == "C":
                C_count = C_count + 1
            elif motif[i] == "G":
                G_count = G_count + 1
            elif motif[i] == "T":
                T_count = T_count + 1
            else:
                print("error: unexpected character")
                return None
        A_list.append(A_count)
        C_list.append(C_count)
        G_list.append(G_count)
        T_list.append(T_count)
    keys=["A","C","G","T"]
    ans=[A_list,C_list,G_list,T_list]
    return dict(zip(keys,ans))

def score(motifs):
    motifs_count=count_matrix(motifs)
    c_score=0
    for i in range(len(motifs[0])):
        high_score = 0
        for key,val in motifs_count.items():
            if val[i]> high_score:
                high_score=val[i]
            else:
                high_score=high_score
        c_score=c_score+high_score
    return c_score

def consensus(motifs):
    motifs_count=count_matrix(motifs)
    consens = ""
    for i in range(len(motifs[0])):
        high_score = 0
        char=""
        for key,val in motifs_count.items():
            if val[i]> high_score:
                high_score=val[i]
                char=key
            else:
                high_score=high_score
        consens=consens+char
    return consens
def hamming_distance(a,b):
    dist=0
    for i in range(len(a)):
        if a[i]!=b[i]:
            dist=dist+1
    return dist

def total_distance(pattern,sequences):
    total=0
    for seq in sequences:
        min=len(pattern)
        dif= len(seq)-len(pattern)
        for i in range(dif+1):
            win_len = len(seq) - dif + i
            window=seq[i:win_len]
            if hamming_distance(pattern,window) < min:
                min = hamming_distance(pattern,window)
        total=total+min
    return total
def ppm_count(motifs):
    counts=count_matrix(motifs)
    motif_num = len(motifs)
    for key,val in counts.items():
        ppm_list = []
        for i in range(len(val)):
            propa=((val[i]+1)/(4+motif_num))
            ppm_list.append(propa)
        counts[key]=ppm_list
    return counts

class MotifProfile :
    def __init__(self, motifs, pseudocount=1):
        self.motifs=motifs
        self.pseudocount=pseudocount
        self.l=len(motifs[0])
        self.ppm=ppm_count(motifs)

    def lmer_probability(self,lmer):
        propa=1
        for i in range(len(lmer)):
            for key,val in self.ppm.items():
                if lmer[i] == key :
                    propa=propa*val[i]
        return propa

    def most_probable_lmer(self, sequence):
        if len(sequence) < self.l :
            print("error canot compute most propable lmer , lmer is shorter than motif")
            return None
        dif = len(sequence)-len(self.motifs[0])
        best_propa=0
        best_lmer=""
        for j in range(dif+1):
            win_len = len(self.motifs[0])
            window = sequence[j:win_len+j]
            if self.lmer_probability(window)>best_propa:
                best_propa = self.lmer_probability(window)
                best_lmer=window
        return best_lmer
    def consensus(self):
        cons = consensus(self.motifs)
        return cons


rng = random.Random(1)                 # a random number generator with seed 1
i = rng.randint(0, 50)                 # random integer, 0 <= i <= 50 (both ends included!)
lmer = rng.choice(["ACG", "CGT", "GTA"])   # one random item of a list

from Bio import SeqIO

sequences = [str(record.seq) for record in SeqIO.parse("planted_motif.fasta", "fasta")]
class Motfinder():
    def __init__(self,sequences,l,seed=1):
        self.sequences=sequences
        self.l =l
        self.seed=seed
