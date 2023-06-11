from itertools import combinations
import random
from datetime import date
from tkinter import *


def min_seq_length(sequences):  #
    return min([len(x) for x in sequences])

######################### Computational Function
def random_singlept_crossover(sequences):      # cross over between 2 sequences at random point
    min_length = min_seq_length(sequences)
    random_point = random.randint(1,min_length)
    new_string1 = sequences[0][:random_point] + sequences[1][:-random_point]
    new_string2 = sequences[1][:random_point] + sequences[0][:-random_point]
    return [new_string1, new_string2]

def defined_singlept_crossover(sequences,position_dict,lengths):   # cross over between 2 sequences at random point
    min_length = min_seq_length(sequences);
    positions1 = position_dict[sequences[0]]; positions2 = position_dict[sequence[1]]
    for pos1, pos2,length in zip(positions1, positions2,lengths):
        random_point = random.randint(1, min_length)
        if ((random_point > pos1 and (random_point - pos1) < length) and
            (random_point > pos2 and (random_point - pos2) < length)):
            random_point = random_point + length;
            if random_point > min_length:
                return sequences
            new_string1 = sequences[0][:random_point] + sequences[1][:-random_point]
            new_string2 = sequences[1][:random_point] + sequences[0][:-random_point]
            return [new_string1, new_string2]
        else:
            return sequences

def crossover(sequences, cross_percent=1, iteration=1, max_cap=1e4,defined=0):
    # cross_percent: number of population crossover per iteration #max_cap: maximum number of sequence allowed to generate
    if defined == 1:
        motifs = list(ast.literal_eval(motif_entry.get())); mot_lengths = [len(i) for i in motifs];
        if motifs == '':
            tkinter.messagebox.showinfo('Warm Warning', 'No, sorry. You can\'t.')
            return None
        else:
            motif_pos_dict = {}
            for seq in sequences:
                motif_pos_dict[seq] = findmotif(seq,motifs) ## list of position of the inputted motifs
    new_sequences = []; generation_pool = sequences; #initialize the first pool to cross over
    for i in range(iteration):
        pairpool = list(combinations(generation_pool,2))
        pairpool_selected = [x for x in pairpool if random.random()< cross_percent] # choose a portion from the combination pool
        generation_pool=[]
        for pair in pairpool_selected:
            if len(generation_pool) >= max_cap: break ### condition to stop adding sequences
            if defined == 1:
                newpair = defined_singlept_crossover(pair,motif_pos_dict, mot_lengths)  ;### pair : ('acgt','acgt); position:[[3,2],[1,4]]
            else:
                newpair = random_singlept_crossover(pair)  ; #print(newpair) #### newpair: ['ACGGT','ACGT]
            generation_pool.append(newpair)
        generation_pool= [seq for seqpair in generation_pool for seq in seqpair]
        new_sequences.append(generation_pool)
        temp_seqs = [seq for sublist in new_sequences for seq in sublist]
        if len(temp_seqs)>max_cap:
            return temp_seqs
    return temp_seqs

def pointmutation(sequences,mutation_rate, skewrate): #mutation rate is number in a million # skew (A,C,G,T) rate
    if skewrate =='': skewrate = (1,1,1,1)
    tot = sum(skewrate); dicerate = list(np.cumsum(skewrate)/tot-0.25); print(dicerate)
    bases = ['A','C','G','T']
    for seq_num in range(len(sequences)):
        for base in range(len(sequences[seq_num])):
            dice = random.random();
            if dice <= mutation_rate/1e6:
                dice2 = random.random()
                for i in range(len(dicerate)):
                    if dice2 > dicerate[i]:
                        seq_base = list(sequences[seq_num])
                        seq_base[base] = bases[i]
                        sequences[seq_num]=''.join(seq_base)
    return sequences

### Input Parameters

seqs = ['CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA','GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG',\
     'TAGTACCGAGACCGAAAGAAGTATACAGGCGT','TAGATCAAGTTTCAGGTGCACGTCGGTGAACC',\
     'AATCCACCAGCTCCACGTGCAATGTTGGCCTA']

cross_percent=1;
iteration=1;
max_cap = 1e4;
defined=0;     ### Defined motif

## Output Parameters

output = crossover(seqs,cross_percent, iteration, max_cap,defined)
file = open('ISMSequences_'+str(date.today())+'.txt','w')
for item in output:
    file.write('%s\n' %item)