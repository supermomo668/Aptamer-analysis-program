import random, csv, datetime,os, ast, numpy as np, time, math as m, ast, \
    matplotlib as mpl

import tkinter as tk
from tkinter import *
from tkinter.filedialog import askopenfilename
from tkinter import messagebox, font
from itertools import combinations
from scipy.cluster import hierarchy
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
import seaborn
import matplotlib.pyplot as plt
from matplotlib import transforms
import matplotlib.patheffects
from matplotlib.font_manager import FontProperties
plt.style.use('seaborn-ticks')

rootwindow = Tk()
class App:
    def __init__(self, master):
        master.title("The Aptamer App: Holy Moti"); master.geometry('400x300')
        master.iconbitmap("MM_icon.ico")
        self.frame = Frame(master, width=350, height=220).grid(rowspan=7, columnspan=5, sticky=N)
        Label(master, text='Welcome to the duly decorated application\n\nWhat are you doing here?').grid(row=1,column=3)
        button1 = Button(self.frame, text="ISM", command=self.ISMwindow, fg='red')
        button1.grid(row=2, column=3)
        button2 = Button(self.frame, text = 'Motif Search', command = self.MotifSearchWindow,fg='blue')
        button2.grid(row = 3, column = 3)
        self.readylabel = Label(self.frame, text="")
        self.readylabel.grid(row=7, column=3)

    def ISMwindow(self):
        ISM()
        self.readylabel.config(text="Launching... to Aptamer In Silico Maturation")
    def MotifSearchWindow(self):
        MotifSearch()
        self.readylabel.config(text="Launching... to Aptamer Motif Search")
    def AptamerAI(self):
        AptamerAI()
        self.readylabel.config(text = "Still a construction zone... wear a helmet.")

class ISM:
    def __init__(self):
        master_ISM = Tk(); master_ISM.title("ISM Generation" ); master_ISM.geometry('590x720')
        frame = Frame(master_ISM, highlightbackground="grey", highlightcolor="blue", highlightthickness=1,
                      width = 580, height = 710).grid(rowspan = 18, columnspan = 8, pady = 4, padx = 4,sticky= NW)
        Label(master_ISM, text = 'Your input file: \nNote: Line-delimited sequences').grid(row = 1, column = 3,sticky= NW)
        button1 = Button(master_ISM, text = 'Upload file here (.txt/.csv)',command = self.file_upload); button1.grid(row=1, column =4,sticky= NW);
        Label(master_ISM, text = 'Desired output file name (w/o extention):  \nDefault: Name_Date.txt'
                                     ).grid(row = 2, column = 3)

        self.entry_fileout = Entry(master_ISM, width = 30); self.entry_fileout.bind("<Return>")
        self.entry_fileout.grid(row = 2, column = 4)
        Label(master_ISM, text = 'In Silico Maturation', font=("Courier",15,"bold italic")).grid(row = 3, column=3)
        button_fold = Button(master_ISM, text = 'Generate ISM Sequences',fg='Violet', command=self.ISMGenerate)
        button_fold.grid(row = 11, column = 4)


        ## Parameters
        Label(master_ISM, text = 'Crossover Parameters: ').grid(row = 4, column = 2)
        Label(master_ISM, text='Cross Percentage (0-1)\nDefault:1').grid(row=5, column=3)
        self.entry_CrossPercent = Entry(master_ISM); self.entry_CrossPercent.bind("<Return>")
        self.entry_CrossPercent.grid(row=5, column=4); self.entry_CrossPercent.insert(END,1)
        Label(master_ISM, text='No. of Generation\nDefault:1').grid(row=6, column=3)
        self.entry_NumGeneration = Entry(master_ISM); self.entry_NumGeneration.bind("<Return>")
        self.entry_NumGeneration.grid(row=6, column=4); self.entry_NumGeneration.insert(END,1)
        Label(master_ISM, text='Maximum No.Sequences\nDefault:1000').grid(row=7, column=3)
        self.entry_MaxCap = Entry(master_ISM,text = '1000'); self.entry_MaxCap.bind("<Return>")
        self.entry_MaxCap.grid(row=7, column=4); self.entry_MaxCap.insert(END,1000)

        Label(master_ISM, text='Mutation Parameters: ').grid(row=8, column=2)
        Label(master_ISM, text='Mutation Rate(/1000)\nDefault: 20').grid(row=9, column=3)
        self.entry_MutRate = Entry(master_ISM); self.entry_MutRate.bind("<Return>")
        self.entry_MutRate.grid(row=9, column=4); self.entry_MutRate.insert(END,20)
        Label(master_ISM, text='Mutation Skew (A,C,G,T)\nDefault:(1,1,1,1)').grid(row=10, column=3)
        self.entry_MutSkew = Entry(master_ISM); self.entry_MutSkew.bind("<Return>")
        self.entry_MutSkew.grid(row=10, column=4); self.entry_MutSkew.insert(END,'[1,1,1,1]')

        label = Label(master_ISM, text='Status: ', fg='Red'); label.grid(row=11, column=2)
        self.label_status = Label(master_ISM, text="nothing much"); self.label_status.grid(row=11, column=3)

        self.textbox1 = Text(master_ISM, height=13, width=64); self.textbox1.insert(END, 'File Info: ')
        self.textbox1.grid(row = 12,column = 0,columnspan=8, rowspan=5, padx = 5)
        ##

    def ISMGenerate(self):
        if self.type == '.csv':
            self.textbox1.insert(END, 'Not ready for this yet')
        elif self.type == '.txt':
            name_fileout = self.entry_fileout.get()
            if name_fileout == '':
                name_fileout = 'MyResults_'+str(datetime.date.today())
            output_directory = os.path.dirname(self.filedir)+ '/'+name_fileout+'.txt'

            file_out = open(output_directory, 'w')

            crossover_output_sequences = \
                self.crossover(Tolist(self.filedir), self.SwiftProcessEntry(self.entry_CrossPercent),
                self.SwiftProcessEntry(self.entry_NumGeneration),
                               self.SwiftProcessEntry(self.entry_MaxCap))
            mutation_output_sequences = \
                self.pointmutation(crossover_output_sequences, self.SwiftProcessEntry(self.entry_MutRate),
                                   list(ast.literal_eval(self.entry_MutSkew.get())))
            output_sequences = mutation_output_sequences; numcount = len(output_sequences)

            [file_out.write(i+'\n') for i in output_sequences]
            self.textbox1.insert(END, '\n.\n.\nHere you go! \n' + 'Your output file directory is: '
                                 + output_directory + '\n with '+str(numcount)+ ' sequences.')

    def file_upload(self):
        self.filedir = askopenfilename(); self.label_status.config(text = 'Importing...')
        if self.filedir.endswith('.csv'):
            self.label_status.config(text = "will be converted"); self.type = '.csv'
        elif self.filedir.endswith('.txt'):
            self.label_status.config(text = "all good"); self.type = '.txt'
        else:
            messagebox.showinfo("Wait", "Something is wrong")
            return
        self.seqs = Tolist(self.filedir)
        self.textbox1.insert(END, '\nYour input file is located at:'+self.filedir+'\nNumber of sequences: '+
                             str(len(self.seqs))+ '\nFirst sequence: '+self.seqs[0]+'\nLast Sequence: '+self.seqs[-1])

    def min_seq_length(self,sequences):  #
        return min([len(x) for x in sequences])

    ######################### Computational Function
    def random_singlept_crossover(self,sequences):  # cross over between 2 sequences at random point
        min_length = self.min_seq_length(sequences)
        random_point = random.randint(1, min_length)
        new_string1 = sequences[0][:random_point] + sequences[1][:-random_point]
        new_string2 = sequences[1][:random_point] + sequences[0][:-random_point]
        return [new_string1, new_string2]

    def defined_singlept_crossover(self,sequences, position_dict, lengths):  # cross over between 2 sequences at random point
        min_length = self.min_seq_length(sequences)
        positions1 = position_dict[sequences[0]]
        positions2 = position_dict[sequences[1]]
        for pos1, pos2, length in zip(positions1, positions2, lengths):
            random_point = random.randint(1, min_length)
            if ((random_point > pos1 and (random_point - pos1) < length) and
                    (random_point > pos2 and (random_point - pos2) < length)):
                random_point = random_point + length
                if random_point > min_length:
                    return sequences
                new_string1 = sequences[0][:random_point] + sequences[1][:-random_point]
                new_string2 = sequences[1][:random_point] + sequences[0][:-random_point]
                return [new_string1, new_string2]
            else:
                return sequences

    # Cross over a list of sequence at random position. With an option to cross-over at a certain motif position
    def crossover(self,sequences, cross_percent=1, iteration=1, max_cap=1e4, defined=0, motif_pos_dict=None, mot_lengths=None):
        # cross_percent: number of population crossover per iteration #max_cap: maximum number of sequence allowed to generate
        iteration=int(iteration); new_sequences = []
        generation_pool = sequences;  # initialize the first pool to cross over
        for i in range(iteration):
            pairpool = list(combinations(generation_pool, 2))
            pairpool_selected = [x for x in pairpool if random.random() < cross_percent]
            # choose a portion from the combination pool
            generation_pool = []
            for pair in pairpool_selected:
                if len(generation_pool) >= max_cap: break
                if defined == 1:
                    newpair = self.defined_singlept_crossover(pair, motif_pos_dict, mot_lengths)
                    ### pair : ('acgt','acgt); position:[[3,2],[1,4]]
                else:
                    newpair = self.random_singlept_crossover(pair);  # print(newpair) #### newpair: ['ACGGT','ACGT']
                generation_pool.append(newpair)
            generation_pool = [seq for seqpair in generation_pool for seq in seqpair]
            new_sequences.append(generation_pool)
            temp_seqs = [seq for sublist in new_sequences for seq in sublist]
            if len(temp_seqs) >= max_cap:
                return temp_seqs[:max_cap]
        return temp_seqs

    def pointmutation(self, sequences, mutation_rate, skewrate):  # mutation rate is number in a million # skew (A,C,G,T) rate
        if skewrate == '': skewrate = (1, 1, 1, 1)
        tot = sum(skewrate)
        dicerate = list(np.cumsum(skewrate) / tot - 0.25)
        bases = ['A', 'C', 'G', 'T']
        for seq_num in range(len(sequences)):
            for base in range(len(sequences[seq_num])):
                dice = random.random()
                if dice <= mutation_rate / 1e3:
                    dice2 = random.random()
                    for i in range(len(dicerate)):
                        if dice2 > dicerate[i]:
                            seq_base = list(sequences[seq_num])
                            seq_base[base] = bases[i]
                            sequences[seq_num] = ''.join(seq_base)
        return sequences
    def SwiftProcessEntry(self,entry):
        if entry.get() == '':
            return
        else:
            return float(entry.get())

class MotifSearch:
    def __init__(self):
        root = Tk(); root.geometry("950x700"); root.title('Holy Moti')
        root.iconbitmap("MM_icon.ico")

        label1 = Label(root, text="Sequence: (Format: Python 'List')"); label1.grid(row=0)
        self.Dna_entry = Entry(root); self.Dna_entry.bind("<Return>"); self.Dna_entry.grid(row=0, column=1)

        Label(root, text="Motif Length Up to: (Format: Integer < Kmax)").grid(row=1)
        self.k_entry = Entry(root); self.k_entry.bind("<Return>"); self.k_entry.grid(row=1, column=1)

        Label(root, text="No. of Reiteration (Default: 100): (Format: Integer)").grid(row=2)
        self.N_entry = Entry(root); self.N_entry.bind("<Return>"); self.N_entry.grid(row=2, column=1)

        others_label = Label(root, text='Other Functions'); others_label.grid(row=2, column=3);  ##  Make underlined text
        f = font.Font(others_label, others_label.cget("font")); f.configure(underline=True, slant='italic'); others_label.configure(font=f);
        Button(root,text='GQ').grid(row = 3, column = 3)

        upload_button_fastq = Button(root, text='Upload Fastq(.fastq)', command= self.file_upload)
        upload_button_fastq.grid(row=0, column=2)
        upload_button_csv = Button(root, text='Upload CSV(.csv/.txt)', fg='navy blue', command= self.file_upload)
        upload_button_csv.grid(row=1,  column=2)
        ###################### Motifs Button
        frame1 = Frame(root); frame1.grid(row=4, column=0)
        compute_button1 = Button(frame1, text='Randomized\nMotif Search', width=15, height=2, fg='blue',
                                 command=lambda: self.Status_search(1)); compute_button1.pack(side=LEFT)
        frame2 = Frame(root); frame2.grid(row=4, column=1)
        compute_button2 = Button(frame2, text='Gibb\'s Sampling\nMotif Search', width=15, height=2, fg='green',
                                 command=lambda: self.Status_search(2)); compute_button2.pack()
        frame3 = Frame(root); frame3.grid(row=4, column=2)
        compute_button3 = Button(frame3, text='Greedy\nMotif Search', width=15, height=2, fg='purple',
                                 command=lambda: self.Status_search(3)); compute_button3.pack()

        frame5 = Frame(root, width=20, height=4).grid(row=5, column=2)
        reset_button = Button(root, text='Reset Space &  Import', width=17, height=1, fg='black', pady=3,
                              command=self.reset); reset_button.grid(row=5, column=2, sticky='S')
        exit_button = Button(root, text='Exit', width=10, height=1, fg='red',
                             command=root.destroy); exit_button.grid(row=5, column=1, sticky='S')
        graph_button = Button(root, text='Show Graph', width=10, height=1, fg='black', pady=10, anchor='s',
                              command=self.open_graph)
        graph_button.grid(row=5, column=0, sticky='S')

        self.label_status = Label(root, text = 'Status: Search Awaits', height=2)
        self.label_status.grid(row=6, column =1)

        Dna = 'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',\
               'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
        self.T = Text(root, height=30, width=117); self.T.insert(END, 'Input Example: \n' + str(Dna))
        self.T.grid(row=7, column=0, columnspan=4, rowspan=2, padx = 2, pady = 2)

    ### Variables
        self.seq_list = Variable()
        self.BestMotifs = Variable()
        self.scores = Variable()
        self.entropy = Variable()

    def file_upload(self):  ############# Import csv or fastq
        Tk().withdraw(); filename = askopenfilename()
        self.T.insert(END, '\n\nImporting...')
        self.T.insert(END, '\nImported file:\n ' + filename + '\n')
        if filename.endswith('.fastq'): T.insert(END, '\nThis function is not ready ...\'\n')
        elif filename.endswith('.csv'):
            my_seq = self.csv2list(filename)
            self.T.insert(END, 'The file you uploaded contains: ' + str(len(my_seq)) + ' processed sequences.\n')
            self.T.insert(END,'\nFirst sequence: ' + my_seq[0] + '\nLast Sequence: ' +
                          my_seq[-1]+'\nLength of first sequence: '+str(len(my_seq[0]))+'\n')
        elif filename.endswith('.txt'):
            my_seq = Tolist(filename)
            self.T.insert(END, 'The file you uploaded contains: ' + str(len(my_seq)) + ' processed sequences.\n')
            self.T.insert(END,'\nFirst sequence: ' + my_seq[0] + '\nLast Sequence: ' +
                          my_seq[-1] + '\nLength of the first sequence: ' + str(len(my_seq[0]))+'\n')
        else:
            self.T.insert(END, 'Help!! Something has gone horribly wrong :< ')
        return self.seq_list.set(my_seq)

    def upload_list(self):
        return seq_list.get()

    def entry_processor(self,Dna_entry, k_entry, N_entry):  ### take .entry as input
        if self.Dna_entry.get():
            Dna_disp = self.Dna_entry.get()
            self.seq = list(ast.literal_eval(Dna_disp)) # turn a 'list string' into real list
        else:
            self.seq = None
        if self.k_entry.get():
            self.k = int(self.k_entry.get())
        else:
            self.k = None
        if not self.N_entry.get():
            N = 100
        else:
            N = int(self.N_entry.get())
        return [self.seq, self.k, N]


    def fixup_thelist(self, mylist):
        min = len(mylist[0])
        lengthlist = [min]
        cleaned_list = []
        for x in (x for x in mylist if 'N' not in x):
            cleaned_list.append(x[:x.find(' ')])
            lengthlist.append(len(x))
            if lengthlist[-1] <= lengthlist[-2]:
                min = lengthlist[-1]
        return [i[0:min - 1] for i in cleaned_list]

    def csv2list(self, filename):
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            seq_list = list(reader)
        my_seq = []
        for i in range(len(seq_list)):
            my_seq.append(seq_list[i][0])
        return self.fixup_thelist(my_seq)

    def display_results(self, seq, k, delta_t, k_bestmotifs, k_scores, k_entropy, searchname):
        self.T.insert(END, '\nYou ran: ' + str(searchname) + ' (k = ' + str(len(seq[0])) + ')\n| k = ' + str(k) +
                 '| Number of Sequence: ' + str(len(seq)) + '| Runtime: ' + delta_t + 's\n\n')
        self.T.insert(END, 'Motifs:\n ' + str(k_bestmotifs) + '\n\nScores:\n' + str(k_scores) + '\n\nEntropy:\n ' + str(
            k_entropy) + '\n')
        self.scores.set(k_scores); self.entropy.set(k_entropy)
        self.label_status.config(text = 'Search Awaits')
        return self.scores, self.entropy, k_bestmotifs

    def upload_choicegate(self,upload):
        t0 = time.time();  # print('Upload is: ' +str(upload));
        if not self.k_entry.get() or (not self.Dna_entry.get() and upload == 0) or (
                upload == 0 and int(self.k_entry.get()) > len(ast.literal_eval(self.Dna_entry.get())[0])):
            messagebox.showinfo('Opps', 'You did something bad ...')
            return None
        else:
            if upload == 1:
                _, self.k, N = self.entry_processor(None, self.k_entry, self.N_entry)
                self.seq = self.seq_list.get()
            else:
                self.seq, self.k, N = self.entry_processor(self.Dna_entry, self.k_entry, self.N_entry)
            return self.seq, self.k, N, t0

    def display_motif_Random(self,upload):
        if self.upload_choicegate(upload) == None: return None
        seq, k, N, t0 = self.upload_choicegate(upload)
        searchname = 'Randomized Motif Search'
        k_bestmotifs, k_scores, k_entropy = self.RandomizedMotifSearch(seq, k, N)
        delta_t = str(time.time() - t0)
        self.display_results(seq, k, delta_t, k_bestmotifs, k_scores, k_entropy, searchname)
        return seq, k_scores, k_entropy, k_bestmotifs

    def display_motif_Gibbs(self, upload):
        if self.upload_choicegate(upload) == None: return None
        seq, k, N, t0 = self.upload_choicegate(upload)
        searchname = 'Gibb\'s Randomized Sampler'
        k_bestmotifs, k_scores, k_entropy = self.GibbsSampler(seq, k, N)
        delta_t = str(time.time() - t0)
        self.display_results(seq, k, delta_t, k_bestmotifs, k_scores, k_entropy, searchname)
        return seq, k_scores, k_entropy, k_bestmotifs

    def display_motif_Greedy(self, upload):
        if self.upload_choicegate(upload) == None: return None
        seq, k, N, t0 = self.upload_choicegate(upload)
        searchname = 'Greedy Motif Search'
        k_bestmotifs, k_scores, k_entropy = self.GreedyMotifSearchWithPseudocounts(seq, k)
        delta_t = str(time.time() - t0)
        self.display_results(seq, k, delta_t, k_bestmotifs, k_scores, k_entropy, searchname)
        return seq, k_scores, k_entropy, k_bestmotifs

    ################################################################## Search Main Functions
    def RandomizedMotifSearch(self, Dna, k_iter, N):
        k_bestmotifs = []
        k_scores = []
        k_entropy = []
        for k in range(2, k_iter + 1):
            t = len(Dna)
            M = self.RandomMotifs(Dna, k, t)
            BestMotifs = M
            for i in range(N):
                Motifprofile, entropy = self.ProfileWithPseudocounts(M)
                # print(Profile)
                M = self.Motifs(Motifprofile, Dna, k)
                this_score = self.Score(M)
                if this_score < self.Score(BestMotifs):
                    BestMotifs = M
                    BestP, Bestentropy = Motifprofile, entropy
            k_bestmotifs.append(BestMotifs)
            k_scores.append(this_score / k)
            k_entropy.append(Bestentropy / k)
        return k_bestmotifs, k_scores, k_entropy

    ################ Gibb's
    def GibbsSampler(self, Dna, k_iter, N):
        BestMotifs = []; k_bestmotifs = []; k_scores = []; k_entropy = []
        t = len(Dna) 
        for k in range(2, k_iter + 1):
            RandMotifs = self.RandomMotifs(Dna, k, t)
            BestMotifs = RandMotifs
            for i in range(1, N):
                tsided_dice = random.randint(0, t - 1); removed_seq = Dna[tsided_dice]
                remained_seqs = Dna[0:tsided_dice] + Dna[tsided_dice + 1:]
                motif_removed = RandMotifs[tsided_dice]
                motif_selected = RandMotifs[0:tsided_dice] + RandMotifs[tsided_dice + 1:]
                motifprofile, entropy = self.ProfileWithPseudocounts(motif_selected)
                generated_motif = self.ProfileGeneratedString(removed_seq, motifprofile, k)
                motif_selected.append(generated_motif)
                Motifs = motif_selected
                this_score = self.Score(Motifs)
                if this_score < self.Score(BestMotifs):
                    BestMotifs = Motifs; Bestentropy = entropy
            k_bestmotifs.append(BestMotifs)
            k_scores.append(this_score / k)
            k_entropy.append(Bestentropy / k)

        return k_bestmotifs, k_scores, k_entropy

    ################ Greedy
    def GreedyMotifSearchWithPseudocounts(self, Dna, k_iter):
        k_bestmotifs = []
        k_scores = []
        k_entropy = []
        for k in range(2, k_iter + 1):
            t = len(Dna)
            BestMotifs = []
            loop_scores = []
            for i in range(0, len(Dna)):
                BestMotifs.append(Dna[i][0:k])
            for i in range(len(Dna[0]) - k + 1):
                Motifs = []
                Motifs.append(Dna[0][i:i + k])

                for j in range(1, len(Dna)):
                    P, entropy = self.ProfileWithPseudocounts(Motifs[0:j])
                    Motifs.append(self.ProfileMostProbableKmer(Dna[j], k, P))
                    this_score = self.Score(Motifs)
                if this_score < self.Score(BestMotifs):
                    BestMotifs = Motifs
                    loop_scores.append(this_score)
                    BestP, entropy = P, entropy
            k_bestmotifs.append(BestMotifs)
            k_scores.append(loop_scores[-1] / k)
            k_entropy.append(entropy / k)
        return k_bestmotifs, k_scores, k_entropy

    ## Small Functions
    def Score(self,Motifs):
        count, k = CountWithPseudocounts(Motifs, 1)
        median_motif = Consensus(Motifs)
        score = 0
        # print(count)
        for i in range(k):
            for j in 'ACGT':
                if j != median_motif[i]:
                    score += count[j][i]
        return score

    def plainScore(self, Motifs):
        count, k = CountWithPseudocounts(Motifs, 0)
        median_motif = Consensus(Motifs)
        score = 0
        # print(count)
        for i in range(k):
            for j in 'ACGT':
                if j != median_motif[i]:
                    score += count[j][i]
        return score

    def Pr(self,Text, Profile):
        prob = 1
        for i in range(len(Text)):
            prob *= self.Profile[Text[i]][i]
        return prob

    def ProfileWithPseudocounts(self,Motifs):
        count, k = CountWithPseudocounts(Motifs, 1)
        num_seq = len(Motifs) + 4
        totentropy = 0
        probcount = count
        for i in range(k):
            for j in 'ACGT':
                probcount[j][i] = count[j][i] / num_seq
                f = probcount[j][i]
                if f != 0:
                    totentropy += -f * m.log2(f)
        return probcount, totentropy

    def Entropy4Logo(self,Motifs):
        count, k = CountWithPseudocounts(Motifs, 0)
        num_seq = len(Motifs) + 4
        entropy_logo = []
        for i in range(k):
            entropy_logo.append([])
            for j in 'ACGT':
                f = float(count[j][i] / num_seq)
                if f > 0:
                    entropy_logo[i].append((j, -f * m.log(f)))
                else:
                    entropy_logo[i].append((j, 0))
        return entropy_logo

    ################################################################## Randomized
    def Motifs(self, Profile, Dna, k):
        motifs = []
        for i in range(len(Dna)):
            mot = self.ProfileMostProbableKmer(Dna[i], k, Profile)
            motifs.append(mot)
        return motifs

    def RandomMotifs(self, Dna, k, t):
        return [text[random.randint(0, len(text) - k):][0:k] for text in Dna]

    ################################################################## Gibb's
    def Normalize(self, Probabilities):
        return {key: Probabilities[key] / sum(Probabilities.values()) for key in Probabilities}

    def WeightedDie(self, Probabilities):
        num = random.uniform(0, 1)
        cumprob = 0
        for i in Probabilities:
            cumprob += Probabilities[i]
            if num < cumprob:
                return i

    def Pr(self,Text, Profile):
        prob = 1
        for i in range(len(Text)):
            prob *= Profile[Text[i]][i]
        return prob

    def ProfileGeneratedString(self,Text, Profile, k):
        probabilities = {}
        for i in range(0, len(Text) - k + 1):
            probabilities[Text[i:i + k]] = self.Pr(Text[i:i + k], Profile)
        return self.WeightedDie(self.Normalize(probabilities))

    ################################################################### Greedy
    def ProfileMostProbableKmer(self, Text, k, Profile):
        probs = []
        topprobs_seq = []
        for i in range(len(Text) - k + 1):
            probs.append(self.Pr(Text[i:i + k], Profile))
        # print(probs)
        for j in range(len(probs)): 
            if probs[j] == max(probs):
                topprobs_seq = Text[j:j + k]
                # print(max(probs))
            # topprob_index = probs.index(topprob)
            if max(probs) == 0:
                topprobs_seq = Text[0:k]
        return topprobs_seq

    #### GUI
    def reset(self):
        self.T.delete('1.0', END)
        return self.label_status.config(text='cleared'), self.seq_list.set('')

    def Status_search(self, number):
        self.label_status.configure(text='Searching for the bottom of your patience...')
        upload_status = 0
        if self.seq_list.get() != '': upload_status = 1;  # print('sequence is: '+ str(seq_list.get()))
        print('Upload is: ' + str(upload_status))
        if number == 1:
            if self.display_motif_Random(upload_status) == None: return None
            _, scores, entropy, motifs = self.display_motif_Random(upload_status)
            self.BestMotifs.set(motifs); self.scores.set(scores); self.entropy.set(entropy)
        if number == 2:
            if self.display_motif_Gibbs(upload_status) == None: return None
            _, scores, entropy, motifs = self.display_motif_Gibbs(upload_status)
            self.BestMotifs.set(motifs); self.scores.set(scores); self.entropy.set(entropy)
        if number == 3:
            if self.display_motif_Greedy(upload_status) == None: return None
            _, scores, entropy, motifs = self.display_motif_Greedy(upload_status)
            self.BestMotifs.set(motifs); self.scores.set(scores); self.entropy.set(entropy)
        return self.label_status.config(text='Status: Done! I''m ready.')

    def Status_running(self):
        status_var = 'Searching for the bottom of your patience...'
        return status_var

    def exit_graph(self, window):
        window.destroy()
        T.insert(END, '\nGraph Closed \n\n')

    def display_sequence(self, seq_T, entry):
        k_input = entry.get()
        int_k = int(k_input)
        if not k_input or int(k_input) > int(k_entry.get()) or int(k_input) < 2:
            messagebox.showinfo('Wait?!', 'Sorry, you cannot do that :(')
        else:
            motifs_k = list(BestMotifs.get())
            motif = motifs_k[int_k - 2]
            [seq_T.insert(END, 'The Motifs with length k = ' + k_input + ' are: \n' + str(
                motif) + '\n' + 'Consensus Motif: ' + str(Consensus(motif)) + '\n\n')]
            return motifs_k.set(motifs_k[int_k - 2])

    def open_graph(self):
        GraphDisplay(self.BestMotifs.get(), self.k, self.scores.get(), self.entropy.get())

def Consensus(Motifs):
    count, k = CountWithPseudocounts(Motifs, 1)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

def CountWithPseudocounts(Motifs, psuedo):  # Count Parallel Motifs Content (get Median)
    count = {}  # initializing the count dictionary
    k = len(Motifs[0])
    t = len(Motifs)
    for symbol in "ACGT":
        count[symbol] = []
        if psuedo == 1:
            count[symbol] = list(np.ones(k, dtype='int'))
        elif psuedo == 0:
            count[symbol] = list(np.zeros(k, dtype='int'))
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count, k

###################################### Dendrogram with Levenstein Distance of motifs
class GraphDisplay:
    def __init__(self, Motifs, k, scores, entropy):
        ### Parameters
        self.k_motif = Variable();self.Motifs = Motifs; self.k = k
        if not scores or not entropy:
            self.label_status.config(text='Status: Except it is not working. See why.')
        else:
            window = Tk(); window.title("Graphs"); window.geometry('820x740')
            # data
            num = len(scores); x = np.linspace(2, num + 1, num)

            fig = plt.figure(figsize=(8, 5)); plt.subplot(2, 1, 1)
            plt.plot(x, scores, 'b-*'); plt.xticks(np.arange(2, k + 1, step=1))
            plt.ylabel('Motif Scores (/base)'); plt.title('Motif Search Results')

            plt.subplot(2, 1, 2); plt.plot(x, entropy, 'r-*'); plt.xticks(np.arange(2, k + 1, step=1))
            plt.xlabel('K of k-mers'); plt.ylabel('Motif Entropy (/base)')

            canvas = FigureCanvasTkAgg(fig, master=window)
            canvas.draw()
            canvas.get_tk_widget().pack(side= tk.TOP, fill = tk.BOTH, expand=1)

            ## Textbox
            window_frame2 = Frame(window, width=650, height=120);window_frame2.pack(side=BOTTOM, expand=1)
            self.seq_T = Text(window_frame2, height=8, width=98); self.seq_T.pack(expand=1)
            self.seq_T.insert(END, 'The motifs resulted are as follow...\n\n')

            window_frame3 = Frame(window, width=300, height=20); window_frame3.pack(side=BOTTOM)
            # Entry
            window_frame1 = Frame(window, width=800, height=20); window_frame1.pack()
            label_seq = Label(window_frame1, text= u"Enter the \'K\' of K-mers to be examined: (<K max)")
            label_seq.pack(side=LEFT)
            self.k_input = Entry(window_frame1); self.k_input.bind("<Return>"); self.k_input.pack(side=LEFT)
            enter_button = Button(window_frame1, text='Enter', command=self.display_sequence)
            enter_button.pack()
            cluster_button = Button(window_frame3, text="Display Motif Dendrogram",
                                    command=self.dendrogram)
            cluster_button.pack(side=RIGHT)
            save_motif_button = Button(window_frame1, text="Save Motifs to File", command=self.savemotif)
            window.mainloop()

    def display_sequence(self):
        inp_k = int(self.k_input.get())
        if not self.k_input or inp_k > self.k or inp_k < 2:
            messagebox.showinfo('Wait?!', 'Sorry, you cannot do that :(')
        else:
            motifs_k = list(self.Motifs)
            motif = motifs_k[inp_k - 2]
            [self.seq_T.insert(END, 'The Motifs with length k = ' + str(inp_k) + ' are: \n' + str(motif) + '\n' +
                          'Consensus Motif: ' + str(Consensus(motif)) + '\n\n')]
            self.k_motif.set(motif)

    def save_motif(self):
        if not self.Motifs:
            messagebox('Hang on', 'You have not produce the motif for a desired k yet')
        with open("Motifs_"+str(datetime.date())) as file:
            file.write('\n'.join(self.Motifs))
    ##
    def dendrogram(self):
        Dendrogram(self.k_motif.get(), self.k_input.get())

    def exit_graph(self,window):
        window.destroy()
        self.T.insert(END, '\nGraph Closed \n\n')

## Produce Dendrogram Graphs  
class Dendrogram:
    def __init__(self, motif, k_input):
        if k_input == '':
            messagebox.showinfo('Actually...', 'You are missing something...')
            return None
        k_s = int(k_input)
        if k_s < 3:
            messagebox.showinfo('Actually...', 'Nobody does that. Really.')
            return None
        else:
            window2 = Tk(); window2.title("Cluster"); window2.geometry('1220x540')

            dist_matrix = self.levenstein_distances(motif, (1, 1, 1))
            fig = plt.figure(figsize=(12, 5), dpi=100, facecolor='w', edgecolor='k')
            Z = hierarchy.linkage(dist_matrix)
            plt.xlabel('(Unique) Motif Candidates ');plt.ylabel('Levinstein Edit Distances (No. of Edits)')
            plt.title('Motif Candidates Cluster')
            hierarchy.dendrogram(Z, leaf_rotation=14, leaf_font_size=7, labels=motif, color_threshold='default')
            fig_x, fig_y = 1, 1
            canvas1 = FigureCanvasTkAgg(fig, master=window2)
            canvas1.draw()
            canvas1.get_tk_widget().pack(side= tk.TOP, fill = tk.BOTH, expand=1)
            # Draw Motif Logo
            self.draw_logo2(self.Entropy4Logo(motif))
            window2.mainloop()

    # Apply Iterative Levenstein to produce distance matrix for creating cluster dendrogram
    def levenstein_distances(self, consensus_motifs, cost):  ### matrix of motif
        motif_dict = {}
        for i in consensus_motifs:
            try:
                motif_dict[i] += 1
            except KeyError:
                motif_dict[i] = 1
        motif_list = list(motif_dict.keys())
        #### Matrix for levenstein_distances
        cols = len(motif_list); rows = cols
        distance_matrix = np.zeros((rows, cols))
        for i in range(rows):
            distance_matrix[i][i] = 0
            for ii in range(i + 1, cols):
                distance_matrix[i][ii] = self.iterative_levenshtein(motif_list[i], motif_list[ii], costs=cost)
                ################# flip matrix
                distance_matrix[ii][i] = distance_matrix[i][ii]
        return distance_matrix

    def iterative_levenshtein(self, s, t, costs=(2, 2, 2)):
        rows = len(s) + 1; cols = len(t) + 1
        deletes, inserts, substitutes = costs
        dist = [[0 for x in range(cols)] for x in range(rows)]
        ### Setup 1st row and 1st column
        for row in range(1, rows):
            dist[row][0] = row * deletes
        for col in range(1, cols):
            dist[0][col] = col * inserts

        for col in range(1, cols):
            for row in range(1, rows):
                if s[row - 1] == t[col - 1]:
                    cost = 0
                else:
                    cost = substitutes
                dist[row][col] = min(dist[row - 1][col] + deletes, dist[row][col - 1] + inserts,
                                     dist[row - 1][col - 1] + cost)  # substitution
        return dist[row][col]

    def draw_logo2(self, all_scores, fontfamily='Arial', size=80):
        COLOR_SCHEME = {'G': 'orange', 'A': 'red', 'C': 'blue', 'T': 'darkgreen'}
        BASES = list(COLOR_SCHEME.keys())

        if fontfamily == 'xkcd':
            plt.xkcd()
        else:
            mpl.rcParams['font.family'] = fontfamily
        window = Tk(); window.title("Cluster"); window.geometry('1120x300')
        
        fig, ax = plt.subplots(figsize=(len(all_scores), 2.5))
        canvas2 = FigureCanvasTkAgg(fig, master=window)
        font = FontProperties(); font.set_size(size); font.set_weight('bold')

        # font.set_family(fontfamily)

        ax.set_xticks(range(1, len(all_scores) + 1)); ax.set_yticks(range(0, 3))
        ax.set_xticklabels(range(1, len(all_scores) + 1), rotation=90)
        ax.set_yticklabels(np.arange(0, 3, 1))
        seaborn.despine(ax=ax, trim=True)

        trans_offset = transforms.offset_copy(ax.transData, fig=fig, x=1, y=0, units='dots')

        for index, scores in enumerate(all_scores):
            yshift = 0
            for base, score in scores:
                score *= 2
                txt = ax.text(index + 1, 0, base, transform=trans_offset, fontsize=80, color=COLOR_SCHEME[base],
                              ha='center', fontproperties=font,

                              )
                txt.set_path_effects([Scale(1.0, score)])
                fig.canvas.draw()
                canvas2.draw()
                window_ext = txt.get_window_extent(txt._renderer)
                yshift = window_ext.height * score
                trans_offset = transforms.offset_copy(txt._transform, fig=fig, y=yshift, units='points')
            trans_offset = transforms.offset_copy(ax.transData, fig=fig, x=1, y=0, units='points')  
        
        canvas2.get_tk_widget().pack(side= tk.TOP, fill = tk.BOTH, expand=1)
        window.mainloop()

    def Entropy4Logo(self, Motifs):
        count, k = CountWithPseudocounts(Motifs, 0)
        num_seq = len(Motifs) + 4
        entropy_logo = []
        for i in range(k):
            entropy_logo.append([])
            for j in 'ACGT':
                f = float(count[j][i] / num_seq)
                if f > 0:
                    entropy_logo[i].append((j, -f * m.log(f)))
                else:
                    entropy_logo[i].append((j, 0))
        return entropy_logo

## scaling for motif diagram
class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

def Tolist(filename):
    with open(filename, 'r') as f:
        reader = csv.reader(f); seq_list = list(reader)
    my_seq = []
    for i in range(len(seq_list)):
        my_seq.append(seq_list[i][0])
    return my_seq

App(rootwindow)
rootwindow.mainloop()