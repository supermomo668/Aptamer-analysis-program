import ast
import csv
import datetime
import math as m
import os
import random
import time
import tkinter as tk
from collections import Counter
from itertools import combinations
from tkinter import *
from tkinter import font, messagebox
from tkinter.filedialog import askopenfilename
from tkinter.font import ITALIC
from tkinter.ttk import Checkbutton, Radiobutton

import matplotlib as mpl
import matplotlib.patheffects
import matplotlib.pyplot as plt
import numpy as np
import seaborn
from matplotlib import transforms
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2Tk)
from matplotlib.font_manager import FontProperties
from scipy.cluster import hierarchy

plt.style.use('seaborn-ticks')

rootwindow = Tk()
class App:
    def __init__(self, master):
        master.title("The Aptamer App: Holy Moti"); master.geometry('400x300')
        master.iconbitmap("MM_icon.ico")
        self.frame = Frame(master, width=350, height=220).grid(rowspan=7, columnspan=5, sticky=N)
        Label(master, text='Welcome to the duly decorated application\n\nWhat are you doing here?').grid(row=1,column=3)
        button2 = Button(self.frame, text = 'Motif Search/Discovery', command = self.MotifSearchWindow,fg='blue')
        button2.grid(row = 2, column = 3)
        button1 = Button(self.frame, text="ISM", command=self.ISMwindow, fg='red')
        button1.grid(row=3, column=3)
        
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
                      width = 580, height = 710).grid(rowspan = 19, columnspan = 8, pady = 4, padx = 4,sticky= NW)
        Label(master_ISM, text = 'Your input file: \nNote: Line-delimited sequences').grid(row = 1, column = 3,sticky= NW)
        button1 = Button(master_ISM, text = 'Upload file here (.txt/.csv)',command = self.file_upload); button1.grid(row=1, column =4,sticky= NW)
        Label(master_ISM, text = 'Desired output file name (w/o extention):  \nDefault: Name_Date.txt').grid(row = 2, column = 3)

        self.entry_fileout = Entry(master_ISM, width = 30); self.entry_fileout.bind("<Return>")
        self.entry_fileout.grid(row = 2, column = 4)
        Label(master_ISM, text = 'In Silico Maturation', font=("Courier",15,"bold italic underline")).grid(row = 3, column=3)
        button_fold = Button(master_ISM, text = 'Generate ISM Sequences',fg='Violet', command=self.ISMGenerate)
        button_fold.grid(row = 12, column = 4)

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
        Label(master_ISM, text="Motif Preservation").grid(row=11, column=2)
        self.button_motif_upload = Button(master_ISM, text = "Upload Motifs (.txt)",
                                state=DISABLED, command=lambda: self.ISM_motif_upload)
        self.button_motif_upload.grid(row=11, column = 3)
        self.use_motif_ism = BooleanVar()
        button_motif_check = Checkbutton(master_ISM, text="w/w/o Motif Preservation",
                                        variable=self.use_motif_ism, command=self.motif_ism_toggle)
        button_motif_check.grid(row=11, column = 4)
        label = Label(master_ISM, text='Status: ', fg='Red'); label.grid(row=12, column=2)
        self.label_status = Label(master_ISM, text="pending"); self.label_status.grid(row=12, column=3)

        self.textbox1 = Text(master_ISM, height=13, width=64); self.textbox1.insert(END, 'File Info: ')
        self.textbox1.grid(row = 13, column = 0, columnspan=8, rowspan=5, padx = 5)
        
    ## ISM Main Generating Button Command
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
                self.crossover(Tolist(self.filedir), self.Entry_processString(self.entry_CrossPercent),
                self.Entry_processString(self.entry_NumGeneration), self.Entry_processString(self.entry_MaxCap))
            mutation_output_sequences = \
                self.pointmutation(crossover_output_sequences, self.Entry_processString(self.entry_MutRate),
                                   list(ast.literal_eval(self.entry_MutSkew.get())))
            output_sequences = mutation_output_sequences; numcount = len(output_sequences)

            [file_out.write(i+'\n') for i in output_sequences]
            self.textbox1.insert(END, '\n.\n.\nHere you go! \n' + 'Your output file directory is: '
                                 + output_directory + '\n with '+str(numcount)+ ' sequences.')

    # Handle Upload of Sequences to undergo ISM
    def file_upload(self):
        self.filedir = askopenfilename(); self.label_status.config(text = 'Importing...')
        if self.filedir.endswith('.csv'):
            self.label_status.config(text = "will be converted"); self.type = '.csv'
        elif self.filedir.endswith('.txt'):
            self.label_status.config(text = "all good"); self.type = '.txt'
        else:
            messagebox.showinfo("Wait", "Something is wrong. Did you upload?")
            return
        self.seqs = Tolist(self.filedir)
        self.textbox1.insert(END, '\nYour input file is located at:'+self.filedir+'\nNumber of sequences: '+
                             str(len(self.seqs))+ '\nFirst sequence: '+self.seqs[0]+'\nLast Sequence: '+self.seqs[-1])
    # Motif ISM toggle command
    def motif_ism_toggle(self):
        if not self.use_motif_ism.get():    # Toggle
            self.button_motif_upload["state"] = 'normal'
            self.label_status.config(text="Motif preservation mode On")
            self.use_motif_ism.set(True)
        else:
            self.button_motif_upload["state"] = 'disabled'
            self.label_status.config(text="Motif preservation mode off")
            self.use_motif_ism.set(False)

    # Handle motif upload
    def ISM_motif_upload(self):
        if self.use_motif_ism:
            self.motif_file = askopenfilename(filetypes=['.txt'])
            self.label_status.config(text = 'Motifs are uploaded. ISM will now preserve motifs')
            if not self.motif_file:
                messagebox.showinfo("Wait", "Something is wrong. Did you upload?")
                return
            self.ISM_motifs = Tolist(self.motif_file)
    
            button_motif_upload["state"] = "normal"

    ######################### ISM Computational Function
    # Input: 2 sequence
    # Output: new sequence in which a random point of crossover is introduced
    def random_singlept_crossover(self,sequences):  # cross over between 2 sequences at random point
        min_length = min([len(x) for x in sequences])   # take shorter sequence of both to simulate actual cross over
        random_point = random.randint(1, min_length)
        new_string1 = sequences[0][:random_point] + sequences[1][:-random_point]
        new_string2 = sequences[1][:random_point] + sequences[0][:-random_point]
        return [new_string1, new_string2]

    # Input: 2 strings of equal length
    # Output: hamming distnace between the two
    def HammingDistance(self, motif, text):
        assert( len(motif) == len(text))
        dist = 0
        for i in range(len(motif)):
            if motif[i] == text[i]:
                dist+=1
        return dist    

    # Input: a string and a text
    # Output: bool array of string overlapping in text if it has sufficiently small distance threshold 
    def ScanHamming(self, motif, text, tol = 0.8):
        minnum_overlap = np.ceil(len(motif)*tol)
        pos_motif = np.full(len(text), False, dtype = bool)    # Initialize boolean array
        overlap_count = self.HammingDistance(motif, text[:len(motif)])
        for i in range(0, len(text)-len(motif)):    # Moving window of comparison
            if motif[0] == text[i]:     # if first letter matches in the current window
                if motif[0] != text[i+1]:   # Subtract overlap count 
                    overlap_count -= 1
            else:                       # if first letter do not match
                if motif[0] == text[i+1]:
                    overlap_count += 1
            if motif[-1] == text[len(motif)+i-1]:   # if last letter match
                if motif[-1] != text[len(motif)+i]:
                    overlap_count -= 1
            else:                                   # if last letter do not match
                if motif[0] == text[i+1]:             
                    overlap_count += 1
            if overlap_count >= minnum_overlap:
                pos_motif[i:i+len(motif)] = True    # All bases spanning across the motif is blocked out (marked TRUE)
        return pos_motif                            # return boolean array when TRUE is where motif is presented    

    # Input: 2 sequences, ism_motifs
    # Output: 2 new sequences after random crossing-over with perservation of motif
    def defined_singlept_crossover(self,sequences, ism_motifs ,toler = 0.8):  # cross over between 2 sequences at random point
        min_length = self.min_seq_length(sequences)     # take shorter sequence of both to simulate actual cross over
        final_motif_pos = np.full(min_length, False, dtype=bool)
        for motif in ism_motifs:
            motif_pos1 = self.ScanHamming(motif, sequences[0])
            motif_pos2 = self.ScanHamming(motif, sequences[1])
            final_motif_pos *= motif_pos1[:min_length] * motif_pos2[:min_length]  # Union so final motif pos contain only those index which does not have motif
        possible_index = np.where(final_motif_pos==False)[0]
        if len(possible_index) < 1:     # No possible crossover point
            return sequences
        else:
            random_point = rrandom.choice(possible_index)  # choose indexes which motif will NOT be affected
            new_string1 = sequences[0][:random_point] + sequences[1][:-random_point]
            new_string2 = sequences[1][:random_point] + sequences[0][:-random_point]
            return [new_string1, new_string2]

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

    # Input: A list of sequence at random position. With an option to cross-over 
    # at a certain motif positions and more options on cross-over
    # output: New Sequences
    def crossover(self,sequences, cross_percent=1, iteration=1, max_cap=1e4, mot_lengths=None):
        # cross_percent: number of population crossover per iteration #max_cap: maximum number of sequence allowed to generate
        iteration=int(iteration); new_sequences = []
        generation_pool = sequences;  # initialize the first pool to cross over
        for i in range(iteration):  # How many iteration  to perform random selection & cross over
            pairpool = list(combinations(generation_pool, 2))   # Create random combinations of pair
            pairpool_selected = [x for x in pairpool if random.random() < cross_percent] # choose a portion from the combination pool
            generation_pool = []
            for pair in pairpool_selected:  # cross over that randomly chosen pair
                if len(generation_pool) >= max_cap: break
                if self.ISM_motifs:
                    newpair = self.defined_singlept_crossover(pair, self.ISM_moifs)
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
    
    # Process entry strings into float
    def Entry_processString(self,entry):
        if entry.get() == '':
            return
        else:
            return float(entry.get())

########################## GUI Class
class MotifSearch:
    def __init__(self):
        root = Tk(); root.geometry("950x700"); root.title('Holy Moti')
        root.iconbitmap("MM_icon.ico")

        label1 = Label(root, text="Sequence: (Format: Python 'List')"); label1.grid(row=0)
        frame_Dna_entry = Frame(root); frame_Dna_entry.grid(row=0, column=1)
        Label(frame_Dna_entry, text="upload").pack(side=LEFT)
        self.Check_upload = BooleanVar(); self.Check_upload.set(False)
        self.Checkbutton_upload = Checkbutton(frame_Dna_entry, variable=self.Check_upload, command=self.check_trigger_upload)
        self.Checkbutton_upload.pack(side=LEFT, padx=10)
        self.Dna_entry = Entry(frame_Dna_entry); self.Dna_entry.bind("<Return>"); self.Dna_entry.pack(side=LEFT)

        Label(root, text="Motif Length Up to/single k (dedicated): (Format: Integer < Kmax)").grid(row=1)
        self.k_entry = Entry(root, state=DISABLED); self.k_entry.bind("<Return>"); self.k_entry.grid(row=1, column=1)

        Label(root, text="No. of Reiteration (Default: 100): (Format: Integer)").grid(row=2)
        frame_N_entry = Frame(root); frame_N_entry.grid(row=2, column=1)
        Label(frame_N_entry, text="Auto").pack(side=LEFT)
        self.Check_N_auto = BooleanVar(); self.Check_N_auto.set(0)
        self.Checkbutton_autoN = Checkbutton(frame_N_entry, variable=self.Check_N_auto, command=self.check_trigger_autoN)
        self.Checkbutton_autoN.pack(side=LEFT, padx=10) 
        self.N_entry = Entry(frame_N_entry, state=DISABLED); self.N_entry.bind("<Return>"); self.N_entry.pack(side=LEFT)

        others_label = Label(root, text='Other Functions'); others_label.grid(row=2, column=3);  ##  Make underlined text
        f = font.Font(others_label, others_label.cget("font")); f.configure(underline=True, slant='italic'); others_label.configure(font=f);
        Button(root,text='G-Quadruplex', state=DISABLED).grid(row = 3, column = 3)

        self.upload_button_fastq = Button(root, text='Upload Fastq(.fastq)', command= self.file_upload, state=DISABLED)
        self.upload_button_fastq.grid(row=0, column=2)
        self.upload_button_csv = Button(root, text='Upload CSV(.csv/.txt)', fg='navy blue', command= self.file_upload, state=DISABLED)
        self.upload_button_csv.grid(row=1,  column=2)
        ###################### Motifs Button
        radio_button_frame = Frame(root); radio_button_frame.grid(row=4, column=0, columnspan= 2)
        Label(radio_button_frame, text="Please select search types", font=font.Font(underline=True)).pack(side=TOP)
        Label(radio_button_frame, text="Definitive    |               Discovery ").pack(side=TOP, padx= 30, anchor=W)
        self.search_type = IntVar(); self.search_type.set(1)
        search_algos = ["Alignment\nRandomized Search",
                        "Randomized\nMotif Search",
                        "Gibb\'s Sampling\nMotif Search"]
        Rbutton_Random = Radiobutton(radio_button_frame, text=search_algos[0], variable= self.search_type, 
                        command=lambda:self.configure_search_entry(1), value=1); Rbutton_Random.pack(side=LEFT, padx = 15, anchor=W)    
        Rbutton_Gibbs = Radiobutton(radio_button_frame, text=search_algos[1], variable= self.search_type, 
                        command=lambda:self.configure_search_entry(2), value=2); Rbutton_Gibbs.pack(side=LEFT, padx = 15, anchor=W)
        Rbutton_Gibbs = Radiobutton(radio_button_frame, text=search_algos[2], variable= self.search_type, 
                        command=lambda:self.configure_search_entry(3), value=3); Rbutton_Gibbs.pack(side=LEFT, padx = 15, anchor=W)

        search_button = Button(root, text = "Search!", width = 15, height = 1, fg="Red", font= font.Font(size=12, weight='bold'),
                                command=self.run_search); search_button.grid(row=4, column = 2, stick='S')
        button_export_motif = Button(root, text= "Export Motifs", width = 8, font=font.Font(weight='bold'), state=DISABLED,
             command=self.export_align_motifs); button_export_motif.grid(row=4, column = 3)

        frame5 = Frame(root, width=20, height=4).grid(row=5, column=2)
        reset_button = Button(root, text='Reset Space &  Import', width=17, height=1, fg='black', pady=3,
                              command=self.reset); reset_button.grid(row=5, column=2, sticky='S')
        exit_button = Button(root, text='Exit', width=10, height=1, 
                             command=root.destroy); exit_button.grid(row=5, column=1, sticky='S')
        self.graph_button = Button(root, text='Show Graph', width=10, height=1, fg='black', pady=10, anchor='s',
                              command=self.open_graph, state=DISABLED)
        self.graph_button.grid(row=5, column=0, sticky='S')

        self.label_status = Label(root, text = 'Status: Search Awaits', height=2)
        self.label_status.grid(row=6, column =1)

        Dna = 'CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA', 'GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG', 'TAGTACCGAGACCGAAAGAAGTATACAGGCGT',\
               'TAGATCAAGTTTCAGGTGCACGTCGGTGAACC', 'AATCCACCAGCTCCACGTGCAATGTTGGCCTA'
        self.T = Text(root, height=25, width=117); self.T.insert(END, 'Input Example: \n' + str(Dna))
        self.T.grid(row=7, column=0, columnspan=4, rowspan=2, padx = 2, pady = 2)

    ### Variables
        self.seq_list = Variable()
        self.BestMotifs = Variable()
        self.scores = Variable()
        self.entropy = Variable()
    
    ############ GUI trigger functions
    # Checkbox triggers
    def check_trigger_upload(self):
        if self.Check_upload.get():
            self.upload_button_fastq['state']= "disabled"
            self.upload_button_csv['state']= "disabled"
            self.Dna_entry['state'] = "normal"
            self.Check_upload.set(False)
        else: 
            self.upload_button_fastq['state']= "normal"
            self.upload_button_csv['state']= "normal"
            self.Dna_entry['state'] = "disabled"
            self.Check_upload.set(True)

    def check_trigger_autoN(self):
        if self.Check_N_auto.get():
            self.N_entry['state'] = "disabled"
            self.N_entry['text'] = "automatic"
            self.Check_N_auto.set(False)
        else:
            self.N_entry['state'] = "normal"
            self.Check_N_auto.set(True)
        
    def configure_search_entry(self, algo):
        if algo == 3:
            self.k_entry['state'] = 'disabled'
            self.N_entry['state'] = 'disabled'
            self.search_type.set(3)
        else:
            self.k_entry['state'] = 'normal'
            self.N_entry['state'] = 'normal'
            if algo == 1:
                self.search_type.set(1)
            else:
                self.search_type.set(2)

    #Import csv or fastq
    def file_upload(self):  
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
            self.T.insert(END, 'Help!! Something has gone horribly wrong :< Please check your upload.')
        return self.seq_list.set(my_seq)

    def upload_list(self):
        return seq_list.get()

    def entry_processor(self, Dna_entry, k_entry, N_entry):  ### take .entry as input
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

    # Input: filename of file with sequences
    # Output: list of sequences in the file
    def csv2list(self, filename):
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            seq_list = list(reader)
        my_seq = []
        for i in range(len(seq_list)):
            my_seq.append(seq_list[i][0])
        return self.fixup_thelist(my_seq)

    # Input: search results 
    # Output: GUI labels of the modified field
    def display_results(self, seq, k, delta_t, k_bestmotifs, k_scores, k_entropy, searchname):
        if searchname == 'Randomized Motif Search' or 'Gibb\'s Randomized Sampler':
            self.T.insert(END, '\n\nYou ran: ' + str(searchname) + ' (k = ' + str(len(seq[0])) + ')\n| k = ' + str(k) +
                    '| Number of Sequence: ' + str(len(seq)) + '| Runtime: ' + delta_t + 's\n\n')
            self.T.insert(END, 'Motifs:\n ' + str(k_bestmotifs) + '\n\nScores:\n' + str(k_scores) + '\n\nEntropy:\n ' + str(
                k_entropy))
            self.scores.set(k_scores); self.entropy.set(k_entropy)
            self.label_status.config(text = 'Search Awaits')
            self.graph_button.config(state=ACTIVE)
        else:   # 
            # k_scores: display as total number of motifs found         # k_entropy: display as number of unique motifs
            self.T.insert(END, '\n\nYou ran: ' + str(searchname) + '| k = ' + str(k) + '\n'+
                    '| Number of Sequence: ' + str(len(seq)) + '| Runtime: ' + delta_t + 's\n\n')  
            self.T.insert(END, 'Motifs:\n ' + str(k_bestmotifs) + '\n\nTotal Number of Motif Moeity Found:\n' + str(k_scores) + 
                        '\n\nMotif Diversity:\n ' + str(k_entropy))                         
            self.label_status.config(text = 'Search Awaits')
            self.graph_button.config(state=DISABLED)
        #return self.scores, self.entropy, k_bestmotifs

    # Input: upload status
    # Output: handle user input and update variables
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

    ## Output: Display on GUI the search results on the text field & return results
    def display_motif(self,upload):
        if self.upload_choicegate(upload) == None: return None
        seq, k, N, t0 = self.upload_choicegate(upload)
        if self.search_type.get() == 3:
            searchname = 'Alignment Motif search'
            self.all_motifs = self.AlignmentMotifSearch(seq)
            delta_t = str(np.round(time.time() - t0, 2))
            self.display_results(seq, k,  delta_t, all_motifs, 
                    sum(all_motifs.values()), len(all_motifs.keys()), searchname) 
            self.export_align_motifs['state'] = "normal"
            return self.all_motifs
        else:
            if self.search_type.get() == 1:
                searchname = 'Randomized Motif Search'
                k_bestmotifs, k_scores, k_entropy = self.RandomizedMotifSearch(seq, k, N)
            elif self.search_type.get() == 2:
                searchname = 'Gibb\'s Randomized Sampler'
                k_bestmotifs, k_scores, k_entropy = self.GibbsSampler(seq, k, N)
            delta_t = str(np.round(time.time() - t0, 2))
            self.display_results(seq, k, delta_t, k_bestmotifs, k_scores, k_entropy, searchname)
        
        return seq, k_scores, k_entropy, k_bestmotifs

    # Output aligmnet motif search results
    def export_align_motifs(self):
        save_dir = askopenfilename(filetypes=[(".txt")])
        with open(save_dir) as file:
            file.write('\n'.join(self.all_motifs.keys()))

    ################################################################## Search Main Functions
    # function to calculate proper N value (auto/input)
    def get_auto_N(self):
        if self.Check_N_auto:
            n_times = int(1000*(num_dna/5)*(10**(avg_k/100)))
        else:
            n_times = int(self.N_entry.get())
            if n_times < 100:
                messagebox.showinfo("Hang on", "Please increase your N to over 100. Higher N yields better results.")
        return n_times
    # Randomized Motif Search iterating over different k length of motif
    def RandomizedMotifSearch(self, Dna, k_iter, N):
        k_bestmotifs = []
        k_scores = []
        k_entropy = []
        N = self.get_auto_N()
        for k in range(2, k_iter + 1):
            BestMotifs, this_score, Bestentropy = self.SingleRandomMotifSearch(Dna, k, N)
            k_bestmotifs.append(BestMotifs)
            k_scores.append(this_score / k)
            k_entropy.append(Bestentropy / k)
        return k_bestmotifs, k_scores, k_entropy

    # Single randomized Motif Search (exahustive as random motifs are reinitilialized on times 
    # based on length of the sequence)
    def SingleRandomMotifSearch(self, Dna, k, N):
        t = len(Dna)
        M = self.RandomMotifs(Dna, k, t)
        BestMotifs = M
        for i in range(N):
            Motifprofile, entropy = ProfileWithPseudocounts(M)
            M = self.Motifs(Motifprofile, Dna, k)
            this_score = self.Score(M)
            if this_score < self.Score(BestMotifs):
                BestMotifs = M
                BestP, Bestentropy  = Motifprofile , entropy
        return BestMotifs, this_score, Bestentropy
    
    #  Alignment Motif Search
    def AlignmentMotifSearch(self, Dna):
        num_dna = len(Dna); avg_k = [len(i) for i in random.choice(Dna, k = 4)]
        # Estimate number of times to run random alignment to make runtime acceptable ~ minute
        n_times = self.get_auto_N()
        self.T.insert(END, "\nCompleted random motif search "+ str(n_times)+ " times.")
        all_motifs_dict = []    # list of dictionary of motifs & counts
        for i in range(n_times):
            random.seed(version = 2)    # seed with time 
            dna_pair = random.choice(Dna, k=2)
            motifs_dict = self.LocalAlignment(dna_pair[0], dna_pair[1])
            all_motifs_dict.append(motifs_dict)
        ## Combine all motifs dictionary
        motifs = Counter()
        for d in all_motifs_dict: 
            motifs.update(d)
        motifs = {k:v for k, v in sorted(motifs.items(), reverse = True)}
        return motifs

    ################ Gibb's Search
    def GibbsSampler(self, Dna, k_iter, N):
        BestMotifs = []; k_bestmotifs = []; k_scores = []; k_entropy = []
        t = len(Dna) 
        N = self.get_auto_N()
        for k in range(2, k_iter + 1):
            RandMotifs = self.RandomMotifs(Dna, k, t)
            BestMotifs = RandMotifs
            for i in range(1, N):
                tsided_dice = random.randint(0, t - 1); removed_seq = Dna[tsided_dice]
                remained_seqs = Dna[0:tsided_dice] + Dna[tsided_dice + 1:]
                motif_removed = RandMotifs[tsided_dice]
                motif_selected = RandMotifs[0:tsided_dice] + RandMotifs[tsided_dice + 1:]
                motifprofile, entropy = ProfileWithPseudocounts(motif_selected)
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

    ###########################################3 Small Functions
    # Input: List of Motifs
    # Output: Distance score between median consensus motif of a profile and count
    def Score(self,Motifs, psuedo=1):
        count, k = CountWithPseudocounts(Motifs, psuedo)
        median_motif = Consensus(Motifs)
        score = 0
        for i in range(k):
            for j in 'ACGT':
                if j != median_motif[i]:
                    score += count[j][i]
        return score

    # Input: a sequence, and motif profile
    # Output: probability of sequence occurring in the profile
    def Pr(self,Text, Profile):
        prob = 1
        for i in range(len(Text)):
            prob *= Profile[Text[i]][i]
        return prob

    # Input: List of Motifs
    # Output: The entropy of motifs
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
    # Use ProfileMostProbablyKmer on each sequence
    def Motifs(self, Profile, Dna, k):
        motifs = []
        for i in range(len(Dna)):
            mot = self.ProfileMostProbableKmer(Dna[i], k, Profile)
            motifs.append(mot)
        return motifs

    # Input: list of dna
    # Output: A random motif of length k among dna
    def RandomMotifs(self, Dna, k, t):
        random.seed(a = None, version=2)
        return [text[random.randint(0, len(text) - k):][0:k] for text in Dna]

    ################################################################## Gibb's
    # Ensure sum of probabilites equal to 1
    def Normalize(self, Probabilities):
        return {key: Probabilities[key] / sum(Probabilities.values()) for key in Probabilities}

    # Input: a list of Probabilities 
    # Output: random number in Probabilities with higher chance of selecting the larger the value 
    def WeightedDie(self, Probabilities):
        num = random.uniform(0, 1)
        cumprob = 0
        for i in Probabilities:
            cumprob += Probabilities[i]
            if num < cumprob:
                return i

    # Input: Dna sequence, profile and k of k-mer
    # Generate a weighted random string based on the profile
    def ProfileGeneratedString(self,Text, Profile, k):
        probabilities = {}
        for i in range(0, len(Text) - k + 1):
            probabilities[Text[i:i + k]] = self.Pr(Text[i:i + k], Profile)
        return self.WeightedDie(self.Normalize(probabilities))

    # Input: Dna, k, target profile
    # Output: a k-mer in Dna that is most likely to appear based on profile
    def ProfileMostProbableKmer(self, Text, k, Profile):
        probs = []
        topprobs_seq = []
        for i in range(len(Text) - k + 1):
            probs.append(self.Pr(Text[i:i + k], Profile))
        for j in range(len(probs)): 
            if probs[j] == max(probs):
                topprobs_seq = Text[j:j + k]
            if max(probs) == 0:     # Null case
                topprobs_seq = Text[0:k]
        return topprobs_seq

    ############################################### GUI functions
    # GUI results display clear
    def reset(self):
        self.T.delete('1.0', END)
        return self.label_status.config(text='cleared'), self.seq_list.set('')

    # GUI update status
    def run_search(self):
        self.label_status.configure(text='Searching for the bottom of your patience...')
        upload_status = 0
        if self.seq_list.get() != '': upload_status = 1;  # print('sequence is: '+ str(seq_list.get()))
        print('Upload is: ' + str(upload_status))
        # number denote the type of search algorithm
        if self.search_type.get() == 0 or 1:
            _, scores, entropy, motifs = self.display_motif(upload_status)            
            self.BestMotifs.set(motifs); self.scores.set(scores); self.entropy.set(entropy)
        elif self.search_type.get() == 2:
            self.AlignmentBestMotifs = self.display_motif(upload_status)
        return self.label_status.config(text='Status: Done! I''m ready.')

    # Set GUI Label Status
    def Status_running(self):
        status_var = 'Searching for the bottom of your patience...'
        return status_var

    # GUI label update: graph close
    def exit_graph(self, window):
        window.destroy()
        T.insert(END, '\nGraph Closed \n\n')

    # GUI display results
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

    # Open new graph window
    def open_graph(self):
        GraphDisplay(self.BestMotifs.get(), self.k, self.scores.get(), self.entropy.get())

#############################3 Alignment Algorithm
#### Local Alignment: Global alignment with 0 weight vertices from source to all nodes
def DNAScoring(GC_bonus):
    bases = ["A", "C", "G", "T"]
    score_mat = pd.DataFrame(data=np.ones((4,4)), columns=bases, index = bases); 
    score_mat["C"] = 1+ GC_bonus; score_mat["G"] = 1 + GC_bonus
    score_mat.loc["C"] = 1 + GC_bonus; score_mat.loc["G"] = 1 + GC_bonus
    return score_mat

class LocalAlignment:
    def __init__(self, ScoringMatrix, indel_penalty):
        self.ScoreMat = ScoringMatrix
        self.indel_penalty = indel_penalty
        print('Size of Scoring Matrix: ', self.ScoreMat.shape)
        print('Types : ', list(self.ScoreMat.columns))
        
    # Form a v x w matrix of scoring map between string1, string2
    def LCSLocal(self, string1, string2):
        self.string1 = string1; self.string2 = string2
        k1 = len(string1); k2 = len(string2)
        self.node_matrix = np.zeros((k1+1, k2+1)); self.backtrack_matrix = np.zeros((k1+1, k2+1))
        print('Backtrack Annotations: {Diagnoal:', 0, 'Down:', 1, 'Right:', 2, 'End(Zero):', 3,'}')
        
        self.backtrack_matrix[0,1:] = np.linspace(2, 2, k2); self.backtrack_matrix[1:,0] = np.linspace(1, 1, k1)
        self.node_matrix[0,1:] = np.linspace(self.indel_penalty, self.indel_penalty*k2, k2)
        self.node_matrix[1:,0] = np.linspace(self.indel_penalty, self.indel_penalty*k1, k1)
        self.max_score = 0; 
        for i in range(0, k1):           ## Along the Row
            for j in range(0, k2):           ## along the column
                diag_vertice = self.node_matrix[i,j]+self.ScoreMat[string1[i]][string2[j]]
                down_vertice = self.node_matrix[i,j+1] + self.indel_penalty
                right_vertice = self.node_matrix[i+1,j] + self.indel_penalty
                all_vertices = [diag_vertice, down_vertice, right_vertice,  0]
                pointer = np.argmax(all_vertices);  
                self.node_matrix[i+1][j+1] = [diag_vertice,down_vertice,right_vertice, 0][pointer]
                
                #### Find Max Cell + Score (Local only)
                current_nodevalue = self.node_matrix[i+1][j+1]    
                if current_nodevalue >= self.max_score:
                    self.max_score = current_nodevalue
                    self.max_node = (i+1,j+1)
                #print('Max Node: ',self.max_node)
                #### Make Backtrack Matrix
                if (i <= k1-1) and (j <= k2-1):
                    self.backtrack_matrix[i+1,j+1] = pointer
        return self.node_matrix, self.backtrack_matrix, self.max_score
    
    # Tolerance for when the the alignment does not increase: 0 , no gap allowed
    def BacktrackLocalAlign(self, max_node = None, tolerance=0):   
        if not max_node:
            i, j = self.max_node
        else:
            i, j = max_node
        align_tracks = set()   # Set object of nodes travelled
        current_tolerance = tolerance
        score = 0; align1 = ''; align2 = ''
        while (i >= 1) or (j >= 1) or tolerance > 0:
            current_node = (i,j)
            align_tracks.add(current_node)
            base1 = self.string1[i-1]; base2 = self.string2[j-1]
            if self.backtrack_matrix[i][j] == 0:   # Diagonal
                align1 = base1 + align1; align2 = base2 + align2
                score += self.ScoreMat[base1][base2]
                i -= 1; j -= 1; 
            elif self.backtrack_matrix[i][j] == 1:  # Up
                align1 = base1 + align1; align2 = '-' + align2
                score += self.indel_penalty
                i -= 1
            elif self.backtrack_matrix[i][j] == 2:   # Left
                align1 = '-' + align1; align2 = base2 + align2
                score += self.indel_penalty
                j -= 1
            elif self.backtrack_matrix[i][j] == 3:   # End immediately
                i=0; j=0; 
                score += 0
            if self.node_matrix[i-1][j-1] >= self.node_matrix[current_node[0]-1][current_node[1]-1]:
                save_last_align1 = align1; save_last_align2 = align2
                current_tolerance -= 1
            else:
                current_tolerance = tolerance 
        return align1, align2, score, align_tracks
                
        # For Finding multiple alignments along 2 string (including suboptimal)
        def AllAlign(self, match_score_lim = 3):
            match_lim = match_lim_base + np.floor(0.2*len(self.string1))
            max_pos = {}   # position for maximum scores to start backtracking at
            all_motifs = {}  # record all alignment motifs
            # Go through node score matrix
            for i in range(len(self.string1), -1, -1):
                for j in range(len(self.string2), -1, -1):
                    if self.node_matrix[i][j] >= match_score_lim:
                        try:
                            max_pos[self.node_matrix[i][j]] = [(i, j)]
                        except KeyError:
                            max_pos[self.node_matrix[i][j]].append((i,j))

            # Sort dictionary by value (ascending)
            max_pos = {k:v for k, v in sorted(max_pos.items(), reverse = True)}
            all_align_tracks = set()
            for score, max_nodes in max_pos.items():
                for max_node in max_nodes:
                    if max_node not in all_align_tracks:
                        align1, _, _, align_track = self.BacktrackLocalAlign(max_node=max_node)
                        try:
                            all_motifs[align1] += 1
                        except KeyError:
                            all_motifs[align1] = 1
                        all_align_tracks = all_align_tracks.add(align_track)        
            return all_motifs

#################################################### Common Functions
# Find the most general consensu among a list of motifs
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

# Input: list of motifs
# Output: Profile of motifs (in probabilies), total entropy
def ProfileWithPseudocounts(Motifs, psuedo = 1):
    count, k = CountWithPseudocounts(Motifs, psuedo)
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

# Input: List of Motifs , bool of pseduo count
# Output: dictionary of lists of counts| columns: bases (key) ; row: base count of sequences (value)
def CountWithPseudocounts(Motifs, psuedo):  
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
            self.save_motif_button = Button(window_frame3, text="Save Motifs to File", 
                                            state=DISABLED, command=self.save_motif)
            self.save_motif_button.pack(side=BOTTOM, pady=3, padx = 10)
            self.first_display = True
            window.mainloop()

    # Displaying motifs search information
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
            self.save_motif_button["state"] = "normal"
            if self.first_display:
                self.seq_T.insert(END, "You may save the above result with \"Save Motifs to File\" button.")
                self.first_display = False
            

    # Allow motifs to be saved to a file
    def save_motif(self):
        if not self.Motifs:
            messagebox('Hang on', 'You have not produce the motif for a desired k yet')
            return
        save_dir = askopenfilename(filetypes=[(".txt")])
        with open(save_dir) as file:
            file.write('\n'.join(self.Motifs))
    
    # Pass motif information to Dendrogram class
    def dendrogram(self):
        self.dend_data = Dendrogram(self.k_motif.get(), self.k_input.get())

    # Exit graph window & update GUI label
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
            messagebox.showinfo('Warning', 'Your k input is too small.')
            return None
        else:
            window2 = Tk(); window2.title("Cluster"); window2.geometry('1220x540')

            dist_matrix = self.levenstein_distances(motif, (1, 1, 1))
            fig = plt.figure(figsize=(12, 5), dpi=100, facecolor='w', edgecolor='k')
            Z = hierarchy.linkage(dist_matrix)
            plt.xlabel('(Unique) Motif Candidates ');plt.ylabel('Levinstein Edit Distances (No. of Edits)')
            plt.title('Motif Candidates Cluster')
            self.cluster_Data = hierarchy.dendrogram(Z, leaf_rotation=14, leaf_font_size=7, labels=motif,\
                 color_threshold='default')
            fig_x, fig_y = 1, 1
            canvas1 = FigureCanvasTkAgg(fig, master=window2)
            canvas1.draw()
            canvas1.get_tk_widget().pack(side= tk.TOP, fill = tk.BOTH, expand=1)
            # Draw Motif Logo
            self.draw_logo2(self.Entropy4Logo(motif))
            window2.mainloop()
    
    def data(self):
        return self.cluster_Data

    # Output: Levenstein distance matrix for creating cluster dendrogram
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

    #  Output: Levenshtein distance between 2 sequence
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

    # plot motif profile logo
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

    # Calculate entropy for Logo generation
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

## scaling class for motif diagram
class Scale(matplotlib.patheffects.RendererBase):
    def __init__(self, sx, sy=None):
        self._sx = sx
        self._sy = sy

    def draw_path(self, renderer, gc, tpath, affine, rgbFace):
        affine = affine.identity().scale(self._sx, self._sy)+affine
        renderer.draw_path(gc, tpath, affine, rgbFace)

## Read in line delimited text file
# Output: List of sequences
def Tolist(filename):
    with open(filename, 'r') as f:
        reader = csv.reader(f); seq_list = list(reader)
    my_seq = []
    for i in range(len(seq_list)):
        my_seq.append(seq_list[i][0])
    return my_seq

App(rootwindow)
rootwindow.mainloop()
