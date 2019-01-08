from itertools import combinations
import random, csv
from datetime import date
from tkinter import *
from tkinter.filedialog import askopenfilename
rootwindow = Tk()
class App:
    def __init__(self, master):
        master.title("Aptamer ISM");
        master.geometry('400x300');
        self.frame = Frame(master, width=350, height=220).grid(rowspan=7, columnspan=5, sticky=N)
        Label(master, text='Welcome to the duly decorated application\n\nWhat are you doing here?').grid(row=1,column=3)
        button1 = Button(self.frame, text="ISM", command=self.ISMwindow, fg='red');
        button1.grid(row=2, column=3)

        self.not_ready = "Surprise! It is not ready yet"
        self.readylabel = Label(self.frame, text="");
        self.readylabel.grid(row=7, column=3)

    def ISMwindow(self):
        ISM()
        self.readylabel.config(text="Launching... to Aptamer wonderland")

class ISM:
    def __init__(self):
        master_ISM = Tk();
        master_ISM.title("ISM Generation" ); master_ISM.geometry('590x720');
        frame = Frame(master_ISM, highlightbackground="grey", highlightcolor="blue", highlightthickness=1,
                      width = 580, height = 710).grid(rowspan = 18, columnspan = 8, pady = 4, padx = 4,sticky= NW)
        Label(master_ISM, text = 'Your input file: \nNote: Line-delimited sequences').grid(row = 1, column = 3,sticky= NW);
        button1 = Button(master_ISM, text = 'Upload file here (.txt/.csv)',command = self.file_upload); button1.grid(row=1, column =4,sticky= NW);
        Label(master_ISM, text = 'Desired output file name:  \nNote: Enter Only filename'
                                     '(without extension)').grid(row = 2, column = 3);

        self.entry_fileout = Entry(master_ISM, width = 30); self.entry_fileout.bind("<Return>");
        self.entry_fileout.grid(row = 2, column = 4);

        button_fold = Button(master_ISM, text = 'Generate ISM Sequences',fg='Violet', command=self.ISMGenerate);
        button_fold.grid(row = 3, column = 4);

        label = Label(master_ISM,text = 'Status: ',fg = 'Red'); label.grid(row = 4,column = 2);
        ## Parameters
        self.label_status = Label(master_ISM,text = "nothing much"); self.label_status.grid(row = 4, column = 3);
        Label(master_ISM, text = 'Crossover Parameters: ').grid(row = 5, column = 2)
        Label(master_ISM, text='Cross Percentage (0-1)\nDefault:1').grid(row=6, column=3);
        self.entry_CrossPercent = Entry(master_ISM); self.entry_CrossPercent.bind("<Return>");
        self.entry_CrossPercent.grid(row=6, column=4);
        Label(master_ISM, text='No. of Generation\nDefault:1').grid(row=7, column=3);
        self.entry_NumGeneration = Entry(master_ISM); self.entry_NumGeneration.bind("<Return>");
        self.entry_NumGeneration.grid(row=7, column=4);
        Label(master_ISM, text='Maximum No.Sequences\nDefault:1000').grid(row=8, column=3);
        self.entry_MaxCap = Entry(master_ISM); self.entry_MaxCap.bind("<Return>");
        self.entry_MaxCap.grid(row=8, column=4);

        Label(master_ISM, text='Mutation Parameters: ').grid(row=9, column=2)
        Label(master_ISM, text='Mutation Rate(/1000)\nDefault: 20').grid(row=10, column=3);
        self.entry_MutRate = Entry(master_ISM); self.entry_MutRate.bind("<Return>");
        self.entry_MutRate.grid(row=10, column=4);
        Label(master_ISM, text='Mutation Skew (A,C,G,T)\nDefault:(1,1,1,1)').grid(row=11, column=3);
        self.entry_MutSkew = Entry(master_ISM); self.entry_MutSkew.bind("<Return>");
        self.entry_MutSkew.grid(row=11, column=4);

        self.textbox1 = Text(master_ISM, height=13, width=64); self.textbox1.insert(END, 'File Info: ');
        self.textbox1.grid(row = 12,column = 0,columnspan=8, rowspan=5, padx = 5)

    def ISMGenerate(self):
        if self.type == '.csv':
            self.textbox1.insert(END, 'Not ready for this yet')
        elif self.type == '.txt':
            name_fileout = self.entry_fileout.get();
            if name_fileout == '':
                name_fileout = 'MyResults_'+str(datetime.date.today())
            output_directory = os.path.dirname(self.filedir)+ '/'+name_fileout+'.txt';

            file_out = open(output_directory, 'w')
            [file_out.write(i) for i in self.pointmutation(
                self.crossover(Tolist(self.filedir),SwiftProcessEntry(entry_CrossPercent),
                               SwiftProcessEntry(entry_NumGeneration),SwiftProcessEntry(entry_MaxCap)),
                SwiftProcessEntry(entry_MutRate),SwiftProcessEntry(entry_MutSkew))]
            self.textbox1.insert(END, '\n.\n.\nHere you go! \n' + 'Your output file directory is: ' + output_directory)


    def file_upload(self):
        self.filedir = askopenfilename(); self.label_status.config(text = 'Importing...');
        if self.filedir.endswith('.csv'):
            self.label_status.config("will be converted"); self.type = '.csv'
        elif self.filedir.endswith('.txt'):
            self.label_status.config(text = "all good"); self.type = '.txt'
        else:
            tkinter.messagebox.showinfo("Wait", "Something is wrong")
            return
        self.seqs = Tolist(self.filedir);
        self.textbox1.insert(END, '\nYour input file is located at:'+self.filedir+'\nNumber of sequences: '+
                             str(len(self.seqs))+ '\nFirst sequence: '+self.seqs[0]+'\nLast Sequence: '+self.seqs[-1])

    def min_seq_length(self,sequences):  #
        return min([len(x) for x in sequences])

    ######################### Computational Function
    def random_singlept_crossover(self,sequences):  # cross over between 2 sequences at random point
        min_length = min_seq_length(sequences)
        random_point = random.randint(1, min_length)
        new_string1 = sequences[0][:random_point] + sequences[1][:-random_point]
        new_string2 = sequences[1][:random_point] + sequences[0][:-random_point]
        return [new_string1, new_string2]

    def defined_singlept_crossover(self,sequences, position_dict, lengths):  # cross over between 2 sequences at random point
        min_length = min_seq_length(sequences);
        positions1 = position_dict[sequences[0]];
        positions2 = position_dict[sequence[1]]
        for pos1, pos2, length in zip(positions1, positions2, lengths):
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

    def crossover(self,sequences, cross_percent=1, iteration=1, max_cap=1e4):
        # cross_percent: number of population crossover per iteration #max_cap: maximum number of sequence allowed to generate

        new_sequences = [];
        generation_pool = sequences;  # initialize the first pool to cross over
        for i in range(iteration):
            pairpool = list(combinations(generation_pool, 2))
            pairpool_selected = [x for x in pairpool if
                                 random.random() < cross_percent]  # choose a portion from the combination pool
            generation_pool = []
            for pair in pairpool_selected:
                if len(generation_pool) >= max_cap: break  ### condition to stop adding sequences
                if defined == 1:
                    newpair = defined_singlept_crossover(pair, motif_pos_dict,
                                                         mot_lengths);  ### pair : ('acgt','acgt); position:[[3,2],[1,4]]
                else:
                    newpair = random_singlept_crossover(pair);  # print(newpair) #### newpair: ['ACGGT','ACGT]
                generation_pool.append(newpair)
            generation_pool = [seq for seqpair in generation_pool for seq in seqpair]
            new_sequences.append(generation_pool)
            temp_seqs = [seq for sublist in new_sequences for seq in sublist]
            if len(temp_seqs) > max_cap:
                return temp_seqs
        return temp_seqs

    def pointmutation(self, sequences, mutation_rate, skewrate):  # mutation rate is number in a million # skew (A,C,G,T) rate
        if skewrate == '': skewrate = (1, 1, 1, 1)
        tot = sum(skewrate);
        dicerate = list(np.cumsum(skewrate) / tot - 0.25);
        print(dicerate)
        bases = ['A', 'C', 'G', 'T']
        for seq_num in range(len(sequences)):
            for base in range(len(sequences[seq_num])):
                dice = random.random();
                if dice <= mutation_rate / 1e6:
                    dice2 = random.random()
                    for i in range(len(dicerate)):
                        if dice2 > dicerate[i]:
                            seq_base = list(sequences[seq_num])
                            seq_base[base] = bases[i]
                            sequences[seq_num] = ''.join(seq_base)
        return sequences
    def SwiftProcessEntry(self,entry):
        if entry.get() == '':
            return None
        else:
            return entry.get()

def Tolist(filename):
    with open(filename, 'r') as f:
        reader = csv.reader(f); seq_list = list(reader)
    my_seq = [];
    for i in range(len(seq_list)):
        my_seq.append(seq_list[i][0])
    return my_seq

App(rootwindow)
rootwindow.mainloop()