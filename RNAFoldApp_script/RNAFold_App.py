# -*- coding: utf-8 -*-
import os, subprocess, datetime, csv, pathlib
from tkinter import *
from tkinter import messagebox
from tkinter.filedialog import askopenfilename

rootwindow = Tk()
current_folder = pathlib.Path(__file__).parent
RNAFold_path = current_folder / 'Vienna_RNA' / 'RNAfold.exe'
class App:
    def __init__(self, master):

        master.title("Aptamer Stuff"); master.geometry('400x300');
        master.iconbitmap("Space-Alien.ico")
        self.frame = Frame(master, width=350, height=200).grid(rowspan = 7, columnspan = 5,sticky = N)
        Label(master, text = 'Welcome to this duly decorated application\n\nWhat are you doing here?').grid(row = 1, column = 3)
        button0 = Button(self.frame, text = "Fastq Extracter", command = self.SeqExtract, fg='blue'); button0.grid(row = 2, column = 3)
        button1 = Button(self.frame, text = "RNA Fold", command = self.RNAfold_window, fg ='red' ); button1.grid(row = 3, column = 3)
        button2 = Button(self.frame, text = "RNA Interaction", command = self.RNAinteract_window); button2.grid(row  =4, column = 3)
        button3 = Button(self.frame, text="", command=self.RNAother_window); button3.grid(row=5, column=3)
        self.readylabel = Label(self.frame, text=""); self.readylabel.grid(row=7, column=3)
    def RNAfold_window(self):
        GUI_RNAfold();
        rootwindow.lower()
        self.readylabel.config(text = "Launching... to Aptamer wonderland")
    def RNAinteract_window(self):
        self.readylabel.config(text = "Surprise! It is not ready yet")
    def RNAother_window(self):
        self.readylabel.config(text= "You can see, I am a just very suspicious button")
    def SeqExtract(self):
        SequenceExtract();
        rootwindow.lower()
        self.readylabel.config(text = "Extract sequences ...")

from Bio import SeqIO
class SequenceExtract:
    def __init__(self):
        master_rnafold = Tk(); master_rnafold.title("Sequence Extracter"); master_rnafold.geometry('580x660'); master_rnafold.lift()
        frame = Frame(master_rnafold, highlightbackground="grey", highlightcolor="blue", highlightthickness=1,
                      width=550, height=640).grid(rowspan=15, columnspan=8, sticky=N, pady=10, padx=10)
        Label(master_rnafold, text='Your input file: \nRequirement: .fastq file').grid(row=1, column=3);
        button1 = Button(master_rnafold, text='Upload file here (.fastq)', command=self.file_upload); button1.grid(row=1, column=4);
        Label(master_rnafold, text='Desired output file name:  \nRequirement: Filename without extension'
                                   '\n (Optional)').grid(row=2, column=3);

        self.entry_fileout = Entry(master_rnafold, width=30); self.entry_fileout.bind("<Return>");
        self.entry_fileout.grid(row=2, column=4);

        button_fold = Button(master_rnafold, text='Extract sequences', fg='red', command=self.ExtractSeq);
        button_fold.grid(row=3, column=4);

        label = Label(master_rnafold, text='Status: '); label.grid(row=4, column=2); self.label_status = Label(master_rnafold, text="nothing much");
        self.label_status.grid(row=4, column=3);

        self.textbox1 = Text(master_rnafold, height=30, width=60); self.textbox1.insert(END, 'File Info: ');
        self.textbox1.grid(row=5, column=0, columnspan=8, rowspan=3, padx=5)

    def ExtractSeq(self):
        name_fileout = self.entry_fileout.get();
        if name_fileout == '':
            name_fileout = 'MyResults_' + str(datetime.date.today())
        name_fileout = name_fileout + '.txt'
        seq_list = []
        for seq_record in SeqIO.parse(self.filedir, "fastq"):
            seq_list.append(seq_record.seq)
        seqlength = len(seq_list)
        output_directory = os.path.dirname(self.filedir) + r'/' + 'FastqExtract_Results'
        try:
            os.makedirs(output_directory)
        except FileExistsError:
            pass
        with open(output_directory + r'/'+ name_fileout, 'w') as file_out:
            for sequence in seq_list:
                file_out.write("%s\n" % sequence)
        if seq_list:
            process_status = 'Sucess!'
        else:
            self.textbox1.insert(END,'\nSomething is wrong.\n')
            return
        self.textbox1.insert(END, '\nThe process is ...' + process_status + '\n.\n' +
                             'Your output file directory is: ' + output_directory + ' with ' + str(seqlength) +
                             ' sequences.\n First Sequence: '+seq_list[0] + '\nLast Sequence: '+seq_list[-1] +
                             '\n(Length: '+str(len(seq_list[-1]))+' )')

    def file_upload(self):
        self.filedir = askopenfilename(); self.label_status.config(text = 'Importing...');
        if self.filedir.endswith('.fastq'):
            self.label_status.config(text = "will be converted"); self.type = '.fastq'
        else:
            messagebox.showinfo("Wait", "Something is wrong")
            return
        self.textbox1.insert(END, '\nYour input file is located at:'+self.filedir+'\n')
        return self.filedir

class GUI_RNAfold:
    def __init__(self):
        master_rnafold = Tk(); master_rnafold.lift();
        master_rnafold.title("RNA Folder"); master_rnafold.geometry('580x660');
        Frame(master_rnafold, highlightbackground="grey", highlightcolor="blue", highlightthickness=1,
                      width = 550, height = 640).grid(rowspan = 15, columnspan = 8, sticky = N,pady = 10, padx = 10)
        Label(master_rnafold, text = 'Your input file: \nRequirement: Line-delimited sequences').grid(row = 1, column = 3);
        button1 = Button(master_rnafold, text = 'Upload file here (.txt/.csv)',command = self.file_upload); button1.grid(row=1, column =4);
        Label(master_rnafold, text = 'Desired output file name:  \nRequirement: Filename without extension'
                                     '\n (Optional)').grid(row = 2, column = 3);

        self.entry_fileout = Entry(master_rnafold, width = 30); self.entry_fileout.bind("<Return>");
        self.entry_fileout.grid(row = 2, column = 4);

        button_fold = Button(master_rnafold, text = 'Give them the folding',fg='red', command=self.folding);
        button_fold.grid(row = 3, column = 4);

        label = Label(master_rnafold,text = 'Status: '); label.grid(row = 4,column = 2);
        self.label_status = Label(master_rnafold,text = "nothing much"); self.label_status.grid(row = 4, column = 3);

        self.textbox1 = Text(master_rnafold, height=30, width=60); self.textbox1.insert(END, 'File Info: ');
        self.textbox1.grid(row = 5,column = 0,columnspan=8, rowspan=3, padx = 5)
        self.filetype = ''

    def folding(self):
        if self.filetype == '.csv':
            self.textbox1.insert(END, 'Not ready for this yet')
        elif self.filetype == '.txt':
            ## IMPORTANT: File Path of RNAFold
            name_fileout = self.entry_fileout.get();
            if name_fileout == '':
                name_fileout = 'MyResults_'+str(datetime.date.today())
            ### Output Directory
            output_directory = os.path.dirname(self.filedir) + r'/' + 'RNAFold_Results';

            ### Make 'Results' directory
            try:
                os.makedirs(output_directory)
            except FileExistsError:
                pass
            #command_line = str(RNAFold_path) + ' ' + '-i ' + self.filedir
            command_line= [str(RNAFold_path),'-i',self.filedir]
            self.textbox1.insert(END, '.\nYour command line is:\n' + str(command_line)+'\n')
            self.textbox1.insert(END, 'Output Directory:\n' + str(output_directory)+'\n')
            try:
                output_file = output_directory + r'/' +name_fileout
                with open(output_file, 'w') as f:
                    subprocess.check_call(command_line, shell=True, stdout=f, stderr=f, stdin=f)
                with open(output_file, 'r') as f:
                    num_line = [i for i,l in enumerate(f)][-1]/2 + 1
                self.textbox1.insert(END, '\n\nThe process is ...Success!\n.\n.\nHere you go! \n' +
                    'Your output file directory is: ' + str(output_directory) + ' with ' + str(int(num_line)) + ' sequences.')
            except:
                messagebox.showinfo('Hang on!','Something has gone wrong with your input file.')
        else:
            self.textbox1.insert(END, 'Not ready for this yet')

    def file_upload(self):
        self.filedir = askopenfilename(); self.label_status.config(text = 'Importing...');
        if self.filedir.endswith('.csv'):
            self.label_status.config(text = "will be converted"); self.filetype = '.csv'
        elif self.filedir.endswith('.txt'):
            self.label_status.config(text = "all good"); self.filetype = '.txt'
        else:
            messagebox.showinfo("Wait", "Something is wrong")
            return
        self.seqs = Tolist(self.filedir);
        if len(self.seqs) < 1:
            messagebox.showinfo("Wait","Your sequences are spoiled.")
            return
        self.textbox1.insert(END, '\nYour input file is located at:'+self.filedir+'\nNumber of sequences: '+
                             str(len(self.seqs))+ '\nFirst sequence: '+self.seqs[0]+'\nLast Sequence: '+self.seqs[-1]+'\n')



def Tolist(filename):
    with open(filename, 'r') as f:
        reader = csv.reader(f); seq_list = list(reader)
    my_seq = [];
    for i in range(len(seq_list)):
        my_seq.append(seq_list[i][0])
    return my_seq


App(rootwindow)
rootwindow.mainloop()