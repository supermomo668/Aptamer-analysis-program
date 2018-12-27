import os, subprocess, datetime
from tkinter import *
from tkinter.messagebox import showinfo
from tkinter.filedialog import askopenfilename
##check directory
dir = r"C:\Users\Mo\Documents\Aptagen\Aptani\aptani\Vienna_RNA"; os.chdir(dir)
## Run command: os.system(r"C:\Users\Mo\Documents\Aptagen\Aptani\aptani\Vienna_RNA\RNAfold.exe")

#FNULL = open(os.devnull, 'w');    #use this if you want to suppress output to stdout from the subprocess


rootwindow = Tk()

class App:
    def __init__(self, master):

        master.title("Aptamer Stuff"); master.geometry('400x300');
        self.frame = Frame(master, width=350, height=200).grid(rowspan = 7, columnspan = 5,sticky = N)
        Label(master, text = 'Welcome to the duly decorated application\n\nWhat are you doing here?').grid(row = 1, column = 3)
        button1 = Button(self.frame, text = "RNA Fold", command = self.RNAfold_window, fg ='red' ); button1.grid(row = 2, column = 3)
        button2 = Button(self.frame, text = "RNA Interaction", command = self.RNAinteract_window); button2.grid(row  =3, column = 3)
        button3 = Button(self.frame, text = ""); button3.grid(row  =4, column = 3)
        button4 = Button(self.frame, text=""); button4.grid(row=5, column=3)
        self.not_ready = "Surprise! It is not ready yet"
        self.readylabel = Label(self.frame, text=""); self.readylabel.grid(row=7, column=3)
    def RNAfold_window(self):
        GUI_RNAfold()
        self.readylabel.config(text = "Launching... to Aptamer wonderland")
    def RNAinteract_window(self):
        self.readylabel.config(text = self.not_ready)

class GUI_RNAfold:
    def __init__(self):
        master_rnafold = Tk();
        master_rnafold.title("RNA Folder" ); master_rnafold.geometry('580x660');
        frame = Frame(master_rnafold, highlightbackground="grey", highlightcolor="blue", highlightthickness=1,
                      width = 550, height = 640).grid(rowspan = 15, columnspan = 8, sticky = N,pady = 10, padx = 10)
        Label(master_rnafold, text = 'Your input file: \nRequirement: Line-delimited sequences').grid(row = 1, column = 3);
        button1 = Button(master_rnafold, text = 'Upload file here (.txt/.csv)',command = self.file_upload); button1.grid(row=1, column =4);
        Label(master_rnafold, text = 'Desired output file name:  \nRequirement: Please enter filename'
                                     '(without extension)').grid(row = 2, column = 3);

        self.entry_fileout = Entry(master_rnafold, width = 30); self.entry_fileout.bind("<Return>");
        self.entry_fileout.grid(row = 2, column = 4);

        button_fold = Button(master_rnafold, text = 'Give them the folding',fg='red', command=self.folding);
        button_fold.grid(row = 3, column = 4);

        label = Label(master_rnafold,text = 'Status: '); label.grid(row = 4,column = 2);
        self.label_status = Label(master_rnafold,text = "nothing much"); self.label_status.grid(row = 4, column = 3);

        self.textbox1 = Text(master_rnafold, height=30, width=60); self.textbox1.insert(END, 'File Info: ');
        self.textbox1.grid(row = 5,column = 0,columnspan=8, rowspan=3, padx = 5)

    def folding(self):
        if self.type == '.csv':
            self.textbox1.insert(END, 'Not ready for this yet')
        elif self.type == '.txt':
            executable = r"C:\Users\Mo\Documents\Aptagen\Aptani\aptani\Vienna_RNA\RNAfold.exe"
            name_fileout = self.entry_fileout.get();
            if name_fileout == '':
                name_fileout = 'MyResults_'+str(datetime.date.today())
            output_directory = os.path.dirname(self.filedir)+ '/'+name_fileout+'.txt';

            command_line = executable + ' ' + '-i ' + self.filedir
            file_out = open(output_directory, 'w')
            subprocess.call(command_line, shell=False,stdout= file_out)
            self.textbox1.insert(END, '\n.\n.\nHere you go! \n' + 'Your output file directory is: ' + output_directory)
                                 #+ ' with the length of: ' + len(Tolist(output_directory)))

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


import csv
def Tolist(filename):
    with open(filename, 'r') as f:
        reader = csv.reader(f); seq_list = list(reader)
    my_seq = [];
    for i in range(len(seq_list)):
        my_seq.append(seq_list[i][0])
    return my_seq

App(rootwindow)
rootwindow.mainloop()