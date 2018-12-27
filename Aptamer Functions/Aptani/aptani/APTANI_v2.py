#!/usr/bin/env python
# -*- coding: utf-8 -*-

def nucleotides():
    global A
    global C
    global T
    global G
    with open('count.csv','r') as t:
        file = t.readlines()
        lista = [item.rstrip('\n') for item in file]
        l1 = []
        number = []
        sequence = []
        A = 0
        C = 0
        T = 0
        G = 0
        for elem in lista:
            l1.extend(elem.strip().split(':'))
        for i in range(len(l1)):
            if i % 2 == 0:
                sequence.append(l1[i])
            else:
                number.append(l1[i])
        for e in range(len(sequence)):
            a = 0
            t = 0
            c = 0
            g = 0
            for i in range(len(sequence[e])):
                if sequence[e][i] == 'A':
                    a += 1
                elif sequence[e][i] == 'T':
                    t += 1
                elif sequence[e][i] == 'G':
                    g += 1
                elif sequence[e][i] == 'C':
                    c += 1
            A = A + (a * float(number[e]))
            C = C + (c * float(number[e]))
            T = T + (t * float(number[e]))
            G = G + (g * float(number[e]))
        pool = A + C + T + G 
        A = (A/pool)*100
        T = (T/pool)*100
        G = (G/pool)*100
        C = (C/pool)*100

def serafini(filename,frequency,cutoff):
    words = []
    fd = open('count.csv','w')
    old_stdout = sys.stdout
    sys.stdout = fd
    with open(filename,'r') as f:
        for line in f:
            match = re.search('{}'.format(args.left.replace('U','T')), line)
            stop = re.search('{}'.format(args.right.replace('U','T')), line)
            if match:
                start = (match.end())
                if stop:
                    end = (stop.start())
                    words.append(line[start:end])
                else:
                    pass
            else:
                pass
        word_freq={}
        tot = len(words)
        for c in words:
            if len(c) > cutoff:
                word_freq[c]=word_freq.get(c,0)+1
            else:
                pass
        lista=sorted(word_freq.items(),key=lambda x:x[1],reverse=True)

        for (word,freq) in lista:                                   # sorting the dictionary
            print('{0}:{1}'.format(word,(freq/tot)))

    fd.close()
    fd = open('subopt_input.csv','w')
    old_stdout = sys.stdout
    sys.stdout = fd
    seq = []
    value = []
    with open('count.csv', 'r') as t:
        for line in t:
            seq.append(line.strip().split(':')[0])
            value.append(line.strip().split(':')[1])
        for i in range(len(seq)):
            if float(value[i]) > frequency:
                print('{}'.format(args.struct_left.replace('U','T'))+seq[i]+'{}'.format(args.struct_right.replace('U','T')))
            else:
                break
    fd.close()

def fastaer():
    with open('subopt_input.csv','r') as A:
        apt = A.readlines()
        fd = open('subopt_fasted.csv','w')
        old_stdout = sys.stdout
        sys.stdout = fd
        for line in apt:
            print('>' + str(line.strip()))
            print(line.strip())
        fd.close()

def cluster():
    global clus
    fastaer()
    os.system('clustalo -i subopt_fasted.csv --clustering-out=clustered.csv --output-order=tree-order & wait') #calls Clustal Omega to create clusters from sequences
    os.system("awk '{print $9,$2}' clustered.csv | sort -n > data.csv & wait") #from Clustal Omega output select only cluster number and seq ID columns
    with open('data.csv','r') as data:
        data = data.readlines()
        clus = [item.rstrip('\n').replace(':','').split(' ')[-1] for item in data]

def selecting(frequency):                # module that selects sequence based on frequency cutoff input    
    fd = open('subopt_input.csv','w')
    old_stdout = sys.stdout
    sys.stdout = fd                 
    seq = []
    value = []
    with open('count.csv', 'r') as t:
        for line in t:
            seq.append(line.strip().split(':')[0])
            value.append(line.strip().split(':')[1])
        for i in range(len(seq)):
            if float(value[i]) > frequency:
                print(seq[i])
            else:
                break
    fd.close() 

def word_count(filename,num):
    global tot
    fd = open('count.csv','w')
    old_stdout = sys.stdout
    sys.stdout = fd
    with open(filename,'r') as f:
        words = []
        for line in f:
            match=re.search(r'[ATGC]{%s}'%num,line)                     #select only DNA sequences 99 bp long
            if match and len(line.strip('\n')) == num:
                words.append(line.strip())  

        word_freq={}
        tot = len(words)
        for c in words:
            word_freq[c]=word_freq.get(c,0)+1
        lista=sorted(word_freq.items(),key=lambda x:x[1],reverse=True)
     
        for (word,freq) in lista:                                   # sorting the dictionary
            print('{0}:{1}'.format(word,(freq/tot))) 
    fd.close()

def tag(file):
    for line in file:
        line = line.strip()
        if args.variable:
            lstart.append(len(args.struct_left))
            rstart.append(len(line)-len(args.struct_right))
        else:
            match = re.search('{}'.format(args.left), line)
            stop = re.search('{}'.format(args.right), line)
            if (line[0] == "A") or (line[0] == "C") or (line[0] == "T") or (line[0] == "G") or (line[0] == "U"):
                if match:
                    lstart.append(match.end())
                    if stop:
                            rstart.append(stop.start())
                    else:
                            del lstart[-1]
                            lstart.append(300)
                            rstart.append(0)
                else:
                    lstart.append(300)
                    rstart.append(0)
            else:
                continue
    return lstart
    return rstart

def hairpin(file):
    with open(file, "r") as t:                                                                          
        count = -1
        hp_patt = re.compile(r"\([^()]+\)")
        ext = t.readlines()                                                                             
        file = [item.rstrip('\n') for item in ext]                                                       
        tag(file)
        for line in file:                                                                               
                if(line[0] == "A") or (line[0] == "C") or (line[0] == "T") or (line[0] == "G") or (line[0] == "U"):        
                        count = count + 1
                        tot_hairpins.append(0)
                else:
                        line = line.strip()
                        for m in hp_patt.finditer(line):                                                
                                if (m.start() > lstart[count]) and (m.end() < rstart[count]):           
                                        hairpins.append(aptamer[count][m.start():m.end()])              
                                        tot_hairpins[count] = tot_hairpins[count] + 1                   
                                else:                                                                   
                                        continue
    return tot_hairpins
    return hairpins

def find_left(file):
    obr = 0
    cbr = 0
    check = -1
    left = re.compile(r"(\(\()")
    right = re.compile(r"(\)\.+\))")
    abort = re.compile(r"(\.+\))")
    count = -1
    global bulge_left_little
    global bulge_right_wide
    global tot_bulge_ll
    with open(file, "r") as t:
        ext = t.readlines()
        file = [item.rstrip('\n') for item in ext]
        for line in file:
                if (line[0] == "A") or (line[0] == "C") or (line[0] == "T") or (line[0] == "G") or (line[0] == "U"):
                    count = count + 1
                    tot_bulge_ll.append(0)
                else:
                        line = line.strip()
                        for m in left.finditer(line):
                            if m.start() > lstart[count]: 
                                bulge_left_little.append(aptamer[count][m.start():m.end()])
                                l = line[m.end():]
                                swap = m.end()
                                abt = re.match(abort, l)
                                if not abt:
                                        for i in l:
                                                check = check + 1
                                                if i == "(":
                                                        obr = obr + 1
                                                elif i == ")":
                                                        cbr = cbr + 1
                                                        if (cbr == obr and obr != 0):    
                                                                op = l[check:]
                                                                swap = swap + check
                                                                sec = re.match(right, op)
                                                                if sec:
                                                                    if sec.end() + swap < rstart[count]:
                                                                        bulge_right_wide.append(aptamer[count][sec.start()+swap:sec.end()+swap])
                                                                        tot_bulge_ll[count] = tot_bulge_ll[count] + 1
                                                                        check = 0
                                                                        break
                                                                    else:
                                                                        del bulge_left_little[-1]
                                                                        break
                                                                else:
                                                                        del bulge_left_little[-1]
                                                                        check = 0
                                                                        break
                                                        else:
                                                                continue
                                                else:
                                                        continue
                                else:
                                    del bulge_left_little[-1]
                            else:
                                continue
    return bulge_left_little
    return bulge_right_wide
    return tot_bulge_ll

def find_right(file):
    obr = 0
    cbr = 0
    check = -1
    left = re.compile(r"(\(\.+\()")
    right = re.compile(r"(\)\))")
    abort = re.compile(r"(\.+\))")
    count = -1
    global bulge_right_little
    global bulge_left_wide
    global tot_bulge_rl
    with open(file, "r") as t:
        ext = t.readlines()
        file = [item.rstrip('\n') for item in ext]
        for line in file:
                if (line[0] == "A") or (line[0] == "C") or (line[0] == "T") or (line[0] == "G") or (line[0] == "U"):
                    count = count + 1
                    tot_bulge_rl.append(0)
                else:
                        line = line.strip()
                        for m in left.finditer(line):
                            if m.start() > lstart[count]: 
                                bulge_left_wide.append(aptamer[count][m.start():m.end()])
                                l = line[m.end():]
                                swap = m.end()
                                abt = re.match(abort, l)
                                if not abt:
                                        for i in l:
                                                check = check + 1
                                                if i == "(":
                                                        obr = obr + 1
                                                elif i == ")":
                                                        cbr = cbr + 1
                                                        if (cbr == obr and obr != 0):    
                                                                op = l[check:]
                                                                swap = swap + check
                                                                sec = re.match(right, op)
                                                                if sec:
                                                                    if sec.end() + swap < rstart[count]:
                                                                        bulge_right_little.append(aptamer[count][sec.start()+swap:sec.end()+swap])
                                                                        tot_bulge_rl[count] = tot_bulge_rl[count] + 1
                                                                        check = 0
                                                                        break
                                                                    else:
                                                                        del bulge_left_wide[-1]
                                                                        break
                                                                else:
                                                                        del bulge_left_wide[-1]
                                                                        check = 0
                                                                        break
                                                        else:
                                                                continue
                                                else:
                                                        continue
                                else:
                                    del bulge_left_wide[-1]
                            else:
                                continue
    return bulge_left_wide
    return bulge_right_little
    return tot_bulge_rl

def find_int(file):
    obr = 0
    cbr = 0
    check = -1
    left = re.compile(r"(\(\.+\()")
    right = re.compile(r"(\)\.+\))")
    abort = re.compile(r"(\.+\))")
    count = -1
    global interior_left
    global interior_right
    global tot_interior
    with open(file, "r") as t:
        ext = t.readlines()
        file = [item.rstrip('\n') for item in ext]
        for line in file:
                if (line[0] == "A") or (line[0] == "C") or (line[0] == "T") or (line[0] == "G") or (line[0] == "U"):
                    count = count + 1
                    tot_interior.append(0)
                else:
                        line = line.strip()
                        for m in left.finditer(line):
                            if m.start() > lstart[count]: 
                                interior_left.append(aptamer[count][m.start():m.end()])
                                l = line[m.end():]
                                swap = m.end()
                                abt = re.match(abort, l)
                                if not abt:
                                        for i in l:
                                                check = check + 1
                                                if i == "(":
                                                        obr = obr + 1
                                                elif i == ")":
                                                        cbr = cbr + 1
                                                        if (cbr == obr and obr != 0):    
                                                                op = l[check:]
                                                                swap = swap + check
                                                                sec = re.match(right, op)
                                                                if sec:
                                                                    if sec.end() + swap < rstart[count]:
                                                                        interior_right.append(aptamer[count][sec.start()+swap:sec.end()+swap])
                                                                        tot_interior[count] = tot_interior[count] + 1
                                                                        check = 0
                                                                        break
                                                                    else:
                                                                        del interior_left[-1]
                                                                        check = 0
                                                                        break
                                                                else:
                                                                    del interior_left[-1]
                                                                    check = 0
                                                                    break
                                                        else:
                                                                continue
                                                else:
                                                        continue
                                else:
                                    del interior_left[-1]
                            else:
                                continue
    return interior_left
    return interior_right
    return tot_interior

def write_HP(filename, listname):
    count = 0
    uniq = []
    fd = open(filename,'w')
    old_stdout = sys.stdout
    sys.stdout = fd
    for i in listname:
        if i not in uniq:
            uniq.append(i)
            count = count + 1
            print('>seq' + str(count))
            print(i)
    fd.close()

def write(filename, listname):
    count = 0
    fd = open(filename,'w')
    old_stdout = sys.stdout
    sys.stdout = fd
    for i in listname:
         count = count + 1
         print('>seq' + str(count))
         print(i)
    fd.close()

def uniques(leftfile, rightfile, leftout, rightout):
    with open(leftfile,'r') as a, open(rightfile,'r') as b:
        left = a.readlines()
        right = b.readlines()
        left = [item.rstrip('\n') for item in left]
        right = [item.rstrip('\n') for item in right]
        uniq = []
        uniq2 = []

        fd = open(leftout,'w')
        old_stdout = sys.stdout
        sys.stdout = fd
    
        for i in range(len(left)):
            if left[i][0] != '>':
                if (left[i] not in uniq) and (right[i] not in uniq2):
                    uniq.append(left[i])
                    uniq2.append(right[i])
                    print(left[i-1])
                    if len(left[i]) > 2:
                        print(left[i][1:-1])
                    else:
                        print(left[i])

                if (left[i] in uniq) and (right[i] not in uniq2):
                    uniq2.append(right[i])
                    print(left[i-1])
                    if len(left[i]) > 2:
                        print(left[i][1:-1])
                    else:
                        print(left[i])

                if (left[i] not in uniq) and (right[i] in uniq2):
                    uniq.append(left[i])
                    print(left[i-1])
                    if len(left[i]) > 2:
                        print(left[i][1:-1])
                    else:
                        print(left[i])

        fd.close()

        uniq = []
        uniq2 = []
        fd = open(rightout,'w')
        old_stdout = sys.stdout
        sys.stdout = fd

        for i in range(len(left)):
            if left[i][0] != '>':
                if (left[i] not in uniq) and (right[i] not in uniq2):
                    uniq.append(left[i])
                    uniq2.append(right[i])
                    print(left[i-1])
                    if len(right[i]) > 2:
                        print(right[i][1:-1])
                    else:
                        print(right[i])

                if (left[i] in uniq) and (right[i] not in uniq2):
                    uniq2.append(right[i])
                    print(left[i-1])
                    if len(right[i]) > 2:
                        print(right[i][1:-1])
                    else:
                        print(right[i])

                if (left[i] not in uniq) and (right[i] in uniq2):
                    uniq.append(left[i])
                    print(left[i-1])
                    if len(right[i]) > 2:
                        print(right[i][1:-1])
                    else:
                        print(right[i])

        fd.close()

def trim(file,out):                             #toglie il primo e l'ultimo nucleotide dalla sequenza prima dell'allineamento
    with open(file,'r') as a:
        left = a.readlines()
        left = [item.rstrip('\n') for item in left]

        fd = open(out,'w')
        old_stdout = sys.stdout
        sys.stdout = fd
    
        for i in range(len(left)):
            if left[i][0] != '>':
                print(left[i-1])
                print(left[i][1:-1])

        fd.close()

def backup(originator,destination):
    try:
        with open(originator,'r') as a, open(destination,'r') as b:
            origin = a.readlines()
            processed = b.readlines()
            origin = [item.rstrip('\n') for item in origin]
            processed = [item.rstrip('\n') for item in processed]
            new = []
            zipped = []
            complete = []

            for i in range(len(processed)):
                for j in range(len(origin)):
                    if (processed[i][0] == '>') and (processed[i] == origin[j]):
                        new.append(processed[i])
                        zipped.extend((origin[j+1][0],processed[i+1],origin[j+1][-1]))

            zipped = [a+b+c for a,b,c in zip(zipped[::3],zipped[1::3],zipped[2::3])]

            for i in range(len(new)):
                complete.append(new[i])
                complete.append(zipped[i])

        fd = open(destination,'w')
        old_stdout = sys.stdout
        sys.stdout = fd

        for i in range(len(complete)):
            print(complete[i])

        fd.close()
    except FileNotFoundError:
        pass

def post_proc(filename, lista):
    global Q
    try:    
        with open(filename,'r') as t:
            for line in t:
                if line[0] != ">":
                    lista.append(line.strip())

            Q = len(lista)
            return lista
    except IOError:
        pass

def single(listname, letter, probable_string, nuc):
    global Q
    global perc_name
    for i in range(len(listname)):
        if listname[i] != 0:
            try:
                result = 100*(listname[i]/Q)
                if (result > perc_name[i]) and (result > nuc):
                    probable_string[i] = letter
                    perc_name[i] = result
                else:
                    continue
            except ZeroDivisionError:
                break
        else:
            continue

def annihilate(lista,name):
    global A,T,C,G
    global K
    global perc_name
    lenght = 0
    fd = open('final_' + name + '.csv' ,'a')
    old_stdout = sys.stdout
    sys.stdout = fd
 
    for i in range(len(lista)):
        if len(lista[i]) > lenght:
            lenght = len(lista[i])
    
    perc_name = [0]*lenght
    a_stri = [0]*lenght
    c_stri = [0]*lenght
    g_stri = [0]*lenght
    t_stri = [0]*lenght
    x_stri = [0]*lenght
    probable_string = [0]*(lenght + 1)

    for i in range(len(lista)):                
        for j in range(len(lista[i])):           
            k = lista[i][j]
            if k == "A":
                a_stri[j] = (a_stri[j] + 1)
            elif k == "C":
                c_stri[j] = (c_stri[j] + 1)
            elif k == "T":
                t_stri[j] = (t_stri[j] + 1)
            elif k == "G":
                g_stri[j] = (g_stri[j] + 1)
            elif k == "-":
                x_stri[j] = (x_stri[j] + 1)
            elif k == "N":
                x_stri[j] = (x_stri[j] + 1)
            elif k == "U":
                t_stri[j] = (t_stri[j] + 1)
            else:
                continue

    single(a_stri, "A", probable_string, A)
    single(g_stri, "G", probable_string, G)
    single(c_stri, "C", probable_string, C)
    single(t_stri, "T", probable_string, T)
    single(x_stri, "-", probable_string, 25)

    new_string = ', '.join(str(x) for x in probable_string)
    new_string = new_string.replace(", ", "")
    top_score = new_string.replace("0", "-")
    K = len(top_score)-1    
    dots = re.compile(r"--")
    vis = re.sub(dots, '', top_score)
    print(vis)
    fd.close()

def coupling(lista1,lista2,name):
    global A,T,C,G
    global K
    global perc_name
    global Score
    lenght = 0
    for i in range(len(lista1)):
        if len(lista1[i]) > lenght:
            lenght = len(lista1[i])
        
    perc_name = [0]*lenght
    a_stri = [0]*lenght
    c_stri = [0]*lenght
    g_stri = [0]*lenght
    t_stri = [0]*lenght
    x_stri = [0]*lenght
    probable_string = [0]*lenght

    for i in range(len(lista1)):                
        for j in range(len(lista1[i])):           
            k = lista1[i][j]
            if k == "A":
                a_stri[j] = (a_stri[j] + 1)
            elif k == "C":
                c_stri[j] = (c_stri[j] + 1)
            elif k == "T":
                t_stri[j] = (t_stri[j] + 1)
            elif k == "G":
                g_stri[j] = (g_stri[j] + 1)
            elif k == "-":
                x_stri[j] = (x_stri[j] + 1)
            elif k == "N":
                x_stri[j] = (x_stri[j] + 1)
            elif k == "U":
                t_stri[j] = (t_stri[j] + 1)
            else:
                continue

    single(a_stri, "A", probable_string, A)
    single(g_stri, "G", probable_string, G)
    single(c_stri, "C", probable_string, C)
    single(t_stri, "T", probable_string, T)
    single(x_stri, "-", probable_string, 25)   

    new_string = ', '.join(str(x) for x in probable_string)
    new_string = new_string.replace(", ", "")
    top_left_score = new_string.replace("0", "-")
    k_left = len(top_left_score)   

    lenght = 0
    for i in range(len(lista2)):
        if len(lista2[i]) > lenght:
            lenght = len(lista2[i])
        
    perc_name = [0]*lenght
    a_stri = [0]*lenght
    c_stri = [0]*lenght
    g_stri = [0]*lenght
    t_stri = [0]*lenght
    x_stri = [0]*lenght
    probable_string = [0]*lenght

    for i in range(len(lista2)):                
        for j in range(len(lista2[i])):           
            k = lista2[i][j]
            if k == "A":
                a_stri[j] = (a_stri[j] + 1)
            elif k == "C":
                c_stri[j] = (c_stri[j] + 1)
            elif k == "T":
                t_stri[j] = (t_stri[j] + 1)
            elif k == "G":
                g_stri[j] = (g_stri[j] + 1)
            elif k == "-":
                x_stri[j] = (x_stri[j] + 1)
            elif k == "N":
                x_stri[j] = (x_stri[j] + 1)
            elif k == "U":
                t_stri[j] = (t_stri[j] + 1)
            else:
                continue

    single(a_stri, "A", probable_string, A)
    single(g_stri, "G", probable_string, G)
    single(c_stri, "C", probable_string, C)
    single(t_stri, "T", probable_string, T)
    single(x_stri, "-", probable_string, 25)

    new_string = ', '.join(str(x) for x in probable_string)
    new_string = new_string.replace(", ", "")
    top_right_score = new_string.replace("0", "-")
    k_right = len(top_right_score)    

    dots = re.compile(r"--")
    right = re.sub(dots, '', top_right_score)
    left = re.sub(dots, '', top_left_score)
    
    K = k_right + k_left
    fd = open('final_' + name + '.csv' ,'a')
    old_stdout = sys.stdout
    sys.stdout = fd
    print(left,'||',right)
    fd.close()

def end(filename,name,total,loop):   #filename: input file   name: output file      total: list of loop for aptamers   loop: global list of loop
  global clus
  with open(filename,'r') as t, open('count.csv','r') as freq:                  
    t = t.readlines()
    freq = freq.readlines()
    value = []
    for line in freq:
        value.append(line.strip().split(':')[1])

    ciao = [item.rstrip('\n').replace('T','U') for item in t]

    fd = open(name + '.csv','a')   
    old_stdut = sys.stdout
    sys.stdout = fd

    for item in ciao:
        count = 0
        inizio  = 0
        times = 0
        uniq = []
        fine = total[0]                  #
        while count != (len(total)-1):
            for l in loop[inizio:fine]:     #    
                punteggio = 0
                if len(item) == len(l):
                    if l not in uniq:
                        uniq.append(l)
                        for i in range(len(item)):
                            if item[i] == '-':              
                                punteggio += 0.5
                            elif item[i] == l[i]:
                                punteggio += 1
                            elif item[i] != l[i]:
                                punteggio = punteggio - 1
                        punteggio = punteggio/len(item)
                        obj = str(l)
                        times = (loop.count(obj)/len(loop))*100
                        if punteggio > 0:
                            print(item,'\t',punteggio,'\t',l,'\t',times,'\t',aptamer[count][:-17],'\t',value[count],'\t',clus[count])

            inizio = fine
            count = count + 1
            fine = fine + total[count]
            times = 0
            
    fd.close()
    
def couple_end(filename,name,total,left,right):
    global clus
    with open(filename,'r') as t, open('count.csv','r') as freq:                          #input per funzione: file da elaborare
        t = t.readlines()
        freq = freq.readlines()
        value = []
        for line in freq:
            value.append(line.strip().split(':')[1])

        ciao = [item.rstrip('\n').replace(' || ',' ').replace(' \t ',' ').replace('T','U').rsplit() for item in t]  #splittare per interior e bulge
        fd = open(name + '.csv','a')                             #input per funzione: file da scrivere
        old_stdut = sys.stdout
        sys.stdout = fd

        for item in ciao:
            count = 0
            inizio  = 0
            times_left = 0
            times_right = 0
            uniq = []
            fine = total[0]                                  #input per funzione: numero loop per aptameri
            while count != (len(total)-1):
                for l,k in zip(left[inizio:fine],right[inizio:fine]):                     #input per funzione: lista totale loop
                    punteggio = 0
                    try:
                        if (len(item[0]) == len(l)) and (len(item[1]) == len(k)):
                            if l and k not in uniq:
                                uniq.append(l)
                                uniq.append(k)
                                for i in range(len(item[0])):
                                    if item[0][i] == '-':              
                                        punteggio += 0.5
                                    elif item[0][i] == l[i]:
                                        punteggio += 1
                                    elif item[0][i] != l[i]:
                                        punteggio = punteggio - 1
                                obj = str(l)
                                times_left =(left.count(obj)/len(left))*100
                                for i in range(len(item[1])):
                                    if item[1][i] == '-':              
                                        punteggio += 0.5
                                    elif item[1][i] == k[i]:
                                        punteggio += 1
                                    elif item[1][i] != k[i]:
                                        punteggio = punteggio - 1
                                obj = str(k)
                                times_right = (right.count(obj)/len(right))*100
                                
                        punteggio = punteggio/(len(item[0]) + len(item[1]))
                    except IndexError:
                        pass
                    try:
                        if punteggio > 0:
                            print(item[0],'||',item[1],'\t',punteggio,'\t',l,'||',k,'\t',times_left,'\t',times_right,'\t',aptamer[count][:-14],'\t',value[count],'\t',clus[count])
                    except IndexError:
                        pass                
                inizio = fine
                count += 1
                fine = fine + total[count]
                    
        fd.close()

def sorter(filename,output):
    with open(filename,'r') as file:
        test = file.readlines()
        test = [item.rstrip('\n') for item in test]
        test.sort()
        uniq = []

        fd = open(output + '_data.csv','a')                             #input per funzione: file da scrivere
        old_stdut = sys.stdout
        sys.stdout = fd

        try:
            if len(test[0].split(' \t ')) == 8:
                print('Consensus Sequence','\t','Identity Score','\t','Retrieved Structure','\t','Population Left Strand','\t','Population Right Strand','\t','Aptamer Sequence','\t','Aptamer Frequency','\t','Cluster')
            else:
                print('Consensus Sequence','\t','Identity Score','\t','Retrieved Structure','\t','Population','\t','Aptamer Sequence','\t','Aptamer Frequency','\t','Cluster')
            for i in test:
                if i not in uniq:
                    uniq.append(i)
                    print(i)       #da controllare se printa anche parentesi o no
                else:
                    continue
        except IndexError:
            print('Sorry, no matches were found for this structure during your analysis.')
        
        fd.close()

import sys,math,random,re,subprocess,string,locale,os,time,argparse

parser = argparse.ArgumentParser(
    prog='python3.3 APTANI.py',
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=("""
 .----------------.  .----------------.  .----------------.  .----------------.  .-----------------. .----------------. 
| .--------------. || .--------------. || .--------------. || .--------------. || .--------------. || .--------------. |
| |      __      | || |   ______     | || |  _________   | || |      __      | || | ____  _____  | || |     _____    | |
| |     /  \     | || |  |_   __ \   | || | |  _   _  |  | || |     /  \     | || ||_   \|_   _| | || |    |_   _|   | |
| |    / /\ \    | || |    | |__) |  | || | |_/ | | \_|  | || |    / /\ \    | || |  |   \ | |   | || |      | |     | |
| |   / ____ \   | || |    |  ___/   | || |     | |      | || |   / ____ \   | || |  | |\ \| |   | || |      | |     | |
| | _/ /    \ \_ | || |   _| |_      | || |    _| |_     | || | _/ /    \ \_ | || | _| |_\   |_  | || |     _| |_    | |
| ||____|  |____|| || |  |_____|     | || |   |_____|    | || ||____|  |____|| || ||_____|\____| | || |    |_____|   | |
| |              | || |              | || |              | || |              | || |              | || |              | |
| '--------------' || '--------------' || '--------------' || '--------------' || '--------------' || '--------------' |
 '----------------'  '----------------'  '----------------'  '----------------'  '----------------'  '----------------' 
Aptani processes a SELEX or Cell-SELEX output and retrieves
the most probable aptamer-related binding motifs"""))
parser.add_argument("-v","--version", action="store_true",help="Print program version.")
parser.add_argument("Input File", help="SELEX/Cell-SELEX output data to be analyzed.  Must be last input in command line. Format accepted is FastQ")
parser.add_argument("-c","--cycles", type=int,help="Number of cycles of APTANI. Default is 50", default=50)
parser.add_argument("-p","--percentage", type=int,help="Percentage of aptamers analyzed in APTANI. Default is 20 percent", default=20)
parser.add_argument("-e","--energy", type=int, help="Energy Range in RNAsubopt calculation. Default is 1", default=1)
parser.add_argument("-f","--frequency",type=float, help="Frequency cutoff for aptamers processing. Default is 0.00001", default=0.00001)
parser.add_argument("-n","--number",type=int,help="Lenght of the Aptamer sequence in SELEX/Cell-SELEX output. Default is 99", default=99)
parser.add_argument("-d","--delete",action="store_false",help="Deletes Count.csv file. Default is TRUE. ATTENTION: this file may occupy several GB", default='TRUE')
parser.add_argument("--tree", action="store_true",help="Creates a PDF tree based upon clusters for the selected sequences")
parser.add_argument("--width",type=int,help="Number of pixels of width of the Cluster image PDF. Default is 10000",default=10000)
parser.add_argument("--height",type=int,help="Number of pixels of height of the Cluster image PDF. Default is 10000",default=10000)
parser.add_argument("--left",type=str, help="Left Tag to be searched. Input as RNA nt [AUCG]. Default is AUGCGG", default= 'AUGCGG')
parser.add_argument("-L","--struct-left", type=str, help="Left tag to be attached during secondary structure calculations. Default is AUGCGG", default= 'AUGCGG')
parser.add_argument("-R","--struct-right", type=str, help="Right tag to be attached during secondary structure calculations. Default is is CAGACG", default='CAGACG')
parser.add_argument("--right",type=str, help="Right Tag to be searched. Input as RNA nt [AUCG]. Default is CAGACG", default='CAGACG')
parser.add_argument("--max-hmm-iterations",type=int, help="Number of maximum HMM iterations. Clustal Omega parameter. Default is 1",default=1)
parser.add_argument("--cluster-size",type=int, help="Number of maximum population for each cluster. Clustal Omega parameter. Default is 100",default=100)
parser.add_argument("-x","--variable",action="store_true",help="Searches and elaborates ONLY the variable region of aptamers comprised between the defined tags")
parser.add_argument("--cutoff",type=int,help="Lenght cutoff for the searched variable region between left and right tags. Default is 30",default=30)
args = parser.parse_args()
if args.version:
    print("APTANI V1.0 developed by Jimmy Caroli\nCenter for Genome Research (CGR)\nUniversity of Modena and Reggio Emilia 2014")
    exit()

tot = 0
HMM = args.max_hmm_iterations
size = args.cluster_size
W = args.width
H = args.height
frequency = args.frequency
cycle_input = args.cycles
cycles = range(1,cycle_input+1)
filename = sys.argv[-1]
cutoff = args.cutoff

if args.variable:
    serafini(filename,frequency,cutoff)
    args.left = args.struct_left
    args.right = args.struct_right
else:
    word_count(filename,args.number)
    selecting(frequency)

nucleotides()
if args.tree:
    fastaer()
    os.system('clustalo -i subopt_fasted.csv --max-hmm-iterations={0} --cluster-size={1} --clustering-out=cluster_out.csv --guidetree-out=tree.csv & wait'.format(HMM,size))
    os.system('java -jar figtree.jar -graphic PDF -width {0} -height {1} tree.csv tree.pdf &'.format(W,H))
os.system('./RNAsubopt -e {0} < subopt_input.csv > subopt_output.csv & wait'.format(args.energy))
K = 0 					
Q = 0
aptamer = []
with open('subopt_output.csv', 'r') as file:      #'subopt_output.csv'
    for line in file:
        if (line[0] == "A") or (line[0] == "C") or (line[0] == "T") or (line[0] == "G") or (line[0] == "U"):
            aptamer.append(line)
        else:
            continue

lstart = []
rstart = []
tot_hairpins = []
hairpins = []
bulge_left_little = []
bulge_left_wide = []
tot_bulge_ll = []
bulge_right_little = []
bulge_right_wide = []
tot_bulge_rl = []
interior_left = []
interior_right = []
tot_interior = []

percentage = args.percentage/100
rand_apt = (len(aptamer)*percentage)*cycle_input
rand_apt = int(rand_apt)
times = range(1,rand_apt+1) 
        
hairpin('subopt_output.csv')           #'subopt_output.csv'
find_left('subopt_output.csv')         #'subopt_output.csv'
find_right('subopt_output.csv')         #'subopt_output.csv'
find_int('subopt_output.csv')         #'subopt_output.csv'
for i in cycles:
    rand_hp = []
    rand_bulge_ll = []
    rand_bulge_lw = []
    rand_bulge_rw = []
    rand_bulge_rl = []            # ogni volta azzera le liste degli item
    rand_int_left = []
    rand_int_right = []           # presi a random dal ciclo successivo
    for i in times:                            # per (times) volte
            n = random.randrange(0,len(aptamer))   # scegli un numero n random
            a = sum(tot_hairpins[0:n-1])           # e poi da quell'n
            a1 = tot_hairpins[n-1]                 # risali agli item nelle diverse 
            for t in hairpins[a:a+a1]:# liste e li aggiungi alle  
                    rand_hp.append(t)                  # liste inizializzate pre ciclo 
            b = sum(tot_bulge_ll[0:n-1])
            b1 = tot_bulge_ll[n-1]
            for t in bulge_left_little[b:b+b1]:
                    rand_bulge_ll.append(t)
            for t in bulge_right_wide[b:b+b1]:
                    rand_bulge_rw.append(t)
            c = sum(tot_bulge_rl[0:n-1])
            c1 = tot_bulge_rl[n-1]
            for t in bulge_right_little[c:c+c1]:        #controllare perchè non sono sicuro
                    rand_bulge_rl.append(t)
            for t in bulge_left_wide[c:c+c1]:
                    rand_bulge_lw.append(t)
            d = sum(tot_interior[0:n-1])
            d1 = tot_interior[n-1]
            for t in interior_left[d:d+d1]:
                    rand_int_left.append(t)
            for t in interior_right[d:d+d1]:
                    rand_int_right.append(t)
                    
    write_HP('hairpin.fasta',rand_hp)         # crea i file .fasta per  
    write('bulge_ll.fasta',rand_bulge_ll)        # ogni set generato a random
    write('bulge_lw.fasta',rand_bulge_lw)
    write('bulge_rl.fasta',rand_bulge_rl)
    write('bulge_rw.fasta',rand_bulge_rw)
    write('interior_left.fasta',rand_int_left)
    write('interior_right.fasta',rand_int_right)

    trim('hairpin.fasta','hp.fasta')                                                                #trim del primo e dell'ultimo nucleotide a singolo input
    uniques('interior_left.fasta','interior_right.fasta','int_left.fasta','int_right.fasta')        #trim del primo e dell'ultimo nucleotide a doppio input
    uniques('bulge_ll.fasta','bulge_rw.fasta','little_left_bulge.csv','wide_right_bulge.fasta')
    uniques('bulge_lw.fasta','bulge_rl.fasta','wide_left_bulge.fasta','little_right_bulge.csv')

    os.system('clustalo -i hp.fasta -o output_hp.csv --force')                          #allineamenti con clustal omega
    os.system('clustalo -i int_left.fasta -o output_interior_left.csv --force')
    os.system('clustalo -i int_right.fasta -o output_interior_right.csv --force')
    os.system('clustalo -i wide_left_bulge.fasta -o output_wide_left.csv --force')
    os.system('clustalo -i wide_right_bulge.fasta -o output_wide_right.csv --force')
                                                                                
    backup('hairpin.fasta','output_hp.csv')                         #ripristino primo e ultimo nucleotide trimmato precedentemente
    backup('bulge_lw.fasta','output_wide_left.csv')                 #primo file --> pre-allineato
    backup('bulge_rw.fasta','output_wide_right.csv')                #secondo file --> output-allineato
    backup('interior_left.fasta','output_interior_left.csv')
    backup('interior_right.fasta','output_interior_right.csv')
     
    int_left_proc = []                                          # interior left
    post_proc('output_interior_left.csv', int_left_proc)        #
    int_right_proc = []                                         # interior right
    post_proc('output_interior_right.csv', int_right_proc)      #
    coupling(int_left_proc,int_right_proc,'interior')           # coupling
    
    hp_proc = []                                                # hairpins
    post_proc('output_hp.csv', hp_proc)                         #
    annihilate(hp_proc,'hairpin')                               # processing
    
    bulge_ll_proc = []                                          # bulge left little
    post_proc('little_left_bulge.csv', bulge_ll_proc)           #
    bulge_rw_proc = []                                          # bulge right wide
    post_proc('output_wide_right.csv', bulge_rw_proc)           #
    coupling(bulge_ll_proc,bulge_rw_proc,'bulge_left_little')   # coupling
    
    bulge_lw_proc = []                                          # bulge left wide
    post_proc('output_wide_left.csv', bulge_lw_proc)            #
    bulge_rl_proc = []                                          # bulge right little
    post_proc('little_right_bulge.csv', bulge_rl_proc)          #
    coupling(bulge_lw_proc,bulge_rl_proc,'bulge_left_wide')     # coupling

clus = []
cluster()

if len(tot_hairpins) != 0:
    end('final_hairpin.csv','output_hairpins',tot_hairpins,hairpins)
    sorter('output_hairpins.csv','Hairpins')
else:
    pass
couple_end('final_interior.csv','output_intra_strand',tot_interior,interior_left,interior_right)
sorter('output_intra_strand.csv','Intra_Strand')
couple_end('final_bulge_left_little.csv','output_right_bulges', tot_bulge_ll, bulge_left_little, bulge_right_wide)
sorter('output_right_bulges.csv','Right_Bulges')
couple_end('final_bulge_left_wide.csv','output_left_bulges', tot_bulge_rl, bulge_left_wide, bulge_right_little)
sorter('output_left_bulges.csv','Left_Bulges')

os.system('rm *.fasta final_* little_* output_*')# subopt_* clustered.csv data.csv')

if args.delete:
    os.system('rm -r count.csv')
else:
    pass
