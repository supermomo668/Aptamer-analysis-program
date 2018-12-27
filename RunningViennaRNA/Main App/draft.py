import csv, os
def csv2list(filename):
    with open(filename, 'r') as f:
        reader = csv.reader(f); seq_list = list(reader)
    my_seq = [];
    for i in range(len(seq_list)):
        my_seq.append(seq_list[i][0])
    return my_seq

print(csv2list("asdf.txt"))
print(os.getcwd())