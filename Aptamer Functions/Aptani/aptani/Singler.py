import sys

file = open('Aptamers_list.csv','r')
file = file.readlines()
Unique = []
Stripped = []
for item in file:
 if item not in Unique:
  Unique.append(item)

for item in Unique:
 Stripped.append(item.strip(' ').strip('\n'))

fd = open('Unique.csv','w')
old_stdout = sys.stdout
sys.stdout = fd
for item in Stripped:
 print(item)
fd.close()

