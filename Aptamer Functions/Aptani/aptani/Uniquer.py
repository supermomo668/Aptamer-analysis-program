import sys,math,random,re,subprocess,string,locale,os,time,argparse

parser = argparse.ArgumentParser(
    prog='python3.3 uniquer.py',
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
the most probable aptamer-related binding motifs.
This tool aims at finding unique data from APTANI output."""))
parser.add_argument("Input File", help="APTANI output")
parser.add_argument("-p","--hairpins",action="store_true",help="Elaborates HAIRPINS data")
parser.add_argument("-s","--splitter",type=str, help="The row split element of the input file. Deafult ", default= '\t')

args = parser.parse_args()

filename = sys.argv[-1]

Unique = []
#read input file
file = open(filename,'r')
file = file.readlines()
Output = []
Clear = []
#step that splits the csv line text and clears from spaces and newlines symbols
for item in file:
  pippo = item.replace(' ','').replace('\n','').split(args.splitter)
  if args.hairpins:
    pippo.append(pippo[2]+'_'+pippo[4])
  else:
    pippo.append(pippo[2]+'_'+pippo[5])
  Clear.append(pippo)
#selects the Unique sets
for item in Clear:
  if item[-1] not in Unique:
    Unique.append(item[-1])
    Output.append(item)


fd = open('Elaborated_'+filename,'a')   
old_stdut = sys.stdout
sys.stdout = fd

for item in Output:
  print(', '.join(item))

fd.close


  


