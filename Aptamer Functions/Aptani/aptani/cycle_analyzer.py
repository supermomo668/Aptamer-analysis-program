import sys,argparse

parser = argparse.ArgumentParser(
    prog='python3 cycle_analyzer.py',
    formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("-i","--input-files",help="Input files to be investigated. Program accepts any number of input files.",nargs='+')
parser.add_argument("-r","--reference",help="Reference file for sequences investigation.")
parser.add_argument("-t","--top",type=int,help="Number of investigated sequences. Default is 10.",default=10)
parser.add_argument("-o","--output",help="Output File.")
args = parser.parse_args()
ref_seq = []
ref_freq = []

with open(args.reference,'r') as file:    #create the reference lists for investigated sequences
    reference = file.readlines()	  #and relative frequencies (ref_seq and ref_freq)
    for i in range(args.top):
        ref_seq.append(reference[i].rstrip('\n').split(':')[0])
        ref_freq.append(reference[i].rstrip('\n').split(':')[1])

Dict = {} 				  #create the dictionary to store data
for i in ref_seq:
    Dict['%s' % i] = []

for item in args.input_files:		  #read every input file and temporary stores
    with open(item,'r') as file:	  #sequences and frequencies
        tested = file.readlines()
        tested_seq = [item.rstrip('\n').split(':')[0] for item in tested]
        tested_freq = [item.rstrip('\n').split(':')[1] for item in tested]
        for i in ref_seq:		  #matches stored sequences with reference sequences
            if i in tested_seq:
                Dict[i].append(tested_freq[tested_seq.index(i)]) #if match retrieves relative frequency and store in dictionary
            else:
                Dict[i].append('Not Found')			 #if not match stores "Not Found" in dictionary
counter = 0

for i in ref_seq:			#append to dictionary the reference data frequence
    Dict[i].append(ref_freq[counter])	
    counter +=1

fd = open(args.output,'w')	#creates output files
old_stdout = sys.stdout
sys.stdout = fd

print('Sequence'+'\t'+'\t'.join(str(x) for x in args.input_files)+'\t'+args.reference)
for key,values in Dict.items():	#print dictionary in tab delimited format into output file
    print(str(key)+'\t'+str(values).replace("['","").replace("', '","\t").replace("']",""))

fd.close()
