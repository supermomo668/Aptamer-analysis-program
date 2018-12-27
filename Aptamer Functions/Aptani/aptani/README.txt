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

APTANI 1.0

Aptani processes a SELEX or Cell-SELEX output and retrieves the most probable aptamer-related binding motifs.

Written and Developed by Jimmy Caroli, University of Modena and Reggio Emilia, Center for Genome Research (CGR), MO.

APTANI is a novel bioinformatic tool able to analyze SELEX/Cell-SELEX data, retrieve the most represented
secondary structure motifs and select those aptamers that shows the most interesting binding potential.

APTANI has been developed in Python 3.3, and requires to run:
- python 3.3 or later release
- sudo privileges
- RNAsubopt as additional program to calculate secondary structure, included in the package
- Clustal Omega to perform MSA (Multi Sequence Alignment), not included in the package
- FigTree to generate the image tree of clustered data, included in the package

*** INSTALLATION & USAGE ***

APTANI does not need any installation process.
Download the package and extract in a folder. This folder will contain every file necessary for the execution of APTANI.
The basic run of APTANI is the following:

python3.3 APTANI.py -c (CYCLES) -r (RANDOMS) --left (LEFT TAG) --right (RIGHT TAG) Input File

We recommend to copy the input file in the same folder where you run APTANI.
APTANI will write output files in the same folder it is run.

*** OPTIONS & COMMANDS ***

APTANI features different command line commands as input. Although they are cited in the help of the program,
here we list and explain them.

-c, --cycles CYCLES: 		Input required. This parameter set the number of cycles to perform in a run of APTANI.
		     		No default parameter for this input. We suggest to set from 50 to 100 cycles for each run.

-r, --randoms RANDOMS:		Input required. This parameter set the number of random sequences picked at each cycle
				of APTANI. No default paramater for this input. Dimension for this parameter may vary
				according to the investigated number of sequences. We suggest to set it not under the 20%
				of the total investigated number of sequences.

-e, --energy ENERGY:		This parameter set the energy range for the RNAsubopt calculation step. Default is 1 (Kcal).
				Increasing this parameter may lead to a more extensive analysis of the secondary structure
				conformation. On the other side, this will generate large dimensions file and slow down
				considerably the calculation process. We suggest to not pass the threshold of 3 (Kcal).

-f, --frequency FREQUENCY:	This parameter set the frequency cutoff for the sequence analysis in APTANI.
				Default for this parameter is 10e-7 (0.0000001). This paramater may change according to
				the quality of your input files and the analyzed cycle of SELEX defined as input file.

-n, --number NUMBER:		This parameter set the lenght of the aptamer sequence to be searhed in the first step
				of the process. Default for this parameter is 99. Lenght of aptamer sequences derived
				from SELEX/Cell-SELEX experiments may vary according to the experiment itself.
				This lenght takes into account both left and right tag lenght, so pay attention at this.

-d, --delete:			Command line parameter. When written as a command, decline the deletion from the disk
				of both files Frequency.csv and Count.csv. Frequency.csv contains all the sequences
				that have passed the frequency cutoff, while Count.csv contains all the sequences
				in the SELEX/Cell-SELEX input file, with the calculated frequency. The default
				configuration is to delete both these files.

--tree:				This command will allow the creation of a clustered tree of the processed
				sequences in APTANI. The sequences that will be clustered are the one who have passed
				the frequency cutoff parameter. The tree will be generated as a .PDF image. 
				Image size and computational time for this function are related to the number of
				sequences to be clustered. We recommend to NOT run this command when the total number of
				sequences analyzed exceeds three to five thousands (3.000 - 5.000).

--width, WIDTH:			This parameter set the number of pixels of width in the generated tree .PDF image.
				The default value is 10.000. Modifications in this parameter may generate larger file
				as output for the clustered tree process.

--height HEIGHT: 		This parameter set the number of pixels of height in the generated tree .PDF image.
				The default value is 10.000. Modifications in this parameter may generate larger file
				as output for the clustered tree process.

--left LEFT:			Parameter that set the left tag to be searched when analyzing the aptamer sequences.
				Aptamer variable sequences generated through evolutionary selection are usually
				wrapped between two tags, left and right. Both tags needs to be set as a input 
				parameter in order to achieve flexibility in data analysis. Each SELEX/Cell-SELEX
				experiment may have different left adn/or right tags. Default is AUGCGG. Input the tag
				in RNA code (U instead of T).

--right RIGHT:			Parameter that set the right tag to be searched when analyzing the aptamer sequences.
				Aptamer variable sequences generated through evolutionary selection are usually
				wrapped between two tags, left and right. Both tags needs to be set as a input 
				parameter in order to achieve flexibility in data analysis. Each SELEX/Cell-SELEX
				experiment may have different left adn/or right tags. Default is CAGACG. Input the tag
				in RNA code (U instead of T).

--max-hmm-iterations:		Number of maximum HMM iterations. This is Clustal Omega parameter, see the help of
				the program for further informations about this. Default is 1.

--cluster-size CLUSTER_SIZE:	Number of maximum individuals for each generated cluster. This is Clustal Omega parameter, 
				see the help of the program for further informations about this. Default is 100.

--variable:			Parameter that allows to investigate ONLY the
                                variable region between the left and right tags of the retrieved aptamer. Left and right tag will be
                                added to the retrieved variable region in order to calculate the secondary structure of the nucleotide string.

--cutoff CUTOFF:		Parameter the define the lenght of the variable region searched in the first step of the program.
                                Use this parameter in your command line if you
				select to investigate only the variable region of the aptamer. Default value
				is set to 30 for this parameter.

*** OUTPUT & ANALYSIS ***

APTANI will write four different files in .csv format.
These files will be:
 
- Hairpins_data.csv
- Intra_Strand_data.csv
- Left_Bulges_data.csv
- Right_Bulges_data.csv

which will contain the data generated for each of the four different loop structures. Data sheet will feature
the query search, the motif retrieved, the identity score, the population percentage for that motif in that loop structure pool 
and the aptamer sequence in which the motif has been found. We suggest to set a threshold for the population percentage
and then sort for identity score, in order to select first the most reprensented motifs and then select those who
where more similar to the query sequence.
The same aptamer sequence may be present in more than one loop structure data sheet. Those sequences are the most interesting.

*** REFERENCES & BUG ***

For any problem or if you find a bug in the code, feel free to mail your questions at:
jimmy.caroli@unimore.it

Enjoy APTANI
