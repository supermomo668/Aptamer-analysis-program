Welcome to Aptamer App - plus version!

The fully GUI-fied interface aims to take you through finding important motif information in your Aptamer sequences with ease. The GUI should be able to easily guide you through the use of application and prevent user confusion.

Aptamer relies on certain sequence structure that appears as motif among a pool of selected Aptamer to bind to their targets. Such motifs are anticipated to perform certain Bio/physical interaction with the target. Identification of such aptameric motifs help greatly in the design of Aptamers.
Subsequently, it could allow In Silico Maturation of Aptamers to randomly generate more Aptameric candidates that could potentially exceed the originally selected. Increasing diversity of Aptamers are often helpful and could greatly benefit further selection/optimization. See Ikebukuro's lab on Aptamer's ISM technology.

The Aptamer app aims at a 2 part delivery:
1. Motif discovery/identification
	a. Discovery
		i. Gibb's Motifs Sampling: using Markov chain's idea of increasing chances towards more frequent motifs appearance, the algorithm slowly narrow down to likely motifs candidate.
		ii. Randomized Motifs Search: Sample the sequences at random points many times until motifs create a higher scoring motif profile.
	
	b. Identification
		i. Alignment Motif Search: Using a custom local alignment search, the algorithm collects motif identified between random pairs in a pool of sequence and store motifs into a dictionary and score them by frequency of occurring. 

2. ISM (In silico Maturation)
	a. Custom Inputs:
		- Crossovers between sequence
		- Point mutation of sequence 
		- Number of rounds 
		* Motif Protection: Upload the motifs you identified to perform ISM while protecting these motifs that is highly anticipated to be biophyscially functional. Let's leave them alone.

Using the application:
Use any package installer from Anaconda/pip to install all application dependencies in your environment with requirements.txt.

Run/execute AptamerApp-plus.py within the installed environment.


