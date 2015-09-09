# GprotPRED
a tool for detection and classification of G-proteins

GprotPRED is a tool for accurate detection and classification of G-proteins. 
It uses profile Hidden Markov Models (pHMMs) for the four known heterotrimeric Galpha protein families, 
the Gbeta and the Ggamma subunit in order to classify a set of protein sequences 
into the appropriate G-protein family.
This is a standalone version of GprotPRED tool for local execution.

DEPENDENCIES

-hmmer software package (https://www.ebi.ac.uk/Tools/hmmer/software)<br>
-perl


USAGE

	$ GprotPRED.pl -i <input_fasta_file> -o <output_directory> -nofasta <no_fasta_output> -p <selected_profiles>

		i | input_fasta_file		Input file that contains the query protein sequences in fasta format. Required

		o | output_directory		The directory where the result output files will be stored. Required

	nofasta | no_fasta_output	Use this option if you don't want fasta output. Optional, default: fasta files of the predicted proteins are generated

		p | selected_profiles		Use selected profiles only, i.e. -p Gs,Gio. Optional, default: All profiles are used (Galpha,Gs,Gio,Gq11,G1213,Gbeta,Ggamma)
	

EXAMPLES OF USAGE

 	$ GprotPRED.pl -i inFasta.fa -o ./results -p Gs,Gio

Input file: inFasta.fa<br>
Output dir: ./results<br>
Profiles to be run: Gs and Gio only<br>

After the run a new directory has been created inside ./results.<br>
The directory's name is the time-stamp of the run,<br>
i.e. 20150910_135842, which means 10th of September 2015, 13:58:42<br>

Inside ./results/10092015_135842 the below files are stored:<br>
summary.txt - text summary of the results<br>
Gs.res - hmmsearch output for Gs profile<br>
Gio.res - hmmsearch output for Gio profile<br>
fasta_output/Gs.fa - file that contains the proteins that match Gs profile in fasta format<br>
fasta_output/Gio.fa - file that contains the proteins that match Gio profile in fasta format<br>

	$ GprotPRED.pl -i inFasta.fa -o ./results

Input file: inFasta.fa<br>
Output dir: ./results<br>
All profiles will be run by default<br>

After the run a new directory has been created inside ./results.<br>
The directory's name is the time-stamp of the run,<br>
i.e. 20150910_135842, which means 10th of September 2015, 13:58:42<br>

Inside ./results/10092015_135842 the below files are stored:<br>
summary.txt - text summary of the results<br>
Galpha.res - hmmsearch output for Galpha profile<br>
Gs.res - hmmsearch output for Gs profile<br>
Gio.res - hmmsearch output for Gio profile<br>
Gq11.res - hmmsearch output for Gq11 profile<br>
G1213.res - hmmsearch output for G1213 profile<br>
Gbeta.res - hmmsearch output for Gbeta profile<br>
Ggamma.res - hmmsearch output for Ggamma profile<br>
fasta_output/Galpha.fa - file that contains the proteins that match Galpha profile in fasta format<br>
fasta_output/Gs.fa - file that contains the proteins that match Gs profile in fasta format<br>
fasta_output/Gio.fa - file that contains the proteins that match Gio profile in fasta format<br>
fasta_output/G1213.fa - file that contains the proteins that match G1213 profile in fasta format<br>
fasta_output/Gq11.fa - file that contains the proteins that match Gq11 profile in fasta format<br>
fasta_output/Gbeta.fa - file that contains the proteins that match Gbeta profile in fasta format<br>
fasta_output/Ggamma.fa - file that contains the proteins that match Ggamma profile in fasta format

