This code is written to search for a motif in a protein or in multiple proteins. Rather than a discrete motif sequence the input format of the motif matches a sequence logo format allowing for a more accurate description of sequence similarity. More information on this and any relevant trouble shooting can be found in the attached word document 'Instructions_for_motif_finder.docx'.

# Packages to install prior to use 
	bio
	pandas 

# Input files
Two input files are needed for this code to work 
	protein_codes.csv
	motif.txt
Example data for these can be found in the repository and must be saved to the same directory as the code. Protein codes can be inputted via accession code (leave sequence field empty) or by sequence (input '#' followed by the chosen name for your protein in access_code field). For more information on the modification and formatting of this data please read the instructions document 'Instructions_for_motif_finder.docx' found in the repository.


# Results output 
Results will print on the screen and also save to two text files: one containing the incomplete found motifs (found_sequences.txt), and one containing just the complete found motifs (results.txt) saved to the same directory as the code. Note that these text files will NOT overwrite with each run and instead append any new results to the end of the textfile. 
