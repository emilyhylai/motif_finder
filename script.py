########## importing packages ##########
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd


########## reading in protein sequences ##########

# function to return fasta sequence from accession code
def fetch_sequence(accession_code):
    Entrez.email = "ug19932@bristol.ac.uk"
    handle = Entrez.efetch(db="protein", id=accession_code, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    sequence = record.seq
    return sequence

# reading in list of accession codes
protein_df = pd.read_csv('protein_codes.csv').dropna(how='all') 

# adding sequence to protein_df
temp_seq_list = []

for i in range(protein_df['access_code'].count()):
    if protein_df.loc[i,'access_code'].count('#') == 0:
        temp_seq_list.append(fetch_sequence(protein_df.loc[i,'access_code']))
    else:
        bio_sequence = Seq(protein_df.loc[i,'sequence'])       
        temp_seq_list.append(bio_sequence)      

protein_df['sequence'] = temp_seq_list


########## reading in motif sequences ###########

# reading motif text file and putting into a list 
motif_file = open('motif.txt', 'r')
motif_list = motif_file.read().split('\n')
motif_file.close()

# makes a list of possible amino acids in each position 
motif_pos_list = []
for i_position in range(len(motif_list[0])):
    temp_string = []
    for i_motifs in range(len(motif_list)):
        temp_string.append(motif_list[i_motifs][i_position])
    motif_pos_list.append(''.join(temp_string))

# makes a dataframe with a count of each amino acid in each position     
motif_df = pd.DataFrame(columns=['position','A','R','N','D','C','Q','G','E','H','I','L','K','M','F','P','S','T','W','Y','V'])
all_amino_acids = ['A','R','N','D','C','Q','G','E','H','I','L','K','M','F','P','S','T','W','Y','V']

for i_position in range(len(motif_pos_list)):
    new_row = {'position':i_position,
               'A':motif_pos_list[i_position].count('A'),
               'R':motif_pos_list[i_position].count('R'),
               'N':motif_pos_list[i_position].count('N'),
               'D':motif_pos_list[i_position].count('D'),
               'C':motif_pos_list[i_position].count('C'),
               'Q':motif_pos_list[i_position].count('Q'),
               'G':motif_pos_list[i_position].count('G'),
               'E':motif_pos_list[i_position].count('E'),
               'H':motif_pos_list[i_position].count('H'),
               'I':motif_pos_list[i_position].count('I'),
               'L':motif_pos_list[i_position].count('L'),
               'K':motif_pos_list[i_position].count('K'),
               'M':motif_pos_list[i_position].count('M'),
               'F':motif_pos_list[i_position].count('F'),
               'P':motif_pos_list[i_position].count('P'),
               'S':motif_pos_list[i_position].count('S'),
               'T':motif_pos_list[i_position].count('T'),
               'W':motif_pos_list[i_position].count('W'),
               'Y':motif_pos_list[i_position].count('Y'),
               'V':motif_pos_list[i_position].count('V')}
    motif_df.loc[len(motif_df)] = new_row

# normalises the counts to the max amino acids 
motif_df.set_index('position',inplace=True, drop=True)    
max_values = motif_df.max(axis=1)
motif_df = motif_df.divide(max_values, axis="rows")


########## searching through the sequences and writing results to a text file ##########

# variables to hold full motif data 
motif_df_full = motif_df
full_motif_length = len(motif_df_full)
print('####### motif scores #######')
results_motif_df = motif_df.round(2)
results_motif_df = results_motif_df.replace(0,'-')
print(results_motif_df.T)
print()

# initialising results dataframe 
results_df = pd.DataFrame(columns=['accession code', 'found position', 'found motif','score'])

# printing title for results 
print('####### All found sequences #######') 
print()

# looping through each protein in the df
for i_protein in range(protein_df['access_code'].count()):
    sequence = protein_df.loc[i_protein,'sequence']
    motif_df = motif_df_full # reset to full motif (the df is cut short later to find every permutation 
    df = pd.DataFrame(columns=['length','location','score','rel score', 'seq']) # defining df for incomplete sequences 
    
    # loop to remove starting amino acid in motif to check for each permutation  
    for i_remove in range(len(motif_df)): 
        
        # makes 2D dictionary for motif with each possible amino acid in each position and their respective scores 
        motif_dict = {}
        for i_position in range(len(motif_df)):
            for i_aminoacids in range(len(motif_df.columns)):
                aminoacid = motif_df.columns[i_aminoacids]
                percentage = motif_df.loc[i_position,aminoacid]
                if motif_df.loc[i_position,aminoacid] != 0:
                    if i_position in motif_dict:
                        motif_dict[i_position].update({aminoacid:percentage})
                    else: 
                        motif_dict[i_position] = {aminoacid:percentage}       
        motif_dict = {k: dict(sorted(v.items(), key=lambda x: x[1], reverse=True)) for k, v in motif_dict.items()}           
        
    
        # makes 2D list of possible amino acids in each position no scores 
        aa_2d_list = [] # list of amino acids 
        for i_motif_pos in range(len(motif_dict)):
            aa_2d_list.append([key[0] for key in motif_dict[i_motif_pos].keys()])
            
        motif_length = len(aa_2d_list)

        # searching for the first amino acid in motif and making a list with every found location
        first_aa_location = []   
        for i_aa in range(len(aa_2d_list[0])): # looping through each possiblle amino acid in first position 
            seq_pointer = 0 
            location = 0 
            aa = aa_2d_list[0][i_aa] 
            temp_first_aa_location = [] # 2D list holding location of each possible amino acid in first position of motif 
            while location != -1:
                location = sequence.find(aa, seq_pointer)
                seq_pointer = location + 1
                if location != -1:
                    temp_first_aa_location.append(location)
            first_aa_location.append(temp_first_aa_location)

        # looping through each possible amino acid in first position 
        for i_first_aa in range(len(first_aa_location)): 
            aa = aa_2d_list[0][i_first_aa] 
            # looping through each found location 
            for i_location in range(len(first_aa_location[i_first_aa])): 
                score = motif_dict[0][aa] #setting inital score 
                length_count = 1
                seq_pointer = first_aa_location[i_first_aa][i_location] 
                found_seq = sequence[seq_pointer]
                isfound = True 
                for motif_pos in range(motif_length): # for the length of the motif 
                    if motif_pos+1 == motif_length: # if reached the end of the motif set final score 
                        score = score / full_motif_length
                        rel_score = (score*full_motif_length)/length_count
                        if length_count > 2:
                            new_row = {'length': length_count,'location': first_aa_location[i_first_aa][i_location],  'score': score, 'rel score': rel_score, 'seq': found_seq}
                            df.loc[len(df)] = new_row
                    elif seq_pointer+1 == len(sequence):
                        isfound = False
                    elif isfound == True:
                        seq_pointer += 1
                        isfound = False
                        for next_pos in range(len(aa_2d_list[motif_pos+1])): #for the length of the next position in the motif 
                            if isfound == False: 
                                if  sequence[seq_pointer] == aa_2d_list[motif_pos+1][next_pos]:                         
                                    score = score + motif_dict[motif_pos+1][aa_2d_list[motif_pos+1][next_pos]]
                                    length_count += 1
                                    isfound = True
                                    found_seq = ''.join([found_seq, aa_2d_list[motif_pos+1][next_pos]])

                    else:
                        score = score  
        motif_df = motif_df.drop(index=0) 
        motif_df = motif_df.reset_index(drop=True)
        
########## writing results to textfiles and printing results ##########

    # makes textfile with incomplete found sequences  
    
    df = df.sort_values(by='location') ############ CAN CHANGE SORTING PREFERENCE FOR TEXTFILE 
    df_string = df.to_string(index = False)
    f = open("found_sequences.txt", "a")
    f.write(df_string)
    f.write('\n')
    f.write('\n')
    f.close()
    # printing found sequences
    
    df = df.sort_values(by='length', ascending=False) ########### CAN CHANGE SORTING PREFERENCE FOR PRINTING
    print(protein_df['access_code'][i_protein])
    print(df.to_string(index=False)) 
    print()
    
    # making dataframe and textfile of results (if whole motif is found) 
    found_pos = []
    found_score = []
    found_motif = []
    for i_found in range(len(df)):
        if df.loc[i_found,'length'] == len(motif_df_full):
            found_pos.append(df.loc[i_found, 'location'])
            found_motif.append(df.loc[i_found, 'seq'])
            found_score.append(df.loc[i_found, 'score'])           

    new_row = {'accession code': protein_df['access_code'][i_protein], 
               'found position': found_pos, 
               'found motif': found_motif,
               'score': found_score}
    results_df.loc[len(results_df)] = new_row
    
results_df_string = results_df.to_string(index = False)
f = open("results.txt", "w")
f.write(results_df_string)
f.close()

# printing results 
print('####### Complete motifs found #######')
print()
print(results_df.to_string(index=False))
print()