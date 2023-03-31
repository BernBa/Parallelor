# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 10:40:03 2022

@author: bauerber
"""

# =============================================================================
# 
# IMPORTANT THINGS TO CHECK
# lower case all the glosses!
# make sure the glosses are sorted according to ID
#
# =============================================================================

import pandas as pd
from unidecode import unidecode # To convert the Greek characters into Latin
import linecache # To read the lines of the .txt-file we import linecache
from textwrap import TextWrapper
from termcolor import colored

# To find consecutive stars in the file
def count_consec_stars(starlist):
    consec = [0]
    for x, y in zip(starlist, starlist[1:]):
        if x == y - 1:
            consec[-1] += 1
            if x not in consecutive_stars:
                consecutive_stars.append(x)
            consecutive_stars.append(y)                   
        else:
            consec.append(1)
    return (consecutive_stars)

# To turn a list into a string
def listToString(s):
    parallel_gloss = ""
    for element in s:
        parallel_gloss += element
    return parallel_gloss

# To delete duplicates in lists
def delete_duplicates(x):
    return list(dict.fromkeys(x))

# The number of manuscripts
# number_of_manuscripts = int(input('How many manuscripts do you have? '))
number_of_manuscripts = 2

# Number of the loop/run
run = 0
new_ms_list = []

# Start the the loop
again = 'y'
while again != 'n':
    # print('This is run number: ', run)
    # Get a list of all the ms abbreviations
    if run == 0:
        ms_list = []
        for manuscript in range(number_of_manuscripts):
            ms = input('Please input the abbreviation of manuscript ' + str(manuscript) + ': ')
            ms_list.append(ms)
    else:
        ms_list = new_ms_list
    
    
    # Establish a dictionary containing "Abbreviation of the manuscript" and a dataframe for each manuscript
    d = {}
    d2 = {}

    for ms in ms_list:
        inputfile = ms + '.csv'
        d[ms] = pd.read_csv(inputfile)
        d2[ms] = d[ms]
        d2[ms] = d2[ms].sort_values('ID')
        d[ms]['Gloss'] = d[ms]['Gloss'].replace({'á': 'a','é': 'e','í': 'i','ó': 'o','ú': 'u', r'[^\w\s]+': '', ' ': ''}, regex=True) #maybe include ' ': 'w'
        d[ms]['Gloss'] = d[ms]['Gloss'].apply(unidecode)
        d[ms] = d[ms].sort_values('ID')
        
    
    # Establishing and printing the output in the FASTA format
    for key in d:
        ms = key
        text = d[key]['Gloss'].tolist()
        text = 'w'.join([str(a) for a in text]) # glosses are divided by 'W' because it does not occur otherwise
        tw = TextWrapper()
        tw.width = 80
        text = ('\n'.join(tw.wrap(text)))
        print('>' + ms + '\n' + text)
    
    print('\n', 'Please copy and paste this text into the text into the correspoding field on the website (https://www.genome.jp/tools-bin/clustalw): ')
    print('\n', '********************************************************************'
          '\n', 'Use the following settings:', '\n'
          # 'K-tuple 1, Window size 10, Gap 3', '\n', 
          # 'Diagonals 5, Scoring ABSOLUTE', '\n',
          'Pairwise Alignment: SLOW/ACCURATE', '\n',
          'Enter your sequences [...]: PROTEIN', '\n',
          'Gap Open 1, Gap Extension 0.05', '\n',
          'Gap Open 1, Gap Extension 0.05', '\n'
          'Weight and Hydrophilic NO')
    
    txt_name = input("When you saved the the result as .txt in the home folder, please input the filename (without the .txt): ")    
    txt_name = txt_name + '.txt'
    
    
    # Input manuscript text
    fasta_lines = []
    
    with open(txt_name, 'r') as fp:
        x = fp.readlines()
        total_lines = len(x)
           
    for ln in range(1, total_lines):
        fasta_lines.append(ln)
    
    
    # Lists for the indices of the first manuscript and the stars
    space_between = (number_of_manuscripts + 2) # the extra line is for the starline
    ms1_lines = fasta_lines[::space_between]
    star_lines = fasta_lines[(number_of_manuscripts - 1)::space_between]
    
    
    # Establish the text of the first manuscript
    ms1_text = []
    # n = 0
    
    for pos in range(1, total_lines):
        if pos in ms1_lines:
            # if n > 0:
            #     pos = pos+1
            ms1_text.append(linecache.getline(txt_name, pos)[16:]) 
            # n =+ 1
            
    ms1_text = listToString(ms1_text).replace("\n", "")
    ms1_text = ms1_text.lower()
    
    
    # Establish the text of the star line
    star_text = []
    n = 0
    
    for pos in range(1, total_lines):
        if pos in star_lines:
            if n > 0:
                pos = pos+n
            star_text.append(linecache.getline(txt_name, pos)[16:]) 
        n =+ 1
    
    star_text = listToString(star_text).replace("\n", "")
    
    
    # Empty lists for getting the star characters in the first manuscript
    starlist = []
    consecutive_stars = []
    consecutive_chars = []
    
    # Establish the list of consecutive stars in the stars text
    # And insert 'w' when consecutive numbers != plus 1 (i.e. 9, 10, 14,...)
    starlist = [pos for pos, char in enumerate(star_text) if char == '*']
    consec_list = count_consec_stars(starlist)
    sliced_consec_list = []
    for s1, s2 in zip(consec_list, consec_list[1:]):
        if s2 == s1 + 1:
            if s1 not in sliced_consec_list:
                sliced_consec_list.append(s1)
            sliced_consec_list.append(s2) 
        else:
            sliced_consec_list.append('w')
    
    # The sliced_consec_list is sliced into sublists at the occurring w
    size = len(sliced_consec_list)
    idx_list = [idx + 1 for idx, val in
                enumerate(sliced_consec_list) if val == ('w' or ' ')]
    
    
    sliced_list = [sliced_consec_list[i: j] for i, j in 
                   zip([0] + idx_list, idx_list + ([size] if idx_list[-1] != size else []))]
    
    # =============================================================================
    # Establish the list of the matching characters in the first manuscript and 
    # turn it into a string (consecutive_chars_string) and list (consecutive_chars_list)
    # =============================================================================
    
    for char in sliced_list:
        for ch in char:
            if ch == "w":
                consecutive_chars.append(' ')
            else:
                consecutive_chars.append(ms1_text[ch])
    
    consecutive_chars_string = listToString(consecutive_chars).replace('w', ' ')
    consecutive_chars_list = consecutive_chars_string.split(' ')
    consecutive_chars_list = list(filter((' ').__ne__, consecutive_chars_list))
    consecutive_chars_list = list(filter(None, consecutive_chars_list))
    consecutive_chars_list = [x for x in consecutive_chars_list if len(x)>2] # Removing those instance in which the 'w' caused the second *
    
    
    # =============================================================================
    # This section maps the consecutive characters onto the glosses
    # =============================================================================
    
    n = - number_of_manuscripts
    prev_ID = 0
    prev_gloss = 0
    parallel_list = []
    
    print(consecutive_chars_list)
    
    for consecutive_chars in consecutive_chars_list:
        p = 0
        print('\n', '**************************************************************')
        print(colored(consecutive_chars, 'cyan'))
        for key in d:
            ID = d[key].index[d[key]['Gloss'].str.contains(consecutive_chars) == True].tolist()
            ID.sort()
            if key == ms_list[p]:
                if parallel_list != []:
                    if ID[0] >= int((parallel_list[n])[1]):
                        the_gloss = ID[0]
                    else:
                        higher_ID = [value for value in ID if value >= int((parallel_list[n])[1])]
                        higher_ID.sort()
                        ID = higher_ID
            else:
                the_gloss = ID[0]
                
            it = 1
            for i in ID:
                print('\n', key, colored(i, 'yellow'), d2[key]['Gloss'].loc[i], d2[key]['ID'].loc[i])
                if it == 12:
                    break
                it += 1
            print('__________________________________________________________')    
            if p == (number_of_manuscripts - 1):
                for key in d:
                    gloss_number = input('Please input the number of the gloss of ' + key + ' (if none fits press enter): ')
                    if gloss_number == '':
                        n -= number_of_manuscripts
                        break
                    else:
                        parallel_list.append((key, gloss_number))
            p += 1
            n += 1
    
    
    # =============================================================================
    # The output is an ordered list of the found parallel glosses
    # =============================================================================
    
    parallel_list = delete_duplicates(parallel_list)
    outputlist = []
    delete_index_list = []
    for parallel in parallel_list:
        df = d[parallel[0]]
        outputlist.append(df['ID'].loc[int(parallel[1])])
        output_item = df['ID'].loc[int(parallel[1])]
        delete_this = df.index[df['ID'] == output_item].tolist()
        delete_index_list.append(delete_this[0])
    
    
    # print(delete_index_list)
    
    count = 0
    count_odd = 0
    count_even = 0
    for item in range(0, len(delete_index_list)):
        delete_index = int(delete_index_list[count])    
        if item % 2:
            ms = ms_list[1]
            d2[ms] = d2[ms].drop(d2[ms].index[(delete_index - count_even)])
            count_even += 1
        else:
            ms = ms_list[0]
            d2[ms] = d2[ms].drop(d2[ms].index[(delete_index - count_odd)])
            count_odd += 1
        count += 1
    
    run += 1
    new_ms_list = []
    for ms in ms_list:
        d2[ms].to_csv((ms + str(run) + '.csv'), encoding='utf-8')
        new_ms_list.append(ms + str(run))
        
        
    with open('parallels.txt', 'a+') as f:
        for gloss in outputlist:
            f.write(f"{gloss}\n")
        f.write(f"Run {run}\n")    
    again = input('Do you want to try another run (y/n)? ')    
      
    if again == 'n':
        break
        # delete_index = int(delete_index_list[0])
        # print(delete_index)
        # df = df.drop([df.index[delete_index]])
        # d[parallel[0]] = df
    
    # outputlist.sort()
            
    # for ms in ms_list:
    #     d[ms].to_csv((ms + '1.csv'), encoding='utf-8')
        
    # print(outputlist)