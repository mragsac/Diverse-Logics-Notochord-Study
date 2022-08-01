# This script takes in FASTQ files that have been collapsed
# on unique sequences, then conducts preliminary filtering
# based on an input dictionary file

# michelle ragsac (mragsac@eng.ucsd.edu)

##### * ##### * ##### * ##### * ##### 

from datetime import datetime
import argparse
import time

import pandas as pd

##### * ##### * ##### * ##### * ##### 
# DEFINE ARGUMENTS
##### * ##### * ##### * ##### * ##### 

parser = argparse.ArgumentParser()

# Required arguments for the script
parser.add_argument("infile", help="location of collapsed fastq file to perform filters on", type=str)
parser.add_argument("dictionary", help="location of enhancer-barcode tag dictionary to reference", type=str)

# Optional arguments for the script
parser.add_argument("-v", "--verbosity", help="increase output verbosity", action="store_true")
parser.add_argument("--linker", help="define sequence separating the barcode tag (linker) to search for",
                    default="ACTAGTGAGTTTGTC")
parser.add_argument("--location", help="define expected location of the linker sequence",
                    default=30, type=int)
parser.add_argument("--window", help="define bp window to look for linker sequence",
                    default=5, type=int)
parser.add_argument("--minimum", help="define the minimum read count for a sequence",
                    default=25, type=int)
parser.add_argument("--barcode", help="define the minimum number of barcodes for an enhancer sequence",
                    default=2, type=int)

args = parser.parse_args()

##### * ##### * ##### * ##### * ##### 
# HELPER METHODS
##### * ##### * ##### * ##### * ##### 

def current_time():
    """Helper Method:
    Prints out the current time for verbose statements
    """
    return(time.strftime('%Y-%m-%d %H:%M:%S',time.localtime()))

def hamming_distance(seq1,seq2):
    """Helper Method:
    Determines the Hamming Distance (# mismatches) between two sequences
    """
    return(sum([1 for i,j in zip(seq1,seq2) if i != j]))

def find_linker(sequence,hamming_threshold=2):
    """Given a sequence, split the sequence according to the
    linker sequence that is present (with mismatches allowed according
    to the hamming threshold)
    """
    # Attempt looking for an exact match of the linker, 
    # then splitting the sequence according to this exact match
    exact_split = sequence.split(linker_sequence)
    if len(exact_split) == 2:
        return(exact_split)
    
    # If we can't find an exact match, search for the linker sequence
    # and split with hamming distance of less than some threshold
    potential_linker_sequences = [
        sequence[i:i+len(linker_sequence)] for i in range(
            linker_location-linker_window,linker_location+linker_window)]
    
    for subseq in potential_linker_sequences:
        if hamming_distance(subseq,linker_sequence) <= hamming_threshold:
            approximate_split = sequence.split(subseq)
            return(approximate_split)
    
    # If we can't find the linker sequence, return the original sequence
    return([sequence])

##### * ##### * ##### * ##### * ##### 

START_TIME = datetime.now()

print(f"\n|| Running '{__file__}' ...")
if args.verbosity:
    print("")
    print(f"{current_time()} =========================================================================================")
    print(f"{current_time()} --- Starting ANALYSIS STEP 1: Preliminary Filtering ---")
    print(f"{current_time()} =========================================================================================")
    print(f"{current_time()} STEP 1: Reading in the input COLLAPSED FASTQ file (.tsv) ...", end=" ")

# Read in the collapsed FASTQ file
FILENAME = args.infile
df = pd.read_csv(FILENAME, sep='\t', header=None, names=['Raw_Sequence','Raw_Counts'])

if args.verbosity:
    print("Done!")
    print(f"{current_time()} FILENAME: {FILENAME}")

num_sequences = df.shape[0]
num_rawcounts = df['Raw_Counts'].sum()

if args.verbosity:
    print(f"{current_time()} =========================================================================================")
    print(f"{current_time()} --- Starting File Information ---")
    print(f"{current_time()} =========================================================================================")
    print(f"{current_time()} Number of Unique Read Sequences:\t{num_sequences}")
    print(f"{current_time()} Number of Total Reads:\t\t{num_rawcounts}")
    print(f"{current_time()} =========================================================================================")

##### * ##### * ##### * ##### * #####

if args.verbosity:
    print(f"{current_time()} STEP 2: Reading in the input DICTIONARY file (.tsv) ...", end=" ")

# Import the dictionary
DICTIONARY_FILENAME = args.dictionary

if args.verbosity:
    print("Done!")
    print(f"{current_time()} FILENAME: {DICTIONARY_FILENAME}")
    print(f"{current_time()} STEP 2: Generating barcode2enhancer DICTIONARY ...", end=" ")

# Consolidate the barcodes on the line and add them to the barcode2enhancer dictionary
barcode2enhancer = {}
isHeader = True
with open(DICTIONARY_FILENAME,'r') as file:
    for line in file:
        if isHeader == True:
            isHeader = False
            continue
        enhancer,barcodes = line.strip().split()
        barcodes = barcodes.split(",")
        
        for bc in barcodes:
            if bc not in barcode2enhancer:
                barcode2enhancer[bc] = {}
            barcode2enhancer[bc][enhancer] = 1

if args.verbosity:
    print("Done!")
    print(f"{current_time()} =========================================================================================")
    print(f"{current_time()} --- Dictionary Statistics ---")
    print(f"{current_time()} =========================================================================================")
    print(f"{current_time()} Number of Enhancers:\t\t{len(set([enh for bc in barcode2enhancer for enh in barcode2enhancer[bc]]))}")
    print(f"{current_time()} Number of Barcodes:\t\t\t{len(barcode2enhancer)}")
    print(f"{current_time()} Number of Unique Barcodes:\t\t{sum([1 for bc in barcode2enhancer if len(barcode2enhancer[bc]) == 1])}")
    print(f"{current_time()} Number of Mult-Mat Barcodes:\t{sum([1 for bc in barcode2enhancer if len(barcode2enhancer[bc]) > 1])}")
    print(f"{current_time()} =========================================================================================")
            
##### * ##### * ##### * ##### * #####

if args.verbosity:
    print(f"{current_time()} STEP 3: Identifying and Filtering Out Sequences with N's Present ...", end=" ")

# Filter out reads with N's present 
with pd.option_context('mode.chained_assignment', None):
    df['hasN'] = df.apply(lambda r: True if 'N' in r['Raw_Sequence'] else False, axis=1) 
    df_filtN = df.loc[df['hasN'] == False]

if args.verbosity:
    print("Done!")
    print(f"{current_time()} =========================================================================================")

##### * ##### * ##### * ##### * #####

# Define linker information
linker_sequence = args.linker
linker_location = args.location
linker_window = args.window

if args.verbosity:
    print(f"{current_time()} --- Linker Information ---")
    print(f"{current_time()} Linker Sequence:\t\t\t{linker_sequence}")
    print(f"{current_time()} Linker Location:\t\t\t{linker_location}th bp")
    print(f"{current_time()} Linker Flexible Window:\t\t{linker_window} bp")
    print(f"{current_time()} =========================================================================================")
    print(f"{current_time()} STEP 4: Identifying and Filtering Out Sequences missing the LINKER ...", end=" ")

# Filter out reads missing a linker
with pd.option_context('mode.chained_assignment', None):
    df_filtN['Split_Sequence'] = df_filtN.apply(lambda r: find_linker(r['Raw_Sequence']), axis=1)
    df_filtN['Barcode_Length'] = df_filtN.apply(lambda r: len(r['Split_Sequence'][0]), axis=1) # assume barcode before first split
    df_filtN['Number_Of_Splits'] = df_filtN.apply(lambda r: len(r['Split_Sequence']), axis=1)

    df_filtN_filtL = df_filtN.loc[df_filtN['Number_Of_Splits'] == 2]
    df_filtN_filtL_filtB = df_filtN_filtL.loc[(df_filtN_filtL['Barcode_Length'] == 29) | (df_filtN_filtL['Barcode_Length'] == 30)]

if args.verbosity:
    print("Done!")
    print(f"{current_time()} =========================================================================================")

##### * ##### * ##### * ##### * #####

if args.verbosity:
    print(f"{current_time()} STEP 5: Determining Total Number of Read Counts per BARCODE ...", end=" ")

# Gather the number of total barcode counts 
with pd.option_context('mode.chained_assignment', None):
    df_filtN_filtL_filtB['Barcode'] = df_filtN_filtL_filtB.apply(lambda r: r['Split_Sequence'][0], axis=1)
    df_bcToTotalCounts = df_filtN_filtL_filtB.groupby('Barcode').apply(lambda r: r['Raw_Counts'].sum()).reset_index()
    df_bcToTotalCounts.columns = ['Barcode','Total_Raw_Counts']

if args.verbosity:
    print("Done!")
    print(f"{current_time()} =========================================================================================")

##### * ##### * ##### * ##### * #####

if args.verbosity:
    print(f"{current_time()} STEP 6: Identifying and Filtering Out Sequences with BARCODE-ENHANCER Associations -NOT- in the DICTIONARY ...", end=" ")

# Determine the enhancer associated with each barcode
with pd.option_context('mode.chained_assignment', None):
    df_bcToTotalCounts['Enhancer'] = df_bcToTotalCounts.apply(
        lambda r: barcode2enhancer[r['Barcode']] if r['Barcode'] in barcode2enhancer else 'Not in Dictionary', axis=1)
    df_bcToTotalCounts_filtD = df_bcToTotalCounts.loc[df_bcToTotalCounts['Enhancer'] != 'Not in Dictionary']

if args.verbosity:
    print("Done!")
    print(f"{current_time()} STEP 6: Identifying and Filtering Out Sequences with Multiple-Matched BARCODE-ENHANCER Associations ...", end=" ")

with pd.option_context('mode.chained_assignment', None):
    df_bcToTotalCounts_filtD['Number_of_Associated_Enhancers'] = df_bcToTotalCounts_filtD.apply(lambda r: len(r['Enhancer']), axis=1)
    df_bcToTotalCounts_filtD_filtE = df_bcToTotalCounts_filtD.loc[df_bcToTotalCounts_filtD['Number_of_Associated_Enhancers'] < 2]
    df_bcToTotalCounts_filtD_filtE['Enhancer'] = df_bcToTotalCounts_filtD_filtE.apply(lambda r: list(r['Enhancer'])[0], axis=1)

if args.verbosity:
    print("Done!")
    print(f"{current_time()} =========================================================================================")

##### * ##### * ##### * ##### * #####

if args.verbosity:
    print(f"{current_time()} STEP 7: Identifying and Filtering Out ENHANCERS <2 BARCODES ...", end=" ")

# Filter out enhancers that do not have at least 2 barcodes
enh_num_bc = pd.DataFrame(df_bcToTotalCounts_filtD_filtE.groupby('Enhancer')['Barcode'].nunique())
enh_num_bc = enh_num_bc.loc[enh_num_bc['Barcode'] >= args.barcode]
df_bcToTotalCounts_filtD_filtE_filt2bc = df_bcToTotalCounts_filtD_filtE.loc[df_bcToTotalCounts_filtD_filtE['Enhancer'].isin(enh_num_bc.index)]

if args.verbosity:
    print("Done!")
    print(f"{current_time()} =========================================================================================")

##### * ##### * ##### * ##### * #####

if args.verbosity:
    print(f"{current_time()} STEP 8: Identifying and Filtering Out BARCODE-ENHANCER associations with <25 Total Raw Reads ...", end=" ")

# Filter out sequences that do not make the minimum read cutoff 

cutoff_read_count = args.minimum

df_bcToTotalCounts_filtD_filtE_filt2bc_filtC = df_bcToTotalCounts_filtD_filtE_filt2bc.loc[
    df_bcToTotalCounts_filtD_filtE_filt2bc['Total_Raw_Counts'] >= cutoff_read_count]

if args.verbosity:
    print("Done!")
    print(f"{current_time()} =========================================================================================")

##### * ##### * ##### * ##### * #####


if args.verbosity:
    print(f"{current_time()} STEP 9: Exporting output files ...", end=" ")

# Export the dataframe with the processed information
OUTPUT_FILENAME = FILENAME.split("/")[-1].split(".txt")[0] + "_filtered-collapsed.txt"
df_bcToTotalCounts_filtD_filtE_filt2bc_filtC.to_csv(OUTPUT_FILENAME,sep="\t")

if args.verbosity:
    print("Done!")
    print("\n")
    
##### * ##### * ##### * ##### * #####

END_TIME = datetime.now()

# Print all of the summary information for the filtering
print(f"|| =============================================================================")
print(f"||                 --- Preliminary Filtering Summary Counts ---")
print(f"|| =============================================================================")

num_sequences = df.shape[0]
num_rawcounts = df['Raw_Counts'].sum()
print(f"||                      *** STARTING NUMBER OF SEQUENCES ***")
print(f"|| -----------------------------------------------------------------------------")
print(f"|| Number of Unique Read Sequences:\t{num_sequences}")
print(f"|| Number of Total Reads:\t\t{num_rawcounts}")

num_sequences = df_filtN.shape[0]
num_rawcounts = df_filtN['Raw_Counts'].sum()
print(f"|| -----------------------------------------------------------------------------")
print(f"||                *** FILTER: Remove Sequences with N's Present ***")
print(f"|| -----------------------------------------------------------------------------")
print(f"|| Number of Unique Read Sequences:\t{num_sequences}")
print(f"|| Number of Total Reads:\t\t{num_rawcounts}")

num_sequences = df_filtN_filtL_filtB.shape[0]
num_rawcounts = df_filtN_filtL_filtB['Raw_Counts'].sum()
print(f"|| -----------------------------------------------------------------------------")
print(f"||               *** FILTER: Remove Sequences Missing the LINKER ***")
print(f"|| -----------------------------------------------------------------------------")
print(f"|| Number of Unique Read Sequences:\t{num_sequences}")
print(f"|| Number of Total Reads:\t\t{num_rawcounts}")

num_sequences = df_bcToTotalCounts_filtD.shape[0]
num_rawcounts = df_bcToTotalCounts_filtD['Total_Raw_Counts'].sum()
print(f"|| -----------------------------------------------------------------------------")
print(f"||  *** FILTER: Remove BARCODE-ENHANCER Associations -NOT- in the DICTIONARY ***")
print(f"|| -----------------------------------------------------------------------------")
print(f"|| Number of Unique Read Sequences:\t{num_sequences}")
print(f"|| Number of Total Reads:\t\t{num_rawcounts}")

num_sequences = df_bcToTotalCounts_filtD_filtE.shape[0]
num_rawcounts = df_bcToTotalCounts_filtD_filtE['Total_Raw_Counts'].sum()
print(f"|| -----------------------------------------------------------------------------")
print(f"||     *** FILTER: Remove Multiple-Matched BARCODE-ENHANCER Associations ***")
print(f"|| -----------------------------------------------------------------------------")
print(f"|| Number of Unique Read Sequences:\t{num_sequences}")
print(f"|| Number of Total Reads:\t\t{num_rawcounts}")

num_sequences = df_bcToTotalCounts_filtD_filtE_filt2bc.shape[0]
num_rawcounts = df_bcToTotalCounts_filtD_filtE_filt2bc['Total_Raw_Counts'].sum()
print(f"|| -----------------------------------------------------------------------------")
print(f"||              *** FILTER: Remove ENHANCERS with <2 BARCODES ***")
print(f"|| -----------------------------------------------------------------------------")
print(f"|| Number of Unique Read Sequences:\t{num_sequences}")
print(f"|| Number of Total Reads:\t\t{num_rawcounts}")

num_sequences = df_bcToTotalCounts_filtD_filtE_filt2bc_filtC.shape[0]
num_rawcounts = df_bcToTotalCounts_filtD_filtE_filt2bc_filtC['Total_Raw_Counts'].sum()
print(f"|| -----------------------------------------------------------------------------")
print(f"|| *** FILTER: Remove BARCODE-ENHANCER Associations with <25 Total Raw Reads ***")
print(f"|| -----------------------------------------------------------------------------")
print(f"|| Number of Unique Read Sequences:\t{num_sequences}")
print(f"|| Number of Total Reads:\t\t{num_rawcounts}")

print(f"|| =============================================================================")

runtime = END_TIME - START_TIME
print(f"|| Total Runtime: {runtime}\n")

##### * ##### * ##### * ##### * #####
