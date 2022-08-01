# This script takes in a GZIPPED FASTQ file and collapses 
# the file on unique sequences that are present within the file 

# michelle ragsac (mragsac@eng.ucsd.edu)

##### * ##### * ##### * ##### * ##### * #####

import time
import sys
import gzip

##### * ##### * ##### * ##### * #####

# Read in the desired FLASH-combined FASTQ file
print(f"{time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())} Reading in the input GZIPPED FASTQ file (.fasta/.fa format)")
FILENAME = sys.argv[1]
print(f"{time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())} FILENAME: {FILENAME}")

# Go through each line of the file to extract just the sequences
sequence_tallies = {}
passedHeader     = True
counter          = 0
with gzip.open(FILENAME, 'rb') as f:
    for l in f:
        line = str(l)[2:-1]
        if line[0] == "@": # sequence is after this line
            passedHeader = True
            continue
        if line[0] == "+": # quality scores are after this line
            passedHeader = False
            continue
        
        # Grab the sequence beneath the header line and tally
        # the number of times we've seen a certain sequence in the file
        if passedHeader == True:
            sequence = line.strip()
            if sequence not in sequence_tallies:
                sequence_tallies[sequence] = 0
            sequence_tallies[sequence] += 1
            counter += 1
            
            if counter % 100000 == 0:
                print(f"{time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())} Evaluated {counter} sequences")
            
print(f"{time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())} There are {len(sequence_tallies)} unique sequences present in the file!")
print(f"{time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())} There are {counter} total reads in the file!")

##### * ##### * ##### * ##### * ##### * #####

print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} Exporting output files!")

# Consolidate output file
sequence_tallies_file = ""
for sequence in sequence_tallies:
    sequence_tallies_file += f"{sequence}\t{str(sequence_tallies[sequence])}\n"
    
# Output results
OUTPUT_FILENAME = FILENAME.split("/")[-1].split(".fa")[0] + "_collapsed.txt"
open(OUTPUT_FILENAME,"w").write(sequence_tallies_file)

print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())} DONE ; EXITING")
