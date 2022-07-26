# Author: Michelle Franc Ragsac (mragsac@eng.ucsd.edu) 

# This script searches an INPUT genome for canonical ZicL binding sites that are close to ETS core binding
# sites, then checks for the funcitonal orientation of these sites while scoring the binding affinities.

############################################################################################################
# Import necessary packages 
############################################################################################################

from Bio import SeqIO 
import time
import sys
import re


def status():
    """HELPER METHOD: Prints out the local time"""
    return time.strftime('%Y-%m-%d %H:%M:%S',time.localtime())

############################################################################################################
# Initialize the methods and data structures needed to calculate Ets binding affinities 
############################################################################################################

print(f"{status()} Reading in 8mer binding affinity information for Ets from Uniprobe")
ETS_FILENAME = '../ref/Ets1_8mers.txt'

# Store motif <-> intensity information from the Uniprobe file
ets2intensity = {}
maximum_intensity = 0
for line in open(ETS_FILENAME).read().rstrip().split("\n"): 
	a = line.split("\t")
	intensity = float(a[3])

    # Define intensities for a given 8mer sequence and its reverse-complement
	ets2intensity[a[0]] = intensity 
	ets2intensity[a[1]] = intensity
	
	if (a[0] == "CCGGAAGT") or (a[1] == "CCGGAAGT"):
		maximum_intensity = intensity
    
	
def score_ets_site(Ets_regex_match,sequence):
	"""Given a motif detected within a sequence, return information regarding the affinity score"""

	# Grab the dinucleotide flanking regions to the Ets core site
	Ets_Coordinates = Ets.span()
	left_coordinate = Ets_Coordinates[0] - 2
	if Ets_Coordinates[0] < 2: # grab what sequence we can if the motif is at the end
		left_coordinate = 0
	right_coordinate = Ets_Coordinates[1] + 2
	motif = sequence[left_coordinate:right_coordinate]

	# Score the motif using the 8mer information then return the information
	motif_score = 0
	if motif in ets2intensity:
		motif_score = min(1.0,ets2intensity[motif]/maximum_intensity)
	else:
		motif_score = 0
	return(motif,motif_score,left_coordinate,right_coordinate)
 
	
def determine_total_affinity(locations,affinities):
	"""
    Given the location of Ets sites present, return the total affinity of the sequence 
    and the separation between binding sites. If there are overlaps between sites,
    only consider the sites with the higher affinity.
    
    ***This method uses brute force to calculate the total affinity and can 
    ***be optimized in the future for the calculations...
	"""
    
	if sum(affinities) == 0:
		return([[],0])

	# Filters for valid binding sites whose affinities are greater than 0
	affinity_filter = [True if e > 0 else False for e in affinities]
	locations_f, affinities_f = [], []
	for i in range(len(affinity_filter)):
		if affinity_filter[i] == True:
			locations_f.append(locations[i])
			affinities_f.append(affinities[i])

	# Calculates the distances between sites
	site_distances = [(locations_f[i][0]-locations_f[i-1][1],(i-1,i)) \
		for i in range(1,len(locations_f))]

	# Determine if all of the sites present are far enough away from each other,
	# and if so, add all of the affinities together and return
	number_of_valid_sites = sum([1 for site in site_distances if site[0] > 2])
	if number_of_valid_sites == len(site_distances):
		return([(locations_f[i],affinities_f[i]) for i in range(len(locations_f))],
				sum([e for e in affinities_f]))

	# Grab the highest affinity site if there are just two sites competing
	if len(site_distances) == 1: 
		max_affinity_index = affinities_f.index(max(affinities_f))
		return([(locations_f[max_affinity_index],affinities_f[max_affinity_index])],
			max(affinities_f))

	if len(site_distances) == 2: # if we have three sites present
		# Compare if the outer two sites combined have a higher 
		# affinity than the inner site 
		if number_of_valid_sites == 0: 
			outer_site_affinities = affinities_f[0] + affinities_f[2]
			if outer_site_affinities > affinities_f[1]:
				return([(locations_f[0],affinities_f[0]),(locations_f[2],affinities_f[2])],
					outer_site_affinities)
			else:
				return([(locations_f[1],affinities_f[1])],affinities_f[1])
		
		# Evaluate if the first or second pair of sites are invalid,
		# then figure out which of the pair to include with the site
		# that is not in an invalid region
		else: 
			if site_distances[0][0] <= 2: 
				invalid_pair = site_distances[0][1]
				if affinities_f[invalid_pair[0]] > affinities_f[invalid_pair[1]]:
					return([(locations_f[0],affinities_f[0]),(locations_f[2],affinities_f[2])],
						sum([affinities_f[0],affinities_f[2]]))
				else:
					return([(locations_f[1],affinities_f[1]),(locations_f[2],affinities_f[2])],
						sum([affinities_f[1],affinities_f[2]]))
			else:
				invalid_pair = site_distances[1][1]
				if affinities_f[invalid_pair[0]] > affinities_f[invalid_pair[1]]:
					return([(locations_f[0],affinities_f[0]),(locations_f[1],affinities_f[1])],
						sum([affinities_f[0],affinities_f[1]]))
				else:
					return([(locations_f[0],affinities_f[0]),(locations_f[2],affinities_f[2])],
						sum([affinities_f[0],affinities_f[2]]))

	if len(site_distances) == 3: # if we have four sites present
		# Compare the 1st and 3rd Site, 2nd and 4th Site, and 1st and 4th Site,
		# then return the pair that has the highest affinity
		if number_of_valid_sites == 0:
			affinity_1_3 = affinities_f[0] + affinities_f[2] # [SITE - include] invalid [SITE] invalid [SITE - include] invalid [SITE]
			affinity_2_4 = affinities_f[1] + affinities_f[3] # [SITE] invalid [SITE - include] invalid [SITE] invalid [SITE - include]
			affinity_1_4 = affinities_f[0] + affinities_f[3] # [SITE - include] invalid [SITE] invalid [SITE] invalid [SITE - include]
			if max([affinity_1_3,affinity_2_4,affinity_1_4]) == affinity_1_3:
				return([(locations_f[0],affinities_f[0]),(locations_f[2],affinities_f[2])],affinity_1_3)
			elif max([affinity_1_3,affinity_2_4,affinity_1_4]) == affinity_2_4:
				return([(locations_f[1],affinities_f[1]),(locations_f[3],affinities_f[3])],affinity_2_4)
			else:
				return([(locations_f[0],affinities_f[0]),(locations_f[3],affinities_f[3])],affinity_1_4)
		
		# Evaluate the pair of sites that have an invalid separation between
		# them, select the one with the largest affinity, and return all other sites
		elif number_of_valid_sites == 2:
			if site_distances[0][0] <= 2: # [SITE*] invalid [SITE*] [SITE - include] [SITE - include]
				invalid_pair = site_distances[0][1]
				if affinities_f[invalid_pair[0]] > affinities_f[invalid_pair[1]]:
					return([(locations_f[0],affinities_f[0]),(locations_f[2],affinities_f[2]),(locations_f[3],affinities_f[3])],
						sum([affinities_f[0],affinities_f[2],affinities_f[3]]))
				else:
					return([(locations_f[1],affinities_f[1]),(locations_f[2],affinities_f[2]),(locations_f[3],affinities_f[3])],
						sum([affinities_f[1],affinities_f[2],affinities_f[3]]))
			elif site_distances[1][0] <= 2: # [SITE - include] [SITE*] invalid [SITE*] [SITE - include]
				invalid_pair = site_distances[1][1]
				if affinities_f[invalid_pair[0]] > affinities_f[invalid_pair[1]]:
					return([(locations_f[0],affinities_f[0]),(locations_f[1],affinities_f[1]),(locations_f[3],affinities_f[3])],
						sum([affinities_f[0],affinities_f[1],affinities_f[3]]))
				else:
					return([(locations_f[0],affinities_f[0]),(locations_f[2],affinities_f[2]),(locations_f[3],affinities_f[3])],
						sum([affinities_f[0],affinities_f[2],affinities_f[3]]))
			else: # [SITE - include] [SITE - include] [SITE*] invalid [SITE*]
				invalid_pair = site_distances[2][1]
				if affinities_f[invalid_pair[0]] > affinities_f[invalid_pair[1]]:
					return([(locations_f[0],affinities_f[0]),(locations_f[1],affinities_f[1]),(locations_f[2],affinities_f[2])],
						sum([affinities_f[0],affinities_f[1],affinities_f[2]]))
				else:
					return([(locations_f[0],affinities_f[0]),(locations_f[1],affinities_f[1]),(locations_f[3],affinities_f[3])],
						sum([affinities_f[0],affinities_f[1],affinities_f[3]]))
		
		# Evaluate grouping of invalid sites to look for the sites that
		# have the largest affinity with proper spacing, return all other sites
		elif number_of_valid_sites == 1:
			if site_distances[0][0] <= 2 and site_distances[1][0] <= 2: # [SITE*] invalid [SITE*] invalid [SITE*] [SITE - include] 
				outer_site_affinities = affinities_f[0] + affinities_f[2]
				if outer_site_affinities > affinities_f[1]: # [SITE - include] invalid [SITE] invalid [SITE - include] [SITE - include] 
					return([(locations_f[0],affinities_f[0]),(locations_f[2],affinities_f[2]),(locations_f[3],affinities_f[3])],
						outer_site_affinities+affinities_f[3])
				else: # [SITE] invalid [SITE - include] invalid [SITE] [SITE - include] 
					return([(locations_f[1],affinities_f[1]),(locations_f[3],affinities_f[3])],
						sum([affinities_f[1],affinities_f[3]]))
			elif site_distances[1][0] <= 2 and site_distances[2][0] <= 2: # [SITE - include] [SITE*] invalid [SITE*] invalid [SITE*]
				outer_site_affinities = affinities_f[1] + affinities_f[3]
				if outer_site_affinities > affinities_f[2]: # [SITE - include] [SITE - include] invalid [SITE] invalid [SITE - include]
					return([(locations_f[0],affinities_f[0]),(locations_f[1],affinities_f[1]),(locations_f[3],affinities_f[3])],
						outer_site_affinities+affinities_f[0])
				else: # [SITE - include] [SITE] invalid [SITE - include] invalid [SITE]
					return([(locations_f[1],affinities_f[1]),(locations_f[2],affinities_f[2])],
						sum([affinities_f[1],affinities_f[2]]))
			else:
				affinity_1_3 = affinities_f[0] + affinities_f[2] # [SITE] invalid [SITE - include] [SITE] invalid [SITE - include]
				affinity_2_4 = affinities_f[1] + affinities_f[3] # [SITE - include] invalid [SITE] [SITE - include] invalid [SITE]
				if affinity_1_3 > affinity_2_4:
					return([(locations_f[0],affinities_f[0]),(locations_f[2],affinities_f[2])],
						sum([affinities_f[0],affinities_f[2]]))
				else:
					return([(locations_f[1],affinities_f[1]),(locations_f[3],affinities_f[3])],
						sum([affinities_f[1],affinities_f[3]]))
                
# { seqName : { motif : bindingaffinity , motif2 : ...  } }
name2Ets2intensity = {} 

############################################################################################################
# Initialize the genome informaiton that was taken in as input
############################################################################################################

# Read in the desired genome file
print(f"{status()} Reading in the input genome file (.fasta/.fa format)")
FILENAME = sys.argv[1]
chr2seq = SeqIO.to_dict(SeqIO.parse(open(FILENAME,'r'), "fasta")) 

# Defines flanking window to look for nearby ETS sites to a Zic site
WINDOW = int(sys.argv[2])
print(f"{status()} Input Flanking Window Size around Zic: {WINDOW}")

############################################################################################################
# Define sequence motifs to include within the search
############################################################################################################

# Save sequencing adaptor sequences for future reference
prime5 = "CTGGAGTTCAGACGAGGAGAAACCAGCCT" 
prime3 = "GAAAACCATTCTCCTCTTCTAGAGGATCT"

# Save regex pattern for Zic sites with Ben's annotations
ZicBinding = "CAGCTGTG|CACAGCTG" 		# zic1/2/3 + zic1/2/3 reverse-complement 
ZicBinding += "|CCGCAGT|ACTGCGG" 		# zic7/3/1 + zic7/3/1 reverse-complement
ZicBinding += "|CCGCAGTC|GACTGCGG" 		# zic6 + zic6 reverse-complement
ZicBinding += "|CCCGCTGTG|CACAGCGGG" 	# zic1 + zic1 reverse-complement
ZicBinding += "|CCAGCTGTG|CACAGCTGG" 	# zic3 + zic3 reverse-complement
ZicBinding += "|CCGCTGTG|CACAGCGG" 		# zic2/zicC + zic2/zicC reverse-complement
ZicBinding += "|CCCGCAGTC|GACTGCGGG" 	# zic5 + zic5 reverse-complement

# Save regex pattern for Ets binding site cores
EtsBinding = "GGAA|TTCC" 	# common core motif + reverse-complement
EtsBinding += "|GGAT|ATCC" 	# slightly less common core motif + reverse-complement

# Save all of the reverse-complement sites 
ZicReverse = ["CACAGCTG","CACAGCTG","CACAGCTG","ACTGCGG","GACTGCGG",
	"CACAGCGGG","CACAGCTGG","CACAGCGG","GACTGCGGG"]

############################################################################################################
# Define the output files and their associated headers
############################################################################################################

line_out = "{}\n".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
line_out_orientation = "{}\n".format(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))

line_out += "Location\tSequence\tSequence+Adapters\tTotal_ETS_Affinity\tEts_Locations_with_Affinities\n"
line_out_orientation += "Location\tSequence\tDistance_Between_Zic_and_Ets\tZic_Sequence\tEts_Sequence\tComplete_Ets_Motif\tLocation_of_Motif\tEts_Binding_Affinity\n"

############################################################################################################
# Perform the genome search on the data
############################################################################################################

num_enhancers_found = 0
num_enhancers_orientation_found = set()
for Chr in chr2seq:
	# Compile all of the Zic binding sites possible into regex
	ZicLFinder = re.compile(ZicBinding,re.IGNORECASE)

	# Go through each Zic binding site found and evaluate 
	# flanking window for nearby ETS binding sites
	for ZicL in ZicLFinder.finditer(str(chr2seq[Chr].seq)):
		# Gather location of the Zic binding site, and the flanking regions
		Start = ZicL.start() - WINDOW
		End = ZicL.start() + len(ZicL.group()) + WINDOW
		seq = str(chr2seq[Chr].seq)[Start:End]
		location = Chr + ":" + str(Start) + "-" + str(End)

		# Evaluate just the flanking sites around Zic
		temp = " ".join(seq.split(ZicL.group()))
		if len(temp) == 0:
			continue

		# Check for the orientation of the Zic site 
		# (5'->3' Forward = True; Reverse = True)
		# then gather position of the site within the enhancer
		ZicOrientation = True 
		Zic_Pointing_Position = temp.index(" ") + 1
		if ZicL.group() in ZicReverse:
			ZicOrientation = False
			Zic_Pointing_Position = temp.index(" ")

		# Look for ETS binding sites nearby using *just* the core motif
		if len(re.findall(EtsBinding,temp,re.IGNORECASE)) > 1: 
			EtsFinder = re.compile(EtsBinding,re.IGNORECASE) 
			flanking_sequences = temp.split(" ")
			if len(flanking_sequences) == 1:
				continue

			# Stores which side of the sequences is facing the Zic site to look for proper orientation
			if ZicOrientation: 
				facing_sequence = flanking_sequences[1] # forward direction Zic
			else:
				facing_sequence = flanking_sequences[0] # reverse-complement Zic

			Ets_Locations = []
			Ets_Affinities = []

			left_flank = True
			for sequence_to_search in flanking_sequences: 
				for Ets in EtsFinder.finditer(sequence_to_search):
					# Score the ETS site that we've found
					if location not in name2Ets2intensity.keys():
						name2Ets2intensity[location] = {}
					[motif, motif_score, left_coordinate, right_coordinate] = score_ets_site(Ets,sequence_to_search)
					name2Ets2intensity[location][motif] = motif_score

					# Store the location and affinity of the ETS motif we've found
					Ets_Affinities.append(motif_score)
					if left_flank:
						Ets_Locations.append((left_coordinate,right_coordinate))
					if not left_flank:
						Ets_Locations.append((left_coordinate+31,right_coordinate+31))

					##### * ##### * ##### * ##### * ##### * #####
					# Check the orientation of the ETS site then gather the position of the site in the enhancer
					EtsOrientation = True
					Ets_Pointing_Position = int(Ets.end())
					if Ets.group() in ["TTCC","ATCC"]:
						EtsOrientation = False
						Ets_Pointing_Position = int(Ets.start()) 
						if left_flank != True: 
							Ets_Pointing_Position = int(Ets.start()) + 31
					
					# Determine if the ETS site is present in the sequence facing Zic,
					# and that the ETS and Zic sites are facing each other
					if facing_sequence == sequence_to_search and ZicOrientation != EtsOrientation:
						line_out_orientation += location + "\t" + \
							seq + "\t" + str(abs(Zic_Pointing_Position-Ets_Pointing_Position)) + "\t" + \
							ZicL.group() + "\t" + Ets.group() + "\t" + motif + "\t" + \
							"({},{})\t".format(left_coordinate,right_coordinate) + str(motif_score) + "\n"

						num_enhancers_orientation_found.add(location)

				left_flank = False # next sequence is the right-hand sequence with respect to Zic

			##### * ##### * ##### * ##### * ##### * #####
			# Determine the total affinity of the enhancer sequence with overlaps taken into account
			total_ets_affinity = determine_total_affinity(Ets_Locations,Ets_Affinities)
			total_ets_affinity = (0,0)

			# File with the location of the enhancer, the subsequent sequence,
			# sequence with primers, total ETS affinity, and ETS locations with associated affinities
			line_out += location + "\t" + seq + "\t" + prime5 + seq + prime3 + "\t" + \
				str(total_ets_affinity[1]) + "\t" + str(total_ets_affinity[0]) + "\n"
			
			num_enhancers_found += 1

	print(f"{status()} Done evaluating {Chr}")
	
############################################################################################################
# Export all of the output files
############################################################################################################

PROJECT_NAME = sys.argv[3]
print(f"{status()} Exporting Output Files for genome-search_*_{PROJECT_NAME}.tsv")

open(f"genome-search_total-ets-affinity_{PROJECT_NAME}.tsv","w").write(line_out)
print(f"{status()} Finished exporting genome-search_total-ets-affinity_{PROJECT_NAME}.tsv !")

# open(f"genome-search_adapters_{PROJECT_NAME}.tsv","w").write(line_out_with_adapters)
# print(f"{status()} Finished exporting genome-search_adapters_{PROJECT_NAME}.tsv !")

open(f"genome-search_orientation_{PROJECT_NAME}.tsv","w").write(line_out_orientation)
print(f"{status()} Finished exporting genome-search_orientation_{PROJECT_NAME}.tsv !")
	
############################################################################################################
# Share final numbers on the number of the enhancers found within the input genome
############################################################################################################

print(f"{status()} Enhancers Found: {num_enhancers_found}")
print(f"{status()} Enhancers with Correct Orientation Found: {len(num_enhancers_orientation_found)}")
print(f"{status()} DONE")