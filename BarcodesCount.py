#!/bin/bash
import re
import os.path
import json
import subprocess
from subprocess import call

# Dependencies: pear, bowtie2, bartender

# Define directories
fastqDir = "Demultiplexed/" # direstory with fastq files
alnDir = "Aligned/" # directory for aligned reads in sam format and extracted barcodes
BarDir = "Barcodes/" # directory with updated barcodes, counts, aggregate counts and updated counter_star
TrajDir = "Trajectories/" # directory with timecourse data

# Define library prefixes
Fitness_i7 = [['D710', 'D711', 'D712', 'A701', 'A702', 'A703', 'A704', 'A705'],['A701', 'A702', 'A703', 'A704', 'A705'],['D710', 'D711']]
Fitness_i5 = [['A506', 'D506', 'A505', 'A507', 'A508'],['D501'],['A504', 'D501']]

# txt file that corresponds fastq prefix from demultiplexing to meaningful library ID and BC1 (low complexity) sequences expected to be present in the library
# formatted as follows: 3 columns tab-separated and no heading. 1st column is the library prefix, 2nd column library identifier and 3rd column a list of all BC1 in the library
LibraryKey = 'FitnessLibraryIndex'

#software calls
bt2 = 'bowtie2'
pear = 'pear'

# the following dictionary corresponds initial pool to assays. It is specific to the data and needs modification
Pools = {'1432':['c1','c2','c7','c8'],'E4':['c4','c5','c10','c11']}

# calculates hamming distance
# can be used instead of levenshtein distance calculator
def hamdist(str1, str2):
 #Count the # of differences between equal length strings str1 and str2

        diffs = 0
        for ch1, ch2 in zip(str1, str2):
                if ch1 != ch2:
                        diffs += 1
        return diffs

# calculates levenshtein distance
def LevenshteinDistance(s1, s2):
    if len(s1) < len(s2):
        return LevenshteinDistance(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]

# takes information from the library key txt file and puts it in json files
def json_files(input, output, output2):
	F1 = open(input + '.txt','r')
	Lib_names = {}
	Lib_BC1 = {}
	for line in F1:
		L1 = line.split('\t')
		Lib_names[L1[0]] = L1[1]
		Lib_BC1[L1[1]] = []
		bc = L1[2].strip('\"\n').split(',')
		print (bc)
		for i in bc:
			Lib_BC1[L1[1]].append(i)
	with open(output + '.json', 'w') as file2:
		json.dump(Lib_names, file2)
	with open(output2 + '.json', 'w') as file2:
		json.dump(Lib_BC1, file2)

# the following code has been broken down in functions that are called at the end.
# the functions process the data as follows
# 'align' merges paired-end reads with pear, then aligns to barcode locus with bowtie2
# 'barcodes' extracts barcodes and quality scores from the resulting sam files. It stores sequence ID, UMI and sequences and quality scores for BC1 and BC2 in json format
# 'bartender' collects the BC2 (high complexity) from all libraries in a single file. Then it calls the bartender software in that file and BC2 are clustered
# 'updateReadsV2' 1. uses known lists of BC1 (Low complexity) to cluster BC1 using levenshtein (or hamming) distance 2. the bartender output to update the BC2 barcodes 3. updates libaray names based on a key
# 'counts' goes into the updated read files and outputs a file with counts
# 'aggregateCounts' calculates aggregate counts for assays where the same pools were used
# The 'aggregateCounts' output is plotted as histograms and based on that I decide on a filtering step (the minimum number of aggregate counts in order for a barcode to be considered as lineage representative)
# Run 'AggregateCountsHist.R', check the pdf file
# 'filter' uses the threshold in aggregate counts and outputs filtered versions of all libraries
# 'trajectories' assembles all barcode filtered counts per assay and outputs a single file (time trajectories with just raw counts): these are not normalized counts, use as input in fitness estimation algorithm
# 'norm_trajectories' are also the timecourse trajectories but normalized based on the total counts. I plot these to decide on which timepoints to use for fitness estimation
# (fitness estimation algorithm from Venkataram et al Cell 2016)


# merge paired-end reads with pear, then align to barcode locus with bowtie2
def align(dirIN, dirOUT, index_pair):
	# merge reads with pear
	subprocess.call(pear + ' -f ' + dirIN + '/Read_1_' + index_pair + '-read-1.fastq.gz -r ' + dirIN + '/Read_4_' + index_pair + '-read-1.fastq.gz -o ' + dirOUT + '/' + index_pair + ' -j 8 -n 150', shell=True) # -j threads, Merge F and R Fastq files when overlapping and adjust Qscore
	# aligns with bowtie2
	subprocess.call(bt2 + " -U " + dirOUT + '/' + index_pair + ".assembled.fastq -p 8 --reorder --no-hd --n-ceil L,0,0.50 -L 12 --np 0 -x /share/ceph/gil213group/dia218/DiploidsR1/BCseq/bowtie_ref/BC_ref -S " + dirOUT + '/' + index_pair + ".sam --norc", shell=True)

# extracts barcodes and quality scores from the sam files. It stores sequence ID, UMI and sequences and quality scores for BC1 and BC2 in json file
def barcodes(dir, index_pair):
	D = {}
	infile = dir + '/' + index_pair
	handle_sam =  open(infile + '.sam','r+') #Reads the extracted parts of the SAM File and breaks down cig_str to find start & end point of the barcode
	counter_star = 0
	count_diff_size = 0
	count_26 = 0
	while True:
		#Split the SAM file, output of Bowtie2, to break the cigar string and extract the BC
		#	The output here is the seqID BC2 qualBC2  BC1 qualBC1
		read = handle_sam.readline()   #Splits SAM file line by line
		if not read:
			break

		try:

			read_list = read.rstrip().split('\t')
			if len(read) != 1 :

				seq_id,flag,start,cigar_string,sequence,quality = read_list[0],read_list[1],read_list[3],read_list[5],read_list[9],read_list[10]

			pos = 0  # Defines our position in the read
			start_aln = int(start) - 1 # Defines the starting point of the alignment. It should be 1

			bar_start_f = 40  #Defines where the first BC start (this the BC2 or high complexity)
			bar_start_f = (bar_start_f - start_aln)  #Defines where the first BC start (this the BC2 or high complexity) if the alignement does not start at 1
			bar_end_f = bar_start_f + 26  # Defines where the first BC end (based on the size of the BC 26 nucleotides)

			bar_start_r = 155  #Defines where the second BC start (this the BC1 or low complexity), named because this BC is at the 3' end of the read so mostly on reverse reads
			bar_start_r = (bar_start_r - start_aln)  #Defines where the second BC start (this the BC1 or low complexity) if the alignement does not start at 1
			bar_end_r = bar_start_r + 26  # Defines where the second BC end (based on the size of the BC 26 nucleotides)


			if cigar_string == '*':  #If cig_str = * , it adds + 1 to counter_* and continues reading the rest of the lines
				counter_star += 1
				continue

			if cigar_string != '250M':
				cigar_op = re.compile("([M,I,D])").split(cigar_string)
				#print (cigar_op)
				for i in range(0, len(cigar_op) -1, 2):

					if pos < bar_start_f:
						#if we are before the BC

						if cigar_op[i+1] == 'M':
							pos = pos + int(cigar_op[i])  #Counts the number of ' M  ' and moves current position positively accordingly

						elif cigar_op[i+1] == 'I':
							bar_start_f = bar_start_f + int(cigar_op[i])  #Counts number of ' I ' and skips current position positively accordingly
							bar_end_f = bar_start_f + 26

							bar_start_r = bar_start_r + int(cigar_op[i])  #Counts number of ' I ' and skips current position  accordingly
							bar_end_r = bar_start_r + 26

							pos = pos + int(cigar_op[i])  #moves current position positively

						elif cigar_op[i+1] == 'D':  #Counts number of ' D ' and moves current position accordingly
							bar_start_f = bar_start_f - int(cigar_op[i])
							bar_end_f = bar_start_f + 26

							bar_start_r = bar_start_r - int(cigar_op[i])  #Counts number of ' D ' and skips current position  accordingly
							bar_end_r = bar_start_r + 26

							pos = pos - int(cigar_op[i])  #moves current position positively

					elif pos >= bar_start_f and pos < bar_end_f:
						#if we are within the first BC do not change the start of the first BC but the end and the second BC (named bar_r)

						if cigar_op[i+1] == 'I':  #Counts number of ' I ' within the forward BC and moves end point and reverse BC positively accordingly
							bar_end_f = (bar_start_f + 26) + int(cigar_op[i])
							bar_start_r = bar_start_r + int(cigar_op[i])
							bar_end_r = bar_start_r + 26

							pos = pos + int(cigar_op[i])  #moves current position positively

						elif  cigar_op[i+1] == 'D':  #Counts number of ' D ' between position forward BC and moves end point and reverse BC negatively accordingly
							bar_end_f = (bar_start_f + 26) - int(cigar_op[i])
							bar_start_r = bar_start_r - int(cigar_op[i])  #Counts number of ' D ' and skips current position  accordingly
							bar_end_r = bar_start_r + 26

							pos = pos - int(cigar_op[i])  #moves current position positively

						elif cigar_op[i+1] == 'M':
							pos_= pos + int(cigar_op[i])  #Counts the number of ' M  ' and moves current position positively accordingly

					elif pos > bar_end_f and pos <= bar_start_r:
						#if we are after the first BC change the point of the second BC

						if cigar_op[i+1] == 'I':  #Counts number of ' I ' within the forward BC and moves end point and reverse BC positively accordingly
							bar_start_r = bar_start_r + int(cigar_op[i])
							bar_end_r = bar_start_r + 26

						elif  cigar_op[i+1] == 'D':  #Counts number of ' D ' between position forward BC and moves end point and reverse BC negatively accordingly
							bar_start_r = bar_start_r - int(cigar_op[i])  #Counts number of ' D ' and skips current position  accordingly
							bar_end_r = bar_start_r + 26

						elif cigar_op[i+1] == 'M':
							pos = pos + int(cigar_op[i])  #Counts the number of ' M  ' and moves current position positively accordingly

					elif pos >= bar_start_r and pos < bar_end_r:
						#if we are within the second BC change the end point of the second BC

						if cigar_op[i+1] == 'I':  #Counts number of ' I ' within the forward BC and moves end point and reverse BC positively accordingly
							bar_end_r = (bar_start_r + 26) + int(cigar_op[i])

						elif  cigar_op[i+1] == 'D':  #Counts number of ' D ' between position forward BC and moves end point and reverse BC negatively accordingly
							bar_end_r = (bar_start_r + 26) - int(cigar_op[i])

						elif cigar_op[i+1] == 'M':
							pos = pos + int(cigar_op[i])  #Counts the number of ' M  ' and moves current position positively accordingly


			f_sequence = sequence[bar_start_f:bar_end_f]
			f_quality = quality[bar_start_f:bar_end_f]

			r_sequence = sequence[bar_start_r:bar_end_r]
			r_quality = quality[bar_start_r:bar_end_r]

			if len(f_sequence) != 26 or len(r_sequence) != 26:
				count_diff_size += 1
			else:
				count_26 += 1

			D[seq_id] = [sequence[0:8],f_sequence,r_sequence] # UMI is 0, BC2 is 1, BC1 is 2
		except StopIteration:
			break
	handle_sam.close()
	with open(infile + '.json', 'w') as file2:
		json.dump(D, file2)
	print('number of unaligned read with * = ' + str(counter_star))	#Prints out the counter(number) of cigar strings with ' * '
	print('number reads with one BC at least with wrong size =' + str(count_diff_size))
	print('number reads with both BC at good size = ' + str(count_26))

# BC2 clustering with bartender
# collects the BC2 from all libraries in a single file
# then feeds that to bartender
def bartender(dir, Fitness_i7, Fitness_i5):
	csvFile = open(dir + 'BC2_bartenderInput.csv','w') # file that will be used as input to bartender
	correctLength = 0
	wrongLength = 0
	for block in range (len(Fitness_i7)):
		cols = Fitness_i7[block]
		tps = Fitness_i5[block]
		for t in tps:
			for c in cols:
				with open(dir + t+c + '.json') as f:
					reads = json.load(f) # key is read ID, UMI is 0, BC2 is 1, BC1 is 2
				for id in reads:
					BC2 = reads[id][1]
					if len(BC2)==26:
						csvFile.write('{},{}\n'.format(BC2, id))
						correctLength += 1
					else:
						wrongLength += 1
	print ('BC2 with 26 nt = ' + str(correctLength))
	print ('BC2 with wrong length = ' + str(wrongLength))
	csvFile.close()
	subprocess.call('bartender_single_com -f ' + dir + 'BC2_bartenderInput.csv -o ' + dir + 'AllBC2 -c 1 -l 8 -z -1 -d 4', shell=True)

# uses the bartender output to update the barcodes from each library
# store in new directory
# also updates names based on a key
def updateReadsV2(dirIN, dirOUT, Fitness_i7, Fitness_i5, BC1list_dictionary, libraryNames_dictionary):
    # construct helpful dictionaries from bartender output
    BC2_clustered = {} # key = BC2, value = center.ID
    BC2_centers = {} # key = center.ID, value = center.sequence, at the end I will want BC2_centers[BC2_clustered[bc2]]
    clusters = open(dirIN + 'AllBC2_barcode.csv', 'r')
	for line in clusters:
		if line.startswith('Unique.reads'):
			continue
		else:
			L = line.strip('\n')
			L1 = L.split(',')
			BC2, clusterID = L1[0], L1[2]
			BC2_clustered[BC2] = clusterID
	clusters.close()
	centers = open(dirIN + 'AllBC2_cluster.csv', 'r')
	for line in centers:
		if line.startswith('Cluster.ID'):
			continue
		else:
			L1 = line.split(',')
			clusterID, center = L1[0], L1[1]
			BC2_centers[clusterID] = center
	centers.close()

	# loop through libraries
	for block in range (len(Fitness_i7)):
		cols = Fitness_i7[block]
		tps = Fitness_i5[block]
		for t in tps:
			for c in cols:
				index_pair = t+c
				library = libraryNames_dictionary[index_pair]
				BC1list = BC1list_dictionary[library]
				with open(dirIN + index_pair + '.json') as f:
					Reads = json.load(f) # key is read ID, UMI is 0, BC2 is 1, BC1 is 2
				# round 1 replace bc1 with their center, if applicable, or reject
				PassBC1 = {} # here store the reads that will pass BC1 filter (BC1 close to true BC1)
				RejectBC1 = {} # here store the reads whose BC1 was not a match
				for key in Reads:
					id, umi, bc2, bc1 = key, Reads[key][0], Reads[key][1], Reads[key][2]
					if bc1 in BC1list:
						PassBC1[id] = [umi,bc1,bc2]
					else:
						match = 0
						for ExpBC1 in BC1list:
							LD = LevenshteinDistance(bc1, ExpBC1)
							if LD < 3:
								match = 1
								PassBC1[id] = [umi,ExpBC1,bc2]
								break
						if match == 0:
							RejectBC1[id] = [umi,bc1,bc2]
				#print (index_pair, library, len(PassBC1), len(RejectBC1))
				# round 2 replace bc2 with center if applicable (things with 'too' small or 'too' big bc are rejected)
				UpdatedReads = {}
				RejectBC2 = {}
				for key in PassBC1:
					id, umi, bc1, bc2 = key, PassBC1[key][0], PassBC1[key][1], PassBC1[key][2]
					if bc2 in BC2_clustered:
						UpdatedReads[id] = [umi,bc1,BC2_centers[BC2_clustered[bc2]]]
					else:
						match = 0
						for bc in BC2_clustered:
							LD = LevenshteinDistance(bc, bc2)
							if LD < 3:
								match = 1
								UpdatedReads[id] = [umi,bc1,bc]
								break
						if match == 0:
							RejectBC2[id] = [umi,bc1,bc2]
				print (index_pair, library, len(PassBC1), len(RejectBC1), len(UpdatedReads), len(RejectBC2))
				with open(dirOUT + library + '_updated.json', 'w') as f:
					json.dump(UpdatedReads, f)

# it goes into the updated read files and outputs a file with counts.
def counts(dir, library):
	with open(dir + library + '_updated.json') as f:
		Reads = json.load(f)	# clustered[ID] = UMI, BC1, BC2
	BCcount = {} # make an entrance for each total BC to store count total and how many times each UMI appears
	for r in Reads:
		UMI = Reads[r][0]
		BC = Reads[r][1]+'_'+Reads[r][2]
		if BC not in BCcount:
			BCcount[BC] = [0,{}]
		BCcount[BC][0] += 1 # first entry is the total count for the barcode
		if UMI not in BCcount[BC][1]: # second entry is the UMI dictionary and how many times each UMI appears
			BCcount[BC][1][UMI] = 0
		BCcount[BC][1][UMI] +=1
	for bc in BCcount:
		BCcount[bc].append(len(BCcount[bc][1])) # 3rd entry: this is the length of the UMI dictionary and represents the deduplicated BC or number of genomes with barcode
		BCcount[bc].append(BCcount[bc][0] - BCcount[bc][2]) # 4th entry: this the number of duplicate reads for the barcode
	with open(dir + library + '_counts.json', 'w') as f:
		json.dump(BCcount, f)
	print (library, len(BCcount))

# calculates aggregate counts for assays where the same pools were used
# I will plot the output as histograms and based on that decide on a filtering step
# arguments: directory, pool, dictionary with libraries
# outputs aggregate counts to put in r and make histogarms
def aggregateCounts(dir, pool, libraries):
	AggrC = {}
	for l in libraries:
		with open(dir + l + '_counts.json') as f: # key = barcode, value list = total counts, umi count dictionary, dedup counts, duplicate counts
			counts = json.load(f)
		for bc in counts:
			if bc not in AggrC:
				AggrC[bc] = [0,0]
			AggrC[bc][0] += counts[bc][2] # counts
			AggrC[bc][1] += 1 # number of libraries
	with open(dir + pool + '_AggregateCounts.json', 'w') as f:
		json.dump(AggrC, f) # key is bc, value = list of total counts and number of libraries
	F1 = open(dir + pool + '_AggregateCounts.txt','w')
	F1.write('BC\tcounts\tlibraries\n')
	for bc in AggrC:
		F1.write('{}\t{}\t{}\n'.format(bc,AggrC[bc][0],AggrC[bc][1]))
	F1.close()

# based on the plots from previous step, I filter out barcodes with an aggregate count <50
def filter(dir, pool, library, threshold):
# load aggregate counts
	with open(dir + pool + '_AggregateCounts.json') as f:
		AggrC = json.load(f) # key is bc, value = list of total counts and number of libraries
# load counts
	with open(dir + library + '_counts.json') as f: # key = barcode, value list = total counts, umi count dictionary, dedup counts, duplicate counts
		counts = json.load(f)
	Filtered_counts = {}
	for bc in counts:
		if AggrC[bc][0] > threshold:
			Filtered_counts[bc] = counts[bc][2]
	with open(dir + library + '_FilteredCounts.json', 'w') as f:
		json.dump(Filtered_counts, f) # key is bc, value = counts
	#print (pool, library, len(counts), len(Filtered_counts))
	F1 = open (dir + 'Filtering_summary.txt','a+')
	F1.write('{}\t{}\t{}\t{}\n'.format(pool, library, len(counts), len(Filtered_counts)))
	F1.close()

# Assemble all barcode counts per assay and output a single file (time trajectories with just raw counts)
# Not normalized, use as input in fitness estimation algorithm
def trajectories(dirIN, dirOUT,assay,name): # takes as argument the directory and the files associated with a single assay in a form of a dictionary and a name for the assay (output file)
    trajectory = {}
    timepoints = [0,2,4,6]
	for i in timepoints:
		with open(dirIN + assay[i] + '_FilteredCounts.json') as f: # key = barcode, value = counts
			filteredcounts = json.load(f)

		for bc in filteredcounts:
			if bc not in trajectory:
				trajectory[bc] = [0]*len(assay)
			trajectory[bc][i] = filteredcounts[bc]

	with open(dirOUT + name + '_Trajectory.json', 'w') as f:
		json.dump(trajectory, f) # key is bc, values = counts
	F1 = open(dirOUT + name + '_Trajectory.txt', 'w')
	F1.write('{}\t'.format('BC'))
	for g in range(3):
		F1.write('{}\t'.format(timepoints[g]))
	F1.write('{}\n'.format(timepoints[len(timepoints)-1]))
	for bc in trajectory:
		F1.write('{}\t'.format(bc))
		for g in range(len(trajectory[bc])-1):
			F1.write('{}\t'.format(trajectory[bc][g]))
		F1.write('{}\n'.format(trajectory[bc][len(trajectory[bc])-1]))
	F1.close()

# goes into the count files and outputs a single file per assay with normalized counts
# I will plot these to decide on timepoints to use for fitness estimation (fitness estimation algorithm from Venkataram et al Cell 2016)
def norm_trajectories(dir,assay,name): # takes as argument the directory and the files associated with a single assay in a form of a dictionary and a name for the assay (output file)
    trajectory = {}
    T = [0,2,4,6]
	for timepoint in T:
		with open(dir + assay[timepoint] + '_FilteredCounts.json') as f: # key = barcode, value = counts
			filteredcounts = json.load(f)
		with open(dir + assay[timepoint] + '_counts.json') as f: # key = barcode, value = counts (use that to calculate the total - denominator)
			counts = json.load(f)
		#print (counts)
		totalCounts = sum(counts[bc][2] for bc in counts.keys())
		for bc in filteredcounts:
			if bc not in trajectory:
				trajectory[bc] = [0]*4
			trajectory[bc][i] = filteredcounts[bc]*1000.0/totalCounts

	with open(dir + name + '_normTrajectory.json', 'w') as f:
		json.dump(trajectory, f) # key is bc, values = counts
	F1 = open(dir + name + '_normTrajectory.txt', 'w')
	F1.write('{}\t'.format('BC'))
	for g in range(3):
		F1.write('{}\t'.format(T[g]))
	F1.write('{}\n'.format(T[len(T)-1]))
	for bc in trajectory:
		F1.write('{}\t'.format(bc))
		for g in range(len(trajectory[bc])-1):
			F1.write('{}\t'.format(trajectory[bc][g]))
		F1.write('{}\n'.format(trajectory[bc][len(trajectory[bc])-1]))
	F1.close()

# merge paired-end reads, align and extract barcodes
for block in range (len(Fitness_i7)):
	index7 = Fitness_i7[block]
	index5 = Fitness_i5[block]
	for i5 in index5:
		for i7 in index7:
			align(fastqDir,alnDir,i5+i7)
			barcodes(alnDir,i5+i7)

# collects barcodes from all libraries in a single file and clusters
bartender(alnDir,Fitness_i7, Fitness_i5)

# store library information in json format
json_files(LibraryKey, LibraryKey, 'BC1LibraryKey')
# load json files in dictionaries
with open(LibraryKey + '.json') as f1: # LibKey: key = index pair, value = library name
	LibKey = json.load(f1)
with open('BC1LibraryKey.json') as f1: # BC1Key: key = library name, value = BC1 list
	BC1Key = json.load(f1)
# use bartender oputput to update barcodes with their centers across all libraries.
# Also update file names with meaningful library identifiers
updateReadsV2(alnDir, BarDir, Fitness_i7, Fitness_i5, BC1Key, LibKey)

# goes over the files with the updated barcodes by name and outputs barcode counts
for key in BC1Key: # BC1Key: key = library name, value = BC1 list
	counts(BarDir,key)

# derive aggregate counts from all libraries that were based of the same pool
for p in Pools:
	libraries_list = [] # store all libraries derived from same initial pool (that is all timepoints from all assays initiated)
	libraries_list.append('t0_'+p)
	for col in Pools[p]:
		for plate in [3,4]:
			for timepoint in [2,4,6]:
				libraries_list.append('t'+str(timepoint)+'_p'+str(plate)+'_'+col)
	aggregateCounts(BarDir, p, libraries_list)
    # then go over the count files and filter out barcodes with counts less than a threshold decided after plotting the aggregate data
	for lib in libraries_list:
		filter(BarDir, p, lib, 49)

# assemble the counts per timepoint into trajectories
# Use normalized counts to plot and exclude potential low coverage timepoints for fitness estimation
# use raw counts as input for fitness estimation.
for p in Pools:
	for col in Pools[p]:
		for plate in [3,4]:
			timecourse = {}
			timecourse[0] = 't0_'+p
			for timepoint in [2,4,6]:
				timecourse[timepoint] = 't'+str(timepoint)+'_p'+str(plate)+'_'+col
			name = 'p'+str(plate)+'_'+col # that is the name of the library
			trajectories(BarDir, TrajDir, timecourse, name)
            norm_trajectories(TrajDir, timecourse, name)
