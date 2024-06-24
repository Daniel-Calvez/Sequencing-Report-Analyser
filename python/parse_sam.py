#!/usr/bin/env python3
#========================================================#
###             Assess DORADO basecalling              ###

#Date of creation: 16th of January 2024 (Update)
#Author: FAUCHOIS Antoine, Bioinformatician engineer
#Site: AP-HP, La Pitié-Salpêtrière, FRANCE
#========================================================#

#Import modules-------------------------------------------
import multiprocessing as mp
import mmap
import os
import sys
import time
import warnings
##Import argparse
try:
	import argparse
except ImportError as exception:
	sys.exit("ERROR: argparse module is not installed.")
##Import pandas=
try:
	import pandas as pd
except ImportError as exception:
	sys.exit("ERROR: pandas module is not installed")

#sample_name = args.sam.split("/")[-1].split(".sam")[0]

#Define function--------------------------------------------

def get_arguments():
	"""
	This function allows to extract given arguments
	in command line

	PARAMETERS:
	None

	RETURN:
	args: a list of given command line arguments
	"""
	#Define parser
	parser = argparse.ArgumentParser(prog = "parse_sam.py",
	description = "This program allows to extract quality score and mapping status from a sam file.")
	#Add arguments
	parser.add_argument("-i", "--input", type = str, help = "Path to SAM file to parse.")
	parser.add_argument("-c", "--cpu", type = int, help = "Number of CPU to use. Default: 1",
					 default = 1)
	parser.add_argument("-o", "--output", type = str, help = "Output path. It will contain four files.")
	args = parser.parse_args()
	return args

def parse_arguments(args):
	"""
	This function allows to parse given command line arguments

	PARAMTERS:
	args: a list of given command line arguments

	RETURN:
	None
	"""
	#Empty input
	if args.input is None:
		sys.exit("ERROR: input argument is required")
	#Empty output
	if args.output is None:
		sys.exit("ERROR: output argument is required")
	
	#Non-existing input file
	if not os.path.exists(args.input):
		sys.exit("ERROR: input file not found or doesn't exist")
	
	#Non existing output path
	if not os.path.exists(args.output):
		sys.exit("ERROR: output path not found or doesn't exist")

def convert_to_bin(number):
	"""
	This function convert a ten base number in binary number (16 bits formats)

	PARAMETERS
	number an integer

	RETURN:
	binary_number: a 16 bytes converted binary number
	"""
	#Bin function create a list with a unusefull part 0b 
	#rjust to adjust size to 16 bytes
	binary_number = bin(number)[2:].rjust(16,'0')[::-1]
	return binary_number

def detect_unmapped(bin_flag):
	"""
	This function is abble to detect unmapped reads from a sam flag

	PARAMETERS:
	bin_flag: a binary sam flag

	RETURN:
	a boolean giving if a reads was aligned on reference.
	"""
	if bin_flag[2] == "1":
		return True
	else:
		return False

def compute_runing_time(start_time, end_time):
    """
    This function compute time in h:m:s from start and end time
    
    PARAMETERS
    start_time: start time
    end_time: end time

    RETURN
    hour, minute, seconde: a tuple of numeric corresponding to time in h:m:s
    """
    diff_time = round(end_time - start_time)
    hour = diff_time // 3600
    minute = (diff_time - hour * 3600) // 60
    seconde = (diff_time - (hour * 3600) - (minute * 60))
    return hour, minute, seconde

def convert_file_size(file_to_read):
	"""
	This function computes file's size
	
	PARAMETERS
	file_to_read: path to SAM file

	RETURN
	file_size_formated, units: tuple with file size and unit
	"""
	units = [(1e3, "Ko"), (1e6, "Mo"), (1e9, "Go")]
	file_size = os.path.getsize(file_to_read)
	for i in range(len(units)):
		file_size_formated = file_size / units[i][0]
		if file_size_formated < 1e3:
			break
	file_size_formated = round(file_size_formated, 2)
	return file_size_formated, units[i][1]


def correct_cpu(file_to_read, nb_cpu):
	"""
	This function define the number of needed CPU among allocated

	PARAMETERS:
	file_to_read: SAM file to parse
	nb_cpu: number of allocated CPU

	RETURN:
	read_number, nb_cpu: a tuple giving the number of reads and needed CPU number
	"""
	lines_count = 0
	with open(file_to_read, "r") as input_file:
		with mmap.mmap(input_file.fileno(), length = 0, access = mmap.ACCESS_READ) as mmap_input: 
			while mmap_input.readline() :
				lines_count += 1
	read_number = lines_count
	while read_number // nb_cpu == 0:
		nb_cpu -= 1
		if nb_cpu == 0:
			sys.exit("ERROR: fastq file contains some mistake")
	
	print(f"Process will use {nb_cpu} CPU and deal arround {read_number // nb_cpu} reads/CPU")
	return read_number, nb_cpu

def new_line_verify(mmap_object, position):
	"""
	This function verify if a position is after a \\n
	
	PARAMETERS
	mmap_object: mmap object corresponding to the file
	position: position to look
	
	RETURN
	boolean: boolean corresponding to the test
	"""
	if position == 0:
		return True
	else:
		mmap_object.seek(position)
		char = mmap_object.read(1)	
		return char.decode() == '\n'
		
def get_next_line_position(mmap_object, position):
	"""
	This function gets position of the next line
	
	PARAMETERS
	mmap_object: a mmap object
	position: position to look
	
	RETURN
	position: position of the next line
	"""
	mmap_object.seek(position)
	mmap_object.readline()
	return mmap_object.tell()

def parallel_read(file_to_read, read_number, nb_cpu):
	"""
	This function allows a parallel read for a huge file.
	
	PARAMETERS
	file_to_read: path to SAM file to read
	nb_cpu: number of CPU to use
	
	RETURN
	chunk_results: list of all chunk's result
	"""
		
	#Extract size and cpu information
	file_size = os.path.getsize(file_to_read)
	
	#Compute chunk size
	chunk_size = file_size // read_number * (read_number // nb_cpu)

	#Define chunk args list
	chunk_args = []
	
	#Get start and end for each chunk
	print("Defining chunk limits...")
	with open(file_to_read, "r") as file_input:
		with mmap.mmap(file_input.fileno(), length = 0, access = mmap.ACCESS_READ) as mmap_input:
		
			#Define position of the first chunk
			chunk_start = 0
			
			while chunk_start < file_size:
				chunk_end = min(chunk_start + chunk_size, file_size)

				#verify if chunk ends at the end of a line
				while not new_line_verify(mmap_input, chunk_end):
					#withdraw 1 up to reach a "\n"
					chunk_end -= 1
					#Case chunk_end has decrease until chunk_start
					if chunk_start == chunk_end:
						chunk_end = get_next_line_position(mmap_input, chunk_end)
						break	

				args = (file_to_read, chunk_start, chunk_end)
				chunk_args.append(args)
				chunk_start = chunk_end + 1
	
	print("Processing chunks...")			
	with mp.Pool(nb_cpu) as p:
		chunk_results = p.starmap(process_chunk, chunk_args)	
		
	return chunk_results	

def process_chunk(file_to_read, chunk_start, chunk_end):
	"""
	This function allows to process each line from a given chunk
	
	PARAMETERS
	file_to_read: Path to SAM file to read
	chunk_start: Start position to read the file
	chunk_end: End psoition to read the file
	
	RETURN:
	chunk_result: Result's dictionnary
	"""
	with open(file_to_read, 'r') as file_input:
		with mmap.mmap(file_input.fileno(), length = 0, access = mmap.ACCESS_READ) as mmap_input:
		
			#Initialize parameters
			count_read = 0
			count_map = 0
			
			#Initialize dataframe
			read_length_df = pd.DataFrame(columns = ['read_length', 'mapped', 'occurence'])
			map_score_df = pd.DataFrame(columns = ['Map_score', 'occurence'])
			seq_score_df = pd.DataFrame(columns = ['Seq_score', 'occurence'])
			
			#Process lines for all chunk
			while chunk_start < chunk_end:
				#Reach chunk_start position
				mmap_input.seek(chunk_start)
				#Extract line
				line = mmap_input.readline()
				#Process line
				line = line.decode()
				if not line.startswith("@"):
					line_split = line.split("\t")
					if len(line_split) >= 11:
						count_read = process_count_read(count_read)
						count_map = process_count_map(line_split, count_map)
						read_length_df = process_read_length(line_split, read_length_df)
						map_score_df = process_map_score(line_split, map_score_df)
						seq_score_df = process_seq_score(line_split, seq_score_df)
				
				#Update position
				len_line = len(line)
				chunk_start += len_line
	#Group results
	chunk_result = [count_read, count_map, read_length_df, map_score_df, seq_score_df]
	
	return chunk_result
				

#Functions for process a line

def process_count_read(count_read):
	"""
	This function update a counter
	
	PARAMETERS
	count_read: the read counter

	RETURN
	count_read: the read counter updated
	"""
	count_read += 1
	return count_read

def process_count_map(line_split, count_map):
	"""
	This function update mapped read count
	
	PARAMETERS
	line_split: a splitted line for each tabulation
	count_map: a mapped read counter
	
	RETURN
	count_map: the mapped read counter updated
	"""
	#Extract sam flag from SAM line
	flag = int(line_split[1])
	#Convert flag in binary
	binary = convert_to_bin(flag)
	#Decode binary number
	is_unmmaped = detect_unmapped(binary)

	if is_unmmaped == False:
		count_map += 1
	return count_map

def process_read_length(line_split, read_length_df):
	"""
	This function update a data frame which contains read length 
	and occurence from a line
	
	PARAMETERS
	line_split: a splitted line at each tabulation
	read_length_df: a data-frame that contains read length and their occurence
	this df is sended to be updated.
	
	RETURN
	read_length_df: Updated data-frame
	"""
	#Extract read length
	read_length = len(line_split[9])

	#Extract sam flag
	flag = int(line_split[1])
	binary = convert_to_bin(flag)
	is_unmapped = detect_unmapped(binary)
	map_value = 0 if is_unmapped else 1

	#Extract line index that match read length and mapped value
	index = read_length_df[(read_length_df["read_length"] == read_length) & (read_length_df["mapped"] == map_value)].index

	if len(index) > 0:
		index = int(index[0])
		read_length_df.loc[index, ["occurence"]] += 1
	else:
		to_add = [read_length, map_value, 1]
		read_length_df.loc[len(read_length_df)] =  to_add
	return read_length_df
	
def process_map_score(line_split, map_score_df):
	"""
	This function update a data frame which contains mapping quality score and
	occurence from a line
	
	PARAMETERS
	line_split: a splitted line at each tabulation
	map_score_df: a data-frame that contains mapping quality score and their occurence
	
	RETURN
	map_score_df: Updated data-frame
	"""
	flag = int(line_split[1])
	binary = convert_to_bin(flag)
	is_unmmaped = detect_unmapped(binary)
	if is_unmmaped == False:
		#Extract mapping score
		map_score = line_split[4]
		#Extract index
		index = map_score_df[(map_score_df["Map_score"] == map_score)].index
		if len(index) > 0:
			index = int(index[0])
			map_score_df.loc[index, ["occurence"]] += 1
		else:
			to_add = [map_score, 1]
			map_score_df.loc[len(map_score_df)] = to_add
	return map_score_df

def process_seq_score(line_split, seq_score_df):
	"""
	This function update a data frame which contains sequencing quality score 
	and occurence from a line
	
	PARAMETERS
	line_split: a splitted line at each tabulation
	seq_score_df: data-frame that contains sequencing quality score and their occurence
	
	RETURN
	read_length_df: Updated data-frame
	"""
	read_length = len(line_split[9])
	cum_score = 0
	for base in line_split[10]:
		##Convert ASCCII into score
		cum_score += ord(base)
	mean_qual_score = round(cum_score / read_length, 2)
		
	if mean_qual_score in list(seq_score_df["Seq_score"]):
		index = list(seq_score_df.index[seq_score_df["Seq_score"] == mean_qual_score])
		index = int(index[0])
		seq_score_df.loc[index, ["occurence"]] += 1
	else:
		to_add = [mean_qual_score, 1]
		seq_score_df.loc[len(seq_score_df)] =  to_add
	return seq_score_df	


def combine_result(chunk_results, file_to_read, output_dir):
	"""
	This function combine chun results and write it in text file saved in specified output dir path
	
	PARAMETERS
	chunk_results: list of chunk result
	file_to_read: path to file to deal
	output_dir: output direction path
	"""
	print("Combining chunk results...")
	#Get sample name
	sample_name = file_to_read.split("/")[-1].split(".sam")[0]
	
	#Initialize parameters
	count_read = 0
	count_map = 0
	
	#Initialize dataframe
	read_length_df = pd.DataFrame(columns = ['read_length', 'mapped', 'occurence'])
	map_score_df = pd.DataFrame(columns = ['Map_score', 'occurence'])
	seq_score_df = pd.DataFrame(columns = ['Seq_score', 'occurence'])
	
	#Compute parameters and build dataframes
	for i in range(len(chunk_results)):
		count_read += chunk_results[i][0]
		count_map += chunk_results[i][1]
		read_length_df = pd.concat([read_length_df, chunk_results[i][2]], ignore_index = True)
		map_score_df = pd.concat([map_score_df, chunk_results[i][3]], ignore_index = True)
		seq_score_df = pd.concat([seq_score_df, chunk_results[i][4]], ignore_index = True)
	
	#Summarise data frames
	read_length_df = read_length_df.groupby(['read_length', 'mapped']).sum().reset_index()
	map_score_df = map_score_df.groupby('Map_score').sum().reset_index()
	seq_score_df = seq_score_df.groupby('Seq_score').sum().reset_index()
	
	#write results
	print("Writing results...")
	with open(os.path.join(output_dir, sample_name + "_map.txt"), "w") as map_file:
		map_file.write("read\toccurence\n")
		map_file.write(f"Mapped\t{count_map}\n")
		map_file.write(f"No mapped\t{count_read - count_map}\n")
		map_file.write(f"total\t{count_read}\n")
	with open(os.path.join(output_dir, sample_name + "_read_length.txt"), "w") as read_length_file:
		read_length_file.write(read_length_df.to_string(index = False))
	with open(os.path.join(output_dir, sample_name + "_seq_score.txt"), "w") as seq_score_file:
		seq_score_file.write(seq_score_df.to_string(index = False))
	with open(os.path.join(output_dir, sample_name + "_map_score.txt"), "w") as map_score_file:
		map_score_file.write(map_score_df.to_string(index = False))
					

def main():
	"""
	This function allows to execute the main
	program

	PARAMETERS:
	None

	RETURN:
	None
	"""
	#Set to remove warnigs from pandas
	warnings.filterwarnings("ignore")

	#Extract given command line argument
	args = get_arguments()

	#Parse arguments
	parse_arguments(args)

	#Extract file size
	file_size, unit = convert_file_size(args.input)

	print(f"File size: {file_size} {unit}")

	#Record start time
	start_time = time.time()

	#Get read number and actually needed CPU
	read_number, nb_cpu = correct_cpu(args.input, args.cpu)

	#Extract results per chunk
	chunk_results = parallel_read(args.input, read_number, nb_cpu)

	#Combine result
	combine_result(chunk_results, args.input, args.output)

	#record end_time
	end_time = time.time()

	#Compute running time
	hour, minute, seconde = compute_runing_time(start_time, end_time)

	print(f"Run time: {hour}h:{minute}m:{seconde}s")



#Main program---------------------------------------------------------------------------------
	
if __name__ == "__main__":
	main()


