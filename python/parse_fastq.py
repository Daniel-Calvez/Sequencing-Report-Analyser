#!/usr/bin/env python3

#Import modules-----------------------------------
import multiprocessing as mp
import mmap
import os
import sys
import time
import warnings
try:
	import argparse
except ImportError as exception:
	sys.exit("ERROR: Argparse module is not installed")

try:
	import pandas as pd
except ImportError as exception:
	sys.exit("ERROR: Pandas module is not installed")

#delete warnings
warnings.filterwarnings("ignore")

	
#Define function-----------------------------------------		 	

##Function to define command line arguments=============

def argument_parser():
	"""
	This function allows to extract command lines
	arguments

	PARAMETERS:
	None

	RETURN:
	args: a argument command line object
	"""
	parser = argparse.ArgumentParser(prog = "parse_fastq.py",
	description = "This program allows to extract read length and sequencing quality score average for each reads")
	parser.add_argument("-i", "--input", help = "Fastq file to parse. Must be in fastq format.")
	parser.add_argument("-c", "--cpu", help = "Number of CPU to use. Default 1.", default = 1, type = int)
	parser.add_argument("-s", "--sample", help = "Desired sample name. Default will be file name.")
	parser.add_argument("-o", "--output", help = "Output directorie. Four files will be created inside it.")
	args = parser.parse_args()
	return(args)

##Function to parse arguments
def parse_arguments(args):
	"""
	This function allows to parse command line arguments

	PARAMETERS:
	args: an argument object

	RETURN:
	args: an argument object
	"""
	#Looking for empty arguments
	if args.input is None:
		sys.exit("ERROR: You must provided a fastq/fastq.gz input file")
	
	if args.output is None:
		print("WARNING: output directorie is missing. Output file will be created on your current working directorie")
		args.output = "."
	
	#Looking file or folder existence
	if not os.path.exists(args.input):
		sys.exit("ERROR: Input file not found or doesn't exist.")
	
	if not os.path.exists(args.output):
		sys.exit("ERROR: Output directorie doesn't exist or not found")
	
	#Looking file extension
	extension_format = args.input.split(".")[-1]
	if extension_format != "fastq":
		sys.exit("ERROR: Input file is not a fastq file")

	#Monitoring asked CPU
	cpu_count = mp.cpu_count()
	args.cpu = int(args.cpu)
	if args.cpu > cpu_count:
		print(f"WARNING: Too much CPU aksed. Your computer contains only {cpu_count} CPU. CPU parameters was reset to {cpu_count}")
		args.cpu = cpu_count

	#Add sample name if is None
	if args.sample is None:
		args.sample = args.input.split("/")[-1].split(".")[0]
	
	return (args)


##Function to compute runing time=======================

def compute_runing_time(start_time, end_time):
	"""
	This function compute running time in hour, minute, second and milisec
	"""
	diff_time = end_time - start_time
	hour = int(diff_time // 3600)
	minute = int((diff_time - hour * 3600) // 60)
	seconde = int((diff_time - (hour * 3600) - (minute * 60)))
	centisec = int((diff_time - (hour * 3600) - (minute * 60) - seconde)*100)
	return hour, minute, seconde, centisec
    
##Function to compute file size===========================

def convert_file_size(file_to_read):
	"""
	This function computes file's size

	PARAMETERS:
	file_to_read: the fastq file to process

	RETURN
	file_size_formated: file size
	units: the size file unit, either Ko, Mo or Go
	"""
	units = [(1e3, "Ko"), (1e6, "Mo"), (1e9, "Go")]
	file_size = os.path.getsize(file_to_read)
	for i in range(len(units)):
		file_size_formated = file_size / units[i][0]
		if file_size_formated < 1e3:
			break
	file_size_formated = round(file_size_formated, 2)
	return file_size_formated, units[i][1]
	

##Function to calibrate real needed cpu===================

def correct_cpu(file_to_read, nb_cpu):
	"""
	This function define the number of needed CPU among allocated

	PARAMETERS:
	file_to_read: the fastq file to process
	nb_cpu: CPU number to use

	RETURN:
	read_number: the total number of reads to process
	nb_cpu: the corrected CPU number to use 
	"""
	lines_count = 0
	with open(file_to_read, "r") as input_file:
		with mmap.mmap(input_file.fileno(), length = 0, access = mmap.ACCESS_READ) as mmap_input:
			while mmap_input.readline():
				lines_count += 1
	if lines_count % 4 != 0:
		sys.exit("ERROR: fastq file contains missing data")
	else:
		read_number = lines_count // 4
	while read_number // nb_cpu == 0:
		nb_cpu -= 1
		if nb_cpu == 0:
			sys.exit("ERROR: fastq file contains some mistake")
	
	print("")
	print(f"Process will use {nb_cpu} CPU and deal arround {read_number // nb_cpu} reads/CPU")
	return read_number, nb_cpu

##Function to verify lines===============================

def verify_new_read(mmap_object, position):
	"""
	This function verify if a position is at the start of a new read

	PARAMETERS:
	mmap_object: a mmap object from a mmap file reading
	position: the position to check

	RETURN:
	a boolean idnicating whether a new read start
	"""
	if position == 0:
		return True
	else:
		mmap_object.seek(position)
		char = mmap_object.read(2).decode()
		return char == "\n@"

def get_next_read(mmap_object, position):
	"""
	This function give the position of the next read

	PARAMETERS:
	mmap_object: a mmap object from a mmap file reading
	position: the psoition to check
	"""
	for i in range(4):
		mmap_object.readline()
	position = mmap_object.tell()
	return position
	
##Function for parallel reading============================

def parallel_read(file_to_read, read_number, nb_cpu):
	"""
	This function allows a parallel reading for a fasq file

	PARAMETERS:
	file_to_read: the fastq file to process
	read_number: number of reads containing in the fastq file
	nb_cpu: CPU number to use

	RETURN:
	chunk_results: a list of all individual chunk analysis
	"""
	print("")
	print("Defining chunk limits...")
	
	#Compute file size
	file_size = os.path.getsize(file_to_read)
	
	#Compute approximate chunk size
	chunk_size_approx = file_size // read_number * (read_number // nb_cpu)
	
	#Define chunk_args list
	chunk_args = []
	
	#Compute chunk position
	with open(file_to_read, "r") as input_file:
		with mmap.mmap(input_file.fileno(), length = 0, access = mmap.ACCESS_READ) as mmap_input:
			chunk_start = 0
			while chunk_start < file_size:
				chunk_end = min(chunk_start + chunk_size_approx, file_size)
				while not verify_new_read(mmap_input, chunk_end):
					chunk_end -= 1
					if chunk_end == chunk_start:
						chunk_end = get_next_read(mmap_input, chunk_end)
						break		
				args = (file_to_read, chunk_start, chunk_end)
				chunk_args.append(args)
				chunk_start = chunk_end + 1
	print("")			
	print("Analyse chunk...")					
	with mp.Pool(nb_cpu) as p:
		chunk_results = p.starmap(process_chunk, chunk_args)
	return chunk_results

##Function for process a chunk=========================================

def process_chunk(file_to_read, chunk_start, chunk_end):
	"""
	This function processes a chunk

	PARAMETERS:
	file_to_read: the fastq file to process
	chunk_start: chunk tsrat position
	chun_end: chunk end position

	RETURN:
	chunk_result: a list of results extracted from chunk analysis
	"""
	with open(file_to_read, "r") as input_file:
		with mmap.mmap(input_file.fileno(), length = 0, access = mmap.ACCESS_READ) as mmap_input:
			#Initialize dataframes
			read_length_df = pd.DataFrame(columns = ['read_length', 'occurence'])
			seq_score_df = pd.DataFrame(columns = ['Seq_score', 'occurence'])
			gc_content_df = pd.DataFrame(columns = ['GC_content', 'occurence'])
			while chunk_start < chunk_end:
				#Reach chunk_start position
				mmap_input.seek(chunk_start)
				#Extract line
				line = mmap_input.readline().decode()
					
				#Process line
				##Parse sequence
				if line.startswith("@"):
					chunk_start += len(line)
					### Extract sequence
					mmap_input.seek(chunk_start)
					line = mmap_input.readline().decode()
					read_length_df = process_read_length(line.strip(), read_length_df)
					gc_content_df = process_GC_content(line.strip(), gc_content_df)
					##Update chunk_start
					chunk_start += len(line)
					
				##Parse sequencing quality score
				if line == "+\n":
					chunk_start += len(line)
					### Extract sequencing quality score
					mmap_input.seek(chunk_start)
					line = mmap_input.readline().decode()
					seq_score_df = process_seq_score(line.strip(), seq_score_df)
					##Update chunk_start
					chunk_start += len(line)
	chunk_result = [read_length_df, seq_score_df, gc_content_df]
	
	return chunk_result
						
						
					
##Functions for process a line====================================			

def process_read_length(line, read_length_df):
	"""
	This function computes read length

	PARAMETERS:
	line: a fastq line stripped containing the sequence
	read_length_df: a read length data-frame to complete

	RETURN:
	read_length_df: the updated data-frame
	"""
	read_length = len(line)
	if read_length in list(read_length_df["read_length"]):
		index = list(read_length_df.index[read_length_df["read_length"] == read_length])
		index = int(index[0])
		read_length_df.loc[index, ["occurence"]] += 1
	else:
		to_add = [read_length, 1]
		read_length_df.loc[len(read_length_df)] =  to_add
	return read_length_df

def process_seq_score(line, seq_score_df):
	"""
	This function computes sequencing quality score average for a read

	PARAMTERS:
	line: a fastq line stripped containing the sequence
	seq_score_df: a sequencing score data-frame to complete

	RETURN:
	seq_score_df: the updated data-frame
	"""
	read_length = len(line)
	cum_score = 0
	for base in line:
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

def process_GC_content(line, gc_content_df):
	"""
	This function computes GC content from a given read

	PARAMETERS:
	line: a fastq line stripped containing the sequence
	gc_content_df: a GC content data-frame to complete

	RETURN:
	gc_content_df: the updates data-frame
	"""
	#Compute read length
	read_length = len(line)
	gc_count = 0
	for base in line.upper():
		if base == "G" or base == "C":
			gc_count += 1
	gc_content = round(gc_count/read_length, 4)
	if gc_content in list(gc_content_df["GC_content"]):
		index = list(gc_content_df.index[gc_content_df["GC_content"] == gc_content])
		index = int(index[0])
		gc_content_df.loc[index, ["occurence"]] += 1
	else:
		to_add = [gc_content, 1]
		gc_content_df.loc[len(gc_content_df)] =  to_add
	
	return gc_content_df

##Function combine results============================================

def combine_results(chunk_results, sample_name, output_dir, 
					file_size, unit_size, read_number):
	"""
	This function combines chunk results and write their in separate text files

	PARAMETERS:
	chunk_results: a list of all individual chunk result
	sample_name: a sample name
	output_dir: an output path
	file_size: the file size
	unit_size: the size unit
	read_number: the number of read
	"""
	#Initilialize data-frame
	read_length_df = pd.DataFrame(columns = ['read_length', 'occurence'])
	seq_score_df = pd.DataFrame(columns = ['Seq_score', 'occurence'])
	gc_content_df = pd.DataFrame(columns = ['GC_content', 'occurence'])

	#build data frame
	for i in range(len(chunk_results)):
		read_length_df = pd.concat([read_length_df, chunk_results[i][0]], ignore_index = True)
		seq_score_df = pd.concat([seq_score_df, chunk_results[i][1]], ignore_index = True)
		gc_content_df = pd.concat([gc_content_df, chunk_results[i][2]], ignore_index = True)
	
	#Summarise data frames
	read_length_df = read_length_df.groupby("read_length").sum().reset_index()
	seq_score_df = seq_score_df.groupby("Seq_score").sum().reset_index()
	gc_content_df = gc_content_df.groupby("GC_content").sum().reset_index()

	#Write data-frame
	##Write read length
	with open(os.path.join(output_dir, sample_name + "_read_length.txt"), "w") as read_length_file:
		read_length_file.write(read_length_df.to_string(index = False))
		read_length_file.write("\n")

	##Write quality scores
	with open(os.path.join(output_dir, sample_name + "_seq_score.txt"), "w") as seq_score_file:	
		seq_score_file.write(seq_score_df.to_string(index = False))
		seq_score_file.write("\n")
	
	##Write GC content
	with open(os.path.join(output_dir, sample_name + "_GC_content.txt"), "w") as gc_content_file:	
		gc_content_file.write(gc_content_df.to_string(index = False))
		gc_content_file.write("\n")

	##Write metadata informations
	with open(os.path.join(output_dir, sample_name + "_global_info.txt"), "w") as global_info:
		global_info.write(f"File_size:{file_size} {unit_size}\n")
		global_info.write(f"Read number:{read_number}\n")

#Function for main program==============================================

def main():
	"""
	This function allows to run main program

	PARAMETERS:
	None

	RETURN:
	None
	"""
	#Extract command line arguments
	args = argument_parser()

	#Parse given arguments
	args = parse_arguments(args)

	#Get file size
	file_size, unit = convert_file_size(args.input)
	print("")
	print(f"File size: {file_size} {unit}")

	#Record start time
	start_time = time.time()

	#Define the real number of needed CPU
	read_number, nb_cpu = correct_cpu(args.input, args.cpu)

	#Extract list of chunk results
	chunk_results = parallel_read(args.input, read_number, nb_cpu)

	#Combine and write results
	combine_results(chunk_results, args.sample, args.output,
				 file_size, unit, read_number)
	
	#record end time
	end_time = time.time()

	#Compute running time
	hour, minute, seconde, centisec = compute_runing_time(start_time, end_time)

	print("")
	print(f"Run time: {hour}h:{minute}m:{seconde}s::{centisec}cs")

#Main program-------------------------------------------

if __name__ == "__main__":
	main()