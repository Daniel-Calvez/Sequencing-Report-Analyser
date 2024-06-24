#========================================================#
###                     filter VCF                     ###

#Date of creation: 15th of March 2024 (update)
#Author: FAUCHOIS Antoine, Bioinformatician engineer
#Site: AP-HP, La Pitié-Salpêtrière, FRANCE
#========================================================#

#Import modules====================================================================================
import sys
import re
import os
import warnings
#Import argparse
try:
    import argparse
except ImportError:
    sys.exit("ERROR: argparse module is not installed.")
#Import PyVCF3
try:
    import vcf
except ImportError:
    sys.exit("ERROR: PyVCF3 module is not installed.")
#Numpy
try:
    import numpy as np
except ImportError:
    sys.exit("ERROR: Numpy module is not installed")

warnings.filterwarnings("ignore", category = SyntaxWarning)


#Define function===================================================================================

#Function to extract given command line argument---------------------------------------------------
def extract_argument():
    """
    This fucntion allows to extract given command line arguments

    PARAMETERS:
    None

    RETURN:
    args: an args object
    """
    #Define parser
    parser = argparse.ArgumentParser(prog = "parse_sam.py",
    description = "This program allows to extract quality score and mapping status from a sam file.")
    
    #Add arguments
    parser.add_argument("-b", "--bam", type = str, help = "Path to input BAM file")
    parser.add_argument("-s", "--sam", type = str, help = "Path to input SAM file")
    parser.add_argument("-f", "--ref", type = str, help = "Path to reference fasta file")
    parser.add_argument("-v", "--vcf", type = str, help = "Path to VCF file")
    parser.add_argument("-t", "--threshold", type = float, help = "Depth threshold value to filter individual variant", default = 0.2)
    parser.add_argument("-d", "--min_depth", help = "Minimimum global depth value", type = int,
                        default = 50)
    parser.add_argument("-c", "--cpu", type = int, help = "Number of CPU to use. Default: 1",
					 default = 1)
    parser.add_argument("-r", "--remove", action = "store_true", 
                        help = "Remove intermediate files", default = False)
    parser.add_argument("-o", "--output", type = str, help = "Path to output VCF file name.")
    #Collect arguments
    args = parser.parse_args()
    return args

#Function to parse arguments-----------------------------------------------------------------------
def parse_arguments(args):
    """
    This function allows to parse given command line
    arguments

    PARAMATERS:
    args: an argument object

    RETURN:
    None
    """
    #Empty bam and sam argument
    if args.bam is None and args.sam is None:
        sys.exit("ERROR: You need to provide oat least one of the following arguments: bam/sam.")
    
    #Empty ref argument
    if args.ref is None:
        sys.exit("ERROR: reference argument is empty.")

    #Empty vcf argument
    if args.vcf is None:
        sys.exit("ERROR: vcf argument is empty.")
    
    #Empty ouptut argument
    if args.output is None:
        sys.exit("ERROR: output argument is empty.")

    #Incorrect requested CPU
    if args.cpu > os.cpu_count():
        sys.exit(f"ERROR: requested CPU number is greater than available one's {os.cpu_count()}")
    
    #No existing BAM file
    if args.bam is not None:
        if not os.path.exists(args.bam):
            sys.exit("ERROR: BAM file not found or doesn't exist.")
    
    #No existing SAM file
    if args.sam is not None:
        if not os.path.exists(args.sam):
            sys.exit("ERROR: SAM file not found or doesn't exist.")
    
    #No existing path to output file
    if args.output.find("/") >= 0:
        ouptut_path = "/".join(args.output.split("/")[:-1])
        if not os.path.exists(ouptut_path):
            sys.exit("ERROR: output path not found or doesn't exist.")
    
    #Incorrect threshold value
    if not 0 <= args.threshold <= 1:
        sys.exit("ERROR: Threshold value must be between 0 and 1.")

#Functions that use samtools package===============================================================

##Function to convert input sam in bam file--------------------------------------------------------
def convert_sam(in_sam, nb_cpu):
    """
    This function allows to convert a sam file in bam file
    for mpileup phase

    PARAMETERS:
    in_sam: path to sam file to parse
    nb_cpu: CPU number requested for the process

    RETURN:
    file_name: the name of created file
    """
    #Define bam output file name
    if in_sam.find("/") >= 0:
        file_name = "".join(in_sam.split("/")[-1])
    else:
        file_name = in_sam
    file_name = file_name.split(".sam")[0] + ".bam"

    #Run samtools view command
    os.system(f"samtools view -b -@ {nb_cpu} {in_sam} > {file_name}")

    return file_name

##Function to sort a BAM file----------------------------------------------------------------------
def sort_bam(in_bam, nb_cpu):
    """
    This function allows to sort a given bam file

    PARAMETERS:
    in_bam: path to input bam file
    nb_cpu: number of requested CPU

    RETURN:
    file_name: the of created file
    """
    #Define bam output file name
    if in_bam.find("/") >= 0:
        file_name = "".join(in_bam.split("/")[-1])
    else:
        file_name = in_bam

    #Run samtools command
    os.system(f"samtools sort -@ {nb_cpu} {in_bam} -o {file_name}")

    return file_name

##Function to index a BAM file---------------------------------------------------------------------
def index_bam(in_bam, nb_cpu):
    """
    This function allows to index a bam file.

    PARAMETERS:
    in_bam: path to input BAM file
    nb_cpu: number of requested CPU

    RETURN:
    file_name: name of created file
    """
    if in_bam.find("/") >= 0:
        file_name = "".join(in_bam.split("/")[-1])
    else:
        file_name = in_bam
    file_name = file_name + ".bai"

    #run samtools command
    if not os.path.exists(file_name):
        os.system(f"samtools index -b -@ {nb_cpu} {in_bam} -o {file_name}")

    return file_name

##Function to lauch samtools mpileup command-------------------------------------------------------
def pileup_reads(in_bam, in_ref):
    """
    This function allows to perform a reads pileup from SAM/BAM
    input via samtools package

    PARAMETERS:
    in_bam: path to input bam file
    in_ref: path to reference fasta file

    RETURN:
    file_name: the name of created file
    """
    #Define pileup output file name
    if in_bam.find("/") >= 0:
        file_name = "".join(in_bam.split("/")[-1])
    else:
        file_name = in_bam 
    file_name = file_name.split(".bam")[0] + ".pileup"

    outputDirectory= re.sub("/[A-z]*\.[A-z]*$","/", args.output)

    #Lauch samtools mpileup
    ##max detph is set to 0 to probe all reads
    os.system(f"samtools mpileup -d 0 -f {in_ref} {in_bam} -o {outputDirectory}{file_name}")

    return file_name

#Function to sort a VCF file=======================================================================

def sort_vcf(in_vcf, out_vcf):
    """
    This function allows to sort a VCF file by position

    PARAMETERS:
    in_vcf: input VCF file
    out_vcf: out sorted VCF file
    """
    os.system(f"touch {out_vcf}")
    os.system(f'cat {in_vcf} | egrep "^#" >> {out_vcf}')
    os.system(f'sort -k1,1V -k2,2n {in_vcf} | egrep -v "^#" >> {out_vcf}')

#Function to define CNV============================================================================

#Function to define deletion-----------------------------------------------------------------------
def define_deletion(reference, alternative):
    """
    This function allows to define a deletion from reported alt and ref

    PARAMETERS:
    reference: reference base(s)
    alternative: alternative base(s)

    RETURN:
    target_deletion: a the corresponding deletion
    """
    len_alt = len(alternative)
    deletion = reference[len_alt:]
    target_deletion = "-" + str(len(deletion)) + "(" + str(deletion).upper() + "|" + str(deletion).lower() + ")"
    return target_deletion

#Function to define insertion----------------------------------------------------------------------
def define_insertion(alternative):
    """
    This function allows to define a deletion from reported alt and ref

    PARAMETERS:
    alternative: alternative base(s)

    RETURN:
    target_insertion: a the corresponding insertion
    """
    len_alt = len(alternative)
    target_insertion = "+" + str(len_alt) + "(" + str(alternative).upper() + "|" + str(alternative).lower() + ")"
    return target_insertion

#Function to browse a pileup file==================================================================

#Function to get the previous pileup line__________________________________________________________
def read_previous_line(in_pileup):
    """
    This function allows to redefine the text pointer position
    to the previous line

    PARAMETERS:
    in_pileup: a TextIOWrapper object corresponding to the pileup file

    RETURN:
    None
    """
    #Define line return count
    count_line_return = 0
    #Define start pointer position
    start_position = in_pileup.tell()

    #Browse position one by one to seek the previous line
    i = 0
    while count_line_return < 2:
        if start_position - i == 0:
            in_pileup.seek(0)
            break
        else:
            in_pileup.seek(start_position - i)
            value = in_pileup.read(1)
            if value == "\n":
                count_line_return += 1
            i += 1

#Function to browse pileup fle---------------------------------------------------------------------
def browse_pileup(in_pileup, position):
    """
    This function allows to browse a pileup file to extract
    the line whose contains the desired position.

    PARAMETERS:
    in_pileup: a TextIOWrapper object corresponding to the pileup file
    position: a integer giving the desired position

    RETURN:
    line: the pileup whise fit with the desired position
    """
    #Looking for the previous position
    read_previous_line(in_pileup)

    #extract the position of the current line
    line = in_pileup.readline()
    #extract position
    pos_pileup = int(line.split("\t")[1])
    while pos_pileup != position:
        line = in_pileup.readline()
        if len(line) != 0:
            pos_pileup = int(line.split("\t")[1])
            if pos_pileup > position:
                return None
        else:
            return None
    return line

#Function to parse pileup line=====================================================================

#Function to clean unusefull flag------------------------------------------------------------------
def clean_flag(pileup_line):
    """
    This function allows to clean a pileup in order to remove character indicating start,
    quality score and end read position.
    
    PARAMETERS:
    pileup: pileup content from one position
    flag_regex: a regex to use to remove specific flag

    RETURN:
    pileup: pileup cleaned
    """
    start_read = re.compile(r"\^[ -~]{1}")
    end_read = re.compile(r"\$")

    for regex in [start_read, end_read]:
        if regex.search(pileup_line):
            pileup_line = "".join(regex.split(pileup_line))
    return pileup_line

#Function to parse indel---------------------------------------------------------------------------
def parse_insertion(pileup):
    """
    This function allows to parse insertion from a pileup
    line

    PARAMETERS:
    pileup: a pileup line

    RETURN:
    insert_count: the count of detected insertion
    pileup: the cleaned pileup line
    """
    #Extract the different type of detected insertion
    insertion_init = re.compile(r"\+[0-9]{1,2}")
    insert_list = insertion_init.findall(pileup)
    insert_list = list(set(insert_list))
    #Define an counter
    insert_count = 0
    #Make loop for each type of insertion
    if len(insert_list) > 0:
        for insert in insert_list:
            #Create specific regex
            insertion_regex = "[\\.,NATGCatgc\\*#]{1}\\" + insert + "[NATGCatgc\\*]{" + insert[1] + "}"
            #Catch insertions
            insertion_regex = re.compile(insertion_regex)
            insertions = insertion_regex.findall(pileup)
            #Add to counter the insertions
            insert_count += len(insertions)
            #Deduplicate list
            insertions = list(set(insertions))
            #clean pielup with found insertions
            for insertion in insertions:
                pileup = "".join(pileup.split(insertion))
    return insert_count, pileup

def parse_deletion(pileup):
    """
    This function allows to parse deletion from a pileup
    line

    PARAMETERS:
    pileup: a pileup line

    RETURN:
    delet_count: the count of deletion
    pileup: the cleaned pileup line
    """
    #Define a deletion counter
    delet_count = 0
    #Parse individual deletion
    for delet in ["\\*/#", "\\*"]: 
        deletion_regex = re.compile(delet)
        deletions = deletion_regex.findall(pileup)
        delet_count += len(deletions)
        if len(deletions) > 0:
            pileup = "".join(pileup.split(deletions[0]))
    #Parse multiple deletions
    deletion_init = re.compile("-[0-9]{1,2}")
    delet_list = deletion_init.findall(pileup)
    delet_list = list(set(delet_list))
    #Make loop for each type of insertion
    if len(delet_list) > 0:
        for delet in delet_list:
            #Create specific regex
            deletion_regex = "[\\.,NATGCatgc\\*#]{1}" + delet + "[NATGCatgc\\*]{" + delet[1] + "}"
            #Catch insertions
            deletion_regex = re.compile(deletion_regex)
            deletions = deletion_regex.findall(pileup)
            #Add to counter the insertions
            delet_count += len(deletions)
            #Deduplicate list
            deletions = list(set(deletions))
            #clean pielup with found insertions
            for deletion in deletions:
                pileup = "".join(pileup.split(deletion))
    return delet_count, pileup

#Function to parse pileup line---------------------------------------------------------------------
def parse_pileup_line(pileup_line):
    """
    This function allows to extract information from a
    pileup line

    PARAMETERS:
    pileup_line: a pileup line to parse

    RETURN:
    pileup_info a dictionnary contaning all extracted information
    """
    #Create an empty list that will contain all extracted information
    pileup_info = []
    #Split pileup line at each tab
    pileup_line_split = pileup_line.split("\t")
    #Extract position, reference base and depth
    pos = int(pileup_line_split[1])
    base_ref = pileup_line_split[2]
    depth = int(pileup_line_split[3])
    pileup_info.extend([("Pos", pos), ("Ref", base_ref), ("Depth", depth)])

    #Extract pileup alignment information
    align = pileup_line_split[4]

    #Remove flag for start and end aligned reads from pileup
    align = clean_flag(align)

    #Analyse CNV
    ##Insertion
    insert_count, align = parse_insertion(align)
    insert_depth = round(insert_count / depth, ndigits = 3) if depth != 0 else 0
    ##Deletion
    del_count, align = parse_deletion(align)
    del_depth = round(del_count / depth, ndigits = 3) if depth != 0 else 0
    pileup_info.extend([("Ins", insert_count), ("InsD", insert_depth), ("Del", del_count), ("DelD", del_depth)])

    #Analyse SNP
    for base in ["A", "T", "G", "C"]:
        forward_count = align.count(base)
        reverse_count = align.count(base.lower())
        base_depth = round((forward_count + reverse_count) / depth, ndigits = 3) if depth != 0 else 0
        pileup_info.extend([(base + "F", forward_count), (base + "R", reverse_count), (base + "D", base_depth)])
    
    #Analyse match
    forward_match_count = align.count(".")
    reverse_match_count = align.count(",")
    match_depth = round((forward_match_count + reverse_match_count) / depth, ndigits = 3) if depth != 0 else 0
    pileup_info.extend([("MatF", forward_match_count), ("MatR", reverse_match_count), ("MatD", match_depth)])

    #Analyse Null base
    null_base = align.count("<") + align.count(">")
    null_base_depth = round(null_base / depth, ndigits = 3) if depth != 0 else 0
    pileup_info.extend([("N", null_base), ("ND", null_base_depth)])

    #Convert pileup info into a dict
    pileup_info = {key: value for key, value in pileup_info}
    return pileup_info

#Functions to ckeck variant========================================================================

##Function to check CNV----------------------------------------------------------------------------
def check_CNV(cnv, depth, pileup_line, min_depth, depth_threshold):
    """
    This function allows to check a given CNV

    PARAMETERS:
    cnv: the CNV to parse
    depth: the reported total depth at this position
    pileup_line: the pileup line at this position
    min_depth: minimum total depth cutoff
    depth_threshold: relative depth part of the CNV (%)

    RETURN:
    a boolean indicating if the variant pass the control
    """
    #Check global depth
    if depth < min_depth:
        return False
    else:
        #Create regex with CNV
        cnv_regex = re.compile(r"([\.,NATGCatgc\*#]{1}" + cnv + ")")
        #Compute CNV depth
        depth_CNV = len(cnv_regex.findall(pileup_line)) / depth
        if depth_CNV < depth_threshold:
            return False
        else:
            return True

##Function to check SNP----------------------------------------------------------------------------
def check_SNP(depth, snp_depth, min_depth, depth_threshold):
    """
    This function allows to check a given SNP

    PARAMETERS:
    depth: the reported total depth at this position
    snp_depth: the snp depth
    min_depth: minimum total depth cutoff
    depth_threshold: relative depth part of the CNV (%)

    RETURN:
    a boolean indicating if the variant pass the control
    """
    #Check global depth
    if depth < min_depth:
        return False
    else:
        #Compute SNP depth
        if snp_depth < depth_threshold:
            return False
        else:
            return True

#Define the core function that parse a VCF file====================================================
def filter_vcf(vcf_path, pileup_path, variant_threshold, depth_threshold):
    """
    This is the core function of this program.

    PARAMETERS:
    vcf_path: a path to a VCF file
    pileup_path: a path to a pileup file

    RETURN:
    pileup_array: a completed pileup array with variant information
    """
    with open(vcf_path, "r") as in_vcf, open(pileup_path, "r") as in_pileup:
        vcf_reader = vcf.Reader(in_vcf)
        for record in vcf_reader:
            pos = record.POS
            ref = record.REF
            alts = record.ALT
            filter = record.FILTER
            if (filter is None or len(filter) == 0):
                #PASS filter or filter status is Unknow
                for alt in alts:
                    #Extract pileup line and corresponding informations
                    pileup_line = browse_pileup(in_pileup, pos)
                    pileup_dict = parse_pileup_line(pileup_line)

                    if alt is None:
                        print(alt)
                        print("Ok is None")
                        break

                    #Check alternative base
                    if len(alt) < len(ref):
                        #Deletion case (CNV)
                        print(f"MESSAGE: Parsing CNV at position {pos}")
                        deletion = define_deletion(ref, alt)
                        #Check variant
                        variant_pass = check_CNV(deletion, pileup_dict["Depth"], pileup_line, depth_threshold, variant_threshold)
                    elif len(alt) > len(ref):
                        #Insertion case (CNV)
                        print(f"MESSAGE: Parsing CNV at position {pos}")
                        insertion = define_insertion(alt)
                        #Check variant
                        variant_pass = check_CNV(insertion, pileup_dict["Depth"], pileup_line, depth_threshold, variant_threshold)
                    else:
                        #Case for SNP
                        #Check variant
                        if len(str(ref)) > 1 and len(str(alt)) > 1:
                            print(f"MESSAGE: Parsing multiple SNP at position {pos}")
                            #Check when there is a multiple SNP
                            variant_pass_list = []
                            for sub_alt in str(alt):
                                result = check_SNP(pileup_dict["Depth"], pileup_dict[str(sub_alt) + "D"] ,depth_threshold, variant_threshold)
                                variant_pass_list.append(result)
                            #Analyse all flags
                            variant_pass = True if all(variant_pass_list) == True else False
                        else:
                            #Single SNP
                            print(f"MESSAGE: Parsing SNP at position {pos}")
                            variant_pass = check_SNP(pileup_dict["Depth"], pileup_dict[str(alt) + "D"] ,depth_threshold, variant_threshold)
                    if variant_pass:
                        #Complete the array
                        if not 'pileup_array' in locals():
                            pileup_array = np.array(list(pileup_dict.values()))
                        else:
                            pileup_array = np.vstack((pileup_array, list(pileup_dict.values())))
            else:
                print(f"MESSAGE: Skip position {pos} due to unpassed quality test")
    #Create an empty pileup array after the loop if none variant had been retained
    if not 'pileup_array' in locals():
        pileup_array = np.array([])
        print("WARNING: None variant had been retained during the filtration process.")
        print("It may be possible that the sequencing quality is too low.")
    return(pileup_array)

#Function to write results=========================================================================

#Define a function to write the filtered VCF file--------------------------------------------------
def write_vcf(in_vcf, out_vcf, pileup_array):
    """
    This function allows to write a filtered VCF according to
    pileup array

    PARAMETERS:
    in_vcf: input VCF file
    out_vcf: output VCF file
    pileup_array: the pileup array

    RETURN:
    None
    """
    #Case when pileup array is not empty
    if pileup_array.shape[0] != 0:    
        #Extract the first column from pileup_array = position list
        positions = [row[0] for row in pileup_array]
    
        #Loop over files
        with open(in_vcf, "r") as vcf_input, open(out_vcf, "w") as vcf_output:
            for line in vcf_input:
                if line.startswith("#"):
                    vcf_output.write(line)
                else:
                    #extract position
                    position = line.split("\t")[1]
                    if position in positions:
                        vcf_output.write(line)
    else:
        #In this case, only write the header
        with open(in_vcf, "r") as vcf_input, open(out_vcf, "w") as vcf_output:
            for line in vcf_input:
                if line.startswith("#"):
                    vcf_output.write(line)

#Define function to write pileup array-------------------------------------------------------------
def write_pileup_array(pileup_array):
    """
    This function allows to write a pileup array

    PARAMATERS:
    pileup_array: the pileup array to write

    RETURN:
    None

    #List of created attribut
    Pos: Position
    Ref: Reference base
    Depth: Depth
    Ins: Insertion count
    InsD: Insertion depth (%)
    Del: Deletion count
    DelD: Deletion depth
    AF: A base forward count
    AR: A base reverse count
    AD: A base deletion count
    TF: T base forward count
    TR: T base reverse count
    TD: T base depth count
    GF: G base forward count
    GR: G base reverse count
    GD: G base depth
    CF: C base forward count
    CR: C base reverse count
    CD: C base depth
    MatF: Match count forward
    MatR: Match reverse count
    MatD: Match depth
    N: Null base count
    ND: null base depth
    """
    with open("pileup_pattern.txt", "w") as out_pileup:
        #Create header
        out_pileup.write("Pos\tRef\tDepth\tIns\tInsD\tDel\tDelD\tAF\tAR\tAD\tTF\tTR\tTD\tGF\tGR\tGD\tCF\tCR\tCD\tMatF\tMatR\tMatD\tN\tND\n")
        for row in pileup_array:
            out_pileup.write("\t".join(row) + "\n")

#Main propgram
if __name__ == "__main__":

    #extract command line arguments
    args = extract_argument()

    #Prase argument
    parse_arguments(args)

    #Use samtools to prepare BAM file for pileup
    if args.sam is None and args.bam is not None:
        #case where only a BAm file is gave
        ##Sort bam file
        print("MESSAGE: Sorting BAM file", end = "\n")
        bam_sort = sort_bam(args.bam, args.cpu)
        ##Index bam file
        print("MESSAGE: Indexing BAM file", end = "\n")
        bam_index = index_bam(bam_sort, args.cpu)

    elif args.sam is not None and args.bam is not None:
        #case where both are povided
        print("MESSAGE: Both SAM and BAM arguments are provided. Priority to BAM file", 
              end = "\n")
        ##sort bam file
        print("MESSAGE: Sorting BAM file", end = "\n")
        bam_sort = sort_bam(args.bam, args.cpu)
        ##Index Bam file
        print("MESSAGE: Indexing BAM file", end = "\n")
        bam_index = index_bam(bam_sort, args.cpu)

    else:
        #case where only SAM file is provided
        ##convert sam in bam file
        print("MESSAGE: Converting SAM input in BAM format", end = "\n")
        bam = convert_sam(args.sam, args.cpu)
        ##Sort bam file
        print("MESSAGE: Sorting BAM file", end = "\n")
        bam_sort = sort_bam(bam, args.cpu)
        ##Index BAM file
        bam_index = index_bam(bam, args.cpu)

    #Perform pileup
    ##Check if pileup already exist
    outputDirectory= re.sub("/[A-z]*\.[A-z]*$","/", args.output)
    pileup_in = outputDirectory + bam_sort.split(".bam")[0] + ".pileup"
    if not os.path.exists(pileup_in):
        print("MESSAGE: Creating pileup file", end = "\n")
        pileup_reads(bam_sort, args.ref)
    #Sort VCF file
    if os.path.exists("tempo.vcf"):
        os.remove("tempo.vcf")
    sort_vcf(args.vcf, "tempo.vcf")
    pileup_array = filter_vcf("tempo.vcf", pileup_in, args.threshold, args.min_depth)

    write_pileup_array(pileup_array)

    write_vcf("tempo.vcf", args.output, pileup_array)

    #Remove files
    if args.remove:
        os.remove("tempo.vcf")
        os.remove(pileup_in)
        if args.bam is None:
            os.remove(bam_sort)
            os.remove(bam_index)
        else:
            os.remove(bam_index)
    