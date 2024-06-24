library(argparse)

#Argparse parameters setup------------------------------------------------------
#Create the parser along with its arguments
parser <- ArgumentParser()
parser$add_argument("-v", "--verbose", action="store_true", default=FALSE,
                    help="Print extra output [default]")
parser$add_argument("-s", "--source", default ="//", type="character",
                    help="The path of the output file")
parser$add_argument("-n", "--name", default = "Rapport",type="character",
                    help="The name of the pdf file")
parser$add_argument("-p", "--path", default = "//",type="character",
                    help="The rmd file path")
parser$add_argument("-r", "--reference", default = "//",type="character",
                    help="The reference fasta file path")
parser$add_argument("-w", "--webpath", default = "//",type="character",
                    help="The rmd html file path")

#Store the args in seprarate variables
args <- parser$parse_args()
source <- args$source
name <- args$name
rmd_path <- args$path
rmd_path2 <- args$webpath
reference <- args$reference

#If no source was giver
if(source == "//"){
  stop("Source argument (-s) missing")
}

#If no rmd file path was given
if(rmd_path == "//"){
  stop("RMD path argument (-p) missing")
}

#Remove all "/" remaining at the end of path
if(substr(source, nchar(source), nchar(source))=="") source <- substr(source, 1, nchar(source)-1)
if(substr(source, 1, 1)!="/") source <- paste0(getwd(), "/", source,collapse = "")
if(substr(reference, nchar(reference), nchar(reference))=="") reference <- substr(reference, 1, nchar(reference)-1)
if(substr(reference, 1, 1)!="/") reference <- paste0(getwd(), "/", reference,collapse = "")

#Generate the report from parameters
rmarkdown::render(rmd_path,
                    params = list(OutputDirectory = source, Reference = reference),
                    output_file = paste0(name, ".pdf"),
                    output_dir = source)
if(rmd_path2 != "//"){
  rmarkdown::render(rmd_path2,
                  params = list(OutputDirectory = source, Reference = reference),
                  output_file = paste0(name, ".html"),
                  output_dir = source)
}