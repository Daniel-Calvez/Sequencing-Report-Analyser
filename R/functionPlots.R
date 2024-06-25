## Functions to get data -------------------------------------------------------

getAllSamParseTables <- function(Output_Directory){
  #' This function get all data tables produced by the script
  #' @param Output_Directory The directory containing the data from the script
  #' @return A list of all the data tables

  # Get all the data into matching data tables
  alignementMapScore <- read.table(paste0(Output_Directory, "/python_parse_sam_results/alignment_map_score.txt"), header = T)
  alignementMap <- read.table(paste0(Output_Directory, "/python_parse_sam_results/alignment_map.txt"), header =  T, sep = "\t")
  alignementReadLength <- read.table(paste0(Output_Directory, "/python_parse_sam_results/alignment_read_length.txt"), header = T)
  alignementSeqScore <- read.table(paste0(Output_Directory, "/python_parse_sam_results/alignment_seq_score.txt"), header = T)
  
  
  return(list(alignementMap, alignementMapScore, alignementReadLength, alignementSeqScore))
}

getAllFastqParseTables <- function(Output_Directory){
  #' This function get all the fastq analysis results and returns it as tables
  #' @param Output_Directory The main output directory
  #' @return A list of tables
  
  
  #List all files associated to the fastq analysis
  pathFastq <- paste0(Output_Directory, "/python_parse_fastq_results", collapse = "")
  pathList <- list.files(pathFastq, full.names = TRUE)
  #Read only the files that are interesting
  gcTable <- read.table(pathList[[grep("GC_content",pathList)]], header = TRUE)
  lengthTable <- read.table(pathList[[grep("read_length",pathList)]], header = TRUE)
  avgScoreTable <- read.table(pathList[[grep("seq_score",pathList)]], header = TRUE)
  
  #Return all files
  return(list(gcTable, lengthTable, avgScoreTable))
}

## Functions to create plots from ----------------------------------------------

plotsSequencingSample <- function(readQuality, readLength){
  #' This function get the data tables about sequencing quality and length
  #' to generate bar plots from them
  #' @param readQuality dataframe containing the quality of each read
  #' @param readLength dataframe containing the length of each read
  #' @return An arrangment of the two bar plots
  
  # Create the mean value of the quality score of sequences 
  meanValue <- weighted.mean(readQuality$Seq_score, readQuality$occurence)
  medianValue <- median(readQuality$Seq_score)
  
  # Build the first plot
  p1 <-
    ggplot(readQuality) +
    # Create a bar chart
    geom_bar(aes(x = Seq_score, y = occurence), stat = "identity", fill = "skyblue3", position = "dodge") +
    # Add the mean and median
    geom_vline(aes(xintercept = weighted.mean(Seq_score, occurence), colour = "Moyenne"), linetype = "dashed") +
    geom_vline(aes(xintercept = weighted.median(Seq_score, occurence), colour = "Médiane"), linetype = "dashed") +
    labs(colour = "")+
    # Set the origin at 0,0
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
    # Add the graph and axis legend
    labs(x = "Score qualité Phred",
         y = "Occurence",
         title = "Scores qualités des reads", 
         subtitle = paste0("Moyenne: ", trunc(meanValue*10^2)/10^2, "\n",
                           "Mediane: ", trunc(medianValue*10^2)/10^2)) +
    theme_classic()+
    theme(legend.position = "bottom")
  
  
  # Create the mean value of the quality score of sequences 
  meanValue <- weighted.mean(readLength$read_length,w = readLength$occurence)
  medianValue <- weighted.median(readLength$read_length, w =  readLength$occurence)
  
  # Build the second plot
  p2 <-
    ggplot(readLength) +
    aes(x = read_length, y = occurence) +
    # Create a bar chart
    geom_bar(stat = "identity", fill = "skyblue3", position = "dodge") +
    # Add the mean and median
    geom_vline(aes(xintercept = weighted.mean(read_length,w = occurence), colour = "Moyenne"), linetype = "dashed") +
    geom_vline(aes(xintercept = weighted.median(read_length, w = occurence), colour = "Médiane"), linetype = "dashed") +
    labs(colour = "")+
    # Set the origin at 0,0
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
    # Add the graph and axis legend
    labs(x = "Longueur des reads (b)",
         y = "Occurence",
         title = "Longueurs des reads",
         subtitle = paste0("Moyenne: ", trunc(meanValue*10^2)/10^2, "\n",
                           "Mediane: ", trunc(medianValue*10^2)/10^2)) +
    theme_classic()+
    theme(legend.position = "bottom")
  
  return(list(p1, p2))
}

## Functions to create plots from mapped reads ---------------------------------

plotsAlignmentMap <- function(alignmentMap){
  #' This function get the data tables about the read mapping
  #' and returns a chart of it
  #' @param alignmentMap dataframe containing the read mappings
  #' @return An arrangment of the two bar plots
  
  
  # Modification on the data to make it more easily readable
  alignmentMap <-    
    alignmentMap %>% 
    mutate(Group = "a") %>% 
    relocate(Group, .before = read) %>% 
    filter(read != "total") %>% 
    mutate(read = factor(read, levels = c("Mapped", "No mapped"), labels = c("Mappés", "Non mappés"))) %>% 
    mutate(percent = round(occurence / sum(occurence) * 100, digits = 2)) %>% 
    mutate(y_lab =  sum(occurence) - cumsum(occurence) + occurence / 2)
  
  
  ###Data vizualisation
  
  ##Map data
  # Build the plot
  p1 <-
    ggplot(alignmentMap) +
    # Index the columns on the axis
    aes(x = Group, y = occurence, fill = read) +
    geom_bar(stat = "identity") +
    coord_flip() +
    geom_text(aes(y = y_lab, label = paste0(percent, "%"))) +
    # Separate mapped from unmapped values
    scale_fill_manual(values = c("Mappés" = "palegreen4", "Non mappés" = "brown3")) +
    # Legends the graph
    labs(x = "", y = "",
         title = "Reads mappés",
         fill = "") +
    # Visual modifications (so it can look better)
    theme(axis.text =  element_blank(),
          axis.ticks = element_blank(),
          panel.background = element_rect(fill = 'white', color = 'white'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top")
  
  return(p1)
}

plotReadMapQuality <- function(readMapQuality){
  #' This function get the data tables about the reads quality
  #' and returns a chart of it
  #' @param alignmentMap dataframe containing the read mappings
  #' @return An arrangment of the two bar plots
  
  
  # Create the mean and median value of the quality score of reads weighted by their occurence
  meanValue <- weighted.mean(readMapQuality$Map_score,w = readMapQuality$occurence)
  medianValue <- weighted.median(readMapQuality$Map_score,w = readMapQuality$occurence)
  p2 <-
    ggplot(readMapQuality) +
    # Associate axis to the right columns and build a bar chart
    aes(x = Map_score, y =occurence)+
    geom_bar(stat = "identity", fill = "skyblue3", position = "dodge") +
    # Draw a vertical line from mean and median 
    geom_vline(aes(xintercept = weighted.mean(Map_score, occurence), colour = "Moyenne"), linetype = "dashed") +
    geom_vline(aes(xintercept = weighted.median(Map_score, occurence), colour = "Médiane"), linetype = "dashed") +
    labs(colour = "")+
    # Set the origin at 0,0
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
    labs(x = "Score qualité Phred",
         y = "Occurence",
         title = "Scores qualités des reads",
         subtitle = paste0("Moyenne: ", trunc(meanValue*10^2)/10^2, "\n",
                           "Mediane: ", trunc(medianValue*10^2)/10^2)) +
    theme_classic()+
    theme(legend.position = "bottom")

  return(p2)
}

plotReadMapLength <- function(readMapLength){
  #' This function get the data tables about the read mapping length
  #' and returns a chart of it
  #' @param alignmentMap dataframe containing the read mappings
  #' @return An arrangment of the two bar plots

  # Store the weighted mean and median values
  meanValue <- weighted.mean(readMapLength$read_length,w = readMapLength$occurence)
  medianValue <- weighted.median(readMapLength$read_length, w =  readMapLength$occurence)
  
  # Build the plot
  p3 <-
    ggplot(readMapLength) +
    # Define axis and add the values to a bar chart
    aes(x = read_length, y = occurence, fill= as.character(mapped)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.75) +
    scale_fill_manual(values = c("0" = "firebrick3", "1" = "royalblue2"))+
    # Draw a vertical line from mean and median 
    geom_vline(aes(xintercept = weighted.mean(read_length,w = occurence), colour = "Moyenne"), linetype = "dashed") +
    geom_vline(aes(xintercept = weighted.median(read_length, w = occurence), colour = "Médiane"), linetype = "dashed") +
    labs(colour = "")+
    # Set the origin at 0,0
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
    # Finalize graph by adding a legend for it
    labs(x = "Longueurs des reads (b)",
         y = "Occurence",
         fill = "Mapped:",
         title = "Longueurs des reads",
         subtitle = paste0("Moyenne: ", trunc(meanValue*10^2)/10^2, "\n",
                           "Mediane: ", trunc(medianValue*10^2)/10^2)) +
    theme_classic()+
    theme(legend.position = "bottom")

  return(p3)
}

plotReadMapSeqScore <- function(readMapSeqScore){
  #' Creates a bar plot from the mapped quality score of sequences
  #' @param readMapSeqScore A data table containing the quality score along occurences
  #' @return The bar plot created

  # Store the mean and median values
  meanValue <- weighted.mean(readMapSeqScore$Seq_score,w = readMapSeqScore$occurence)
  medianValue <- weighted.median(readMapSeqScore$Seq_score, w =  readMapSeqScore$occurence)
  
  # Build the plot
  p4 <-
    ggplot(readMapSeqScore) +
    
    # Defines the dataframe columns as the X and Y axis and create a bar chart from the data
    aes(x = Seq_score, y = occurence) +
    geom_bar(stat = "identity", fill = "skyblue3", position = "dodge") +
    # Draw a vertical line from mean and median 
    geom_vline(aes(xintercept = weighted.mean(Seq_score,w = occurence), colour = "Moyenne"), linetype = "dashed") +
    geom_vline(aes(xintercept = weighted.median(Seq_score, w = occurence), colour = "Médiane"), linetype = "dashed") +
    labs(colour = "")+
    # Set the y axis at 0
    scale_y_continuous(expand = c(0, 0))+
    # Legend the plot and change visual details
    labs(x = "Score d'alignement des séquences mappées (b)",
         y = "Occurence",
         title = "Score d'alignement des séquences mappées",
         subtitle = paste0("Moyenne: ", trunc(meanValue*10^2)/10^2, "\n",
                           "Mediane: ", trunc(medianValue*10^2)/10^2)) +
    theme_classic()+
    theme(legend.position = "bottom")

  return(p4)
}

## Functions to create plots from coverage data --------------------------------

plotCoverage <- function(coverage){
  #' A function creating a plot showing the coverage of each sequences
  #' @param coverage A dataframe describing the coverage of each sequence
  #' @return The plot from the data frame coverage
  
  # Remove uncovered zones (Analyzed separately)
  coverage <- coverage %>% filter(depth!=0)
  
  # Build the plot
  p <- ggplot(data = coverage) + 
    # Define axis from columns and build a rectangle chart from it
    aes(x = start, y = depth) + 
    geom_rect(mapping = aes(xmin = start, ymin = 0, xmax = stop, ymax = depth), color = "lightblue3", fill = "lightblue3") +
    geom_smooth()+
    # Draw a horizontal line from mean and median 
    geom_hline(aes(yintercept = median(depth), col = "Médiane"), linetype = "dashed") +
    geom_hline(aes(yintercept = mean(depth), col = "Moyenne"), linetype = "dashed") +
    # Set the origin at 0,0
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
    # Details the legend
    labs(title = "Profondeur de séquençage",
         colour ="",
         subtitle = paste0("Moyenne : ",trunc(coverage$depth*10^2)/10^2),
         x = "Position",
         y = "Profondeur (X)") +
    theme_classic()
  return(p)
}

unknownRegions <- function(consensusSequence){
  #' This functions creates a plot that show which region has coverage
  #' @param consensusSequence The consensus sequence
  #' @return A rectangular plot showing the coverage of the reference
  
  # Add the locations of N into a vector
  location <- as.data.frame(str_locate_all(consensusSequence, "n"))$start
  #Initialize a dataframe containing the start and end of N
  unknow_position <- data.frame(start = NA, stop = NA)
  start <- location[1]
  # Parse the N positions
  for (ligne in 1:(length(location) - 1))
  {
    # If the N chains ends add it to the dataframe and reloop
    if (location[ligne] + 1 != location[ligne + 1])
    {
      unknow_position[nrow(unknow_position)+1,] <- c(start, ligne)
      start <- location[ligne+1]
    }
  }
  # If the N chain had no gap
  if(!is.empty(location))unknow_position[nrow(unknow_position)+1,] <-c(start, tail(location, n=1))
  unknow_position <- unknow_position[-1,]
  percent_unknow <- round(str_count(consensusSequence, "n") / nchar(consensusSequence) * 100, digits = 1 )
  
  # Build the plot
  plot <- ggplot() + 
    # Create a rectangle chart
    geom_rect(mapping = aes(xmin = 1, xmax = nchar(consensusSequence),
                            ymin = -1, ymax = 1),
              fill = "white", color = "black") +
    # Add red parts to the chart where there are unknown regions
    geom_rect(data = unknow_position,
              mapping = aes(xmin = start, xmax = stop,
                            ymin = -0.99, ymax = 0.99), fill = "brown3") +
    # Add a legend
    labs(x = "Position",
         title = "Régions inconnues",
         subtitle = paste0("Non déterminé: ", percent_unknow, "%")) +
    # Visual modifications
    theme(axis.text.y =  element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = 'white', color = 'white'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  return(plot)
}

## Functions to create plots and tables from GC content ------------------------

plotGC <- function(gcTable, reference){
  #' Creates a bar plot from the gc scores along with a normal law based on reference
  #' @param gcTable A data table containing the quality score along occurences
  #' @return The bar plot created
  
  normalLawValues <- normalLawGC(reference)
  
  # Store the mean, median and standard deviation
  meanValue <- weighted.mean(gcTable$GC_content,w = gcTable$occurence)
  medianValue <- weighted.median(gcTable$GC_content, w =  gcTable$occurence)
  std_deviation <- sd(gcTable$GC_content)
  # Build the plot
  p <-
    ggplot(gcTable) +
    # Defines the dataframe columns as the X and Y axis and create a bar chart from the data
    aes(x = GC_content, y = occurence) +
    geom_bar(stat = "identity", fill = "skyblue4", position = "dodge") +
    # Draw a vertical line from mean and median 
    geom_vline(aes(xintercept = weighted.mean(GC_content,w = occurence), colour = "Moyenne"), linetype = "dashed") +
    geom_vline(aes(xintercept = weighted.median(GC_content, w = occurence), colour = "Médiane"), linetype = "dashed") +
    
    labs(colour = "")+
    # Set the origin at 0,0 and turn y into a logarithmic scale
    scale_x_continuous(expand = c(0, 0)) + scale_y_log10(expand= c(0,0))+
    # Legend the plot and change visual details
    labs(x = "Taux de GC",
         y = "log10 Occurence",
         title = "Score d'alignement de séquences",
         subtitle = paste0("Moyenne: ", trunc(meanValue*10^2)/10^2, "\n",
                           "Mediane: ", trunc(medianValue*10^2)/10^2,"\n")) +
    theme_classic()+
    theme(legend.position = "bottom")
  
  #Create the table to compare mean, median and standard deviation
  dataComparisonTable <- data.frame(matrix(nrow = 0, ncol = 2))
  dataComparisonTable[1,] <- c(trunc(meanValue*10^2)/10^2, trunc(normalLawValues[1]*10^2)/10^2)
  dataComparisonTable[2,] <- c(trunc(medianValue*10^2)/10^2, trunc(normalLawValues[3]*10^2)/10^2)
  dataComparisonTable[3,] <- c(trunc(std_deviation*10^2)/10^2, trunc(normalLawValues[2]*10^2)/10^2)
  colnames(dataComparisonTable) <- c("Reads","Reference")
  rownames(dataComparisonTable) <- c("Moyenne", "Mediane", "Ecart-type")
  
  #Create the Kolmogorov-Smirnov test and student test
  ksTest <- ks.test(x=rep(gcTable$GC_content, gcTable$occurence), "pnorm", normalLawValues[1], normalLawValues[2])
  studentTest <- t.test(rep(gcTable$GC_content, gcTable$occurence), mu = normalLawValues[1])
  
  #Create a table with the Kolmogorov-Smirnov test results
  tableKsTest <- data.frame(matrix(nrow = 1, ncol = 2))
  tableKsTest[1,] <- c(ksTest$statistic, ksTest$p.value)
  colnames(tableKsTest) <- c("D", "p.value")
  
  #Show test results from Student (separating the confidence interval from the other)
  studentTestStats <- data.frame(matrix(nrow = 1, ncol = 3))
  studentTestStats[1,] <- c(studentTest$statistic, studentTest$parameter, studentTest$p.value)
  colnames(studentTestStats) <- c("t", "degrés de libertés","p.value")
  
  studentTestConfidence <- data.frame(matrix(nrow = 1, ncol = 3))
  studentTestConfidence[1,] <- c(studentTest$conf.int, studentTest$null.value)
  colnames(studentTestConfidence) <- c("min", "max", "Moyenne de référence")
  
  return(list(p, dataComparisonTable, tableKsTest, studentTestStats, studentTestConfidence))
}

normalLawGC <- function(reference){
  #' This function takes a reference sequence to return the mean and standard
  #' deviation of the GC scores
  #' @param reference The path to the reference file
  #' @return A vector of number containing the mean and standard deviation
  
  # Read the file and get the mean
  refSequence <- read.fasta(reference, as.string = TRUE, seqonly = TRUE)[[1]]
  mean <- getMeanGC(refSequence)
  # Split the sequence every 30 nucleotides and get the gc scores of each subsequences
  splittedList <- str_split_1(gsub("(.{30})", "\\1 ", refSequence), " ")
  analyzedList <- unlist(lapply(splittedList, getMeanGC))
  # Use the subsequences to compute the standard deviation
  std_deviation <- sqrt(mean(analyzedList^2)-mean(analyzedList)^2)
  median <- median(analyzedList)
  return(c(mean, std_deviation,median))
}

getMeanGC <- function(sequence){
  #' This function get the meanGC of a sequence
  #' @param sequence The gene sequence
  #' @return The gc score of the sequence
  moyenne <- (str_count(sequence, "G") + str_count(sequence, "C"))/str_length(sequence)
  return(moyenne)
}

# Function to create plots and tables from variant calling ---------------------
vcfAnalysis <- function(variants, reference){
  
  # Creation of the table containing all the raw readable data
  info <- unlist(str_split(unique(variants$V9), ":"))
  infoTable <- data.frame(matrix(nrow = 0, ncol = length(info) +4))
  colnames(infoTable) <- c("Gene","position", "reference", "mutation", info)
  infoList <- str_split(variants$V10, ":")
  for(i in 1:length(infoList)){
    infoTable[nrow(infoTable)+1,] <- c(variants$V1[i], variants$V2[i], variants$V4[i], variants$V5[i],infoList[[i]])
  }
  
  
  # Create the table containing all mutations with their position
  mutationTable <- variants[,-c(1,3,6,7,8,9,10)]
  colnames(mutationTable) <- c("Position", "Reference", "Mutated")
  referenceLength <- str_length(read.fasta(reference,as.string = TRUE, seqonly = TRUE))
  
  # Build the plot
  plot <- ggplot() + 
    # Create a rectangle chart
    geom_rect(mapping = aes(xmin = 1, xmax = referenceLength,
                            ymin = -1, ymax = 1),
              fill = "white", color = "black") +
    # Add red segment where there is a mutation
    geom_segment(aes(x = mutationTable$Position, xend = mutationTable$Position, y = -1, yend = 1), col = "red3", alpha=0.65)+
    # Add a legend
    labs(x = "Position",
         title = "Position des mutations") +
    # Visual modifications
    theme(axis.text.y =  element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = 'white', color = 'white'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  
  
  return(list(infoTable,plot))
}