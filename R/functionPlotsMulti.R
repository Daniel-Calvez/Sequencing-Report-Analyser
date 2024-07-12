addPatientID <- function(dataTables, patientName){
  #' This function add a patient ID to each data table
  #' @param dataTables A list of data frames
  #' @param patientName the patient name/id
  #' @return The modified data frames
  
  for(i in 1:length(dataTables)){
    dataTables[[i]] <- dataTables[[i]] %>% mutate(patient = patientName)
  }
  
  return(dataTables)
}



combineFastqTables <- function(listDataTable){
  #' This function aggregate every datatables from fasta information table list 
  #' and do the sum of matching occurences
  #' @param listDataTable the fasta information table list
  #' @return The list of aggregated tables as a data frame
  resTable <- list(data.frame(GC_content = integer(), occurence = integer(), patient = character()),
                   data.frame(read_length = integer(), occurence = integer(), patient = character()),
                   data.frame(Seq_score = integer(), occurence = integer(), patient = character())) 
  
  
  for(dataTable1 in listDataTable){
  # Add the second table to the first
    resTable[[1]] <- bind_rows(resTable[[1]], dataTable1[[1]])
    resTable[[2]] <- bind_rows(resTable[[2]], dataTable1[[2]])
    resTable[[3]] <- bind_rows(resTable[[3]], dataTable1[[3]])

    # Aggregate the occurences grouping by the other column
    resTable[[1]] <- resTable[[1]] %>% group_by(GC_content, patient) %>% summarise(occurence = sum(occurence))
    resTable[[1]] <- resTable[[1]][,colnames(dataTable1[[1]])]
    
    resTable[[2]] <- resTable[[2]] %>% group_by(read_length, patient) %>% summarise(occurence = sum(occurence))
    resTable[[2]] <- resTable[[2]][,colnames(dataTable1[[2]])]

    resTable[[3]] <- resTable[[3]] %>% group_by(Seq_score, patient) %>% summarise(occurence = sum(occurence))
    resTable[[3]] <- resTable[[3]][,colnames(dataTable1[[3]])]
  }
  return(resTable)
}

combineSamTables <- function(listDataTable){
  #' This function aggregate every datatables from SAM information table list 
  #' and do the sum of matching occurences
  #' @param listDataTable the sam information table list
  #' @return The list of the aggregated data tables as a data frame
  
  resTable <- list(data.frame(read = character(), occurence = integer(), patient = character()),
                  data.frame(Map_score = integer(), occurence = integer(), patient = character()),
                  data.frame(read_length = integer(), mapped = integer(), occurence = integer(), patient = character()),
                  data.frame(Seq_score = integer(), occurence = integer(), patient = character()))

  
  for(dataTable1 in listDataTable){
    # Add the second table to the first
    resTable[[1]] <- bind_rows(resTable[[1]], dataTable1[[1]])
    resTable[[2]] <- bind_rows(resTable[[2]], dataTable1[[2]])
    resTable[[3]] <- bind_rows(resTable[[3]], dataTable1[[3]])
    resTable[[4]] <- bind_rows(resTable[[4]], dataTable1[[4]]) 
    
    # Aggregate the occurences grouping by the other column
    resTable[[1]] <- resTable[[1]] %>% group_by(read, patient) %>% summarise(occurence = sum(occurence))
    resTable[[1]] <- resTable[[1]][,colnames(dataTable1[[1]])]
    
    resTable[[2]] <- resTable[[2]] %>% group_by(Map_score, patient) %>% summarise(occurence = sum(occurence))
    resTable[[2]] <- resTable[[2]][,colnames(dataTable1[[2]])]
    
    resTable[[3]] <- resTable[[3]] %>% group_by(read_length, mapped, patient) %>% summarise(occurence = sum(occurence))
    resTable[[3]] <- resTable[[3]][,colnames(dataTable1[[3]])]
    
    resTable[[4]] <- resTable[[4]] %>% group_by(Seq_score, patient) %>% summarise(occurence = sum(occurence))
    resTable[[4]] <- resTable[[4]][,colnames(dataTable1[[4]])]
  }

  return(resTable)
}

combineVariantTables <- function(listDataTable){
  #' This function aggregate every datatables from VCF table list 
  #' @param listDataTable the VCF table list
  #' @return The list of the aggregated data tables as a data frame

  resList <- data.frame(V1 = character(), V2 = integer(), V3 = character(),
                        V4 = character(), V5 = character(), V6 = integer(), 
                        V7 = character(), V8 = character(), V9 = character(), 
                        V10 = character(), patient = character())
  
  # Add the tables together and sort them
  for(dataTable in listDataTable){
    resList <- bind_rows(resList, dataTable)
  }
  resList <- resList[order(resList$patient,resList$V2),]
  return(resList)
}

combineCoverageTable <- function(listDataTable){
  #' This function aggregate every datatables from coverage table list 
  #' @param listDataTable the Coverage table list
  #' @return The list of the aggregated data tables as a data frame
  
  resTable <- data.frame(sample = character(), start = integer(), 
                        stop = integer(), depth = integer(), patient = character())
  globalTable <- data.frame(sample = character(), start = integer(), 
                            stop = integer(), depth = integer())
  # Add the tables together and sort them
  for(dataTable in listDataTable){
    total <- max(dataTable$stop)
    resTable <- bind_rows(resTable, dataTable) %>% mutate(percentage = (stop-start)/total)
    globalTable <- globalTable %>% bind_rows(dataTable 
                                  %>% group_by(sample, start, stop) 
                                  %>% summarise(depth = sum(depth))
                                  %>% mutate(percentage = (stop-start)/total))
  }
  return(list(resTable, globalTable))
}

fastqQualityPlot <- function(fastqGCTab, fastqLengthTab,fastqScoreTab){
  #' This generates sequence quality plots for multiple patients
  #' @param fastqGCTab Information table about GC score per patient
  #' @param fastqLengthTab Information table about read length per patient
  #' @param fastqScoreTab Information table about read quality per patient
  #' @return Plots about sequence quality
  
  # This generates the mandatory tables for creating the plots
  
  
  # This is a lollipop plot with the mean and median quality score sequence 
  # of each patient
  p1 <- ggplot(fastqGCTab, aes(x=Moyenne, y=patient))+
    # Generate the point of the "lollipop"
    geom_point()+
    # Generate the segment of the "lollipop"
    geom_segment(aes(x=Moyenne, y=patient, 
                     xend = 0, yend = patient))+
    # Create the axis legend and the graph title
    labs(x = "Taux de GC moyen",
         y = "Patient",
         title = "Taux de GC par patient")
  
  p2 <- ggplot(fastqLengthTab, aes(x=Moyenne, y=patient))+
    # Generate the point of the "lollipop"
    geom_point()+
    # Generate the segment of the "lollipop"
    geom_segment(aes(x=Moyenne, y=patient, 
                     xend = 0, yend = patient))+
    
    # Create the axis legend and the graph title
    labs(x = "Longueur moyenne",
         y = "Patient",
         title = "Longueur des reads par patient")
  
  p3 <- ggplot(fastqScoreTab, aes(x=Moyenne, y=patient))+
    # Generate the point of the "lollipop"
    geom_point()+
    # Generate the segment of the "lollipop"
    geom_segment(aes(x=Moyenne, y=patient, 
                     xend = 0, yend = patient))+
    
    # Create the axis legend and the graph title
    labs(x = "Score qualité moyen phred",
         y = "Patient",
         title = "Score qualité des reads par patient")
  
  
  return(list(p1, p2, p3))
}

samMultiPlots <- function(samMapTable, samMapScoreTab, samScoreTab){
  #' This generates plots from alignment to compare it with multiple patients
  #' @param samMapTable Information table about the percentage of mapped reads
  #' @param samMapScoreTab Information table about the alignment score
  #' @param samScoreTab Information table about the quality score of mapped reads
  #' @return A list of plots

  
  p1 <- ggplot(samMapTable) +
    # Index the columns on the axis
    aes(x = patient, y = percent, fill = read) +
    geom_bar(stat = "identity") +
    # Separate mapped from unmapped values
    scale_fill_manual(values = c("Mapped" = "palegreen4", "No mapped" = "brown3")) +
    # 
    geom_text(aes(y = (y_lab+0.1)*(dev.size("px")/6), label = paste0(percent, "%")), position = "identity") +
    coord_flip() +
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
  
  p2 <- ggplot(samMapScoreTab, aes(x=Moyenne, y=patient))+
    geom_point()+
    geom_segment(aes(x=Moyenne, y=patient, 
                     xend = 0, yend = patient))+
    labs(x = "Score d'alignement moyen",
         y = "Patient",
         title = "Score d'alignement par patient")
  
  p3 <- ggplot(samScoreTab, aes(x=Moyenne, y=patient))+
    geom_point()+
    geom_segment(aes(x=Moyenne, y=patient, 
                     xend = 0, yend = patient))+
    labs(x = "Score qualité moyen phred des reads mappés",
         y = "Patient",
         title = "Score qualité des reads mappés par patient")
  
  return(list(p1,p2,p3))
}

coverageMultiPlot <- function(patientTable){
  #' This function generates a plot showing the coverage of each patient
  #' @param patientTable The coverage table separated by patient
  #' @return A plot of all the data
  
  # Remove uncovered zones (Analyzed separately)
  patientTable <- patientTable %>% filter(depth!=0)
  
  # Build the plot
  p <- ggplot(data = patientTable) + 
    # Define axis from columns and build a rectangle chart from it
    aes(x = start, y = depth) + 
    geom_rect(mapping = aes(xmin = start, ymin = 0, xmax = stop, ymax = depth, fill = patient)) +
    # Set the origin at 0,0
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0))+
    # Details the legend
    labs(title = "Profondeur de séquençage",
         colour ="",
         subtitle = paste0("Moyenne : ",trunc(patientTable$depth*10^2)/10^2),
         x = "Position",
         y = "Profondeur (X)") +
    theme_classic()
  return(p)
}

unknownRegionsMulti <- function(coverage, percentCovered){
  #' This function creates a plot from regions without coverage
  #' @param coverage The table containing the no coverage zones 
  #' @param percentCovered The percentage of coverage in the genome
  #' @return The plot of the no coverage zones

  if(is.null(coverage$patient)){
    values <- aes(xmin = start, xmax = stop,
                  ymin = -0.99, ymax = 0.99)
  }
  else values <- aes(xmin = start, xmax = stop,
                     ymin = -0.99, ymax = 0.99, 
                     fill = patient, alpha = 0.2)
  
  plot <- ggplot() + 
    # Create a rectangle chart
    geom_rect(mapping = aes(xmin = 1, xmax = max(coverage$stop),
                            ymin = -1, ymax = 1),
              fill = "white", color = "black") +
    # Add red parts to the chart where there are unknown regions
    geom_rect(data = coverage,
              mapping = values) +
    # Add a legend
    labs(x = "Position",
         title = "Régions inconnues",
         subtitle = paste0("Non déterminé moyen: ", percentCovered, "%")
         ) +
    # Visual modifications
    theme(axis.text.y =  element_blank(),
          axis.ticks.y = element_blank(),
          panel.background = element_rect(fill = 'white', color = 'white'),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  return(plot)
  
}

vcfMultiAnalysis <- function(variantTable, reference){
  #' This function creates a table for multi patient analysis and a plot for all patient
  #' @param variantTable The table of each variant associated with patients
  #' @param reference The reference sequence of the vcf table
  #' @return A Dataframe containing VCF information and a plot to visualize it

  
  
  # Creation of the table containing all the raw readable data
  info <- unlist(str_split(unique(variantTable$V9), ":"))
  infoTable <- data.frame(matrix(nrow = 0, ncol = length(info) +5))
  colnames(infoTable) <- c("Gene","position", "reference", "mutation", info, "patient")
  infoList <- str_split(variantTable$V10, ":")
  for(i in 1:length(infoList)){
    infoTable[nrow(infoTable)+1,] <- c(variantTable$V1[i], variantTable$V2[i], 
                                       variantTable$V4[i], variantTable$V5[i],
                                       infoList[[i]], variantTable$patient[i])
  }
  
  
  # Create the table containing all mutations with their position
  mutationTable <- variantTable[,-c(1,3,6,7,8,9,10)]
  colnames(mutationTable) <- c("Position", "Reference", "Mutated", "patient")
  referenceLength <- str_length(read.fasta(reference,as.string = TRUE, seqonly = TRUE))
  
  # Build the plot
  plot <- ggplot() + 
    # Create a rectangle chart
    geom_rect(mapping = aes(xmin = 1, xmax = referenceLength,
                            ymin = -1, ymax = 1),
              fill = "white", color = "black") +
    # Add red segment where there is a mutation
    geom_segment(aes(x = mutationTable$Position, xend = mutationTable$Position, y = -1, yend = 1, col = mutationTable$patient), alpha=0.65)+
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

vcfSingleAnalysis <- function(variantTable, reference){
  #' This function creates a table for multi patient analysis and a plot for all patient
  #' @param variantTable The table of each variant associated with patients
  #' @param reference The reference sequence of the vcf table
  #' @return A Dataframe containing VCF information and a plot to visualize it
  
  # Creation of the table containing all the raw readable data
  info <- unlist(str_split(unique(variantTable$V9), ":"))
  infoTable <- data.frame(matrix(nrow = 0, ncol = length(info) +5))
  colnames(infoTable) <- c("Gene","position", "reference", "mutation", info, "patient")
  infoList <- str_split(variantTable$V10, ":")
  for(i in 1:length(infoList)){
    infoTable[nrow(infoTable)+1,] <- c(variantTable$V1[i], variantTable$V2[i], 
                                       variantTable$V4[i], variantTable$V5[i],
                                       infoList[[i]], variantTable$patient[i])
  }
  
  
  # Create the table containing all mutations with their position
  mutationTable <- variantTable[,-c(1,3,6,7,8,9,10)]
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