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

sequenceQualityPlots <- function(fastqData){
  #' This generates sequence quality plots for multiple patients
  #' @param fastqData
  #' @return Plots about sequence quality
  
  # This is a lollipop plot with the mean and median quality score sequence 
  # Of each patient
  p1 <- ggplot(fastqData[[1]], aes(x=fastqData[[1]]$patient, y=weighted.mean(fastqData[[1]]$GC_content, fastqData[[1]]$occurence)))+
    geom_point()+
    geom_segment(aes(x=fastqData[[1]]$patient, y=weighted.mean(fastqData[[1]]$GC_content, fastqData[[1]]$occurence), xend = 0, yend = weighted.mean(fastqData[[1]]$GC_content, fastqData[[1]]$occurence)))
  
}

multiReadMapPlot <- function(data){
  i <- 0
  plot <- ggplot()
  for(elem in data){
    
    plot <- plot + geom_bar(data = elem[[1]][[1]], aes(x = GC_content, colour = i))
    i <- i+1
  }
  return(plot)
}

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