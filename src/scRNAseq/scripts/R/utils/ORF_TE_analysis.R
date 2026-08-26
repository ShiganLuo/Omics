# Functions
facet_grid_blank = function(){ theme(strip.background = element_blank())}


ORF_analysis_of_TE = function(te,                                  # Name of TE sequence in fasta file
                              return_plot = TRUE,                  # Whether to return plot or dataframe
                              length_ORF = 5,                      # Minimal length of ORF to consider
                              show_stop = TRUE,                    # Whether to show stop codon in final plot
                              show_start = TRUE,                   # Whether to show Methionine in final plot
                              size_peptide_text = 2.5,             # Text size of peptide to show in final plot
                              start_codon = "ATG",                 # regex with start codon to consider
                              fasta_sequences_db = seqs_all,       # fasta file in DNAstringset format
                              peptides_data,                       # dataframe with data on the peptides and from which TEs they came from
                              te_col = "Element_strand",           # name of the column in peptides_data dataframe that contains the name of the TE sequence that code for the peptides
                              seq_col = "Sequence",                # name of the column in peptides_data dataframe that contains the AA sequence of the peptides
                              repel_force_peptide_text = 40        # ggrepel force parameter to avoid peptide text overlapping
                              ){
  
  
  #############################################################################################################################
  ################################################### Retrieving sequences ####################################################
  #############################################################################################################################
  
  # Among sequences in fasta, find fasta sequence with corresponding name
  # Example :
  # >L1PA
  # ATGCATACATGACATG
  # 
  # te should be "L1PA"
  seqs = fasta_sequences_db[te]
  length_ORF_id = paste0("ORF", length_ORF)
  list.to.check = names(seqs)
  
  # Retrieving 3 frames for each strand
  # Forward strand
  seqs_orf1 =  DNAStringSet((as.character(seqs[list.to.check] )))
  seqs_orf2 = DNAStringSet(str_sub(as.character(seqs[list.to.check]), start = 2))
  seqs_orf3 = DNAStringSet(str_sub(as.character(seqs[list.to.check]), start = 3))
  
  # Reverse strand
  seqs_orf4 = reverseComplement(DNAStringSet(as.character(seqs[list.to.check])))
  seqs_orf5 = reverseComplement(DNAStringSet(str_sub(as.character(seqs[list.to.check]), end = -2)))
  seqs_orf6 = reverseComplement(DNAStringSet(str_sub(as.character(seqs[list.to.check]), end = -3)))
  
  # Storing data in list
  list_seqs = list(seqs_orf1, seqs_orf2, seqs_orf3, seqs_orf4, seqs_orf5, seqs_orf6)
  list_sequences = list()
  for (i in 1:6) {
    names(list_seqs[[i]]) <- list.to.check
    list_sequences[[i]] <- list_seqs[[i]][[te]]
  }
  
  # Formating as a DNAStringSet object 
  list_sequences_DNA = Biostrings::DNAStringSet(list_sequences )
  names(list_sequences_DNA) <- c(paste0("ORF", 1:3), paste0("ORF-", 1:3))
  rm(list_seqs, list_sequences)
  
  # Translating DNA sequences to AA sequences
  # Forward strand
  seqs_orf1 =  translate(seqs_orf1, if.fuzzy.codon = "solve")
  seqs_orf2 = translate(seqs_orf2, if.fuzzy.codon = "solve")
  seqs_orf3 = translate(seqs_orf3, if.fuzzy.codon = "solve")
  # Reverse strand
  seqs_orf4 = translate(seqs_orf4, if.fuzzy.codon = "solve")
  seqs_orf5 = translate(seqs_orf5, if.fuzzy.codon = "solve")
  seqs_orf6 = translate(seqs_orf6, if.fuzzy.codon = "solve")
  
  
  # Storing AA data in list
  list_seqs = list(seqs_orf1, seqs_orf2, seqs_orf3, seqs_orf4, seqs_orf5, seqs_orf6)
  list_sequences = list()
  for (i in 1:6) {
    names(list_seqs[[i]]) <- list.to.check
    list_sequences[[i]] <- list_seqs[[i]][[te]]
  }
  
  # Formating as a AAStringSet object 
  list_sequences_AA = Biostrings::AAStringSet(list_sequences )
  names(list_sequences_AA) <- c(paste0("ORF", 1:3), paste0("ORF-", 1:3))
  rm(list_seqs, list_sequences)
  
  ##############################################################################################################################
  ###################################################  Matching stop codons ####################################################
  ##############################################################################################################################
  
  # Defining AA motif that will match stop codon
  stop_set = AAStringSet(x = c("*"))
  
  # List to store results
  list_stop_results = list()
  
  # For each AA sequences among all 6 frames, retrieving coordinates of stop codon
  for (i in 1:6) {
    stop_results = Biostrings::matchPDict(pdict = stop_set, subject =  list_sequences_AA[[i]], max.mismatch  = 0)
    list_stop_results[[i]] <- unlist(start(stop_results))
  }
  names(list_stop_results) <- names(list_sequences_AA)
  
  ##############################################################################################################################
  ###################################################  Matching Methionines ####################################################
  ##############################################################################################################################
  
  # Defining AA motif that will match stop codon
  methionine_set = AAStringSet(x = c("M"))
  
  # List to store results
  list_methionine_results = list()
  
  # For each AA sequences among all 6 frames, retrieving coordinates of methionines
  for (i in 1:6) {
    methionine_results = Biostrings::matchPDict(pdict = methionine_set, subject =  list_sequences_AA[[i]], max.mismatch  = 0)
    list_methionine_results[[i]] <- unlist(start(methionine_results))
  }
  names(list_methionine_results) <- names(list_sequences_AA)
  
  ##############################################################################################################################
  ###################################################  Matching peptides #######################################################
  ##############################################################################################################################
  
  # Retrieving list of peptides to fetch in TEs sequences
  peptide_id = peptides_data %>%
    dplyr::filter(get(te_col) == te) %>%
    dplyr::pull(get(seq_col)) %>%
    unique()
  
  # Creating AAStringset object using peptide sequences
  peptide_aa_set = AAStringSet(x = peptide_id)
  names(peptide_aa_set) = peptide_id
  
  # Lists to store results : start coordinates in DNA sequence, end coordinates in DNA sequence, AA sequence
  list_peptide_DNA_start = list()
  list_peptide_DNA_end = list()
  list_peptide_aa_seq = list()
  
  # For each frame
  for (i in 1:6) {
    # Matching AA sequences with translated AA TE sequence
    peptide_aa_results_df = Biostrings::matchPDict(pdict = peptide_aa_set, subject =  list_sequences_AA[[i]], max.mismatch  = 0)
    peptide_aa_results <- unlist(start(peptide_aa_results_df))
    
    # If peptide found
    if (length(unique(peptide_aa_results)) > 0) {
      # Retrieving and inferring DNA start coordinates
      peptide_nt_results_start = unlist(start(peptide_aa_results_df)) * 3 -2
      # Retrieving and inferring DNA end coordinates
      peptide_nt_results_end = unlist(end(peptide_aa_results_df)) * 3
      # Storing start coordinates, end coordinates and peptide sequence
      list_peptide_DNA_start[[i]] <- peptide_nt_results_start
      list_peptide_DNA_end[[i]] <- peptide_nt_results_end 
      list_peptide_aa_seq[[i]] <- names(peptide_aa_results)
      
    # If no peptide found
    }else{
      list_peptide_DNA_start[[i]] <- ""
      list_peptide_DNA_end[[i]] <- ""
      list_peptide_aa_seq[[i]] <- ""
    }
    
  }
  
  names(list_peptide_DNA_start) <- names(list_sequences_AA)
  names(list_peptide_DNA_end) <- names(list_sequences_AA)
  names(list_peptide_aa_seq) <- names(list_sequences_AA)
  
  ##############################################################################################################################
  ###################################################  Matching ORF ############################################################
  ##############################################################################################################################
  
  # List to store results
  list_orf_results = list()
  
  # For each frame, will find ORF within sequence.
  # Depending on start codon motif submitted.
  # Example :
  # canonical start codon : "ATG"
  # canonical + non canonical start codons : "ATG|CTG|GTG|TTG"
  for (id_orf_loop in 1:6) {
    orf_df <- findORFs(list_sequences_DNA[id_orf_loop], startCodon = start_codon, minimumLength = 5, longestORF = TRUE) 
    
    if (length(orf_df) >= 1 ) {
      names(orf_df) <- names(list_sequences_DNA)[id_orf_loop]
      list_orf_results[[names(list_sequences_DNA)[id_orf_loop]]] <- orf_df[[1]]
    }
  }
  
  # Storing results in IRangesList object
  orf_df_test = IRangesList(list_orf_results)
  
  # Format results and filter ORFs smaller than submitted minimum length
  orf_df  = as.data.frame(orf_df_test) %>%
    dplyr::rename(ORF = group_name) %>%
    dplyr::filter(width >= length_ORF * 3) %>%
    dplyr::select(-group, -width) %>%
    dplyr::mutate(type = length_ORF_id, id = row_number()) %>%
    dplyr::ungroup()
  
  
  ##############################################################################################################################
  ###################################################  Formating results #######################################################
  ##############################################################################################################################
  
  
  # Retrieving peptide start
  df3 = reshape2::melt(list_peptide_DNA_start) %>%
    dplyr::rename(ORF = L1, start = value) %>%
    dplyr::mutate(id = row_number(), type = "Peptide")
  
  # Retrieving peptide end
  df4 = reshape2::melt(list_peptide_DNA_end) %>%
    dplyr::rename(ORF = L1, end = value) %>%
    dplyr::mutate(id = row_number(), type = "Peptide")
  
  # Retrieving peptide sequence
  df3bis = reshape2::melt(list_peptide_aa_seq) %>%
    dplyr::rename(ORF = L1, Sequence = value) %>%
    dplyr::mutate(id = row_number(), type = "Peptide")
  
  # Merging results
  df34 = left_join(x = df3, y = df3bis, by = c("id", "type", "ORF")) %>% 
    left_join(y = df4, by = c("id", "type", "ORF")) %>%
    dplyr::filter(start != "") %>%
    dplyr::group_by(start, end, type, ORF) %>%
    dplyr::top_n(1, wt = id) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(start = as.double(start), end = as.double(end))
  
  # Retrieving stop coordinates
  df5 = reshape2::melt(list_stop_results) %>%
    dplyr::rename(start = value, ORF = L1) %>%
    dplyr::mutate(id = row_number(), type = "Stop")  %>%
    dplyr::mutate(end = start * 3, start = start * 3 - 2)
  
  # Retrieving full sequences coordinates
  df6 = data.frame(start = rep(1, 6), 
                   end = rep(length(list_sequences_DNA$ORF1), 6), 
                   type = rep("Full_seq", 6), id = 1:6, 
                   ORF = c( paste0("ORF", 1:3), paste0("ORF.", 1:3)))
  
  # Retrieving methionine coordinates
  df7 = reshape2::melt(list_methionine_results) %>%
    dplyr::rename(start = value, ORF = L1) %>%
    dplyr::mutate(id = row_number(), type = "Methionine")  %>%
    dplyr::mutate(end = start * 3, start = start * 3 - 2)
  
  
  # If methionine present
  if (nrow(df7) > 0) {
    # Function to create artificial ORF starting with a methionine and finishing at the end of the sequence if no codon stop found
    df8 = bind_rows(list(df5 %>% dplyr::group_by(ORF) %>% dplyr::slice_max(end),
                         df6 %>% dplyr::group_by(ORF) %>% dplyr::slice_max(end) %>% dplyr::mutate(ORF = gsub("\\.", "-", ORF)),
                         df7 %>% dplyr::group_by(ORF) %>% dplyr::slice_max(end))) %>%
      dplyr::select(end, type, ORF) %>%
      tidyr::pivot_wider(names_from = type, values_from = end) %>%
      dplyr::mutate(id = 1:(length(end)), type = length_ORF_id) %>%
      dplyr::filter(Stop < Methionine | is.na(Stop)) %>%
      dplyr::select(-Stop) %>%
      dplyr::rename("end" = Full_seq, "start" = Methionine)
    
    # Merging all results
    df.all = bind_rows(list(df34, df5, df6, orf_df %>% dplyr::filter((start-1) %% 3 == 0), df7, df8))
    
  }else{
    # Merging all results
    df.all = bind_rows(list(df34, df5, df6, orf_df %>% dplyr::filter((start-1) %% 3 == 0), df7))
  }
  
  df.all = df.all %>%
    dplyr::mutate(ORF = factor(gsub("\\.", "-", ORF), levels = c( paste0("ORF", 1:3), paste0("ORF-", 1:3))))  %>%
    dplyr::ungroup() 
  
  ##############################################################################################################################
  ###################################################  Plot ####################################################################
  ##############################################################################################################################
  
  # Color vector
  list_col = c(Stop = "black", Methionine = "forestgreen",  ORF50  = "goldenrod", Peptide = "red")
  names(list_col) = c("Stop", "Methionine", length_ORF_id, "Peptide")
  
  if (return_plot == TRUE) {
    df.seq = df.all %>% dplyr::filter(type == "Full_seq")
    plot = ggplot(df.all) +
      theme_void() +
      facet_grid(ORF~.) +
      # Full sequence illustration
      geom_rect(aes(xmin=start, xmax=end, ymin=0.8, ymax=1.2),
                fill = ifelse(test = df.seq$ORF %in% c("ORF1", "ORF2", "ORF3") , yes = "navyblue", no = "darkslategrey"), 
                color= ifelse(test = df.seq$ORF %in% c("ORF1", "ORF2", "ORF3") , yes = "navyblue", no = "darkslategrey"), 
                alpha=0.7, 
                data = df.seq) +
      # ORF illustration
      geom_rect(aes(xmin=start, xmax=end, ymin=0.5, ymax=1.5, fill=type), 
                color="black", 
                alpha=0.7, 
                data = df.all %>%  dplyr::filter(type == length_ORF_id )) +
      # Stop illustration
      {if(show_stop == TRUE)geom_rect(aes(xmin=start, xmax=end, ymin=0.5, ymax=1.5, fill=type), 
                                      color="black", 
                                      alpha=0.5, 
                                      data = df.all %>% dplyr::filter(type == "Stop"))} +
      # Peptide illustration
      geom_rect(aes(xmin=start, xmax=end, ymin=0, ymax=2, fill=type), 
                color="red", 
                alpha=0.7, 
                data = df.all %>%  dplyr::filter(type == "Peptide")) +
      # Peptide text
      geom_text_repel(data = df.all %>%  dplyr::filter(type == "Peptide"), 
                      aes(label = Sequence, x = start, y = 2.5), 
                      size = size_peptide_text, 
                      direction = "x", 
                      color = "red", 
                      force = repel_force_peptide_text, 
                      segment.size = 0.1)+
      # Methionine illustration
      {if(show_start == TRUE)geom_rect(aes(xmin=start, xmax=end, ymin=0.5, ymax=1.5, fill=type), 
                                       color="forestgreen", 
                                       alpha=0.6, 
                                       data = df.all %>%  dplyr::filter(type == "Methionine"))} +
      
      facet_grid_blank() +
      scale_fill_manual(values = list_col) +
      guides(fill = guide_legend(override.aes=list(color=NA))) +
      labs(title = te, fill = "Type") +
      theme(axis.line.x = element_line(colour = "black", size = 0.5, linetype = 1,  lineend = "butt"), 
            axis.text.x = element_text(vjust = 1, margin = margin(0.1,  unit = "cm")),
            axis.ticks.x  = element_line(colour = "grey50")) +
      coord_cartesian(ylim = c(0, 3))
    return(plot)
  }else{
    return(df.all)
  }
  
  
}


