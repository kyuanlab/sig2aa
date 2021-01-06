#' Create trinucleotide context for a single base substitution
#' 
#' Create a single row data frame containing trinucleotide context based on provided reference 
#' and alternative allele and the position of the mutant allele in the trinucleotide context.
#' @param ref The reference allele, one from the four bases A, T, C, G. 
#' @param alt The alternative allele, one from the four bases A, T, C, G, and different from ref
#' @param pos The of position of the mutant allele in trinucleotide context, only a 
#' @return A 3-column data frame containing 1 row of trinucleotide context based on the provided ref, alt, and pos 
#' @export
#' @example 
#' df <- CreateTrinucleotideContextForSingleBaseSubstitution(ref = "C", alt = T, pos = 1)
#' df <- CreateTrinucleotideContextForSingleBaseSubstitution(ref = "C", alt = T, pos = 2)
#' df <- CreateTrinucleotideContextForSingleBaseSubstitution(ref = "C", alt = T, pos = 3)
#' 
CreateTrinucleotideContextForSingleBaseSubstitution <- function(ref, alt, pos = c(1,2,3)) {
  stopifnot(ref %in% c("A", "T", "C", "G"))
  stopifnot(alt %in% c("A", "T", "C", "G") | !identical(ref, alt) )
  stopifnot(pos %in% c(1, 2, 3) | !length(pos)==1)

  base <- c("A", "T", "C", "G")
  trinucleotideContextRef <- c()
  trinucleotideContextAlt <- c()
  
  if (pos == 2 ) {
    for (i in 1:length(base)) {
      for (j in 1:length(base)) {
        trinucleotideContextRef <- c(trinucleotideContextRef, paste0( base[i], ref, base[j]))
        trinucleotideContextAlt <- c(trinucleotideContextAlt, paste0( base[i], alt, base[j]))
      }
    }
  } else if ( pos == 1 ) {
    for (i in 1:length(base)) {
      for (j in 1:length(base)) {
        trinucleotideContextRef <- c(trinucleotideContextRef, paste0( ref, base[i], base[j]))
        trinucleotideContextAlt <- c(trinucleotideContextAlt, paste0( alt, base[i], base[j]))
      }
    }
  } else {
    for (i in 1:length(base)) {
      for (j in 1:length(base)) {
        trinucleotideContextRef <- c(trinucleotideContextRef, paste0( base[i], base[j], ref))
        trinucleotideContextAlt <- c(trinucleotideContextAlt, paste0( base[i], base[j], alt))
      }
    }
  }
  return(data.frame(ref = trinucleotideContextRef, 
                    alt = trinucleotideContextAlt, 
                    pos = pos))
}

#' Create trinucleotide context and corresponding amino acid annotations 
#' 
#' Create a data frame containing all possible combinations of trinucleotide context and corresponding 
#' amino acid annotations based on provided reference and alternative 
#' allele and the position of the mutant allele in the trinucleotide context.
#' 
#' @param ref The reference allele, one from the four bases A, T, C, G. 
#' @param alt The alternative allele, one from the four bases A, T, C, G, and different from ref
#' @param pos The of position of the mutant allele in trinucleotide context, multiple digits allowed to 
#' create data frames containing different reading frames where the mutant alleles' positions are set 
#' accroding to pos 
#' @return A 5-column data frame containing different trinucleotide context based on the provided ref, alt, 
#' and pos 
#' @export
#' @example 
#' df <- CreateTrinucleotideContextAndAminoAcidAnnotation(ref = "C", alt = T, pos = 1)
#' df <- CreateTrinucleotideContextAndAminoAcidAnnotation(ref = "C", alt = T, pos = 2)
#' df <- CreateTrinucleotideContextAndAminoAcidAnnotation(ref = "C", alt = T, pos = 1:2)
#' df <- CreateTrinucleotideContextAndAminoAcidAnnotation(ref = "C", alt = T, pos = 1:3)
#' 
CreateTrinucleotideContextAndAminoAcidAnnotation <- function (ref, alt, pos = c(1,2,3,1:2,2:3,1:3)) {
  
  stopifnot(ref %in% c("A", "T", "C", "G"))
  stopifnot(alt %in% c("A", "T", "C", "G") | !identical(ref, alt) )
  stopifnot(all(pos %in% c(1, 2, 3)))
 
  trinucleotideContextDF <- 
    do.call (rbind, lapply(pos, 
                           CreateTrinucleotideContextForSingleBaseSubstitution, ref = ref, alt = alt))

  trinucleotideContextDF$ref_aa <- unlist(Map( function(x) Biostrings::GENETIC_CODE[x], trinucleotideContextDF$ref))
  trinucleotideContextDF$alt_aa <- unlist(Map( function(x) Biostrings::GENETIC_CODE[x], trinucleotideContextDF$alt))
  
  return(trinucleotideContextDF)
}

#' Create full trinucleotide context and their amino acid annotations 
#' 
#' Create a data frame containing all possible trinucleotide context according to single or double strand specification.
#' @param strandInfo Strand information, "single" or "double" 
#' @return A data frame containing all possible trinucleotide context according to single or double strand specification.
#' @example 
#' CreateFullTrinucleotideContextAndAminoAcidAnnotation(strandInfo = "single")
#' CreateFullTrinucleotideContextAndAminoAcidAnnotation(strandInfo = "double")
#' 
CreateFullTrinucleotideContextAndAminoAcidAnnotation <- function(strandInfo = "single") {
  
  if (strandInfo == "single") {
    mutationPairs <- data.frame(ref = c("A","A","A","T","T","T","C","C","C","G","G","G"), 
                                alt = c("T","C","G","A","C","G","A","T","G","A","T","C"))
  } else {
    mutationPairs <- data.frame(ref = c("T","T","T","C","C","C"), 
                                alt = c("A","C","G","A","T","G"))
  }
  
  trinucleotideContextDF <- c()
  for (ii in 1:nrow(mutationPairs)) {
    trinucleotideContextDF <- rbind(trinucleotideContextDF, 
                                    CreateTrinucleotideContextAndAminoAcidAnnotation(ref = mutationPairs[ii,]$ref, 
                                                                                     alt = mutationPairs[ii,]$alt,
                                                                                     pos = 2))
  }
  trinucleotideContextDF %>% dplyr::rowwise() %>% 
    dplyr::mutate(type = paste0(substr(ref,start = 2, stop = 2),">",substr(alt,start = 2, stop = 2)))
}

#' Parse Cosmic SBS mutational signatures
#' 
#' A function to parse mutational signatures formatted according to the COSMIC single base substations (SBS) signatures 
#' @param filePath Path to a csv file formatted according the COSMIC single base substations (SBS) signatures.
#' @param type The column to be used to compute frequency related measures.
#' @return A data frame containing single base substation types, trinucleotide context, and frequencies. The data frame will be used to check amino acid hits. 
#' @export
ParseComsicSignatures <- function(filePath, type = "GRCh38") {
  
  # Load COSMIC signatures
  cosmicSig <- read.csv(filePath, stringsAsFactors=F)
  selectedType <- colnames(cosmicSig)[grep(type, colnames(cosmicSig))]
  cosmicSig <- cosmicSig[, c("Type", "Subtype", selectedType)]
  
  # Generate a template
  trinucleotideContextDF <- CreateFullTrinucleotideContextAndAminoAcidAnnotation(strandInfo = "double")
  trinucleotideContextDF <- dplyr::left_join(trinucleotideContextDF, cosmicSig, by = c("type"="Type", "ref"="Subtype"))
  colnames(trinucleotideContextDF)[7] <- "freq" 
  trinucleotideContextDF
}

#' Add a nucleotide to both 5' and 3' directions 
#' 
#' Add a nucleotide to both 5' and 3' directions of the provided sequence and enumerate all possible combinations. 
#' @param context The sequence context, any combination of any length of the four bases A, T, C, G
#' @return A vector of all 16 possible combinations of padded sequences
#' @export 
#' @example 
#' AddNucleotide5pAnd3p("ACG")
#'  [1] "AACGA" "AACGT" "AACGC" "AACGG" "TACGA" "TACGT" "TACGC" "TACGG" "CACGA" "CACGT" "CACGC" "CACGG" "GACGA" "GACGT" "GACGC"
#' [16] "GACGG"
#' AddNucleotide5pAnd3p("AC")
#' [1] "AACA" "AACT" "AACC" "AACG" "TACA" "TACT" "TACC" "TACG" "CACA" "CACT" "CACC" "CACG" "GACA" "GACT" "GACC" "GACG"
AddNucleotide5pAnd3p <- function( context ) {
  
  base <- c("A", "T", "C", "G")
  contextAll <- c()
  for (i in base) {
    for(j in base) {
      contextAll <- c(contextAll, paste0( i, context, j))
    } 
  }
  return( contextAll )
} 

#' Get trinucleotide context reading frames 
#'
#' Create a data frame containing trinucleotide context in 3 different reading frames, ie the mutant allele being the 1st, 2nd, or 3rd of a trinucleotide context  
#' @param x The trinucleotide context without considering reading frame, where the mutant allele is assumed to be in the middle (2nd allele) of three alleles.
#' @return A data frame (9 rows, 3 columns) containing all possible trinucleotide context acrodding to the 3 reading frames. 
#' @export 
#' @example 
#' GetTrinucleotideContextReadingFrames("ACG")
#'       ref pos
#' 1     CGA   1
#' 2     ACG   2
#' 3     AAC   3
#' 4     CGT   1
#' 5     CGC   1
#' 6     CGG   1
#' 7     TAC   3
#' 8     CAC   3
#' 9     GAC   3
GetTrinucleotideContextReadingFrames <- function(x) {
  trinucleotideContextReadingFrames <- x %>% 
    AddNucleotide5pAnd3p %>% # pad a base 5' and 3' of the context
    Biostrings::DNAStringSet()  %>% 
    # Listing all trinucletide context resulting from 3 different reading frames. 
    # Note, the following numbers (start = 3:1, end = 5:3) only works for x being a reading frame 
    # where the mutatnt allele is the 2nd allele.
    lapply(Views, start = 3:1, end = 5:3) %>%   
    lapply(function(x) x %>% lapply(Biostrings::toString) %>% unlist) %>%
    unlist %>% data.frame(ref = ., pos = rep(1:3,16) ) %>% unique # only use unique context
  row.names(trinucleotideContextReadingFrames) <- NULL
  trinucleotideContextReadingFrames
}

#' mutate helper function  
#'
#' A helper function to mutate sequence in dplyr::mutate
#' @param ref
#' @param pos
#' @param alt
#' @return mutant sequence
#' @example 
#' mydata <- data.frame(ref = c("ACA", "ACG"), pos = c(2,2))
#' dplyr::mutate(mydata, alt = .mutate(ref, pos, "T"))
.mutate <- function(ref, pos=2, altAllele) {
  substr(ref, start = pos, stop = pos) <- altAllele
  ref
}


#' Compute context likelihood
#' 
#' A helper function to compute likelihood of each context reading frame in dplyr::mutate
#' @param pos The position of the mutant allele in the trinucleotide context
#' @param ref The reference trinucleotide context
#' @param trinucleotideContextFrequency The frequency of the reference trinucleotide context
#' @param marginalFrequencies The marginal frequencies of the four bases 
#' @param frameProb The probability of the three possible reading frames
#' @return 
#' @example 
#' mydata <- data.frame(ref = "ACG", pos = 2)
#' dplyr::mutate(mydata, likelihood = .likelihood(pos, 
#' ref, 1/16, c("A" = 0.25, "T" = 0.25, "C" = 0.25, "G" = 0.25), c(1/3, 1/3, 1/3)))
.likelihood <- function(pos = c(1,2,3),
                        ref,
                        trinucleotideContextFrequency, 
                        marginalFrequencies, 
                        frameProb = c(1/3, 1/3, 1/3)) {
  
  if (pos == 2) {
    likelihood <- trinucleotideContextFrequency * frameProb[pos]
  } else if (pos == 1) { 
    likelihood <-  trinucleotideContextFrequency * frameProb[pos] *
      marginalFrequencies[substr(ref, start = 3, stop = 3)]
  } else {
    likelihood <-  trinucleotideContextFrequency * frameProb[pos] *
      marginalFrequencies[substr(ref, start = 1, stop = 1)]
  }
  
  likelihood
}

#' Generate all possible amino acid hits and compute their likelihood 
#'
#' Generate a data frame containing all possible amino acid hits and compute their likelihood based on 
#' a single trinucleotide context, its frequency, marginal frequencies of the four bases, probabilities of 3 different reading frames. 
#' @param context A trinucleotide context where the mutation will happen at the 2nd position. 
#' @param altAllele The alternative allele 
#' @param trinucleotideContextFrequency The frequency of the provided trinucleotide context
#' @param marginalFrequencies Marginal frequencies of the four bases 
#' @param frameProb Probabilities of 3 different reading frames. Probabilities are equal by default. 
#' @return A data frame (9 rows, 5 columns) containing all possible amino acid hits and compute their likelihood 
#' @export 
#' @example 
ComputeAminoAcidLikelihood <- function(context, 
                                       altAllele, 
                                       trinucleotideContextFrequency,
                                       marginalFrequencies,
                                       frameProb = c(1/3, 1/3, 1/3)) {
  
  #Get the mutant being 1st, 2nd, and 3rd of a codon
  trinucleotideContextReadingFrames <- 
    GetTrinucleotideContextReadingFrames(context) 
  
  trinucleotideContextReadingFrames <- 
    dplyr::mutate(dplyr::rowwise(trinucleotideContextReadingFrames), 
           alt = .mutate(ref, pos, altAllele),
           ref_aa = Biostrings::GENETIC_CODE[ref],
           alt_aa = Biostrings::GENETIC_CODE[alt])
  
  trinucleotideContextReadingFrames <- 
    dplyr::mutate(dplyr::rowwise(trinucleotideContextReadingFrames), 
           likelihood = .likelihood(pos, 
                                    ref, 
                                    trinucleotideContextFrequency, 
                                    marginalFrequencies, 
                                    frameProb))
  
  attr(trinucleotideContextReadingFrames$ref_aa, "names") <- NULL
  attr(trinucleotideContextReadingFrames$alt_aa, "names") <- NULL
  attr(trinucleotideContextReadingFrames$likelihood, "names") <- NULL
  attr(trinucleotideContextReadingFrames, "groups") <- NULL
  
  return(as.data.frame(trinucleotideContextReadingFrames))
}

#' Convert trinucleotide context to aimno acid context
#'
#' A function to convert trinucleotide context to all possible aimno acid context and compute their likelihood.
#' @param contextDataFrame A data frame containing trinucleotide context and their frequencies 
#' @param marginalFrequencies marginalFrequencies Marginal frequencies of the four bases 
#' @param frameProb Probabilities of 3 different reading frames. Probabilities are equal by default.
#' @return A data frame containing all possible aimno acid replacement and their likelihood
#' @export
#' @example 
ConvertTrinucleotideContextToAimnoAcidContext <- function(contextDataFrame,
                                                          marginalFrequencies,
                                                          frameProb = c(1/3, 1/3, 1/3)) {
    
    cat(paste0("Converting ", nrow(contextDataFrame), " trinucleotide context to aimno acid hits ... \n"))
    trinucleotideContextReadingFrames <- c()
    
    for (i in 1:nrow(contextDataFrame)){
      cat(paste0("\r", i, "/", nrow(contextDataFrame) )  )
      dd <- contextDataFrame[i,]
      trinucleotideContextReadingFrames <- 
        rbind(trinucleotideContextReadingFrames, 
              ComputeAminoAcidLikelihood(context = dd$ref,
                                         altAllele = substr(dd$alt, 2,2),
                                         dd$freq, 
                                         marginalFrequencies,
                                         frameProb = frameProb) 
        )
    } 
    
    cat("\n")
    # Filter out those context that don't generate amino acid change 
    aminoAcidContext <-
      dplyr::filter(trinucleotideContextReadingFrames, ref_aa != alt_aa)
    
    aminoAcidContextUniq <- aminoAcidContext %>% 
      dplyr::group_by(ref_aa, alt_aa) %>% 
      dplyr::summarise(likelihood = sum(likelihood), .groups = 'drop')
    
    aminoAcidContextUniq
}

#' Create amino acid replacement map 
#'
#' A function to generate all possible amino acid replacement with a given editing distance
#' @param numOfChanges Editing distance of the amino acid replacement
#' @return a data frame 
#' @export
#' @example 
#' CreateAminoAcidReplacementMap(numOfChanges = 1)
#' 
CreateAminoAcidReplacementMap <- function(numOfChanges = 1) {
  codonReplacementDistMatrix = as.matrix(Biostrings::stringDist(names(Biostrings::GENETIC_CODE)))
  codonReplacementDistindex = which(codonReplacementDistMatrix==numOfChanges, arr.ind = T)
  codonReplacementDist = data.frame(ref = names(Biostrings::GENETIC_CODE)[codonReplacementDistindex[,1]],
                                    alt = names(Biostrings::GENETIC_CODE)[codonReplacementDistindex[,2]])
  codonReplacementDist = dplyr::mutate(dplyr::rowwise(codonReplacementDist),
                                        ref_aa = Biostrings::GENETIC_CODE[ref],
                                        alt_aa = Biostrings::GENETIC_CODE[alt])
  
  aminoAcidReplacementDist = codonReplacementDist %>% dplyr::filter(ref_aa != alt_aa) 
  
  aminoAcidReplacementDist %>% 
    dplyr::group_by(ref_aa, alt_aa) %>% 
    dplyr::summarise(.groups = 'drop') %>% dplyr::mutate(label = paste0(ref_aa, ">", alt_aa))
}

#' Map amino acid hits 
#' 
#' Map amino acid hits to all possible amino acid replacement with an editing distance of 1.
#' @param aminoAcidHits
#' @return A data frame containing all possible 170 amino acid replacement with an editing distance of 1.
#' @export
MapAminoAcidHits <- function(aminoAcidHits) {
  
  aminoAcidReplacementMap <- CreateAminoAcidReplacementMap(numOfChange = 1)
  aminoAcidReplacementMap <- dplyr::left_join(aminoAcidReplacementMap, 
                                              aminoAcidHits, 
                                              by = c("ref_aa", "alt_aa"))
  
  aminoAcidReplacementMapNorm <- aminoAcidReplacementMap %>% 
    dplyr::rowwise() %>%
    dplyr::mutate(likelihood_norm = if (is.na(likelihood)) {0} else {likelihood}) 
  
  aminoAcidReplacementMapNorm$likelihood_norm <- 
    aminoAcidReplacementMapNorm$likelihood_norm/sum(aminoAcidReplacementMapNorm$likelihood_norm)
  aminoAcidReplacementMapNorm
}

#' Plotting amino acid signatures
#' 
#' Plotting the converted mutational signatures in amino acid replacements. 
#' @param aminoData A data frame containing amino acid replacements and their probabilities
#' @param type A reference or alternative map. 
PlotAminoAcidSignatures <- function(aminoData, type) {
  
  if (type == "ref"){
    aminoData$x_aa <- aminoData$ref_aa
  } else {
    aminoData$x_aa <- aminoData$alt_aa
  }
  
  p1 <-  ggplot(aminoData, aes(x=label, y=likelihood_norm, fill = x_aa) ) +
    geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    guides(fill = guide_legend(ncol = 2))+
    scale_fill_manual(values = glasbey(21))+
    facet_grid(.~x_aa, scales = "free_x", space = "free")+
    ylab("Percentage")+xlab("")+
    theme( axis.ticks = element_blank(),
           panel.background = element_rect(fill = "gray98",colour = "white",
                                           size = 0.5, linetype = "solid"),
           panel.grid.minor.y= element_line(colour = "gray70", linetype = "dotted"),
           panel.grid.major = element_blank(),
           panel.spacing.x = unit(1, "mm"),
           axis.title.y = element_text(size=20),
           axis.title.x = element_text(),
           axis.text.y = element_text(size=10, face = "bold"),
           axis.text.x = element_text(size=10, face = "bold"),
           strip.background.x = element_blank(),
           strip.text.x = element_blank()) 
  
  p2 <- ggplot(aminoData)+
    geom_bar(mapping = aes(x = label, y = 1, fill = x_aa), 
             stat = "identity", 
             width = 1)+
    scale_fill_manual(values = glasbey(21))+
    theme_void()+
    theme(panel.spacing.x = unit(1, "mm"))+
    facet_grid(.~x_aa, scales = "free_x", space = "free")+
    theme(strip.text.x = element_text(size = 18))
  
  p1 <- p1 + theme(legend.position = "none")
  p2 <- p2 + theme(legend.position = "none")
  
  plot <- cowplot::plot_grid(p2, p1, align = "v", ncol = 1, axis = "tb", rel_heights = c(1.5, 15))
  p <- cowplot::plot_grid(plot, nrow = 1, rel_widths = c(10, 1.5))
  p
}

#' Plotting amino acid signatures
#' 
#' Plotting the converted mutational signatures in amino acid replacements. 
#' @param aminoData A data frame containing amino acid replacements and their probabilities
#' @param type A reference or alternative map. 
PlotAminoAcidSignatures1 <- function(aminoData, type) {
  
  if (type == "ref"){
    aminoData$x_aa <- aminoData$ref_aa
    
    p1 <-  ggplot(aminoData, aes(x=label, y=likelihood_norm, fill = x_aa) ) +
      geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      guides(fill = guide_legend(ncol = 2))+
      scale_fill_manual(values = glasbey(21))+
      facet_grid(.~x_aa, scales = "free_x", space = "free")+
      ylab("Percentage")+xlab("")+
      theme( axis.ticks = element_blank(),
             panel.background = element_rect(fill = "gray98",colour = "white",
                                             size = 0.5, linetype = "solid"),
             panel.grid.minor.y= element_line(colour = "gray70", linetype = "dotted"),
             panel.grid.major = element_blank(),
             panel.spacing.x = unit(1, "mm"),
             axis.title.y = element_text(size=20),
             axis.title.x = element_text(),
             axis.text.y = element_text(size=10, face = "bold"),
             axis.text.x = element_text(size=9, face = "bold"),
             strip.background.x = element_blank(),
             strip.text.x = element_blank()) 
    
    p2 <- ggplot(aminoData)+
      geom_bar(mapping = aes(x = label, y = 1, fill = x_aa), 
               stat = "identity", 
               width = 1)+
      scale_fill_manual(values = glasbey(21))+
      theme_void()+
      theme(panel.spacing.x = unit(1, "mm"))+
      facet_grid(.~x_aa, scales = "free_x", space = "free")+
      theme(strip.text.x = element_text(size = 18))
    
    p1 <- p1 + theme(legend.position = "none")
    p2 <- p2 + theme(legend.position = "none")
    
    plot <- cowplot::plot_grid(p2, p1, align = "v", ncol = 1, axis = "tb", rel_heights = c(1.5, 15))
    
  } else {
    aminoData$x_aa <- aminoData$alt_aa
    
    p1 <-  ggplot(aminoData, aes(x=label, y=likelihood_norm, fill = x_aa) ) +
      geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      guides(fill = guide_legend(ncol = 2))+
      scale_fill_manual(values = glasbey(21))+
      facet_grid(.~x_aa, scales = "free_x", space = "free")+
      ylab("Percentage")+xlab("")+
      theme( axis.ticks = element_blank(),
             panel.background = element_rect(fill = "gray98",colour = "white",
                                             size = 0.5, linetype = "solid"),
             panel.grid.minor.y= element_line(colour = "gray70", linetype = "dotted"),
             panel.grid.major = element_blank(),
             panel.spacing.x = unit(1, "mm"),
             axis.title.y = element_text(size=20),
             axis.title.x = element_text(),
             axis.text.y = element_text(size=10, face = "bold"),
             axis.text.x = element_text(size=9, face = "bold"),
             strip.background.x = element_blank(),
             strip.text.x = element_blank()) 
    
    p2 <- ggplot(aminoData)+
      geom_bar(mapping = aes(x = label, y = 1, fill = x_aa), 
               stat = "identity", 
               width = 1)+
      scale_fill_manual(values = glasbey(21))+
      theme_void()+
      theme(panel.spacing.x = unit(1, "mm"))+
      facet_grid(.~x_aa, scales = "free_x", space = "free")+
      theme(strip.text.x = element_text(size = 18))
    
    p1 <- p1 + theme(legend.position = "none")
    p2 <- p2 + theme(legend.position = "none")
    
    plot <- cowplot::plot_grid(p2, p1, align = "v", ncol = 1, axis = "tb", rel_heights = c(1.5, 15))
    
  }
  
  p <- cowplot::plot_grid(plot, nrow = 1, rel_widths = c(10, 1.5))
  p
}

#' Plotting reference and alternative amino acid signatures 
#' 
#' Plotting the converted mutational signatures in amino acid replacements, organised in in both reference and alternative amino acid  
#' @param aminoData A data frame containing amino acid replacements and their probabilities
#' @param filePath Path to save the plot
#' @return plot
#' @export
PlotInAndOutAminoAcidSignatures <- function(aminoAcidSig, filePath)  {
  p1 = PlotAminoAcidSignatures(aminoAcidSig, type = "ref")
  p2 = PlotAminoAcidSignatures(aminoAcidSig, type = "alt")
  plot <- cowplot::plot_grid(p2, p1, align = "v", ncol = 1, rel_widths = c(10, 1.5))
  p <- cowplot::plot_grid(plot, nrow = 1, rel_widths = c(10, 1.5))
  if (!is.null(filePath)) {
    ggsave(filePath, plot = p, width = 25, height = 8, dpi =500)
  }
  p
}


