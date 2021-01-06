rm(list = ls())
options(stringsAsFactors = F)
library("dplyr")
library("tidyr")
library("Biostrings")
library("ggplot2")
library("cowplot")
library("pals")
source("./functions.R")

cosmicSBS2 <- "~/Dropbox/MutationalSignature2AminoAcids/sigProfiler_SBS_signatures_SBS2.csv"
cosmicSBS3 <- "~/Dropbox/MutationalSignature2AminoAcids/sigProfiler_SBS_signatures_SBS3.csv"
cosmicSBS5 <- "~/Dropbox/MutationalSignature2AminoAcids/sigProfiler_SBS_signatures_SBS5.csv"
cosmicSBS6 <- "~/Dropbox/MutationalSignature2AminoAcids/sigProfiler_SBS_signatures_SBS6.csv"
cosmicSBS13 <- "~/Dropbox/MutationalSignature2AminoAcids/sigProfiler_SBS_signatures_SBS13.csv"

trinucleotideContextDF <- ParseComsicSignatures(filePath = cosmicSBS2)

marginalFrequencies <- c("T" = 0.25, "A" = 0.25, "G" = 0.25, "C" = 0.25)
frameProb <- c(1/3, 1/3, 1/3)

aminoAcidHits <-ConvertTrinucleotideContextToAimnoAcidContext(trinucleotideContextDF,
                                                              marginalFrequencies, 
                                                              frameProb)

aminoAcidSig <- MapAminoAcidHits(aminoAcidHits)


PlotInAndOutAminoAcidSignatures(aminoAcidSig, filePath = "./amino_acid_replacement_signature.pdf")
 
 





