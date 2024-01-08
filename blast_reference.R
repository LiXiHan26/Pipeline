suppressMessages(library("here"))
suppressMessages(library("tidyverse"))
suppressMessages(library("magrittr"))
suppressMessages(library("ape"))
suppressMessages(library("optparse"))



#two functions below are quoted from https://raw.githubusercontent.com/boopsboops/UTILITIES/main/RScripts/tab2fas.R
tab2fas <- function(df,seqcol,namecol){
  df <- as.data.frame(df)
  dtmp <- strsplit(as.character(df[[seqcol]]), "")
  names(dtmp) <- as.character(df[[namecol]])
  dat <- ape::as.DNAbin(dtmp)
  return(dat)
}

fas2tab <- function(dnas) {
  if(class(dnas) == "DNAbin") {
    dnas.list <- as.list(dnas)
    dnas.names <- names(dnas.list)
    dnas.char <- sapply(as.character(dnas.list),paste,collapse="")
    dnas.char.lower <- stringr::str_to_lower(dnas.char)
    dnas.clean <- stringr::str_replace_all(dnas.char.lower,"-|n|\\?","")
    dnas.tib <- tidyr::tibble(label=dnas.names,nucleotides=dnas.clean)
    return(dnas.tib)
  } else stop(writeLines("Object must be APE DNAbin format."))
}

#get args
option_list <- list(
  make_option(c("-i","--input"),type = "character"),
  make_option(c("-o","--output"), type = "character")
  )

#set args
opt <- parse_args(OptionParser(option_list = option_list,add_help_option = FALSE))

#test opt
#opt <- NULL
#opt$input <- "18S_output.tsv"
#opt$output <- "18s_blast_reference.fasta"

#input file
inputfile <- read.delim(file = here(opt$input), 
                        header = TRUE, sep = "\t", quote = "\"", dec = ".", fill = TRUE)
colnames(inputfile) <-  c('seqID','taxid','superkingdom','phylum','class','order','family','genus','species','sequence')
#get blast reference, keep seqID and nucleotides 
output_fasta <- tab2fas(df=inputfile, seqcol = "sequence", namecol = "seqID")

#write blast reference
ape::write.FASTA(output_fasta, file = here(opt$output))



