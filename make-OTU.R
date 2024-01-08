suppressMessages(library("ape"))
suppressMessages(library("tidyverse"))
suppressMessages(library("parallel"))
suppressMessages(library("magrittr"))
suppressMessages(library("here"))
suppressMessages(library("optparse"))

#get args
option_list <- list(
  make_option(c("-s","--sintax"),type = "character"),
  make_option(c("-b","--blast"), type = "character"),
  make_option(c("-t","--taxonomy"), type = "character"),
  make_option(c("-u","--otus"), type = "character"),
  make_option(c("-o","--output"), type = "character")
)

#set args
opt <- parse_args(OptionParser(option_list = option_list,add_help_option = FALSE))

# load blast result
local.db.blast.lib1 <- suppressMessages(suppressWarnings(read_tsv(file=here(opt$blast),guess_max=99999999)))
sseqidLocal <- gsub("ref","",local.db.blast.lib1$sseqidLocal)
sseqidLocal <- gsub("dbj","", sseqidLocal)
sseqidLocal <- gsub("gb","", sseqidLocal)
sseqidLocal <- gsub("[|]","",sseqidLocal)
local.db.blast.lib1$sseqidLocal <- sseqidLocal


#set args
opt <- parse_args(OptionParser(option_list = option_list,add_help_option = FALSE))

#add taxonomy
custom.db <- suppressMessages(suppressWarnings(read_tsv(file=here(opt$taxonomy),guess_max=999999)))


# annotate species names
local.db.blast.lib1 %<>% mutate(sciName=pull(custom.db,species)[match(sseqidLocal,pull(custom.db,Accession))],
                                genus=pull(custom.db,genus)[match(sseqidLocal,pull(custom.db,Accession))],
                                family=pull(custom.db,family)[match(sseqidLocal,pull(custom.db,Accession))],
                                order=pull(custom.db,order)[match(sseqidLocal,pull(custom.db,Accession))],
                                class=pull(custom.db,class)[match(sseqidLocal,pull(custom.db,Accession))],
                                phylum=pull(custom.db,phylum)[match(sseqidLocal,pull(custom.db,Accession))],
                                kingdom=pull(custom.db,kingdom)[match(sseqidLocal,pull(custom.db,Accession))])
# chose "best" hit based on bitscore
# also add scinames
local.db.blast.sorted.lib1 <- local.db.blast.lib1 %>% 
  group_by(qseqid) %>%
  arrange(desc(bitscoreLocal),.by_group=TRUE) %>%
  filter(bitscoreLocal==max(bitscoreLocal)) %>%
  arrange(sciName,.by_group=TRUE) %>%
  mutate(sciName=paste(unique(sciName),collapse="; ")) %>%
  slice(1) %>% 
  ungroup()

# read in taxonomy assignment
tax.ass.df.lib1 <- suppressMessages(suppressWarnings(read_tsv(file=here(opt$sintax),col_names=c("ID","taxa_bs","strand","taxa"),guess_max=999999)))

MOTU.lib1 <- read.csv(file=opt$otus, sep="\t")

MOTU.lib1%<>% 
  mutate(qseqid=pull(tax.ass.df.lib1, ID))%>% 
  #left_join(local.db.blast.sorted.lib1, by="qseqid") %>%
  left_join(tax.ass.df.lib1, by="ID") %>% 
  mutate(isFish=if_else(class=="Actinopteri"|class=="Chondrichthyes"|class=="Cephalaspidomorphi"|class=="Elasmobranchii"|class=="Holocephali"|class=="Myxini"|class=="Sarcopterygii"|class=="Teleostei"|class=="Actinopterygii", TRUE, FALSE)) %>%
  #mutate(isIntersted=if_else(class=="Actinopteri"|class=="Chondrichthyes"|class=="Cephalaspidomorphi"|class=="Elasmobranchii"|class=="Holocephali"|class=="Myxini"|class=="Sarcopterygii"|class=="Teleostei"|class=="Actinopterygii"|order=="Carnivora"|order=="Cetacea", TRUE, FALSE))%>%
  write_csv(file=here(opt$output))
