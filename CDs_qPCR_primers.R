# install required packages
install.packages(c("pacman", "devtools", "BiocManager")
# load utilities
devtools::source_gist(id = )devtools::source_gist("7f63547158ecdbacf31b54a58af0d1cc", filename = "Util.R")
# install BioConductor packages (requires R>=3.5)
if(getRversion() >= "3.5.0") {
  BiocManager::install(c("GenomicFeatures", "GenomicRanges","plyranges", "Rsamtools"))
  pacman::p_load(plyranges)
}

# install and load additional packages
pacman::p_load(char = c("seqinr", "tidyverse", "GenomicFeatures", "GenomicRanges", "Rsamtools", "xlsx" )) #  "plyranges", "BiocManager" requires R>=3.5


#### Extract CDs from transcriptome ####
# Read in list of genes
gene_list <- readxl::read_excel("genes for primers.xlsx")
# # Read in transcriptome
# transcriptome <- seqinr::read.fasta("../QTL_analysis/data/CASSAB_K91_Reference_FINAL.fasta")
# # Extract sequences
# qtl_transcripts <- transcriptome[map_lgl(transcriptome, ~grepl(transcripts_to_search, attr(., "Annot")))]
# # Save as fasta
# seqinr::write.fasta(qtl_transcripts, names(qtl_transcripts), 
#                     glue::glue("./qPCR_transcripts.fasta"))


#### Load genome files  ####
# Load genome file
genome_fa <- "../../L_culinaris_genome/Lensculinaris_genome_v1.2.fasta"
if (!file.exists(paste0(genome_fa, ".fai"))) indexFa(genome_fa)
fa <- open(FaFile(genome_fa))
# load gene models file
gff_file <- "../../L_culinaris_genome/Lensculinaris_1.2b_genes_fixed_phase.gff3"
# gff3 <- read_tsv(gff_file, col_names = c("seqid", "source", "type", "start", "end", "score", "strand",
#                                          "phase", "attributes")) 

# gff3 %>% filter(!grepl(paste(discard_pattern, collapse = "|"), attributes))
txdb <- makeTxDbFromGFF(gff_file, format="gff3", dataSource = "v1.2 UoS",
                        organism = "Lens culinaris",
                        taxonomyId = 3864)

# gene_ids <- sub("\\.\\d+$", "", qtl_region_tx$tx_name)
gene_ids <- gene_list$Genes
# gene_pattern <- paste(glue::glue('ID={gene_ids};'), collapse = "|")
# qtl_genes <- gff3 %>% mutate(attributes=gsub('\"', '', attributes, fixed = TRUE)) %>% 
#   separate(attributes, c("ID", "Name", "Description"), sep = ";.+?=") %>%
#   filter(type=="gene", Name %in% gene_ids) #%>% 

#### Extract exons and CDs #####

gene_txs <- transcripts(txdb,filter=list(tx_name=gene_list$Genes),  use.names=TRUE)
gene_cds <- cds(txdb,columns=c("tx_name", "cds_name"),
                filter=list(tx_name=gene_list$Genes),  use.names=TRUE)
txdb[gene_list$Genes]
gene_exons <- exons(txdb,columns=c("exon_name", "exon_id", "tx_name"), filter=list(tx_name=gene_list$Genes),  use.names=TRUE) %>% 
  mutate(exon_num=as.numeric(sub(".+exon(\\d+)", "\\1", exon_name)),
         tx_name=as.character(tx_name))

gene_txs$tx_name
# cds_seqs <- extractTranscriptSeqs(fa,
#                         cdsByOverlaps(txdb, gene_txs))
# tx_seqs <- getSeq(fa, gene_txs) #   cds_seqs[snp_tx$tx_name]
# names(tx_seqs) <- gene_txs$tx_name
exon_seqs <- getSeq(fa, gene_exons)
names(exon_seqs) <- gene_exons$exon_name

# concatenate exons to cds
cds_seqs <- DNAStringSet()
for (tx in gene_txs$tx_name){
  # tx=gene_txs$tx_name[1]
  exon_table <- as_tibble(gene_exons) %>% dplyr::filter(tx_name==tx) %>% 
    arrange(exon_num) 
  cds_seq <- DNAString(x="", start=1, nchar=NA)
  for (j in exon_table$exon_name){
    cds_seq <- xscat(cds_seq, exon_seqs[j])
  }
  names(cds_seq) <- tx
  cds_seqs <- append(cds_seqs, cds_seq)  
}
 

# names(exon_seqs) <- gene_exons$exon_name

#### Primer Design ####
# download a local copy or primer3
# system2("git", args = c("clone", "https://github.com/IdoBar/SSR_pRimer_design.git"))
download.file("https://github.com/IdoBar/SSR_pRimer_design/archive/master.zip", destfile = "primer_design.zip")
unzip("primer_design.zip")
# Load Primer3 functions (make sure primer3_core.exe and settings path are configured)
primer3_home <- "./SSR_pRimer_design-master/primer3/"
# primer3 <- file.path(primer3_home,"bin", "primer3_core.exe")
source(file.path(primer3_home, "callPrimer3.R"))

# Summarise exon junctions in each cds
exon_junctions <- as_tibble(gene_exons) %>% # filter(strand=="+") %>% 
  group_by(tx_name) %>% arrange(exon_num, .by_group = TRUE) %>% 
  mutate(junc_pos=cumsum(width)) %>% dplyr::slice(1:n()-1) %>%
  summarise(junc_list=paste(junc_pos, collapse = " ")) %>% 
  right_join(as_tibble(gene_txs)['tx_name']) %>% 
  mutate(junc_list=if_else(is.na(junc_list), "", junc_list))

# verify that both objects are in the same order
names(cds_seqs)==exon_junctions$tx_name


# Look at the following paper for instructions for primer design: http://onlinelibrary.wiley.com/doi/10.1002/bmb.20461/full
primers <- mapply(.callP3NreadOrg, seq=cds_seqs, name = names(cds_seqs), 
                  junction_list = exon_junctions$junc_list,
                  #sequence_target = paste(blast_GOI$QueryEnd, 1, sep=","),
                  report = filedate(paste("primer3_report", names(cds_seqs), sep="_"), ".txt", 
                                    outdir ="primer3_reports") ,
                  MoreArgs=list(size_range='100-150 80-250', Tm=c(58,60,62), diff_TM=1,
                                primer3_dir = primer3_home, primer_num = 5,
                                settings="primer3_Sajitha_qPCR_settings.txt")
                  ,SIMPLIFY = FALSE,USE.NAMES = TRUE)
# USE.NAMES = FALSE
primers_table <- bind_rows(primers[!is.na(primers)]) %>% mutate(Primer_name=paste(sub("\\.\\d+$", "", Seq_ID), i,  sep = "_"), exon_junction="one of the primers on exon-exon boundary")

# missing_cds <- exon_junctions %>% dplyr::filter(!tx_name %in% unique(primers_table$Seq_ID), junc_list!="") %>% mutate(seq=as.character(cds_seqs[tx_name]), first_junc=as.numeric(sub(" .+", "", junc_list)))
# Add [] at exon junctions to design primers around them
# stringi::stri_sub(missing_cds$seq , missing_cds$first_junc, missing_cds$first_junc-1) <- "[]"
# find for which transcripts we couldn't find primers on the exon-exon boundary
missing_cds <- exon_junctions %>% dplyr::filter(!tx_name %in% unique(primers_table$Seq_ID)) %>% mutate(target_list=gsub("(\\d+)", "\\1,10", junc_list))

# design primers for the "missing" ones
missing_primers <- mapply(.callP3NreadOrg, seq=cds_seqs[missing_cds$tx_name], 
                          name = names(cds_seqs[missing_cds$tx_name]),
                          sequence_target = missing_cds$target_list,
      report = filedate(paste("primer3_report", missing_cds$tx_name, sep="_"), ".txt", 
                                            outdir ="primer3_reports") ,
          MoreArgs=list(size_range='100-150 80-250', Tm=c(58,60,62), diff_TM=1,
                                        primer3_dir = primer3_home, primer_num = 5,
                                        settings="primer3_Sajitha_qPCR_settings.txt")
                          ,SIMPLIFY = FALSE,USE.NAMES = TRUE)

# combine all primers in one table
primers_table <- bind_rows(missing_primers[!is.na(missing_primers)]) %>% mutate(Primer_name=paste(sub("\\.\\d+$", "", Seq_ID), i,  sep = "_"), exon_junction="PCR product spanning exon-exon boundary") %>% bind_rows(primers_table)
primers_table$exon_junction[primers_table$Seq_ID=="Lc13441.1"] <- "No exon-exon boundaries"

#### Save results ####
# Save primer table to file
outfile <- filedate("qPCR_genes_CDs_primers", ".xlsx", dateformat = FALSE)
write.xlsx(primers_table, outfile,row.names = FALSE, sheetName = "primers")

# Add exon information
write.xlsx(as.data.frame(gene_exons), outfile,row.names = FALSE, 
           sheetName = "exon_info", append = TRUE)

# Save transcripts to file
# seqinr::write.fasta(as.list(translate(cds_seqs)), names = names(tx_seqs),
#           file.out = filedate(glue::glue("{stacks_name}_QTL_genes_cds"), 
#                                         ".faa", "QTL_results"))
seqinr::write.fasta(as.list(cds_seqs), names = names(cds_seqs),
                    file.out = filedate(glue::glue("qPCR_genes_CDs"), 
                                        ".fna"))


