meta.data <- read.table("../../copy_number_signatures/britroc_30kb_signature_data_meta.tsv",header = TRUE,sep = "\t")
signatures <- read.table("../../copy_number_signatures/britroc_30kb_signature_data_sig_exposures.tsv",
                         header = TRUE,sep = "\t",row.names = 1,check.names = F)
tissue_mapping <- read.table("data/tissue_mapping_V2.txt",header = TRUE,sep = "\t")
site_table <- read.table(file = "data/site_table.tsv",header = T,sep = "\t",stringsAsFactors = F)

site_table_filt <- site_table[site_table$PATIENT_ID %in% c(81,173,263,264,272),]

combined <- site_table_filt %>% 
    mutate(has_sWGS_Bam = ifelse(SAMPLE_ID %in% samplesheet$SAMPLE_ID,TRUE,FALSE)) %>%
    mutate(has_abs_cnProfile = ifelse(SAMPLE_ID %in% swgs$sample,TRUE,FALSE)) %>%
    mutate(has_Sigs7 = ifelse(SAMPLE_ID %in% colnames(signatures),TRUE,FALSE)) %>%
    select(-c(fk_histological_id,fk_block_id)) %>%
    relocate(SAMPLE_ID,.after = "PATIENT_ID") %>%
    mutate(PATIENT_ID = paste0("BRITROC-",PATIENT_ID)) %>%
    left_join(.,meta.data,by=c("PATIENT_ID","SAMPLE_ID")) %>%
    select(-c("TP53freq")) %>%
    left_join(.,samplesheet,by=c("PATIENT_ID","SAMPLE_ID"))

combined <- combined %>%
            select(-c("group","paired","notes","file"))

swgs_filt <- swgs[swgs$sample %in% combined$SAMPLE_ID,]
signatures_filt <- signatures[,colnames(signatures) %in% combined$SAMPLE_ID]

write.table(combined,"britroc_brainMets_metadata.tsv",quote = F,sep = "\t",
            row.names = F,col.names = T)

write.table(swgs_filt,"britroc_brainMets_absCN_segTable.tsv",quote = F,sep = "\t",
            row.names = F,col.names = T)

write.table(signatures_filt,"britroc_brainMets_ov7_cnSigs.tsv",quote = F,sep = "\t",
            row.names = T,col.names = T)