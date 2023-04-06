library(dplyr)

refit.params <- read.table("metadata_formatting/refitting_parameters_updated.csv",header = T,sep = ",")
meta <- read.table("metadata_formatting/metadata.csv",header = T,sep = ",")

metadata <- meta %>% select(PATIENT_ID,SAMPLE_ID,ALLELE_FREQ,DEPTH) %>%
  group_by(SAMPLE_ID) %>%
  mutate(A_D=ALLELE_FREQ*DEPTH) %>%
  summarise_at(.vars = vars(DEPTH,A_D),sum) %>%
  mutate(TP53freq=A_D/DEPTH) %>%
  select(SAMPLE_ID,TP53freq) %>%
  left_join(x = meta,y = .,by="SAMPLE_ID") %>%
  select(PATIENT_ID,SAMPLE_ID,TP53freq) %>%
  distinct()

smoothed.samples <- as.character(refit.params$name[refit.params$smooth == "Y"])

metadata$smooth <- rep("",times=nrow(metadata))
metadata$smooth[which(metadata$SAMPLE_ID %in% smoothed.samples)] <- "TRUE"
metadata$smooth[which(!metadata$SAMPLE_ID %in% smoothed.samples)] <- "FALSE"

metadata$file <- paste0(metadata$SAMPLE_ID,".bam")

write.table(metadata,file = "sample_sheet.tsv",sep = "\t",row.names = F,col.names = T,quote = F)
