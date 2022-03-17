## Source colour palettes
source("../../colour_palettes.R")

## Extract abs segment calls from QDNAseq segment data
extract_abs_segs <- function(abs_data = x,genes = NULL){
  sample_n <- length(colnames(abs_data))
  to_use <- fData(abs_data)$use
  all_segs <- data.frame()
  gene_bins <- get_gene_bins(abs_data = abs_data,genes = genes)
  sample_segs <- c()
  for(i in 1:sample_n){
    cn_obj <- abs_data[to_use,i]
    #ploidy <- pData(cn_obj)$ploidy
    sample_name <- colnames(cn_obj)
    segments <- assayDataElement(cn_obj,"segmented")
    gene_segs <- c()
    for(j in 1:nrow(gene_bins)){
      bins <- gene_bins[j,]
      start <- bins$start_idx
      end <- bins$end_idx
      seg_mean <- mean(segments[c(start:end)])
      gene_segs <- cbind(gene_segs,seg_mean)
    }
    gene_segs <- as.data.frame(gene_segs)
    colnames(gene_segs) <- gene_bins$gene
    sample_segs <- rbind(sample_segs,gene_segs)
  }
  rownames(sample_segs) <- colnames(abs_data)
  return(as.data.frame(t(sample_segs)))
}

## get copy number alteration calls
get_cna_calls <- function(data=NULL){
  # COSMIC definitions
  #Gain:
  ## average genome ploidy <= 2.7 AND total copy number >= 5
  ## OR average genome ploidy > 2.7 AND total copy number >= 9
  #Loss:
  ## average genome ploidy <= 2.7 AND total copy number = 0
  ## OR average genome ploidy > 2.7 AND total copy number < ( average genome ploidy - 2.7 )
  cna_calls <- do.call(cbind,lapply(colnames(data),FUN = function(x){
    ploidy <- ploidys[x]
    abs_col <- data[,x]
    if(ploidy > 2.7){
      abs_col <- ifelse(abs_col >= 9,"AMP",ifelse(abs_col < (ploidy - 2.7),"DEL","N"))
    } else {
      abs_col <- ifelse(abs_col >= 5,"AMP",ifelse(abs_col <= 0,"DEL","N"))
    }
    abs_col <- data.frame(sample=abs_col)
    colnames(abs_col) <- x
    return(abs_col)
  }))
  rownames(cna_calls) <- rownames(data)
  return(cna_calls)
}

## Get gene CNA rates from CNA calling from abs data
get_cna_rates <- function(data=NULL,genes=NULL){
  if(!is.null(genes)){
    data <- data[rownames(data) %in% genes,]
    genes <- unique(data$Gene)
  } else {
    genes <- unique(data$Gene)
  }
  amp_rate <- count_amp_focal(data = data)
  del_rate <- count_del_focal(data = data)
  cna_rate_table <- rbind(data.frame(gene=names(amp_rate),rate=amp_rate,group=rep("amplification",times=length(amp_rate)),row.names = NULL),
                          data.frame(gene=names(del_rate),rate=del_rate,group=rep("deletion",times=length(del_rate)),row.names = NULL))
  return(cna_rate_table)
}

## Count focal amp
count_amp_focal <- function(data = NULL){
  amp <- c()
  data <- rownames_to_column(data,var = "Gene")
  for(i in as.character(data$Gene)){
    amped <- sum(data[data$Gene == i,][,-1] == "AMP")
    len <- length(data[data$Gene == i,][,-1])
    amp <- append(amp, (amped / len) * 100)
  }
  names(amp) <- data$Gene
  return(amp)
}
## Count focal del
count_del_focal <- function(data = NULL){
  del <- c()
  data <- rownames_to_column(data,var = "Gene")
  for(i in as.character(data$Gene)){
    deled <- sum(data[data$Gene == i,][,-1] == "DEL")
    len <- length(data[data$Gene == i,][,-1])
    del <- append(del, (deled / len) * 100)
  }
  names(del) <- data$Gene
  return(del)
}

## Plot cna rate date - filtered by gene
plot_cna_rates <- function(data=NULL,title="A plot"){
  p <- ggplot(data) +
    geom_col(data = data,aes(x=gene,y=rate,fill=group),position = "dodge") +
    scale_fill_manual(values = colour_palettes$amplification_deletion) +
    ylab(label = "Rate (%)") +
    labs(title = title) +
    theme_bw() + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(vjust = 0.5,hjust = 1,angle = 90))
  return(p)
}

## Perform fishers extact testings between sample groups
fishers_gene_CNA <- function(data = NULL,sub_a = NULL,sub_b = NULL,genes=NULL,p.adj = "bonferroni"){
  data <- rownames_to_column(data,var = "Gene")
  if(!is.null(genes)){
    data <- data[data$Gene %in% genes,]
    genes <- unique(data$Gene)
  } else {
    genes <- unique(data$Gene)
  }
  sub_a_data <- data[,colnames(data) %in% c("Gene",as.character(sub_a))]
  sub_b_data <- data[,colnames(data) %in% c("Gene",as.character(sub_b))]
  
  sub_a_data_amp <- do.call(rbind,
                            lapply(genes,FUN = function(x){
                              sub_g <- sub_a_data[sub_a_data$Gene == as.character(x),]
                              data.frame(Gene=x,
                                         a=sum(sub_g[,-1] == "AMP"),
                                         Total_a=length(sub_g[,-1]) - sum(sub_g[,-1] == "AMP"),
                                         Mean_a=mean(as.numeric(sub_g[which(sub_g[,-1] == "AMP")+1]),na.rm=T),
                                         row.names = NULL)
                            }))
  
  sub_a_data_del <- do.call(rbind,
                            lapply(genes,FUN = function(x){
                              sub_g <- sub_a_data[sub_a_data$Gene == as.character(x),]
                              data.frame(Gene=x,
                                         a=sum(sub_g[,-1] == "DEL"),
                                         Total_a=length(sub_g[,-1]) - sum(sub_g[,-1] == "DEL"),
                                         Mean_a=mean(as.numeric(sub_g[which(sub_g[,-1] == "DEL")+1]),na.rm=T),
                                         row.names = NULL)
                            }))
  
  sub_b_data_amp <- do.call(rbind,
                            lapply(genes,FUN = function(x){
                              sub_g <- sub_b_data[sub_b_data$Gene == as.character(x),]
                              data.frame(Gene=x,
                                         b=sum(sub_g[,-1] == "AMP"),
                                         Total_b=length(sub_g[,-1]) - sum(sub_g[,-1] == "AMP"),
                                         Mean_b=mean(as.numeric(sub_g[which(sub_g[,-1] == "AMP")+1]),na.rm=T),
                                         row.names = NULL)
                            }))
  
  sub_b_data_del <- do.call(rbind,
                            lapply(genes,FUN = function(x){
                              sub_g <- sub_b_data[sub_b_data$Gene == as.character(x),]
                              data.frame(Gene=x,
                                         b=sum(sub_g[,-1] == "DEL"),
                                         Total_b=length(sub_g[,-1]) - sum(sub_g[,-1] == "DEL"),
                                         Mean_b=mean(as.numeric(sub_g[which(sub_g[,-1] == "DEL")+1]),na.rm=T),
                                         row.names = NULL)
                            }))
  
  data_amp <- merge(x = sub_a_data_amp,y = sub_b_data_amp,by="Gene")
  data_del <- merge(x = sub_a_data_del,y = sub_b_data_del,by="Gene")
  
  data_amp$Mean_a[is.na(data_amp$Mean_a)] <- 0
  data_amp$Mean_b[is.na(data_amp$Mean_b)] <- 0
  
  data_del$Mean_a[is.na(data_del$Mean_a)] <- 0
  data_del$Mean_b[is.na(data_del$Mean_b)] <- 0
  
  data_amp$diff <- abs(data_amp$Mean_a - data_amp$Mean_b)
  data_del$diff <- -abs(data_del$Mean_a - data_del$Mean_b)
  
  data_amp$group <- rep("amp",nrow(data_amp))
  data_del$group <- rep("del",nrow(data_del))
  
  amp.p <- c()
  for(i in 1:nrow(data_amp)){
    amp.p[i] <- fisher.test(x = matrix(unlist(data_amp[i,c(2:3,5:6)]),nrow = 2))$p.value
  }
  del.p <- c()
  for(i in 1:nrow(data_del)){
    del.p[i] <- fisher.test(x = matrix(unlist(data_del[i,c(2:3,5:6)]),nrow = 2))$p.value
  }
  
  data_amp$p.value <- amp.p
  data_amp$p.adj <- p.adjust(p = data_amp$p.value,method = p.adj)
  data_amp <- data_amp[order(data_amp$p.value),]
  
  data_del$p.value <- del.p
  data_del$p.adj <- p.adjust(p = data_del$p.value,method = p.adj)
  data_del <- data_del[order(data_del$p.value),]
  
  out <- list(amplification=data_amp,deletion=data_del)
  return(out)
}

############################################

## Get samples function
getSamples <- function(abs_data = NULL,PATIENT_ID = NULL){
  if(is.null(abs_data)){
    stop("No data provided")
  }
  if(is.null(PATIENT_ID)){
    stop("No patient provided")
  }
  phenoDat <- pData(abs_data)
  samples <- phenoDat$name[phenoDat$PATIENT_ID %in% PATIENT_ID]
  return(samples)
}


## Get seg table
getSegTable<-function(x){
  dat<-x
  sn<-Biobase::assayDataElement(dat,"segmented")
  fd <- Biobase::fData(dat)
  fd$use -> use
  fdfiltfull<-fd[use,]
  sn<-sn[use,]
  segTable<-c()
  for(c in unique(fdfiltfull$chromosome))
  {
    snfilt<-sn[fdfiltfull$chromosome==c]
    fdfilt<-fdfiltfull[fdfiltfull$chromosome==c,]
    sn.rle<-rle(snfilt)
    starts <- cumsum(c(1, sn.rle$lengths[-length(sn.rle$lengths)]))
    ends <- cumsum(sn.rle$lengths)
    lapply(1:length(sn.rle$lengths), function(s) {
      from <- fdfilt$start[starts[s]]
      to <- fdfilt$end[ends[s]]
      segValue <- sn.rle$value[s]
      c(fdfilt$chromosome[starts[s]], from, to, segValue)
    }) -> segtmp
    segTableRaw <- data.frame(matrix(unlist(segtmp), ncol=4, byrow=T),stringsAsFactors=F)
    segTable<-rbind(segTable,segTableRaw)
  }
  colnames(segTable) <- c("chromosome", "start", "end", "segVal")
  segTable
}


## Performs wilcox test for copy changes between groups
cna_median_copies <- function(data = NULL,sub_a = NULL,sub_b = NULL,genes=NULL,min_obs=1,p.adj = "bonferroni"){
  data <- rownames_to_column(data,var = "Gene")
  if(!is.null(genes)){
    data <- data[data$Gene %in% genes,]
    genes <- unique(data$Gene)
  } else {
    genes <- unique(data$Gene)
  }
  sub_a_data <- data[,colnames(data) %in% c("Gene",as.character(sub_a))]
  sub_b_data <- data[,colnames(data) %in% c("Gene",as.character(sub_b))]
  
  p.val <- c()
  obs_a <- c()
  obs_b <- c()
  median_a <- c()
  median_b <- c()
  change <- c()
  notes <- c()
  for(i in as.character(genes)){
    a_values <- as.numeric(sub_a_data[sub_a_data$Gene == i,-1])
    b_values <- as.numeric(sub_b_data[sub_b_data$Gene == i,-1])
    obs_a <- append(obs_a,length(a_values[!is.na(a_values)]))
    obs_b <- append(obs_b,length(b_values[!is.na(b_values)]))
    if(length(a_values[!is.na(a_values)]) > min_obs & length(b_values[!is.na(b_values)]) > min_obs){
      median_a <- append(median_a,median(a_values,na.rm = T))
      median_b <- append(median_b,median(b_values,na.rm = T))
      change <- append(change,
                       ifelse(median(a_values,na.rm = T) < median(b_values,na.rm = T),"GAIN","LOSS"))
      stat <- wilcox.test(x = a_values,y = b_values)
      p.val <- append(p.val,stat$p.value)
      
      notes <- append(notes,NA)
    } else {
      p.val <- append(p.val,NA)
      median_a <- append(median_a,NA)
      median_b <- append(median_b,NA)
      change <- append(change,NA)
      notes <- append(notes,"Too few observations")
    }
  }
  amp_copy_data <- data.frame(Gene=genes,
                              obsA=obs_a,
                              obsB=obs_b,
                              medianA=median_a,
                              medianB=median_b,
                              direction=change,
                              pVal=p.val,
                              qVal=p.adjust(p.val,method = p.adj),
                              notes=notes)
  amp_copy_data <- amp_copy_data[order(amp_copy_data$pVal),]
  return(amp_copy_data)
}
## get patient-specific amp changes
get_patient_changes <- function(data=NULL,genes=NULL,patients=NULL){
  #data <- rownames_to_column(data,var = "Gene") %>%
            #mutate(Gene = as.character(Gene))
  if(!is.null(genes)){
    data <- data[rownames(data) %in% genes,]
    genes <- unique(as.character(rownames(data)))
  } else {
    genes <- unique(as.character(rownames(data)))
  }
  amp_changes <- as.data.frame(t(data[rownames(data) %in% genes,
                                           as.character(meta.data$SAMPLE_ID[meta.data$paired == TRUE])]))
  amp_changes <- rownames_to_column(amp_changes,var = "SAMPLE_ID") %>%
    mutate(PATIENT_ID=meta.data$PATIENT_ID[match(SAMPLE_ID,meta.data$SAMPLE_ID)]) %>%
    mutate(group=meta.data$group[match(SAMPLE_ID,meta.data$SAMPLE_ID)]) %>%
    relocate(PATIENT_ID,SAMPLE_ID,group) %>%
    pivot_longer(names_to = "gene",cols = 4:ncol(.))
  amp_changes <- amp_changes[amp_changes$PATIENT_ID %in% patients,]
  return(amp_changes)
}

## Get CN bins for gene bed
get_gene_bins <- function(abs_data = x,genes = y){
  to_use <- fData(abs_data)$use
  cn_obj <- abs_data[to_use,]
  bin_pos <- fData(cn_obj)[,c("chromosome","start","end")]
  gene_bins <- c()
  for(i in 1:nrow(genes)){
    chr <- as.character(genes[i,1])
    start <- as.numeric(genes[i,2])
    end <- as.numeric(genes[i,3])
    g <- as.character(genes[i,4])
    gene_chr_pos <- bin_pos[bin_pos$chromosome == chr,]
    ## Remove genes outside of binned regions
    if(start < gene_chr_pos$start[1] & end < gene_chr_pos$start[1]){
      cat(paste0(g," outside binned genome... skipping\n"))
      next()
    }
    if(start > gene_chr_pos$end[nrow(gene_chr_pos)] & end > gene_chr_pos$end[nrow(gene_chr_pos)]){
      cat(paste0(g," outside binned genome... skipping\n"))
      next()
    }
    min_start <- min(which(min(abs(gene_chr_pos$start - start)) == abs(gene_chr_pos$start - start)))
    min_end <- max(which(min(abs(gene_chr_pos$end - end)) == abs(gene_chr_pos$end - end)))
    # min and max are used for min_start and min_end due to a few (< 5) genes having end and start points equi-distance from bins
    # Additionally, genes spanning the first or last available bin are restricted to only a complete bin as no information is available
    # outside of those bins for cn state
    if(gene_chr_pos$start[min_start] > start & min_start != 1){
      min_start <- min_start - 1
    }
    if(gene_chr_pos$end[min_end] < end & min_end != length(gene_chr_pos$end)){
      min_end <- min_end + 1
    }
    
    index_min <- which(bin_pos$chromosome == chr & bin_pos$start == gene_chr_pos$start[min_start])
    index_max <- which(bin_pos$chromosome == chr & bin_pos$end == gene_chr_pos$end[min_end])
    start.dist <- start - gene_chr_pos$start[min_start]
    end.dist <- gene_chr_pos$end[min_end] - end
    ## Removes rows where either the gene start or end fall outside available bins
    if(end.dist < 0 | start.dist < 0){
      cat(paste0(g," single end out of range... skipping\n"))
      next()
    }
    gene_pos <- data.frame(start_idx=index_min,
                           end_idx=index_max,
                           gene=g,
                           start.dist=start.dist,
                           end.dist=end.dist)
    gene_bins <- rbind(gene_bins,gene_pos)
  }
  return(gene_bins)
}


## TOPGo gene selection function
geneSelection <- function(x){
  x <- x < 0.05
  return(x)
}
## Count CNA events
count_CNA_events <- function(data=NULL){
  cna_calls <- do.call(rbind,lapply(colnames(data),FUN = function(x){
    ploidy <- ploidys[x]
    seg_tab <- getSegTable(data[,x])
    abs_col <- seg_tab$segVal
    if(ploidy > 2.7){
      abs_col <- ifelse(abs_col >= 9,"AMP",ifelse(abs_col < (ploidy - 2.7),"DEL","N"))
    } else {
      abs_col <- ifelse(abs_col >= 5,"AMP",ifelse(abs_col <= 0,"DEL","N"))
    }
    amp_events <- sum(abs_col == "AMP")
    del_events <- sum(abs_col == "DEL")
    all_events <- amp_events + del_events
    abs_col <- data.frame(sample=x,amp_events=amp_events,del_events=del_events,total_events=all_events)
    return(abs_col)
  }))
  return(cna_calls)
}
## Get CN bins for gene bed
get_cytoband_bins <- function(abs_data = x,cyto = y){
  to_use <- fData(abs_data)$use
  cn_obj <- abs_data[to_use,]
  bin_pos <- fData(cn_obj)[,c("chromosome","start","end")]
  cyto_bins <- c()
  for(i in 1:nrow(cyto)){
    chr <- cyto[i,1]
    start <- cyto[i,2]
    end <- cyto[i,3]
    g <- cyto[i,5]
    cyto_chr_pos <- bin_pos[bin_pos$chromosome == chr,]
    ## Due to bin exclusion, some cytobands are not represented by overlapping bins so are dropped
    if(end < min(cyto_chr_pos$start)){
      next
    }
    min_start <- min(which(min(abs(cyto_chr_pos$start - start)) == abs(cyto_chr_pos$start - start)))
    min_end <- max(which(min(abs(cyto_chr_pos$end - end)) == abs(cyto_chr_pos$end - end)))
    # min and max are used for min_start and min_end due to a few (< 5) genes having end and start points equi-distance from bins
    # Additionally, genes spanning the first or last available bin are restricted to only a complete bin as no information is available
    # outside of those bins for cn state
    if(cyto_chr_pos$start[min_start] > start & min_start != 1){
      min_start <- min_start - 1
    }
    if(cyto_chr_pos$end[min_end] < end & min_end != length(cyto_chr_pos$end)){
      min_end <- min_end + 1
    }
    
    index_min <- which(bin_pos$chromosome == chr & bin_pos$start == cyto_chr_pos[min_start,2])
    index_max <- which(bin_pos$chromosome == chr & bin_pos$end == cyto_chr_pos[min_end,3])
    cyto_pos <- data.frame(start_idx=index_min,end_idx=index_max,chr=chr,cytoband=g)
    cyto_bins <- rbind(cyto_bins,cyto_pos)
  }
  return(cyto_bins)
}
## Extract abs segment calls from QDNAseq segment data
extract_broad_CNA <- function(abs_data = x,cyto = NULL,threshold = 0.8){
  sample_n <- length(colnames(abs_data))
  to_use <- fData(abs_data)$use
  all_segs <- data.frame()
  cyto_bins <- get_cytoband_bins(abs_data = abs_data,cyto = cyto)
  # Include all bins in length
  cyto_bins$len <- cyto_bins$end_idx - cyto_bins$start_idx + 1
  # Calculate nubmer of bp per bins in each cytobands
  cyto_bins$bpPerBin <- cytoband$length[match(cyto_bins$cytoband,cytoband$band)] / cyto_bins$len
  # Exclude cytobands with too few bins or bpPerBin above 90% percentile
  cyto_bins <- cyto_bins[cyto_bins$bpPerBin < quantile(cyto_bins$bpPerBin,probs=c(0.9)),]
  sample_segs <- c()
  for(i in 1:sample_n){
    cn_obj <- abs_data[to_use,i]
    #ploidy <- pData(cn_obj)$ploidy
    sample_name <- colnames(cn_obj)
    ploidy <- ploidys[sample_name]
    segments <- as.data.frame(assayDataElement(cn_obj,"segmented"))
    cyto_seg <- c()
    for(j in 1:nrow(cyto_bins)){
      bins <- cyto_bins[j,]
      start <- bins$start_idx
      end <- bins$end_idx
      cyto_segs <- segments[c(start:end),]
      seg_len <- length(cyto_segs)
      if(ploidy > 2.7){
        pct_amp <- sum(cyto_segs >= 9) / seg_len
        pct_del <- sum(cyto_segs < (ploidy - 2.7)) / seg_len
        cyto_seg_value <- ifelse(pct_amp >= threshold,"AMP",ifelse(pct_del >= threshold,"DEL","N"))
      } else {
        pct_amp <- sum(cyto_segs >= 5) / seg_len
        pct_del <- sum(cyto_segs <= 0) / seg_len
        cyto_seg_value <- ifelse(pct_amp >= threshold,"AMP",ifelse(pct_del >= threshold,"DEL","N"))
      }
      cyto_seg <- cbind(cyto_seg,cyto_seg_value)
      
    }
    cyto_seg <- as.data.frame(cyto_seg)
    colnames(cyto_seg) <- paste0(cyto_bins$chr,":",cyto_bins$cytoband)
    sample_segs <- rbind(sample_segs,cyto_seg)
  }
  rownames(sample_segs) <- colnames(abs_data)
  return(as.data.frame(t(sample_segs)))
}
## Count broad amplifcation events
count_amp_broad <- function(data = NULL){
  amp <- c()
  data <- rownames_to_column(data,var = "cytoband")
  for(i in as.character(data$cytoband)){
    amped <- sum(data[data$cytoband == i,][,-1] == "AMP")
    len <- length(data[data$cytoband == i,][,-1])
    amp <- append(amp, (amped / len) * 100)
  }
  amp <- data.frame(cytoband=data$cytoband,rate=as.numeric(amp))
  amp$chr <- str_split(amp$cytoband,pattern = ":",simplify = T,n = 2)[,1]
  amp$chr <- factor(amp$chr,levels = unique(amp$chr))
  return(amp)
}
## Plot broad amplification events
plot_cna_rates_broad <- function(data=NULL,title="A plot"){
  p <-  ggplot(data) +
    geom_col(data = data,aes(x=cytoband,y=rate),position = "dodge") +
    scale_y_continuous(limits = c(0,100)) +
    scale_fill_manual(values = colour_palettes$AMP_DEL) +
    ylab(label = "Rate (%)") +
    labs(title = title) +
    theme_bw() + 
    theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90))
  return(p)
}
## Plot CNV event means
cnv_event_means <- function(data=x){
  plot.data <- data
  pvals <- c()
  n <- c("arx","rlps")
  groups <- c("amp_events","del_events","total_events")
  for(i in groups){
    pvals <- append(pvals,
                    t.test(plot.data$value[plot.data$name == i & plot.data$group == n[1]],plot.data$value[plot.data$name == i & plot.data$group == n[2]])$p.value)
  }
  names(pvals) <- groups
  p <- ggplot(plot.data,aes(x = group,y = value)) +
    geom_point(color=alpha("grey20",0.5),position="jitter") +
    geom_violin(aes(fill=group),alpha=0.6) +
    #geom_density_ridges(aes(fill=group),binwidth=10) +
    labs(title = "CNV event summarys") +
    facet_wrap(. ~ name,nrow = length(groups)) +
    coord_flip() +
    theme_bw()
  
  plim <- layer_scales(p)$y$get_limits()[2] - layer_scales(p)$y$get_limits()[2]*0.05
  p + geom_text(data = data.frame(name=groups,label=paste0("p = ",signif(pvals,digits = 3))),aes(y = plim,x = 2.3,label = label))
}