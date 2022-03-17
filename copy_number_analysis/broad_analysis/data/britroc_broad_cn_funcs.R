## Britroc broad cn analysis functions

get_cytoband_bins <- function(data = x,cyto = y,percentile=0.9){
    to_use <- fData(data)$use
    cn_obj <- data[to_use,]
    bin_pos <- fData(cn_obj)[,c("chromosome","start","end")]
    cyto_bins <- c()
    for(i in 1:nrow(cyto)){
        chr <- as.character(cyto[i,1])
        start <- as.numeric(cyto[i,2])
        end <- as.numeric(cyto[i,3])
        g <- as.character(cyto[i,5])
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
        
        index_min <- which(bin_pos$chromosome == as.character(chr) & bin_pos$start == cyto_chr_pos[min_start,2])
        index_max <- which(bin_pos$chromosome == chr & bin_pos$end == cyto_chr_pos[min_end,3])
        cyto_pos <- data.frame(start_idx=index_min,end_idx=index_max,chr=chr,cytoband=g)
        cyto_bins <- rbind(cyto_bins,cyto_pos)
    }
    # Include all bins in length
    cyto_bins$len <- cyto_bins$end_idx - cyto_bins$start_idx + 1
    # Calculate nubmer of bp per bins in each cytobands
    cyto_bins$bpPerBin <- cyto$length[match(cyto_bins$cytoband,cyto$band)] / cyto_bins$len
    # Exclude cytobands with too few bins or bpPerBin above 90% percentile
    cyto_bins <- cyto_bins[cyto_bins$bpPerBin < quantile(cyto_bins$bpPerBin,probs=c(percentile)),]
    return(cyto_bins)
}

# Extract abs segment calls from QDNAseq segment data
extract_broad_CNA <- function(data = x,cyto = NULL,threshold = 0.8,percentile=0.9,ploidys=NULL){
    sample_n <- length(colnames(data))
    to_use <- fData(data)$use
    data <- data[to_use,]
    cyto_bins <- get_cytoband_bins(data = data,cyto = cyto,percentile=percentile)
    segments <- as.data.frame(assayDataElement(data,"segmented"))
    cyto_seg <- c()
    for(i in 1:nrow(cyto_bins)){
        bins <- cyto_bins[i,]
        start <- bins$start_idx
        end <- bins$end_idx
        segment_data <- segments[c(start:end),]
        called_segs <- call_cna(data = segment_data,ploidys = ploidys)
        amp_rate <- colSums(called_segs == "AMP") / nrow(called_segs)
        del_rate <- colSums(called_segs == "DEL") / nrow(called_segs)
        cyto_seg_value <- ifelse(amp_rate > threshold, "AMP",ifelse(del_rate >= threshold,"DEL","N"))
        cyto_seg <- cbind(cyto_seg,cyto_seg_value)
    }
    cyto_seg <- as.data.frame(cyto_seg)
    colnames(cyto_seg) <- paste0(cyto_bins$chr,":",cyto_bins$cytoband)
    rownames(cyto_seg) <- colnames(data)
    as.data.frame(t(cyto_seg))
}

## Count amp and del events
count_broad <- function(data = NULL,cohort=NULL){
    amp <- data.frame()
    del <- data.frame()
    data <- rownames_to_column(data,var = "cytoband")
    for(i in as.character(data$cytoband)){
        amped <- sum(data[data$cytoband == i,][,-1] == "AMP")
        len <- length(data[data$cytoband == i,][,-1])
        ampneg <- len - amped
        amppct <- (amped / len) * 100
        amp <- rbind(amp,data.frame(cytoband=i,
                                    rate=as.numeric(amppct),
                                    p=amped,
                                    n=ampneg,
                                    group="amp",
                                    stringsAsFactors = F))
    }
    amp$cohort <- rep(cohort,times=nrow(amp))
    amp$chr <- str_split(amp$cytoband,pattern = ":",simplify = T,n = 2)[,1]
    amp$chr <- factor(amp$chr,levels = unique(amp$chr))
    
    for(i in as.character(data$cytoband)){
        deled <- sum(data[data$cytoband == i,][,-1] == "DEL")
        len <- length(data[data$cytoband == i,][,-1])
        delneg <- len - deled
        delpct <- (deled / len) * 100
        del <- rbind(del,data.frame(cytoband=i,
                                    rate=as.numeric(delpct),
                                    p=deled,
                                    n=delneg,
                                    group="del",
                                    stringsAsFactors = F))
    }
    del$cohort <- rep(cohort,times=nrow(del))
    del$chr <- str_split(del$cytoband,pattern = ":",simplify = T,n = 2)[,1]
    del$chr <- factor(del$chr,levels = unique(del$chr))
    return(list(amp=amp,del=del))
}

## Call CNA events
call_cna <- function(data = NULL,ploidys = ploidys){
    cnames <- colnames(data)
    data_c <- do.call(cbind,lapply(cnames,FUN = function(x){
        data_sub <- data[,colnames(data) == x]
        ploidy <- ploidys[x]
        if(ploidy > 2.7){
            data_sub <- case_when(
                data_sub < (ploidy - 2.7) ~ "DEL",
                data_sub >= (ploidy - 2.7) & data_sub < ploidy - 1 ~ "LOSS",
                data_sub >= ploidy - 1 & data_sub < ploidy + 1 ~ "N",
                data_sub >= ploidy + 1 & data_sub < 9 ~ "GAIN",
                data_sub >= 9 ~ "AMP"
            )
        } else {
            data_sub <- case_when(
                data_sub <= 0 ~ "DEL",
                data_sub >  0 & data_sub < ploidy - 1 ~ "LOSS",
                data_sub >= ploidy - 1 & data_sub < ploidy + 1 ~ "N",
                data_sub >= ploidy + 1 & data_sub < 5 ~ "GAIN",
                data_sub >= 5 ~ "AMP"
            )
        }
        data_sub
    }))
    colnames(data_c) <- cnames
    rownames(data_c) <- rownames(data)
    as.data.frame(data_c)
}

## Plot cna_rates
plot_cna_rates_broad <- function(data=NULL,title="A plot"){
    p <-  ggplot(data) +
        geom_col(data = data,aes(x=cytoband,y=rate),position = "dodge") +
        scale_y_continuous(limits = c(0,100)) +
        ylab(label = "Rate (%)") +
        labs(title = title) +
        theme_bw() +
        theme(axis.title.x = element_blank(),axis.text.x = element_text(angle = 90))
    return(p)
}

## Count CNA events
count_CNA_events_broad <- function(data=NULL,ploidys = NULL){
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

## Plot CNV event means
cnv_event_means_broad <- function(data=x,paired=FALSE){
    plot.data <- data
    p <- ggplot(plot.data,aes(x = group,y = value)) +
        geom_point(color=alpha("grey20",0.5),position="jitter") +
        geom_violin(aes(fill=group),alpha=0.6) +
        facet_wrap(. ~ name,nrow = 3,labeller = labeller(name = c(amp_events="amplifications",
                                                                  del_events="deletions",
                                                                  total_events="all"))) +
        scale_fill_manual(values = colour_palettes$diagnosis_relapse) +
        coord_flip() +
        theme_bw()
    
    p
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

## Get bin cna (broad calls)
get_bin_cna <- function(abs_data = NULL,ploidys=NULL){
    segments <- assayDataElement(abs_data,elt = "segmented")
    segments_c <- call_cna(data = segments,ploidys = ploidys)
    freq <- count_cna(data = segments_c)
    freq
}

## Count types of CNA events
count_cna <- function(data = NULL){
    cnames <- colnames(data)
    rname <- rownames(data)
    segments_a <- as.data.frame(t(apply(data,MARGIN = 1,FUN = function(x){
        N <- sum(x == "N") / length(cnames)
        G <- sum(x == "GAIN") / length(cnames)
        A <- sum(x == "AMP") / length(cnames)
        L <- sum(x == "LOSS") / length(cnames)
        D <- sum(x == "DEL") / length(cnames)
        val <- c(neutral=N,gain=G,amplification=A,loss=-L,deletion=-D)
        val
    })))
    rownames(segments_a) <- rname
    segments_a <- segments_a %>%
        rownames_to_column(var = "pos") %>%
        mutate(chr = str_split(pos,pattern = ":",n = 2,simplify = T)[,1]) %>%
        mutate(chr = factor(chr,levels = str_sort(unique(chr),numeric = T))) %>%
        mutate(pos = as.numeric(str_split(string = str_split(pos,pattern = ":",n = 2,simplify = T)[,2],
                                          pattern = "-",simplify = T)[,1])) %>%
        dplyr::select(chr,pos,neutral,gain,amplification,loss,deletion)
    segments_a
}

## calculate fishers exact for CNA event proportions
calc_broad_freq_fishers <- function(a=NULL,b=NULL){
    a_amp <- a$amp
    b_amp <- b$amp
    a_del <- a$del
    b_del <- b$del
    amp_fishers <- c()
    for(i in seq_len(nrow(a_amp))){
        a_amp_l <- a_amp[i,]
        b_amp_l <- b_amp[i,]
        a_P <- a_amp_l$p
        a_N <- a_amp_l$n
        b_P <- b_amp_l$p
        b_N <- b_amp_l$n
        amp_fishers <- append(amp_fishers,fisher.test(matrix(c(a_P,a_N,b_P,b_N),nrow = 2))$p.value)
    }
    amp_out <- cbind(cytoband=a_amp$cytoband,
                     a_amplifcation=a_amp$p,
                     a_nonAmplification=a_amp$n,
                     b_amplification=b_amp$p,
                     b_nonAmplification=b_amp$n,
                     pval=amp_fishers,
                     group=a_amp$group)
    del_fishers <- c()
    for(i in seq_len(nrow(a_del))){
        a_del_l <- a_del[i,]
        b_del_l <- b_del[i,]
        a_P <- a_del_l$p
        a_N <- a_del_l$n
        b_P <- b_del_l$p
        b_N <- b_del_l$n
        del_fishers <- append(del_fishers,fisher.test(matrix(c(a_P,a_N,b_P,b_N),nrow = 2))$p.value)
    }
    del_out <- cbind(cytoband=a_del$cytoband,
                     a_deletions=a_del$p,
                     a_nonDeletions=a_del$n,
                     b_deletions=b_del$p,
                     b_nonDeletions=b_del$n,
                     pval=del_fishers,
                     group=a_del$group)
    
    out <- as.data.frame(rbind(amp_out,del_out))
    out <- out[order(out$pval),]
    out
}

## Get named vector of ploidy values
get_named_ploidy <- function(data = NULL){
    pl <- pData(data)$ploidy
    names <- pData(data)$name
    names(pl) <- names
    pl
}

get_paired_genome_diffs <- function(data,method="mean"){
    samples <- colnames(data)
    sub_meta <- meta.data[meta.data$SAMPLE_ID %in% samples,]
    sample_by_pat_list <- split(sub_meta$SAMPLE_ID,f=as.factor(sub_meta$PATIENT_ID),drop = T)
    
    arxdiffs <- do.call(cbind,lapply(names(sample_by_pat_list),FUN = function(x){
        pat_name <- x
        abs_samples <- sub_meta$SAMPLE_ID[sub_meta$PATIENT_ID %in% x]
        abs_subset <- data[,colnames(data) %in% abs_samples]
        arx_samps <- colnames(abs_subset)[colnames(abs_subset) %in% arx_samples]
        abs_arx <- abs_subset[,colnames(abs_subset) %in% arx_samps]
        abs_arx_bins <- assayDataElement(abs_arx,elt = "segmented")
        if(method == "mean"){
            arx_bin <- rowMeans(abs_arx_bins,na.rm = T)
        } else if(method == "median"){
            arx_bin <- rowMedians(abs_arx_bins,na.rm = T)
        } else {
            stop("unknown method")
        }
        return(arx_bin)
    }))
    
    rlpsdiffs <- do.call(cbind,lapply(names(sample_by_pat_list),FUN = function(x){
        pat_name <- x
        abs_samples <- sub_meta$SAMPLE_ID[sub_meta$PATIENT_ID %in% x]
        abs_subset <- data[,colnames(data) %in% abs_samples]
        rlps_samps <- colnames(abs_subset)[colnames(abs_subset) %in% rlps_samples]
        abs_rlps <- abs_subset[,colnames(abs_subset) %in% rlps_samps]
        abs_rlps_bins <- assayDataElement(abs_rlps,elt = "segmented")
        
        if(method == "mean"){
            rlps_bin <- rowMeans(abs_rlps_bins,na.rm = T)
        } else if(method == "median"){
            rlps_bin <- rowMedians(abs_rlps_bins,na.rm = T)
        } else {
            stop("unknown method")
        }
        return(rlps_bin)
    }))
    arxdiffs <- data.frame(arxdiffs)
    colnames(arxdiffs) <- paste0(names(sample_by_pat_list),"_arx")
    
    rlpsdiffs <- data.frame(rlpsdiffs) 
    colnames(rlpsdiffs) <- paste0(names(sample_by_pat_list),"_rlps")
    
    diffs <- cbind(arxdiffs,rlpsdiffs)
    rownames(diffs) <- rownames(data)
    return(diffs)
}

get_gene_seg <- function(target=NULL,abs_data=NULL){
    to_use <- fData(abs_data)$use
    cn_obj <- abs_data[to_use,]
    bin_pos <- fData(cn_obj)[,c("chromosome","start","end")]
    chr <- as.numeric(str_split(string = target,pattern = ":",simplify = T)[1])
    start <- as.numeric(str_split(string = target,pattern = ":|-",simplify = T)[2])
    end <- as.numeric(str_split(string = target,pattern = ":|-",simplify = T)[3])
    
    gene_chr_pos <- bin_pos[bin_pos$chromosome == chr,]
    min_start <- min(which(min(abs(gene_chr_pos$start - start)) == abs(gene_chr_pos$start - start)))
    min_end <- max(which(min(abs(gene_chr_pos$end - end)) == abs(gene_chr_pos$end - end)))
    if(gene_chr_pos$start[min_start] > start & min_start != 1){
        min_start <- min_start - 1
    }
    if(gene_chr_pos$end[min_end] < end & min_end != length(gene_chr_pos$end)){
        min_end <- min_end + 1
    }
    index_min <- which(bin_pos$chromosome == chr & bin_pos$start == gene_chr_pos[min_start,2])
    index_max <- which(bin_pos$chromosome == chr & bin_pos$end == gene_chr_pos[min_end,3])
    gene_pos <- seq.int(index_min,index_max,1)
    return(gene_pos)
}