# BriTROC archival/relapse copy number signature analysis functions

## Functions from Macintyre et al natgen CNsigs paper (adapted or as-is)
# quantify signatures using lcd function in the YAPSA package
# Required files in data/
# "data/feat_sig_mat.rds"
# "data/hg19.chrom.sizes.txt"
# "data/gap_hg19.txt"

quantifySignatures<-function(sample_by_component,component_by_signature=NULL)
{
  if(is.null(component_by_signature))
  {
    component_by_signature<-readRDS(paste("data/feat_sig_mat.rds",sep="/"))
  }
  signature_by_sample<-YAPSA::LCD(t(sample_by_component),
                                  YAPSA:::normalize_df_per_dim(component_by_signature,2))
  signature_by_sample<-normaliseMatrix(signature_by_sample)
  signature_by_sample
}

# extract copy number features from CN profiles
extractCopynumberFeatures<-function(CN_data, cores = 1)
{
  #get chromosome lengths
  chrlen<-read.table(paste("data/hg19.chrom.sizes.txt",sep="/"),sep="\t",stringsAsFactors = F)[1:24,]
  
  #get centromere locations
  gaps<-read.table(paste("data/gap_hg19.txt",sep="/"),sep="\t",header=F,stringsAsFactors = F)
  centromeres<-gaps[gaps[,8]=="centromere",]
  
  if(cores > 1) {
    require(foreach)
    doMC::registerDoMC(cores)
    
    temp_list = foreach::foreach(i=1:6) %dopar% {
      if(i == 1){
        list(segsize = getSegsize(CN_data) )
      } else if (i == 2) {
        list(bp10MB = getBPnum(CN_data,chrlen) )
      } else if (i == 3) {
        list(osCN = getOscilation(CN_data,chrlen) )
      } else if (i == 4) {
        list(bpchrarm = getCentromereDistCounts(CN_data,centromeres,chrlen) )
      } else if (i == 5) {
        list(changepoint = getChangepointCN(CN_data) )
      } else {
        list(copynumber = getCN(CN_data) )
      }
      
    }
    unlist( temp_list, recursive = FALSE )
  } else {  
    
    segsize<-getSegsize(CN_data)
    bp10MB<-getBPnum(CN_data,chrlen)
    osCN<-getOscilation(CN_data,chrlen)
    bpchrarm<-getCentromereDistCounts(CN_data,centromeres,chrlen)
    changepoint<-getChangepointCN(CN_data)
    copynumber<-getCN(CN_data)
    
    list(segsize=segsize,bp10MB=bp10MB,osCN=osCN,bpchrarm=bpchrarm,changepoint=changepoint,copynumber=copynumber)
  }
  
}

# Generate sample-componenet matrix
generateSampleByComponentMatrix<-function(CN_features, all_components=NULL, cores = 1, rowIter = 1000, subcores = 2)
{
  if(is.null(all_components))
  {
    all_components<-readRDS(paste("data/component_parameters.rds",sep="/"))
  }
  
  if(cores > 1){
    require(foreach)
    
    feats = c( "segsize", "bp10MB", "osCN", "changepoint", "copynumber", "bpchrarm" )
    doMC::registerDoMC(cores)
    
    full_mat = foreach(feat=feats, .combine=cbind) %dopar% {
      calculateSumOfPosteriors(CN_features[[feat]],all_components[[feat]], 
                               feat, rowIter = rowIter, cores = subcores)
    }
  } else {
    full_mat<-cbind(
      calculateSumOfPosteriors(CN_features[["segsize"]],all_components[["segsize"]],"segsize"),
      calculateSumOfPosteriors(CN_features[["bp10MB"]],all_components[["bp10MB"]],"bp10MB"),
      calculateSumOfPosteriors(CN_features[["osCN"]],all_components[["osCN"]],"osCN"),
      calculateSumOfPosteriors(CN_features[["changepoint"]],all_components[["changepoint"]],"changepoint"),
      calculateSumOfPosteriors(CN_features[["copynumber"]],all_components[["copynumber"]],"copynumber"),
      calculateSumOfPosteriors(CN_features[["bpchrarm"]],all_components[["bpchrarm"]],"bpchrarm"))
  }
  
  rownames(full_mat)<-unique(CN_features[["segsize"]][,1])
  full_mat[is.na(full_mat)]<-0
  full_mat
}

getSegsize<-function(abs_profiles)
{
  out<-c()
  samps<-getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }
    else
    {
      segTab<-abs_profiles[[i]]
      colnames(segTab)[4]<-"segVal"
    }
    segTab$segVal[as.numeric(segTab$segVal)<0]<-0
    seglen<-(as.numeric(segTab$end)-as.numeric(segTab$start))
    seglen<-seglen[seglen>0]
    out<-rbind(out,cbind(ID=rep(i,length(seglen)),value=seglen))
  }
  rownames(out)<-NULL
  data.frame(out,stringsAsFactors = F)
}

getBPnum<-function(abs_profiles,chrlen)
{
  out<-c()
  samps<-getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }else
    {
      segTab<-abs_profiles[[i]]
      colnames(segTab)[4]<-"segVal"
    }
    chrs<-unique(segTab$chromosome)
    allBPnum<-c()
    for(c in chrs)
    {
      currseg<-segTab[segTab$chromosome==c,]
      intervals<-seq(1,chrlen[chrlen[,1]==paste0("chr",c),2]+10000000,10000000)
      res <- hist(as.numeric(currseg$end[-nrow(currseg)]),breaks=intervals,plot=FALSE)$counts
      allBPnum<-c(allBPnum,res)
    }
    out<-rbind(out,cbind(ID=rep(i,length(allBPnum)),value=allBPnum))
  }
  rownames(out)<-NULL
  data.frame(out,stringsAsFactors = F)
}

getOscilation<-function(abs_profiles,chrlen)
{
  out<-c()
  samps<-getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }else
    {
      segTab<-abs_profiles[[i]]
      colnames(segTab)[4]<-"segVal"
    }
    chrs<-unique(segTab$chromosome)
    oscCounts<-c()
    for(c in chrs)
    {
      currseg<-segTab[segTab$chromosome==c,"segVal"]
      currseg<-round(as.numeric(currseg))
      if(length(currseg)>3)
      {
        prevval<-currseg[1]
        count=0
        for(j in 3:length(currseg))
        {
          if(currseg[j]==prevval&currseg[j]!=currseg[j-1])
          {
            count<-count+1
          }else{
            oscCounts<-c(oscCounts,count)
            count=0
          }
          prevval<-currseg[j-1]
        }
      }
    }
    out<-rbind(out,cbind(ID=rep(i,length(oscCounts)),value=oscCounts))
    if(length(oscCounts)==0)
    {
      out<-rbind(out,cbind(ID=i,value=0))
    }
  }
  rownames(out)<-NULL
  data.frame(out,stringsAsFactors = F)
}

getCentromereDistCounts<-function(abs_profiles,centromeres,chrlen)
{
  out<-c()
  samps<-getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }else
    {
      segTab<-abs_profiles[[i]]
      colnames(segTab)[4]<-"segVal"
    }
    chrs<-unique(segTab$chromosome)
    all_dists<-c()
    for(c in chrs)
    {
      if(nrow(segTab)>1)
      {
        starts<-as.numeric(segTab[segTab$chromosome==c,2])[-1]
        segstart<-as.numeric(segTab[segTab$chromosome==c,2])[1]
        ends<-as.numeric(segTab[segTab$chromosome==c,3])
        segend<-ends[length(ends)]
        ends<-ends[-length(ends)]
        centstart<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,3])
        centend<-as.numeric(centromeres[substr(centromeres[,2],4,5)==c,4])
        chrend<-chrlen[substr(chrlen[,1],4,5)==c,2]
        ndist<-cbind(rep(NA,length(starts)),rep(NA,length(starts)))
        ndist[starts<=centstart,1]<-(centstart-starts[starts<=centstart])/(centstart-segstart)*-1
        ndist[starts>=centend,1]<-(starts[starts>=centend]-centend)/(segend-centend)
        ndist[ends<=centstart,2]<-(centstart-ends[ends<=centstart])/(centstart-segstart)*-1
        ndist[ends>=centend,2]<-(ends[ends>=centend]-centend)/(segend-centend)
        ndist<-apply(ndist,1,min)
        
        all_dists<-rbind(all_dists,sum(ndist>0))
        all_dists<-rbind(all_dists,sum(ndist<=0))
      }
    }
    if(nrow(all_dists)>0)
    {
      out<-rbind(out,cbind(ID=i,ct1=all_dists[,1]))
    }
  }
  rownames(out)<-NULL
  data.frame(out,stringsAsFactors = F)
}


getChangepointCN<-function(abs_profiles)
{
  out<-c()
  samps<-getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }
    else
    {
      segTab<-abs_profiles[[i]]
      colnames(segTab)[4]<-"segVal"
    }
    segTab$segVal[as.numeric(segTab$segVal)<0]<-0
    chrs<-unique(segTab$chromosome)
    allcp<-c()
    for(c in chrs)
    {
      currseg<-as.numeric(segTab[segTab$chromosome==c,"segVal"])
      allcp<-c(allcp,abs(currseg[-1]-currseg[-length(currseg)]))
    }
    if(length(allcp)==0)
    {
      allcp<-0 #if there are no changepoints
    }
    out<-rbind(out,cbind(ID=rep(i,length(allcp)),value=allcp))
  }
  rownames(out)<-NULL
  data.frame(out,stringsAsFactors = F)
}

getCN<-function(abs_profiles)
{
  out<-c()
  samps<-getSampNames(abs_profiles)
  for(i in samps)
  {
    if(class(abs_profiles)=="QDNAseqCopyNumbers")
    {
      segTab<-getSegTable(abs_profiles[,which(colnames(abs_profiles)==i)])
    }
    else
    {
      segTab<-abs_profiles[[i]]
      colnames(segTab)[4]<-"segVal"
    }
    segTab$segVal[as.numeric(segTab$segVal)<0]<-0
    cn<-as.numeric(segTab$segVal)
    out<-rbind(out,cbind(ID=rep(i,length(cn)),value=cn))
  }
  rownames(out)<-NULL
  data.frame(out,stringsAsFactors = F)
}

getSampNames<-function(abs_profiles)
{
  if(class(abs_profiles)=="QDNAseqCopyNumbers")
  {
    samps<-colnames(abs_profiles)
  }
  else
  {
    samps<-names(abs_profiles)
  }
  samps
}

getSegTable<-function(x)
{
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

calculateSumOfPosteriors<-function(CN_feature,components,name, rowIter = 1000, cores = 1)
{
  
  if(cores > 1){
    require(foreach)
    require(doMC)
    
    len = dim(CN_feature)[1]
    iters = floor( len / rowIter )
    lastiter = iters[length(iters)]
    
    registerDoMC(cores)
    curr_posterior = foreach( i=0:iters, .combine=rbind) %dopar% {
      start = i*rowIter+1
      if(i != lastiter) { end = (i+1)*rowIter } else { end = len }
      flexmix::posterior(components,data.frame(dat=as.numeric(CN_feature[start:end,2])))
    }
  } else {
    curr_posterior<-flexmix::posterior(components,data.frame(dat=as.numeric(CN_feature[,2])))
  }
  
  mat<-cbind(CN_feature,curr_posterior)
  posterior_sum<-c()
  
  ## foreach and parallelising doesn't make the following code faster.
  for(i in unique(mat$ID))
  {
    posterior_sum<-rbind(posterior_sum,colSums(mat[mat$ID==i,c(-1,-2)]))
  }
  params<-flexmix::parameters(components)
  if(!is.null(nrow(params)))
  {
    posterior_sum<-posterior_sum[,order(params[1,])]
  }
  else
  {
    posterior_sum<-posterior_sum[,order(params)]
  }
  colnames(posterior_sum)<-paste0(name,1:ncol(posterior_sum))
  rownames(posterior_sum)<-rownames(unique(mat$ID))
  posterior_sum
}

normaliseMatrix<-function(signature_by_sample,sig_thresh=0.01)
{
  norm_const<-colSums(signature_by_sample)
  sample_by_signature<-apply(signature_by_sample,1,function(x){x/norm_const})
  sample_by_signature<-apply(sample_by_signature,1,lower_norm,sig_thresh)
  signature_by_sample<-t(sample_by_signature)
  norm_const<-apply(signature_by_sample,1,sum)
  sample_by_signature<-apply(signature_by_sample,2,function(x){x/norm_const})
  signature_by_sample<-t(sample_by_signature)
  signature_by_sample
}

lower_norm<-function(x,sig_thresh=0.01)
{
  new_x<-x
  for(i in 1:length(x))
  {
    if(x[i]<sig_thresh)
    {
      new_x[i]<-0
    }
  }
  new_x
}

# Calculate CN stddev across CN states 1-4
# Function written by Dr Dilrini De Silva
cn_stddev <- function(cn_object = cn_object, tolerance = 0.3){
  y <- cn_object
  sdlist <- list()
  # iterate over each sample in absCN object
  for(id in sampleNames(y)){
    x <- y[,id]
    # Grab used reads
    nreads <- pData(x)$used.reads
    # Extract CN and Segment data and form dataframe
    copynumber <- as.numeric(assayDataElement(x, "copynumber"))
    segmented <- as.numeric(assayDataElement(x, "segmented"))
    df <- data.frame(copynumber, segmented, stringsAsFactors=FALSE)
    # Subselect segments which are within tolerance range of the CN states 1-4
    sdvals <- fData(x) %>%
      bind_cols(df) %>%
      filter(use) %>%
      filter(near(segmented, 1, tol=tolerance) |
               near(segmented, 2, tol=tolerance) |
               near(segmented, 3, tol=tolerance) |
               near(segmented, 4, tol=tolerance)) %>% #filter out segments that are not at integer states
      # Assign CN state to each segment over CN states 1-4
      mutate(cn_state=ifelse(near(segmented, 1, tol=tolerance), 1,NA)) %>%
      mutate_at(.vars ="cn_state", funs(ifelse(near(segmented, 2, tol=tolerance),2,.))) %>% 
      mutate_at(.vars ="cn_state", funs(ifelse(near(segmented, 3, tol=tolerance),3,.))) %>%
      mutate_at(.vars ="cn_state", funs(ifelse(near(segmented, 4, tol=tolerance),4,.))) %>%
      group_by(cn_state) %>%
      # Summarise sample information and calculate CN std dev. for segements at each CN state
      summarize(
        name = id,
        rd = nreads,
        #cn_state = cn_state, # included as it is grouping variable
        sd = sd(copynumber)
      )
    # Add object to list
    sdlist[[id]] <- sdvals
    
    
  }
  # bind rows to single dataframe and return
  sddf <- bind_rows(sdlist)
  return(sddf)
}

# Function to plot the feature distribution for each derived feature
plot_feat_dist <- function(featlist = featlist){
  p1 <- ggplot() + geom_density(data = as.data.frame(featlist[[1]]),aes(as.numeric(unlist(featlist[[1]][2])))) + 
    labs(title = names(featlist[1])) +
    theme_bw() + theme(axis.title.x = element_blank())
  p2 <- ggplot() + geom_histogram(data = as.data.frame(featlist[[2]]),aes(as.numeric(unlist(featlist[[2]][2]))),binwidth = 1) + 
    labs(title = names(featlist[2])) + 
    theme_bw() + theme(axis.title.x = element_blank())
  p3 <- ggplot() + geom_histogram(data = as.data.frame(featlist[[3]]),aes(as.numeric(unlist(featlist[[3]][2]))),binwidth = 1) + 
    labs(title = names(featlist[3])) + 
    theme_bw() + theme(axis.title.x = element_blank())
  p4 <- ggplot() + geom_histogram(data = as.data.frame(featlist[[4]]),aes(as.numeric(unlist(featlist[[4]][2]))),binwidth = 1) + 
    labs(title = names(featlist[4])) + 
    theme_bw() + theme(axis.title.x = element_blank())
  p5 <- ggplot() + geom_density(data = as.data.frame(featlist[[5]]),aes(as.numeric(unlist(featlist[[5]][2])))) + 
    labs(title = names(featlist[5])) +
    theme_bw() + theme(axis.title.x = element_blank())
  p6 <- ggplot() + geom_density(data = as.data.frame(featlist[[6]]),aes(as.numeric(unlist(featlist[[6]][2])))) + 
    labs(title = names(featlist[6])) +
    theme_bw() + theme(axis.title.x = element_blank())
  p_list <- list(p1,p2,p3,p4,p5,p6)
  # Plot all as grid of plots
  plot_grid(plotlist = p_list,nrow = 3,ncol = 2)
}

# Plot abs profile
plot_abs <- function(cn_data = x, sample_name = NULL){
  if(is.null(sample_name)){
    stop("No sample specified")
  } else if(!sample_name %in% colnames(cn_data)){
    stop("Sample not present in CN data")
  }
  sub_cn <- cn_data[,colnames(cn_data) == sample_name]
  TP53freq <- pData(sub_cn)$TP53freq
  ploidy <- pData(sub_cn)$ploidy
  purity <- pData(sub_cn)$purity
  if(ploidy>5){
    yrange=15
  } else {
    yrange=10
  }
  # Plot abs fit
  par(mfrow = c(1,1))
  plot(sub_cn,doCalls=FALSE,showSD=TRUE,logTransform=FALSE,ylim=c(0,yrange),ylab="Absolute tumour CN",
       main=paste(sample_name," AF=", round(TP53freq,2), " p=",round(purity,2)," pl=",round(ploidy,2), sep=""))
  abline(h=1:9, col = "blue")
}

# Plot rel exposure plot
plot_rel_expo <- function(sig_quants = NULL){
  if(is.null(sig_quants)){
    stop("No exposure data provided - data should be generated by the function 'quantifySignatures()'")
  }
  melted_CN_data <- reshape2::melt(sig_quants)
    melted_CN_data$Var2 <- factor(melted_CN_data$Var2, levels= patient_sorted_samples)
    ggplot2::ggplot(melted_CN_data,aes(x=Var2,y=value,fill=Var1))+geom_bar(stat="identity")+theme_bw()+
      scale_fill_manual(values = cbPalette,name="Signature") +
      theme_bw() +
      xlab("") +
      ylab("exposure") +
      theme(#axis.text=ggplot2::element_text(size=5),
            axis.text.x = element_text(angle = 90),
            axis.title=ggplot2::element_text(size=5),
            strip.text.x = ggplot2::element_text(size = 7),
            strip.text.y = ggplot2::element_text(size = 7),
            legend.text = ggplot2::element_text(size = 7),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank())
}

# Plot individual rel exposure plots
plot_idv_expo <- function(sig_quants = NULL){
  if(is.null(sig_quants)){
    stop("No exposure data provided - data should be generated by the function 'quantifySignatures()'")
  }
  melted_CN_data <- reshape2::melt(sig_quants)
  melted_CN_data$Var2 <- factor(melted_CN_data$Var2, levels= patient_sorted_samples)
  ggplot2::ggplot(melt(melted_CN_data)) +
    geom_col(aes(Var1,value,fill=Var1)) +
    theme_bw() +
    scale_fill_manual(values = cbPalette,name="Signature") +
    facet_grid(~Var2)
}

# Get patient id for sample names
get_patient <- function(samples = NULL){
  if(is.null(samples)){
    stop("No samples provided - character vector of sample names")
  }
  pat_vec <- c()
  for(i in samples){
    pat_vec <- append(pat_vec,unique(as.character(meta.data$PATIENT_ID[meta.data$SAMPLE_ID %in% i])))
  }
  return(pat_vec)
}

# Misc colour palette - 7 sigs
cbPalette <- c(RColorBrewer::brewer.pal(8,"Dark2"),RColorBrewer::brewer.pal(9,"Set1"),"black")

## additional 

lighten <- function(color, factor = 0.5) {
    if ((factor > 1) | (factor < 0)) stop("factor needs to be within [0,1]")
    col <- col2rgb(color)
    col <- col + (255 - col)*factor
    col <- rgb(t(col), maxColorValue=255)
    col
}

# Generate signature box plot
sig_box_plot <- function(sig_quants = NULL,test="t.test",plot_title = NULL,paired=FALSE){
    if(is.null(sig_quants)){
        stop("No signature exposure quants provided")
    }
    x <- sig_quants
    if(paired){
        sig_box_data <- base::merge(x = rownames_to_column(
            as.data.frame(t(x)),var = "Sample"),
            y = patient.meta[patient.meta$SAMPLE_ID %in% colnames(x),c("PATIENT_ID","SAMPLE_ID","group")],
            by.x = "Sample",
            by.y = "SAMPLE_ID")
        sig_box_data_long <- pivot_longer(data = sig_box_data,cols = c(2:8),names_to = "Signature") %>%
            mutate(group = case_when(group == "arx" ~ "diagnosis",group == "rlps" ~ "relapse")) %>%
            group_by(PATIENT_ID,Signature,group) %>%
            summarise(across(.cols = value,.fns = median))
    } else {
        sig_box_data <- base::merge(x = rownames_to_column(
            as.data.frame(t(x)),var = "Sample"),
            y = patient.meta[patient.meta$SAMPLE_ID %in% colnames(x),c("SAMPLE_ID","group")],
            by.x = "Sample",
            by.y = "SAMPLE_ID")
        sig_box_data_long <- pivot_longer(data = sig_box_data,cols = c(2:8),names_to = "Signature") %>%
            mutate(group = case_when(group == "arx" ~ "diagnosis",group == "rlps" ~ "relapse"))
    }
    if(paired){
        if(test == "t.test"){
            test.p <- c()
            for(i in unique(sig_box_data_long$Signature)){
                p <- t.test(x = sig_box_data_long$value[sig_box_data_long$group == "diagnosis" & sig_box_data_long$Signature == i],
                            y = sig_box_data_long$value[sig_box_data_long$group == "relapse" & sig_box_data_long$Signature == i],paired = TRUE)$p.value
                test.p <- append(test.p,p)
            }
        } else if(test == "wilcox.test"){
            test.p <- c()
            for(i in unique(sig_box_data_long$Signature)){
                p <- wilcox.test(x = sig_box_data_long$value[sig_box_data_long$group == "diagnosis" & sig_box_data_long$Signature == i],
                                 y = sig_box_data_long$value[sig_box_data_long$group == "relapse" & sig_box_data_long$Signature == i],paired = TRUE)$p.value
                test.p <- append(test.p,p)
            }
        } else {
            stop("unknown statistical test provided - t.test or wilcox.test")
        }
    } else {
        if(test == "t.test"){
            test.p <- c()
            for(i in unique(sig_box_data_long$Signature)){
                p <- t.test(x = sig_box_data_long$value[sig_box_data_long$group == "diagnosis" & sig_box_data_long$Signature == i],
                            y = sig_box_data_long$value[sig_box_data_long$group == "relapse" & sig_box_data_long$Signature == i])$p.value
                test.p <- append(test.p,p)
            }
        } else if(test == "wilcox.test"){
            test.p <- c()
            for(i in unique(sig_box_data_long$Signature)){
                p <- wilcox.test(x = sig_box_data_long$value[sig_box_data_long$group == "diagnosis" & sig_box_data_long$Signature == i],
                                 y = sig_box_data_long$value[sig_box_data_long$group == "relapse" & sig_box_data_long$Signature == i])$p.value
                test.p <- append(test.p,p)
            }
        } else {
            stop("unknown statistical test provided - t.test or wilcox.test")
        }
    }
    if(is.null(plot_title)){
        plot_title <- "Signature exposures by signature between groups"
    }
    ggplot(data = sig_box_data_long) +
        geom_jitter(aes(x = Signature,y = value,color=group),
                    alpha=0.3,position = position_jitterdodge(jitter.width = 0.2)) +
        geom_boxplot(aes(x = Signature,y = value,fill=group),outlier.colour = NA) +
        geom_signif(aes(x = Signature,y = value),
                    annotations = ifelse(test.p < 0.01,formatC(test.p,format = "e",digits = 2),signif(test.p,digits = 2)),
                                         xmin = seq.int(0.5,6.5,1)+0.2,
                                         xmax = seq.int(1.5,7.5,1)-0.2,
                                         y_position = 0.9,vjust = 1.4) +
        labs(title = plot_title) +
        scale_y_continuous(name = "Exposure",limits = c(0,1)) +
        scale_fill_manual(values = colour_palettes$diagnosis_relapse) +
        scale_color_manual(values = colour_palettes$diagnosis_relapse) +
        theme_bw()
}

get_seg_counts <- function(abs_data = NULL){
    if(is.null(abs_data)){
        stop("No data provided")
    }
    sample_n <- length(colnames(abs_data))
    to_use <- fData(abs_data)$use
    sample_segs <- c()
    for(i in 1:sample_n){
        cn_obj <- abs_data[to_use,i]
        sample_name <- colnames(cn_obj)
        segments <- assayDataElement(cn_obj,"segmented")
        run.enc <- rle(as.numeric(segments))
        seg_count <- as.numeric(length(run.enc$lengths))
        sample_segs <- rbind(sample_segs,c(sample_name,seg_count))
    }
    sample_segs <- as.data.frame(sample_segs)
    colnames(sample_segs) <- c("sample","segments")
    sample_segs$segments <- as.numeric(sample_segs$segments)
    return(sample_segs)
}