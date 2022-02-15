#' @title Standardize OTU abundance table
#' @description Standardize phyloseq OTU table with with methods from \code{\link[vegan]{decostand}} from vegan package.
#' @param physeq A phyloseq-class object
#' @param method Standardization method
#' @param ... Additional parameters may be passed to vegan \code{\link[vegan]{decostand}} function
#' @return phyloseq object with standardized OTU table.
#' @seealso \code{\link[vegan]{decostand}}, \code{\link{phyloseq_transform_css}}, \code{\link{phyloseq_transform_vst_blind}}, \code{\link{phyloseq_transform_rlog_blind}}, \code{\link{physeq_transform_anderson_log}}
#' 
#' @details
#' Supported methods:
#' \itemize{
#' \item \strong{"total"} - convert data to relative abundances (divide by sample total);
#' \item \strong{"pa"} - convert OTU abundances to presence/absence scale (0/1);
#' \item \strong{"log"} - logarithmic transformation as suggested by Anderson et al. (2006): log (x) + 1 for x > 0, zeros are left as zeros, logarithm base = 2. Please note this is not log(x+1);
#' \item \strong{"hellinger"} - square root of method = "total" (Legendre & Gallagher 2001);
#' \item \strong{"max"} - divide by sample maximum;
#' \item \strong{"frequency"} - divide by sample total and multiply by the number of non-zero items, so that the average of non-zero entries is one;
#' \item \strong{"normalize"} - make sample sum of squares equal to one;
#' \item \strong{"range"} - standardize values into range 0 ... 1. If all values are constant, they will be transformed to 0;
#' \item \strong{"rank"} - replace abundance values by their increasing ranks leaving zeros unchanged. Average ranks are used for tied values;
#' \item \strong{"rrank"} - replace abundance values by relative ranks with maximum 1. Average ranks are used for tied values;
#' \item \strong{"standardize"} - scale OTU abundances within sample to zero mean and unit variance;
#' \item \strong{"wisconsin"} - Wisconsin double standardization where species are first standardized by maxima and then sites by site totals;
#' \item \strong{"chi.square"} - divide by sample sums and square root of OTU sums, and adjust for square root of matrix total (Legendre & Gallagher 2001). When used with the Euclidean distance, the distances should be similar to the Chi-square distance used in correspondence analysis.
#' }
#' 
#' For the implementation of "total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log" methods see \code{\link[vegan]{decostand}}.
#' 
#' @export
#'
#' @examples
#' # Load data
#' data("esophagus")
#'
#' # Total-sum scaling (TSS) normalization
#' phyloseq_standardize_otu_abundance(esophagus, method = "total")
#' # the same as
#' transform_sample_counts(esophagus, function(OTU) OTU/sum(OTU) )
#' identical(
#'   phyloseq_standardize_otu_abundance(esophagus, method = "total"),
#'   transform_sample_counts(esophagus, function(OTU) OTU/sum(OTU)) )
#'
#' # Presence-absence scaling (0/1)
#' phyloseq_standardize_otu_abundance(esophagus, method = "pa")
#'
#' # Logarithmic transformation as in Anderson et al., 2006
#' phyloseq_standardize_otu_abundance(esophagus, method = "log")
#'
#' # Hellinger standardization
#' phyloseq_standardize_otu_abundance(esophagus, method = "hellinger")
#'
phyloseq_standardize_otu_abundance <- function(physeq, method = "total",
                                               rhea_depth = NULL, rhea_round = TRUE, rhea_threshold = NULL,
                                               ...){
    
    ## Method implemente in vegan
    vegan_methods <- c("total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log")
    other_methods <- c("rhea", "wisconsin")
    
    ## Check the orientation of the OTU table
    trows <- phyloseq::taxa_are_rows(physeq)
    if(trows == TRUE){ marg <- 2 } else { marg <- 1 }
    
    ## Extact OTU table
    comm <- as(object = phyloseq::otu_table(physeq), Class = "matrix")
    
    ## Standardize with vegan methods
    if(method %in% vegan_methods){
        
        ## Standardize community table with vegan
        comm_std <- vegan::decostand(comm, method, MARGIN = marg, ...)
        
        ## Transpose table back
        if(method == "chi.square"){ comm_std <- t(comm_std) }
        
    }
    
    ## Normalize counts via division to the sample size and then multiplication by the size of the smaller sample (or other seq depth)
    ## as implemented in "Rhea" package (Lagkouvardos et al., 2017, DOI:10.7717/peerj.2836; Reitmeier et al. 2020, DOI:10.21203/rs.2.21240/v1)
    if(method == "rhea"){
        
        ## If there is no user-provided fixed sample size - use minimum sampling size
        if(is.null(rhea_depth)){
            rhea_depth <- min( phyloseq::sample_sums(physeq) )
        }
        
        ## Convert data to relative abundances
        comm <- vegan::decostand(comm, method = "total", MARGIN = marg)
        
        ## Remove OTUs with low relative abundance
        if(!is.null(rhea_threshold)){
            
            comm[ comm < rhea_threshold ] <- 0
            
            ## Re-normalize by total abundance again (sample sums should be 1)
            comm <- vegan::decostand(comm, method = "total", MARGIN = marg)
        }
        
        ## Multiply by the desired samples size
        comm_std <- comm * rhea_depth
        
        ## Round abundance to the nearest integer below its current value
        # if(rhea_round == TRUE){ comm_std <- floor(comm_std) }
        
        ## Round abundance to the nearest integer while preserving the rounded sum
        if(rhea_round == TRUE){ 
            
            ## Function by josliber
            # https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
            smart_round <- function(x) {
                y <- floor(x)
                indices <- tail(order(x-y), round(sum(x)) - sum(y))
                y[indices] <- y[indices] + 1
                return(y)
            }
            
            comm_std <- as.data.frame(comm_std)
            comm_std <- apply(comm_std, MARGIN = marg, FUN = smart_round)
        }
    } # end of "rhea"
    
    ## Wisconsin Double standardization
    ## Species (MARGIN=2) are first standardized by maxima and then sites by site totals
    if(method == "wisconsin"){
        
        if(marg == 2){ mm <- c(1, 2) }
        if(marg == 1){ mm <- c(2, 1) }
        
        comm_std <- vegan::decostand(comm, "max", MARGIN = mm[1])
        comm_std <- vegan::decostand(comm_std, "tot", MARGIN = mm[2])
        
    }
    
    ## Replace old otu_table with the new one
    phyloseq::otu_table(physeq) <- phyloseq::otu_table(comm_std, taxa_are_rows = trows)
    
    return(physeq)
}

#' @title Remove samples from phyloseq object that have less than n taxa
#'
#' @param physeq A phyloseq-class object
#' @param mintaxa Minimum number of taxa that should be present in a sample (default, 10)
#'
#' @return Trimmed phyloseq object (All samples will have >= N taxa)
phyloseq_richness_filter <- function(physeq, mintaxa = 10){
    
    ## Estimate number of OTUs per sample
    sp <- phyloseq::estimate_richness(physeq, measures = "Observed")
    samples_to_keep <- rownames(sp)[ which(sp$Observed >= mintaxa) ]
    
    if(length(samples_to_keep) == 0){
        stop("All samples will be removed.\n")
    }
    
    if(length(samples_to_keep) == phyloseq::nsamples(physeq)){
        cat("All samples will be preserved\n")
        res <- physeq
    }
    
    if(length(samples_to_keep) < phyloseq::nsamples(physeq)){
        res <- phyloseq::prune_samples(samples = samples_to_keep, x = physeq)
    } else {
        res <- NULL
    }
    
    return(res)
}



#' @title Remove taxa with small mean relative abundance.
#'
#' @param physeq A phyloseq-class object
#' @param frac The minimum cutoff for the relative OTU abundance
#' @details This function searches for taxa with small mean relative abundance and removes them. Result will be returned with original counts in the abundance table.
#' @return Phyloseq object with a subset of taxa.
phyloseq_filter_taxa_rel_abund <- function(physeq, frac = 1e-4){

    ## Transform OTU counts to relative abundance
    rel <- phyloseq::transform_sample_counts(physeq, function(x) x / sum(x) )
    
    ## Filter OTUs
    rel.subs <- phyloseq::filter_taxa(rel, function(x){ mean(x) > frac }, prune = FALSE)
    
    ## Taxa to remove
    tr <- names(rel.subs)[ which(rel.subs == FALSE) ]
    
    if(length(tr) == phyloseq::ntaxa(physeq)){
        ## If all taxa should be removed
        res <- NULL
    } else if(length(tr) == 0){
        ## If there is nothing to remove
        res <- physeq
    } else {
        ## Keep taxa which satisfies the truncation threshold
        res <- phyloseq::prune_taxa(taxa = rel.subs, physeq)
    }
     return(res)
}


#' @title Remove taxa with abundance less then a certain fraction of total abundance.
#'
#' @param physeq A phyloseq-class object
#' @param frac The minimum cutoff for the OTU abundance in the table. This number is a fraction, not a percent.
#' @details
#' If frac = 0.0001, this will retain all OTU's that have at least a 0.01\% total abundance in the OTU table.
#' If you wanted to retain OTUs with at least 1\% total abundance, you must specify, 0.01.
#'
#' @return Phyloseq object with a subset of taxa.
phyloseq_filter_taxa_tot_fraction <- function(physeq, frac = 0.01){

    ## Estimate total abundance of OTUs
    tot <- sum(phyloseq::taxa_sums(physeq))
    
    ## Remove OTUs
    res <- phyloseq::filter_taxa(physeq, function(x){ ( sum(x)/tot ) > frac }, prune = TRUE)
    return(res)
}


#' @title Filter low-prevalence OTUs.
#' @description This function will remove taxa (OTUs) with low prevalence, where prevalence is the fraction of total samples in which an OTU is observed.
#' @param physeq A phyloseq-class object
#' @param prev.trh Prevalence threshold (default, 0.05 = 5\% of samples)
#' @param abund.trh Abundance threshold (default, NULL)
#' @param threshold_condition Indicates type of prevalence and abundance conditions, can be "OR" (default) or "AND"
#' @param abund.type Character string indicating which type of OTU abundance to take into account for filtering ("total", "mean", or "median")
#' @details
#' Abundance threshold defines if the OTU should be preserved if its abundance is larger than threshold (e.g., >= 50 reads).
#' Parameter "threshold_condition" indicates whether OTU should be kept if it occurs in many samples AND/OR it has high abundance.
#' @return  Phyloseq object with a subset of taxa.
phyloseq_filter_prevalence <- function(physeq, prev.trh = 0.05, abund.trh = NULL, threshold_condition = "OR", abund.type = "total"){
    
    ## Threshold validation
    if(prev.trh > 1 | prev.trh < 0){ stop("Prevalence threshold should be non-negative value in the range of [0, 1].\n") }
    if(!is.null(abund.trh)){ 
        if(abund.trh <= 0){ stop("Abundance threshold should be non-negative value larger 0.\n") }
    }
    
    ## Check for the low-prevalence species (compute the total and average prevalences of the features in each phylum)
    prevdf_smr <- function(prevdf){
        ddply(prevdf, "Phylum", function(df1){ data.frame(Average = mean(df1$Prevalence), Total = sum(df1$Prevalence))})
    }

    ## Define prevalence threshold as % of total samples
    ## This function is located in 'phyloseq_prevalence_plot.R' file
    prevalenceThreshold <- prev.trh * phyloseq::nsamples(physeq)
    
    ## Calculate prevalence (number of samples with OTU) and OTU total abundance
    prevdf <- prevalence(physeq)
    
    ## Get the abundance type
    if(abund.type == "total") { prevdf$AbundFilt <- prevdf$TotalAbundance }
    if(abund.type == "mean")  { prevdf$AbundFilt <- prevdf$MeanAbundance }
    if(abund.type == "median"){ prevdf$AbundFilt <- prevdf$MedianAbundance }
    
    ## Which taxa to preserve
    if(is.null(abund.trh)) { tt <- prevdf$Prevalence >= prevalenceThreshold }
    if(!is.null(abund.trh)){
        ## Keep OTU if it either occurs in many samples OR it has high abundance
        if(threshold_condition == "OR"){
            tt <- (prevdf$Prevalence >= prevalenceThreshold | prevdf$AbundFilt >= abund.trh)
        }
        
        ## Keep OTU if it occurs in many samples AND it has high abundance
        if(threshold_condition == "AND"){
            tt <- (prevdf$Prevalence >= prevalenceThreshold & prevdf$AbundFilt >= abund.trh)
        }
    }
    
    ## Extract names for the taxa we want to keep
    keepTaxa <- rownames(prevdf)[tt]
    
    ## Execute prevalence filter
    if (length(keepTaxa) > 0){
        res <- phyloseq::prune_taxa(keepTaxa, physeq)
    } else {
        res <- NULL
    }
    return(res)
}

#' @title Filter rare OTUs based on minimum abundance threshold.
#' @description This function performs sample-wise OTU abundance trimming.
#' @param physeq A phyloseq-class object
#' @param minabund Abundance threshold (default, 10)
#' @param relabund Logical; perform trimming based on relative abundances (default, FALSE)
#' @param rm_zero_OTUs Logical, remove OTUs with zero total abundance
#' @details 
#' OTUs can be considered as rare if they comprise fewer than X (e.g., 10) sequences within a sample. 
#' This function is intented to censore OTU abundance (unsing an arbitrary threshold) on a sample-wise basis. 
#' 
#' Trimming can be performed based on relative abundances of OTUs within a sample (`relabund = TRUE`), but the orginal OTU count will be returned. 
#' For this purpose `minabund` parameter should be provided in a range of (0,1] (e.g., use `minabund = 0.1, relabund = TRUE` to remove OTUs with relative abundance < 10% in each sample).
#' 
#' @return Phyloseq object with a filtered data.
phyloseq_filter_sample_wise_abund_trim <- function(physeq, minabund = 10, relabund = FALSE, rm_zero_OTUs = TRUE){
    
    ## Censore OTU abundance
    if(relabund == FALSE){     # trim based on absolute OTU counts
        
        res <- phyloseq::transform_sample_counts(physeq, function(OTU, ab = minabund){ ifelse(OTU <= ab,  0, OTU) })
        
    } else {                   # trim based on relative abundances within sample, but return original counts
        
        if(!minabund > 0 & minabund <= 1){
            stop("Error: for relative abundance trimmin 'minabund' should be in (0,1] interval.\n")
        }
        
        ## Convert data to relative abundances
        res <- phyloseq_standardize_otu_abundance(physeq, method = "total")
        
        ## Remove relative abundances less than the threshold value
        res <- phyloseq::transform_sample_counts(res, function(OTU, ab = minabund){ ifelse(OTU <= ab,  0, OTU) })
        
        ## Sample sums and data orientation
        smps <- phyloseq::sample_sums(physeq)
        if(phyloseq::taxa_are_rows(physeq) == TRUE){
            mar <- 2
        } else {
            mar <- 1
        }
        
        ## Convert back to counts by multiplying relative abundances by sample sums
        phyloseq::otu_table(res) <- phyloseq::otu_table(
            sweep(x = phyloseq::otu_table(res), MARGIN = mar, STATS = smps, FUN = `*`),
            taxa_are_rows = phyloseq::taxa_are_rows(physeq))
    }
    
    ## Remove zero-OTUs
    if(rm_zero_OTUs == TRUE){
        if (any(taxa_sums(res) > 0)){
            res <- phyloseq::prune_taxa(taxa_sums(res) > 0, res)
        } else {
            res <- NULL
        }
    }
    return(res)
}



#' @title Extract the most abundant taxa.
#' @param physeq A phyloseq-class object
#' @param perc Percentage of the most abundant taxa to retain
#' @param n Number of the most abundant taxa to retain (this argument will override perc argument)
#' @return Phyloseq object with a filtered data.
phyloseq_filter_top_taxa <- function(physeq, perc = 10, n = NULL){
    
    ## Arguments validation
    if(perc <= 0 | perc > 100){ stop("Error: percentage should be in 1-100 range.\n") }
    
    ## Get total abundances for all taxa
    taxx <- sort(phyloseq::taxa_sums(physeq), decreasing = TRUE)
    
    ## Find how many taxa to preserve (if percentage is specified)
    if(is.null(n)){
        n <- phyloseq::ntaxa(physeq) * perc / 100
        n <- floor(n)
    }
    
    ## Extract names for the taxa that should be preserved
    keepTaxa <- names(taxx)[1:n]
    
    ## Extract this taxa
    physeq_pruned <- phyloseq::prune_taxa(keepTaxa, physeq)
    
    return(physeq_pruned)
}

#' @title Filter by sample richness
#' @param physeq A phyloseq-class object
#' @param mintaxa Minimum number of taxa to retain a sample
phyloseq_richness_filter <- function(physeq, mintaxa = 10){
    
    ## Estimate number of OTUs per sample
    sp <- suppressWarnings(phyloseq::estimate_richness(physeq, measures = "Observed"))
    samples_to_keep <- rownames(sp)[ which(sp$Observed >= mintaxa) ]
    
    if(length(samples_to_keep) == 0){
        stop("All samples will be removed.\n")
    }
    
    if(length(samples_to_keep) == phyloseq::nsamples(physeq)){
        cat("All samples will be preserved\n")
        res <- physeq
    }
    
    if(length(samples_to_keep) < phyloseq::nsamples(physeq)){
        res <- phyloseq::prune_samples(samples = samples_to_keep, x = physeq)
    }
    
    return(res)
}

# Input filtering function
filter_all <- function(input, physeq){
    
    # Get initial phyloseq object
    ps0 <- physeq
    
    # Create list for returning summaries
    out <- list(
        starting_taxa = ntaxa(ps0),
        filter_sample_wise_abund_counts = NULL,
        filter_sample_wise_abund_ra = NULL,
        filter_taxa_mean_threshold = NULL,
        filter_taxa_total_threshold = NULL,
        koverA = NULL,
        filter_top_taxa = NULL,
        filter_sample_sums_threshold = NULL,
        filter_richness = NULL,
        ps = NULL
    )
    
    # Per sample counts
    if( input$filter_taxa_sums_threshold > 0 & is(ps0, "phyloseq")){
        ps0 <- phyloseq_filter_sample_wise_abund_trim(
            ps0, 
            minabund = input$filter_taxa_sums_threshold,
            relabund = FALSE,
            rm_zero_OTUs = TRUE)
        out$filter_sample_wise_abund_counts <- ntaxa(ps0)
        
    }
    
    # Per sample RA
    if( input$filter_taxa_ra_threshold > 0 & is(ps0, "phyloseq")){
        ps0 <- phyloseq_filter_sample_wise_abund_trim(
            ps0, 
            minabund = input$filter_taxa_ra_threshold / 100,
            relabund = TRUE,
            rm_zero_OTUs = TRUE)
        out$filter_sample_wise_abund_ra <- ntaxa(ps0)
    }
    
    # Total dataset mean relative abundance.
    if( input$filter_taxa_mean_threshold > 0 & is(ps0, "phyloseq")){
        ps0 <- phyloseq_filter_taxa_rel_abund(
            ps0, 
            frac = input$filter_taxa_mean_threshold / 100)
        out$filter_taxa_mean_threshold <- ntaxa(ps0)
    }
    
    #Remove taxa with abundance less then a certain fraction of total abundance.
    if( input$filter_taxa_total_threshold > 0 & is(ps0, "phyloseq")){
        ps0 <- phyloseq_filter_taxa_tot_fraction(
            ps0, 
            frac = input$filter_taxa_total_threshold / 100)
        out$filter_taxa_total_threshold <- ntaxa(ps0)
        
    }
    
    # K over A filtering 
    if( input$filter_kOverA_sample_threshold > 1 & is(ps0, "phyloseq")){
        ps0 <- phyloseq_filter_prevalence(
            ps0, 
            prev.trh = input$filter_kOverA_sample_threshold / 100,
            abund.trh = input$filter_kOverA_count_threshold / 100,
            threshold_condition = "AND",
            abund.type = "mean")
        out$koverA <- ntaxa(ps0)
    }
    
    # Extract most abundant taxa
    if( input$filter_top_taxa < Inf & is(ps0, "phyloseq")){
        ps0 <- phyloseq_filter_top_taxa(
            ps0, 
            n = input$filter_top_taxa
        )
        out$filter_top_taxa <- ntaxa(ps0)
    }
    
    # Filter sample  total sums
    if( input$filter_sample_sums_threshold > 0 & is(ps0, "phyloseq")){
        # Sample sums filtering
        if(any(sample_sums(ps0) > input$filter_sample_sums_threshold)){
            ps0 <- prune_samples({sample_sums(ps0) > input$filter_sample_sums_threshold}, ps0)
        } else {
            ps0 <- NULL
        }
        out$filter_sample_sums_threshold <- ntaxa(ps0)
    }
    
    # Sample richness filtering
    if( input$filter_sample_taxa_threshold < Inf & is(ps0, "phyloseq")){
        ps0 <- phyloseq_richness_filter(ps0, mintaxa = input$filter_sample_taxa_threshold)
        out$filter_richness <- ntaxa(ps0)
    }
    out$ps <- ps0
    return(out)
}
