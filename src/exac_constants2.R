library(dplyr)
library(magrittr)

source("lib/load_exac_data2.R")

# Load exac data into environment. Steps 0/6 and 1/6
exac <- load_exac_data2()

# Determining sequence contexts
message('Now processing data (2/6): Determining sequence contexts...')

colnames(exac) <- tolower(colnames(exac))
set.seed(99)
exac_sample <- exac[sample(1:nrow(exac), nrow(exac)*.1),] %>% 
  mutate(pos_id = paste(chrom, formatC(pos, width=9, flag='0'), ref, alt, sep='_'),
         indel = nchar(ref) != nchar(alt),
         bases_inserted = nchar(alt) - nchar(ref),
         transition = (ref == 'A' & alt == 'G') | (ref == 'G' & alt == 'A') | (ref == 'C' & alt == 'T') | (ref == 'T' & alt == 'C'),
         # transition = if_else(indel, NA, transition), # seems pointless to set it NA when it's honest FALSE
         cpg = (ref == 'C' & alt == 'T' & substr(context, 2, 3) == 'CG') | (ref == 'G' & alt == 'A' & substr(context, 1, 2) == 'CG'),
         # cpg = if_else(indel, NA, cpg), # seems pointless to set it NA when it's honest FALSE
         alt_is_ancestral = alt == ancestral)

# exac$alt_is_ancestral[exac$ref != exac$ancestral & exac$alt != exac$ancestral]<-NA

# Computing allele frequencies...
message('(3/6): Computing allele frequencies...')
calc_ratio <- function(n, dn){
  frac <- n/dn
  ifelse(is.finite(frac), frac, 0)
}

exac_sample %<>% 
  mutate(af_popmax = calc_ratio(ac_popmax, an_popmax),
         af_global = calc_ratio(ac_adj, an_adj),
         singleton = ac_adj == 1,
         maf_global = pmin(af_global, 1-af_global),
         daf_global = ifelse(alt_is_ancestral, 1-af_global, af_global),
         daf_popmax = ifelse(alt_is_ancestral, 1-af_popmax, af_popmax))

# Calculating bad regions of the genome
message('(4/6): Calculating bad regions of the genome...')
resolution <- 1000
number_of_regions <- 10
intervals <- (0:(250000000/resolution))*resolution

exac_sample %<>% 
  mutate(pos_bin = cut(pos, intervals),
         bin_start = as.numeric(sub("\\((.+),.*", "\\1", pos_bin)))
  

allelic_state <- plyr::count(subset(exac_sample, select=c(chrom, pos)))
multiallelics <- subset(allelic_state, freq > 3, select=c(chrom, pos))
multiallelics$pos_bin <- cut(multiallelics$pos, intervals)
multiallelics$bin_start <- as.numeric( sub("\\((.+),.*", "\\1", multiallelics$pos_bin))
multiallelic_counts <- plyr::count(multiallelics, vars = c('chrom', 'bin_start'))

# why only head?
bad_sectors <- subset(head(multiallelic_counts[order(multiallelic_counts$freq, decreasing = T),], number_of_regions), select=c(chrom, bin_start))
bad_sectors$sector <- paste(bad_sectors$chrom, bad_sectors$bin_start, sep='_')

# Determining what to use
message('(5/6): Determining what to _use_...')
exac$sector <- paste(exac$chrom, exac$bin_start, sep='_')
exac$bad <- exac$sector %in% bad_sectors$sector
exac$use <- exac$an_adj > .80*max(exac$an_adj, na.rm=TRUE) & exac$ac_adj > 0 & exac$filter=='PASS' & !exac$bad
exac$lof_use <- !is.na(exac$lof) & exac$lof == 'HC' & is.na(exac$lof_flags)

