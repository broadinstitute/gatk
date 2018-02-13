# Collects sizes of insertions, deletions and tandem duplications, 
# and return them as objects in memory
#########################################################################
collectSizes <- function(mantaDelSizesFile,
                         mantaInsSizesFile, 
                         mantaDupSizesFile,
                         gatkDelSizesFile,
                         gatkInsSizesFile,
                         gatkDupSizesFile,
                         pacbioDelSizesFile,
                         pacbioInsSizesFile, 
                         mantaScarDelVCF, 
                         gatkScarDelVCF) {
    
    library("plyr")
    library("stringr")
    
    result <- list()
    
    ############
    # clean deletion
    result[["cleanDelSizes_manta"]] <- sort(scan(mantaDelSizesFile))
    result[["cleanDelSizes_gatk"]]  <- sort(scan(gatkDelSizesFile))
    
    ############ 
    # clean insertion and tandem duplications
    cleanInsSizes_manta <- sort(scan(mantaInsSizesFile))
    cleanInsSizes_gatk <- sort(scan(gatkInsSizesFile))

    tandupSizes_manta <- sort(scan(mantaDupSizesFile))
    tandupSizes_gatk <- sort(scan(gatkDupSizesFile))
    
    result[["insertionAndDup_manta"]] <- sort(c(tandupSizes_manta, cleanInsSizes_manta))
    result[["insertionAndDup_gatk"]] <- sort(c(tandupSizes_gatk, cleanInsSizes_gatk))
    
    ############ 
    # scarred deletion
    # each entry an array of 2: deleted bases length and inserted bases length
    mantaScars <- read.table(mantaScarDelVCF)
    names(mantaScars) <- c("CHROM", "POS", "ID", "REF", "ALT", 
                           "QUAL", "FILTER", "INFO", 
                           "FORMAT", "SAMPLE")
    temp <- adply(mantaScars, 1, 
                  function(s) {
                      x <- str_extract(s$INFO, "1M[0-9]+I[0-9]+D")
                      del <- 0
                      ins <- 0
                      if (is.na(x)) { # those without CIGAR but with SVINSSEQ/SVINSLEN
                          del <- as.integer( sub("SVLEN=-", "", str_extract(s$INFO, "SVLEN=-[0-9]+")) )
                          ins <- as.integer( sub("SVINSLEN=", "", str_extract(s$INFO, "SVINSLEN=[0-9]+")))
                      } else {
                          ins <- as.integer( sub("I", "", str_extract(x, "[0-9]+I")) )
                          del <- as.integer( sub("D", "", str_extract(x, "[0-9]+D")) )
                      }
                      c(del, ins)
                  })
    
    result[["scarDel_manta_del"]] <- as.integer(temp$V1)
    result[["scarDel_manta_ins"]] <- as.integer(temp$V2)
    
    gatkScars <- read.table(gatkScarDelVCF)
    names(gatkScars) <- c("CHROM", "POS", "ID", "REF", "ALT", 
                          "QUAL", "FILTER", "INFO")
    
    temp <- adply(gatkScars, 1, 
                  function(s) {
                      del <- as.integer( sub("SVLEN=-", "", str_extract(s$INFO, "SVLEN=-[0-9]+")) )
                      ins <- str_extract(s$INFO, "INSSEQ=[A-Z]+")
                      ins <- as.integer( nchar(ins) - nchar("INSSEQ=") )
                      
                      c(del, ins)
                  })
    
    result[["scarDel_gatk_del"]] <- as.integer(temp$V1)
    result[["scarDel_gatk_ins"]] <- as.integer(temp$V2)
    
    ############ 
    # pacbio call set
    if ( !is.null(pacbioDelSizesFile) & !is.null(pacbioInsSizesFile) ) {
        result[["pacbioDeletionSizes"]] <- sort(scan(pacbioDelSizesFile))
        result[["pacbioInsertionSizes"]] <- sort(scan(pacbioInsSizesFile))
    }

    result
}

save.image("func_collectSizes.Rdata")
