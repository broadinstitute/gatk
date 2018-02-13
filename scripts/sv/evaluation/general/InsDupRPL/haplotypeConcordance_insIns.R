## Script for comparing haplotypes between insertion/insertion records
## for provided overlapping files

# extract information to construct haplotypes
extractInsertedSeqInfo <- function(x) {
    library(stringr)
    insAnnotation <- sub("INSSEQ=", "", str_extract(x$gatkINFO, "INSSEQ=[A-Z]+"))
    
    validationHaptype <- ""
    if ( x$validationID=="." ){ # pacbio doesn't provide a custom ID for its variants
        gatkHaptype <- insAnnotation
        validationHaptype <- toupper(sub("SEQ=", "", str_extract(x$validationINFO, "SEQ=[:alpha:]+")))
    } else {
        gatkHaptype <- paste0(x$gatkREF, insAnnotation)
        validationHaptype <- x$validationALT
    }
    
    info <- list()
    info[["editDist"]] <- as.integer(adist(gatkHaptype, validationHaptype)[1,1])
    info[["gatkAltLen"]] <- nchar(gatkHaptype)
    info[["validationAltLen"]] <- nchar(validationHaptype)
    unlist(info)
}

haplotypeConcordance_insIns <- function(gatkVSvalidationOverlapFile, 
                                        output) {
    
    insertionWithInsertion <- read.table(gatkVSvalidationOverlapFile)

    if (ncol(insertionWithInsertion)==16) {
        names(insertionWithInsertion) <- c("gatkCHR", "gatkPOS", "gatkID", "gatkREF", 
                                           "gatkALT", "gatkQUAL", "gatkFILTER", "gatkINFO", 
                                           "validationCHR", "validationPOS", "validationID", "validationREF", "validationALT", "validationQUAL", 
                                           "validationFILTER", "validationINFO")
    } else {
        names(insertionWithInsertion) <- c("gatkCHR", "gatkPOS", "gatkID", "gatkREF", 
                                           "gatkALT", "gatkQUAL", "gatkFILTER", "gatkINFO", 
                                           "validationCHR", "validationPOS", "validationID", "validationREF", "validationALT", "validationQUAL", 
                                           "validationFILTER", "validationINFO", 
                                           "validationFORMAT", "validationSAMPLE")
    }
    
    library(plyr)
    result <- adply(insertionWithInsertion, 1, extractInsertedSeqInfo)
    write.table(subset(result, select=c(gatkID, editDist, gatkAltLen, validationAltLen)), 
                file = output, 
                row.names = F, sep = "	", quote = F)
    cat("Results of comparing insertion/insertion haplotypes is written to file:\n")
    cat(output, fill = T)
    cat("Please inspect the file for concordance.\n")
}

save(haplotypeConcordance_insIns, file="func_haplotypeConcordnace_insIns.RData")