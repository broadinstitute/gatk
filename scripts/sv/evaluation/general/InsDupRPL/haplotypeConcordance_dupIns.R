## Script for comparing haplotypes between ins/dup records
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

haplotypeConcordance_dupIns <- function(gatkVSvalidationOverlapFile, 
                                        referenceFASTA, 
                                        output) {
    
    tandupWithInsertion <- read.table(gatkVSvalidationOverlapFile,
                                      stringsAsFactors = F, header = F)
    
    if (ncol(tandupWithInsertion)==16) {
        names(tandupWithInsertion) <- c("gatkCHR", "gatkPOS", "gatkID", "gatkREF", 
                                        "gatkALT", "gatkQUAL", "gatkFILTER", "gatkINFO", 
                                        "validationCHR", "validationPOS", "validationID", "validationREF", "validationALT", "validationQUAL", 
                                        "validationFILTER", "validationINFO")
    } else {
        names(tandupWithInsertion) <- c("gatkCHR", "gatkPOS", "gatkID", "gatkREF", 
                                        "gatkALT", "gatkQUAL", "gatkFILTER", "gatkINFO", 
                                        "validationCHR", "validationPOS", "validationID", "validationREF", "validationALT", "validationQUAL", 
                                        "validationFILTER", "validationINFO", 
                                        "validationFORMAT", "validationSAMPLE")
    }
    
    # first check that all Manta insertions here are not upstream of the overlapping duplication call in GATK
    if ( all(tandupWithInsertion$validationPOS >= tandupWithInsertion$gatkPOS) ){
        cat("No validation call set insertions here is upstream of the matching duplication call in GATK.\n")
    } else {
        stop("Some validation call set insertions here are upstream of the matching duplication call in GATK!")
    }
    
    # then check that all GATK variants have duplication numbers on the reference and 
    # on the alt differ by 1
    cat("Check if all GATK duplications here have duplication numbers on the reference and on the alt differ by 1.\n")
    all(apply(tandupWithInsertion, 
              1, 
              function(x) {
                  dupSituation <- str_extract(x[8], "DUP_NUM=[0-9]+,[0-9]+")
                  dupNums <- strsplit(substr(dupSituation, 9, nchar(dupSituation)), 
                                      ",")[[1]]
                  as.integer(dupNums[2]) - as.integer(dupNums[1])
              })==1
    )
    
    # extract duplicated sequence from reference fasta via bedtools
    tandupWithInsertion <- adply(tandupWithInsertion, 
                              1, function(x) {
                                  result <- list()
                                  span <- str_extract(x$gatkINFO, "DUP_REPET_UNIT_REF_SPAN=(chr)?([0-9]{1,2}|X|Y):[0-9]+-[0-9]+")
                                  span <- sub("DUP_REPET_UNIT_REF_SPAN=", "", span)
                                  span <- gsub("(:|-)", "	", span)
                                  result[["gatkDupRefEnd"]] <- strsplit(span, "	")[[1]][3]
                                  unlist(result)
                              }
    )
    write.table(subset(tandupWithInsertion, select=c(gatkCHR, gatkPOS, gatkDupRefEnd)), 
                file = "temp.gatk.bed",
                col.names = F, row.names = F, sep = "	", quote = F)
    system(paste0("bedtools getfasta -fi ", referenceFASTA, " -bed temp.gatk.bed -fo temp.gatk.fasta"))
    dupSeq_gatk <- readLines("temp.gatk.fasta")
    dupSeq_gatk <- dupSeq_gatk[seq(from=2, to=length(dupSeq_gatk), by=2)]
    ret <- file.remove("temp.gatk.bed", "temp.gatk.fasta")
    
    # Haplotype reconstruction
    gatkHaptype <- rep("", nrow(tandupWithInsertion))
    mantaHaptype <- rep("", nrow(tandupWithInsertion))
    for(i in c(1:length(mantaHaptype))){
        
        gatkDupSeq <- str_extract(tandupWithInsertion$gatkINFO[i], "DUPLICATED_SEQUENCE=[A-Z]+")
        gatkDupSeq <- substr(gatkDupSeq, 21, nchar(gatkDupSeq))
        gatkHaptype[i] <- paste0(tandupWithInsertion$gatkREF[i], gatkDupSeq)
        
        gatkInsSeq <- str_extract(tandupWithInsertion$gatkINFO[i], "INSSEQ=[A-Z]+")
        if( ! is.na(gatkInsSeq)) {
            gatkHaptype[i] <- paste0(gatkHaptype[i], sub("INSSEQ=", "", gatkInsSeq))
        }
        
        mantaHaptype[i] <- tandupWithInsertion$mantaALT[i]
    }
    sink(output)
    cat("gatkID	editDist	gatkAltLen	mantaAltLen\n")
    for(i in c(1:length(mantaHaptype))){
        
        dist <- as.integer(adist(gatkHaptype[i], mantaHaptype[i]))
        cat(tandupWithInsertion$gatkID[i], 
            dist, 
            nchar(gatkHaptype[i]), 
            nchar(mantaHaptype[i]), 
            sep="	", fill = T)
    }
    sink()
    cat("Results of comparing duplication/insertion haplotypes is written to file:\n")
    cat(output, fill = T)
    cat("Please inspect the file for concordance.\n")
}

save(haplotypeConcordance_insDup, file="func_haplotypeConcordance_insDup.RData")