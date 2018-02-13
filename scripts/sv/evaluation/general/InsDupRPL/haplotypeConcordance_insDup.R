## Script for comparing haplotypes between ins/dup records
## for provided overlapping files
## note that this is only applicable between GATK and Manta

haplotypeConcordance_insDup <- function(gatkVSmantaOverlapFile, 
                                        referenceFASTA, 
                                        output) {
    
    insertionWithTandup <- read.table(gatkVSmantaOverlapFile,
                                      stringsAsFactors = F, header = F)
    names(insertionWithTandup) <- c("gatkCHR", "gatkPOS", "gatkID", "gatkREF", "gatkALT", "gatkQUAL", 
                                    "gatkFILTER", "gatkINFO", 
                                    "mantaCHR", "mantaPOS", "mantaID", "mantaREF", "mantaALT", "mantaQUAL", 
                                    "mantaFILTER", "mantaINFO", 
                                    "mantaFORMAT", "mantaSAMPLE")
    
    # first check that all GATK variants here are not upstream of the duplication call in Manta
    if ( all(insertionWithTandup$gatkPos >= insertionWithTandup$mantaPos) ){
        cat("No GATK insertions here is upstream of the matching duplication call in Manta.\n")
    } else {
        stop("Some GATK insertions here are upstream of the matching duplication call in Manta!")
    }
    
    ###############
    # two cases are possible in this type of intersection: 
    #    when the Manta and GATK records have the same POS, and 
    #    when the Manta POS is upstream of that of GATK
    # define 
    #    $\alpha$ := refSeq(mantaPOS, gatkPOS-1]
    #    $\beta$  := refSeq(gatkPOS, mantaEND]
    #    $\gamma$ := refSeq(mantaPOS, mantaEND]
    # for the 1st case,  $\alpha$ would be empty while $\beta$ and $\gamma$ point to the same interval
    #    GATK's haplotype is : gatkRefBase + INSSEQ(annotation from gatk) + $\beta$
    #    Manta's haplotype is: mantaRefBase + $\gamma$ + SVINSSEQ(annotation from manta) + $\gamma$
    # for the 2nd case, we have 
    #    GATK's haplotype is : $\alpha$ + gatkRefBase + INSERTED_QUENCE(annotation from gatk) + $\beta$
    #    Manta's haplotype is: mantaRefBase + $\gamma$ + SVINSSEQ(annotation from manta) + $\gamma$
    # it is then apparent then that the the 2nd case covers the 1st case, in terms of 
    #   implementation, one just need to check if $\alpha$ is empty
    
    ###############
    # ask Bedtools for reference bases
    tempBed <- cbind(insertionWithTandup[, c(9, 10)], # mantaChr and mantaPOS
                     apply(insertionWithTandup, 
                           1, function(x) {
                               endPosition <- str_extract(x[16], "END=[0-9]+") # extract mantaEND
                               as.integer(substr(endPosition, 5, nchar(endPosition)))
                           }
                     ))
    
    write.table(tempBed, file = "temp.bed", quote=F, sep="	", row.names = F, col.names = F)
    system(paste0("bedtools getfasta -fi ", referenceFASTA, " -bed temp.bed -fo temp.fasta"))
    refSeq <- readLines("temp.fasta")
    refSeq <- refSeq[seq(from=2, to=length(refSeq), by=2)]
    ret <- file.remove("temp.bed", "temp.fasta")
    
    # extract reference base
    gatkRefBase <- insertionWithTandup$gatkRef
    mantaRefBase <- insertionWithTandup$mantaRef
    
    # do the job
    gatkHaptype <- rep("", nrow(insertionWithTandup))
    mantaHaptype <- rep("", nrow(insertionWithTandup))
    for(i in c(1:nrow(insertionWithTandup))) {
        
        x <- insertionWithTandup$gatkPOS[i] - insertionWithTandup$mantaPOS[i] - 1
        alpha <- substr(refSeq[i], 1, x)
        beta <- substr(refSeq[i], x+2, nchar(refSeq[i]))
        # gamma is refSeq[i]
        
        gatkInsertedSeq <- str_extract(insertionWithTandup$gatkINFO[i], "INSSEQ=[A-Z]+")
        gatkInsertedSeq <- substr(gatkInsertedSeq, 19, nchar(gatkInsertedSeq))
        
        mantaInsertedSeq <- str_extract(insertionWithTandup$mantaINFO[i], "SVINSSEQ=[A-Z]+")
        if(is.na(mantaInsertedSeq)) {
            mantaInsertedSeq <- ""
        } else {
            mantaInsertedSeq <- substr(mantaInsertedSeq, 10, nchar(mantaInsertedSeq))
        }
        
        gatkHaptype[i] <- paste0(alpha, gatkRefBase[i], gatkInsertedSeq, beta)
        mantaHaptype[i] <- paste0(mantaRefBase[i], refSeq[i], mantaInsertedSeq, refSeq[i])
    }
    
    # compute the edit distances between the haplotypes
    sink(output)
    cat("gatkID	editDist	gatkAltLen	mantaAltLen\n")
    for(i in c(1:nrow(insertionWithTandup))){
        cat(insertionWithTandup$gatkID[i], 
            adist(gatkHaptype[i], mantaHaptype[i]), 
             nchar(gatkHaptype[i]), 
            nchar(mantaHaptype[i]), 
            sep="	", fill = T)
    }
    sink()
    cat("Results of comparing insertion/duplication haplotypes is written to file:\n")
    cat(output, fill = T)
    cat("Please inspect the file for concordance.\n")
}
save(haplotypeConcordance_insDup, file="func_haplotypeConcordance_insDup.RData")