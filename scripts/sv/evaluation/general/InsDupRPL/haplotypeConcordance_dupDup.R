## Script for comparing haplotypes between dup/dup records
## for provided overlapping files
## note that this is only applicable between GATK and Manta


haplotypeConcordance_dupDup <- function(gatkVSmantaOverlapFile, 
                                        referenceFASTA, 
                                        output) {
    tandupWithTandup <- read.table(gatkVSmantaOverlapFile, 
                                   stringsAsFactors = F, header = F)
    
    names(tandupWithTandup) <- c("gatkCHR", "gatkPOS", "gatkID", "gatkREF", "gatkALT", "gatkQUAL", 
                                 "gatkFILTER", "gatkINFO", 
                                 "mantaCHR", "mantaPOS", "mantaID", "mantaREF", "mantaALT", "mantaQUAL", 
                                 "mantaFILTER", "mantaINFO", 
                                 "mantaFORMAT", "mantaSAMPLE")
    
    library(plyr)
    
    # first extract the duplicated sequences using bedtools
    tandupWithTandup <- adply(tandupWithTandup, 
                              1, function(x) {
                                  result <- list()
                                  result[["mantaEND"]] <- sub("END=", "", str_extract(x$mantaINFO, "END=[0-9]+"))
                                  span <- str_extract(x$gatkINFO, "DUP_REPET_UNIT_REF_SPAN=(chr)?([0-9]{1,2}|X|Y):[0-9]+-[0-9]+")
                                  span <- sub("DUP_REPET_UNIT_REF_SPAN=", "", span)
                                  span <- gsub("(:|-)", "	", span)
                                  result[["gatkDupRefEnd"]] <- strsplit(span, "	")[[1]][3]
                                  unlist(result)
                              }
    )
    write.table(subset(tandupWithTandup, select=c(mantaCHR, mantaPOS, mantaEND)), 
                file = "temp.manta.bed",
                col.names = F, row.names = F, sep = "	", quote = F)
    write.table(subset(tandupWithTandup, select=c(gatkCHR, gatkPOS, gatkDupRefEnd)), 
                file = "temp.gatk.bed",
                col.names = F, row.names = F, sep = "	", quote = F)
    system(paste0("bedtools getfasta -fi ", referenceFASTA, " -bed temp.manta.bed -fo temp.manta.fasta"))
    system(paste0("bedtools getfasta -fi ", referenceFASTA, " -bed temp.gatk.bed -fo temp.gatk.fasta"))
    refSeq_manta <- readLines("temp.manta.fasta")
    refSeq_manta <- refSeq_manta[seq(from=2, to=length(refSeq_manta), by=2)]
    refSeq_gatk <- readLines("temp.gatk.fasta")
    refSeq_gatk <- refSeq_gatk[seq(from=2, to=length(refSeq_gatk), by=2)]

    ret <- file.remove("temp.manta.bed", "temp.manta.fasta", "temp.gatk.bed", "temp.gatk.fasta")

    # a simplification exists that the duplication calls in GATK are all 1 -> 2 duplications
    cat("Check if all duplication calls in GATK overlapping a uniq Manta duplication call have DUP_NUM=1,2.\n")
    all("DUP_NUM=1,2"==apply(tandupWithTandup,
                             1,
                             function(x) {
                                 str_extract(x[8], "DUP_NUM=[0-9]+,[0-9]+")
                             }))

    # then it is relatively easy to reconstruct the haplotypes
    gatkHaptype <- rep("", nrow(tandupWithTandup))
    mantaHaptype <- rep("", nrow(tandupWithTandup))
    for(i in c(1:nrow(tandupWithTandup))){

        gatkHaptype[i] <- refSeq_gatk[i]
        insSeq <- str_extract(tandupWithTandup$gatkINFO[i], "INSSEQ=[A-Z]+")
        if(!is.na(insSeq)) {
            gatkHaptype[i] <- paste0(gatkHaptype[i], sub("INSSEQ=", "", insSeq))
        }

        mantaHaptype[i] <- refSeq_manta[i]
        insSeq <- str_extract(tandupWithTandup$mantaINFO[i], "SVINSSEQ=[A-Z]+")
        if(!is.na(insSeq)) {
            mantaHaptype[i] <- paste0(mantaHaptype[i],  sub("SVINSSEQ=", "", insSeq))
        }
    }
    # write result to file
    sink(output)
    cat("gatkID	editDist	gatkAltLen	mantaAltLen\n")
    for(i in c(1:length(gatkHaptype))){
        cat(tandupWithTandup$gatkID[i],
            as.integer(adist(gatkHaptype[i], mantaHaptype[i])[1,1]),
            nchar(gatkHaptype[i]),
            nchar(mantaHaptype[i]),
            sep="	", fill = T)
    }
    sink()
    cat("Results of comparing dup/dup is written to file:\n")
    cat(output, fill = T)
    cat("Please inspect the file for concordance.\n")
}
save(haplotypeConcordance_dupDup, file="func_haplotypeConcordance_dupDup.RData")