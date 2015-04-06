args = commandArgs(trailingOnly=TRUE)
opt = list(details.fn=args[1], summary.fn=args[2], output.fn=args[3])

read_metrics_file = function(metrics.fn) {
  contents = read.delim(metrics.fn, comment.char="#", stringsAsFactors=FALSE)
  return(contents)
}

equals_or_is_na = function(x1, x2) {
  if (is.na(x1)) {
    return(is.na(x2))
  } else {
    return(x1 == x2)
  }
}

details = read_metrics_file(opt$details.fn)
summary = read_metrics_file(opt$summary.fn)

pdf(opt$output.fn)
par(mfrow=c(2,2), oma=c(0,0,2,0))

for (i in seq_len(nrow(summary))) {
  cur_summary = summary[i, ]
  cur_sample = cur_summary[1, "SAMPLE"]
  cur_library = cur_summary[1, "LIBRARY"]
  cur_read_group = cur_summary[1, "READ_GROUP"]
  cur_details = details[which((equals_or_is_na(cur_library, details[, "LIBRARY"]) &
    (equals_or_is_na(cur_sample, details[, "SAMPLE"])) &
    (equals_or_is_na(cur_read_group, details[, "READ_GROUP"])))), ]
  
  
  ## Plot conversion rates
  cpg.converted = sum(cur_details$CONVERTED_SITES)
  cpg.seen = sum(cur_details$TOTAL_SITES)
  cpg.conversion = cpg.converted / cpg.seen
  total.conversion = (cpg.converted + cur_summary$NON_CPG_CONVERTED_BASES) / (cpg.seen + cur_summary$NON_CPG_BASES)

  barplot(c("non-CpG"=cur_summary$PCT_NON_CPG_BASES_CONVERTED, "Combined"=total.conversion, "CpG"=cpg.conversion),
          ylim=c(0.95, 1), ylab="% Conversion", xlab="Distribution", main="Bisulfite Conversion Rate",
          col="blue", xpd=FALSE)
  abline(h=0.995, col="grey")
  
  ## Plot histogram of CpG counts by conversion rate
  hist(cur_details$PCT_CONVERTED, 10, xlab="Conversion Rate Of CpGs", ylab="# CpGs",
       main="CpG Conversion Rate Distribution", col="blue")
  
  ## Plot pie chart showing distribution of CpG coverage
  coverage_breaks = c(0, 1, 5, 10, 25, 50, 100, Inf)
  coverage_cut = cut(cur_details$TOTAL_SITES, coverage_breaks)
  cpg_coverage = split(cur_details$TOTAL_SITES, coverage_cut)
  coverages = sapply(cpg_coverage, length)[2:7]
  names(coverages) = paste(">=", c(1, 5, 10, 25, 50, 100), sep="")
  ## If we have 0s all across the pie chart will be effectively meaningless but put in a 100% >= 0 field instead
  ## to avoid an error on pie(). Normally it'd just be a pain to see these, but ...
  if (all(coverages == 0)) {
    coverages = c("No Coverage"=1)
  }
  color_ramp = colorRampPalette(c("white", "#538ED5", "blue"), bias=1, space="Lab")
  colors = color_ramp(length(coverages))[2:length(coverages)]
  pie(coverages, main="Distribution Of CpGs By Coverage", col=colors, clockwise=TRUE)
  
  discards = log10(c("Mismatches"=cur_summary$READS_IGNORED_MISMATCHES, "Size"=cur_summary$READS_IGNORED_SHORT))
  ## Protect against -Inf in the case where we had 0 discards
  discards = ifelse(is.finite(discards), discards, 0)
  barplot(discards, ylab="Number Discarded (log10)", xlab="Reason",
          main="Reads Discarded", col="blue", ylim=c(0, ceiling(max(discards))))

  header_txt = character()
  if (!is.na(cur_sample) && cur_sample != "") {
    header_txt = paste(header_txt, " SAMPLE=", cur_sample, sep="")
  }
  if (!is.na(cur_library) && cur_library != "") {
    header_txt = paste(header_txt, " LIBRARY=", cur_library, sep="")
  }
  if (!is.na(cur_read_group) && cur_read_group != "") {
    header_txt = paste(header_txt, " READ GROUP=", cur_read_group, sep="")
  }
  if (length(header_txt) > 0) {
    mtext(header_txt, outer=TRUE, line=1)
  }
}

dev.off()
  
