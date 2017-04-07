package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;


public class StructuralVariationDiscoveryArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static class FindBreakpointEvidenceSparkArgumentCollection implements Serializable {
        private static final long serialVersionUID = 1L;

        //--------- parameters ----------

        // no-arg constructor for Params object establishes default values
        @VisibleForTesting
        static final FindBreakpointEvidenceSpark.Params defaultParams = new FindBreakpointEvidenceSpark.Params();

        @Argument(doc = "Kmer size.", fullName = "kSize")
        int kSize = defaultParams.kSize;

        @Argument(doc = "maximum kmer DUST score", fullName = "kmerMaxDUSTScore")
        int maxDUSTScore = SVConstants.MAX_DUST_SCORE;

        @Argument(doc = "The minimum mapping quality for reads used to gather evidence of breakpoints.",
                fullName = "minEvidenceMapQ", optional = true)
        int minEvidenceMapQ = defaultParams.minEvidenceMapQ;

        @Argument(doc = "The minimum length of the matched portion of an interesting alignment.  "+
                "Reads that don't match at least this many reference bases won't be used in gathering evidence.",
                fullName = "minEvidenceMatchLength", optional = true)
        int minEvidenceMatchLength = defaultParams.minEvidenceMatchLength;

        @Argument(doc = "Intervals with more than this much coverage are filtered out, because the reads mapped to "+
                "that interval are clearly not exclusively local to the interval.", fullName = "maxIntervalCoverage")
        int maxIntervalCoverage = defaultParams.maxIntervalCoverage;

        @Argument(doc = "Minimum number of reads in cluster to declare an interval of interest.",
                fullName = "minEvidenceCount")
        int minEvidenceCount = defaultParams.minEvidenceCount;

        @Argument(doc = "Minimum number of localizing kmers in a valid interval.", fullName="minKmersPerInterval")
        int minKmersPerInterval = defaultParams.minKmersPerInterval;

        @Argument(doc = "KmerCleaner maximum number of intervals for a localizing kmer.", fullName = "cleanerMaxIntervals")
        int cleanerMaxIntervals = defaultParams.cleanerMaxIntervals;

        @Argument(doc = "KmerCleaner minimum kmer count.", fullName = "cleanerMinKmerCount")
        int cleanerMinKmerCount = defaultParams.cleanerMinKmerCount;

        @Argument(doc = "KmerCleaner maximum kmer count.", fullName = "cleanerMaxKmerCount")
        int cleanerMaxKmerCount = defaultParams.cleanerMaxKmerCount;

        @Argument(doc = "KmerCleaner unique error-free kmers per partition", fullName = "cleanerKmersPerPartitionGuess")
        int cleanerKmersPerPartitionGuess = defaultParams.cleanerKmersPerPartitionGuess;

        @Argument(doc = "Maximum number of templates containing an assembly kmer.", fullName = "maxQNamesPerKmer")
        int maxQNamesPerKmer = defaultParams.maxQNamesPerKmer;

        @Argument(doc = "Guess at number of clean kmers per assembly partition.", fullName = "assemblyKmerMapSize")
        int assemblyKmerMapSize = defaultParams.assemblyKmerMapSize;

        @Argument(doc = "Guess at the ratio of reads in the final assembly to the number reads mapped to the interval.",
                fullName = "assemblyToMappedSizeRatioGuess")
        int assemblyToMappedSizeRatioGuess = defaultParams.assemblyToMappedSizeRatioGuess;

        @Argument(doc = "Maximum FASTQ file size.", fullName = "maxFASTQSize")
        int maxFASTQSize = defaultParams.maxFASTQSize;

        @Argument(doc = "Exclusion interval padding.", fullName = "exclusionIntervalPadding")
        int exclusionIntervalPadding = defaultParams.exclusionIntervalPadding;

        @Argument(doc = "Include read mapping location in FASTQ files.", fullName = "includeMappingLocation")
        boolean includeMappingLocation = true;

        @Argument(doc = "Don't look for extra reads mapped outside the interval.", fullName = "intervalOnlyAssembly")
        boolean intervalOnlyAssembly = false;

        // --------- locations ----------

        @Argument(doc = "bwa-mem index image file", fullName = "alignerIndexImage")
        String alignerIndexImageFile;

        @Argument(doc = "file for read metadata", fullName = "readMetadata", optional = true)
        String metadataFile;

        @Argument(doc = "directory for evidence output", fullName = "breakpointEvidenceDir", optional = true)
        String evidenceDir;

        @Argument(doc = "file for breakpoint intervals output", fullName = "breakpointIntervals", optional = true)
        String intervalFile;

        @Argument(doc = "file for mapped qname intervals output", fullName = "qnameIntervalsMapped", optional = true)
        String qNamesMappedFile;

        @Argument(doc = "file for kmer intervals output", fullName = "kmerIntervals", optional = true)
        String kmerFile;

        @Argument(doc = "file for mapped qname intervals output", fullName = "qnameIntervalsForAssembly", optional = true)
        String qNamesAssemblyFile;

        @Argument(doc = "output dir for assembled fastqs", fullName = "fastqDir", optional = true)
        String fastqDir;

        @Argument(doc = "output dir for assemblies", fullName = "gfaDir", optional = true)
        String gfaDir;

        /**
         * This is a file that calls out the coordinates of intervals in the reference assembly to exclude from
         * consideration when calling putative breakpoints.
         * Each line is a tab-delimited interval with 1-based inclusive coordinates like this:
         *  chr1	124535434	142535434
         */
        @Argument(doc = "file of reference intervals to exclude", fullName = "exclusionIntervals", optional = true)
        String exclusionIntervalsFile;

        /**
         * This is a path to a file of kmers that appear too frequently in the reference to be usable as probes to localize
         * reads.  We don't calculate it here, because it depends only on the reference.
         * The program FindBadGenomicKmersSpark can produce such a list for you.
         */
        @Argument(doc = "file containing ubiquitous kmer list. see FindBadGenomicKmersSpark to generate it.",
                fullName = "kmersToIgnore")
        String kmersToIgnoreFile;

    }

    public static class DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection implements Serializable {
        private static final long serialVersionUID = 1L;

        // todo: document this better
        // Currently the discovery stage requires a reference parameter in 2bit format (to broadcast) and
        // a reference in FASTA format (to get a good sequence dictionary for sorting variants).
        @Argument(doc = "FASTA formatted reference", shortName = "fastaReference",
                fullName = "fastaReference")
        String fastaReference;

        @Argument(doc = "To log simple statistics of the contig alignments or not", shortName = "logContigAlignmentSimpleStats",
                fullName = "logContigAlignmentSimpleStats", optional = true)
        boolean logContigAlignmentSimpleStats = false;

        @Argument(doc = "Minimum flanking alignment length", shortName = "minAlignLength",
                fullName = "minAlignLength", optional = true)
        Integer minAlignLength = SVConstants.DiscoveryStepConstants.DEFAULT_MIN_ALIGNMENT_LENGTH;
    }

}
