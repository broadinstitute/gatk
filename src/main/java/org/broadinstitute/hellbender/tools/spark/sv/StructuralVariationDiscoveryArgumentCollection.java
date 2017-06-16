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
        public static final Params defaultParams = new Params();

        @Argument(doc = "Kmer size.", fullName = "kSize")
        public int kSize = defaultParams.kSize;

        @Argument(doc = "maximum kmer DUST score", fullName = "kmerMaxDUSTScore")
        public int maxDUSTScore = SVConstants.MAX_DUST_SCORE;

        @Argument(doc = "The minimum mapping quality for reads used to gather evidence of breakpoints.",
                fullName = "minEvidenceMapQ", optional = true)
        public int minEvidenceMapQ = defaultParams.minEvidenceMapQ;

        @Argument(doc = "The minimum length of the matched portion of an interesting alignment.  "+
                "Reads that don't match at least this many reference bases won't be used in gathering evidence.",
                fullName = "minEvidenceMatchLength", optional = true)
        public int minEvidenceMatchLength = defaultParams.minEvidenceMatchLength;

        @Argument(doc = "Intervals with more than this much coverage are filtered out, because the reads mapped to "+
                "that interval are clearly not exclusively local to the interval.", fullName = "maxIntervalCoverage")
        public int maxIntervalCoverage = defaultParams.maxIntervalCoverage;

        @Argument(doc = "Minimum weight of the corroborating read evidence to validate some single piece of evidence.",
                fullName = "minEvidenceCount")
        public int minEvidenceWeight = defaultParams.minEvidenceWeight;

        @Argument(doc = "Minimum weight of the evidence that shares a distal target locus to validate the evidence.",
                fullName = "minCoherentEvidenceCount")
        public int minCoherentEvidenceWeight = defaultParams.minCoherentEvidenceWeight;

        @Argument(doc = "Minimum number of localizing kmers in a valid interval.", fullName="minKmersPerInterval")
        public int minKmersPerInterval = defaultParams.minKmersPerInterval;

        @Argument(doc = "KmerCleaner maximum number of intervals for a localizing kmer.", fullName = "cleanerMaxIntervals")
        public int cleanerMaxIntervals = defaultParams.cleanerMaxIntervals;

        @Argument(doc = "KmerCleaner minimum kmer count.", fullName = "cleanerMinKmerCount")
        public int cleanerMinKmerCount = defaultParams.cleanerMinKmerCount;

        @Argument(doc = "KmerCleaner maximum kmer count.", fullName = "cleanerMaxKmerCount")
        public int cleanerMaxKmerCount = defaultParams.cleanerMaxKmerCount;

        @Argument(doc = "KmerCleaner unique error-free kmers per partition", fullName = "cleanerKmersPerPartitionGuess")
        public int cleanerKmersPerPartitionGuess = defaultParams.cleanerKmersPerPartitionGuess;

        @Argument(doc = "Maximum number of templates containing an assembly kmer.", fullName = "maxQNamesPerKmer")
        public int maxQNamesPerKmer = defaultParams.maxQNamesPerKmer;

        @Argument(doc = "Guess at number of clean kmers per assembly partition.", fullName = "assemblyKmerMapSize")
        public int assemblyKmerMapSize = defaultParams.assemblyKmerMapSize;

        @Argument(doc = "Guess at the ratio of reads in the final assembly to the number reads mapped to the interval.",
                fullName = "assemblyToMappedSizeRatioGuess")
        public int assemblyToMappedSizeRatioGuess = defaultParams.assemblyToMappedSizeRatioGuess;

        @Argument(doc = "Maximum FASTQ file size.", fullName = "maxFASTQSize")
        public int maxFASTQSize = defaultParams.maxFASTQSize;

        @Argument(doc = "Exclusion interval padding.", fullName = "exclusionIntervalPadding")
        public int exclusionIntervalPadding = defaultParams.exclusionIntervalPadding;

        @Argument(doc = "Include read mapping location in FASTQ files.", fullName = "includeMappingLocation")
        public boolean includeMappingLocation = true;

        @Argument(doc = "Don't look for extra reads mapped outside the interval.", fullName = "intervalOnlyAssembly")
        public boolean intervalOnlyAssembly = false;

        @VisibleForTesting public static class Params {
            public final int kSize;
            public final int maxDUSTScore;
            public final int minEvidenceMapQ;
            public final int minEvidenceMatchLength;
            public final int maxIntervalCoverage;
            public final int minEvidenceWeight;
            public final int minCoherentEvidenceWeight;
            public final int minKmersPerInterval;
            public final int cleanerMaxIntervals;
            public final int cleanerMinKmerCount;
            public final int cleanerMaxKmerCount;
            public final int cleanerKmersPerPartitionGuess;
            public final int maxQNamesPerKmer;
            public final int assemblyKmerMapSize;
            public final int assemblyToMappedSizeRatioGuess;
            public final int maxFASTQSize;
            public final int exclusionIntervalPadding;

            public Params() {
                kSize = SVConstants.KMER_SIZE;          // kmer size
                maxDUSTScore = SVConstants.MAX_DUST_SCORE;// maximum for DUST-like kmer complexity score
                minEvidenceMapQ = 20;                   // minimum map quality for evidential reads
                minEvidenceMatchLength = 45;            // minimum match length
                maxIntervalCoverage = 1000;             // maximum coverage on breakpoint interval
                minEvidenceWeight = 15;                  // minimum number of evidentiary reads in called cluster
                minCoherentEvidenceWeight = 7;           // minimum number of evidentiary reads in a cluster that all point to the same target locus
                minKmersPerInterval = 20;               // minimum number of good kmers in a valid interval
                cleanerMaxIntervals = 3;                // KmerCleaner maximum number of intervals a localizing kmer can appear in
                cleanerMinKmerCount = 3;                // KmerCleaner min kmer count
                cleanerMaxKmerCount = 125;              // KmerCleaner max kmer count
                cleanerKmersPerPartitionGuess = 600000; // KmerCleaner guess for number of unique error-free kmers per partition
                maxQNamesPerKmer = 500;                 // maximum template names for an assembly kmer
                assemblyKmerMapSize = 250000;           // guess for unique, error-free, scrubbed kmers per assembly partition
                assemblyToMappedSizeRatioGuess = 7;     // guess for ratio of total reads in assembly to evidentiary reads in interval
                maxFASTQSize = 3000000;                 // maximum assembly size (total input bases)
                exclusionIntervalPadding = 0;           // exclusion interval extra padding
            }

            public Params( final int kSize, final int maxDUSTScore, final int minEvidenceMapQ,
                           final int minEvidenceMatchLength, final int maxIntervalCoverage,
                           final int minEvidenceWeight, final int minCoherentEvidenceWeight,
                           final int minKmersPerInterval, final int cleanerMaxIntervals, final int cleanerMinKmerCount,
                           final int cleanerMaxKmerCount, final int cleanerKmersPerPartitionGuess,
                           final int maxQNamesPerKmer, final int assemblyKmerMapSize, final int assemblyToMappedSizeRatioGuess,
                           final int maxFASTQSize, final int exclusionIntervalPadding ) {
                this.kSize = kSize;
                this.maxDUSTScore = maxDUSTScore;
                this.minEvidenceMapQ = minEvidenceMapQ;
                this.minEvidenceMatchLength = minEvidenceMatchLength;
                this.maxIntervalCoverage = maxIntervalCoverage;
                this.minEvidenceWeight = minEvidenceWeight;
                this.minCoherentEvidenceWeight = minCoherentEvidenceWeight;
                this.minKmersPerInterval = minKmersPerInterval;
                this.cleanerMaxIntervals = cleanerMaxIntervals;
                this.cleanerMinKmerCount = cleanerMinKmerCount;
                this.cleanerMaxKmerCount = cleanerMaxKmerCount;
                this.cleanerKmersPerPartitionGuess = cleanerKmersPerPartitionGuess;
                this.maxQNamesPerKmer = maxQNamesPerKmer;
                this.assemblyKmerMapSize = assemblyKmerMapSize;
                this.assemblyToMappedSizeRatioGuess = assemblyToMappedSizeRatioGuess;
                this.maxFASTQSize = maxFASTQSize;
                this.exclusionIntervalPadding = exclusionIntervalPadding;
            }
        }

        // --------- locations ----------

        @Argument(doc = "bwa-mem index image file", fullName = "alignerIndexImage")
        public String alignerIndexImageFile;

        @Argument(doc = "file for read metadata", fullName = "readMetadata", optional = true)
        public String metadataFile;

        @Argument(doc = "directory for evidence output", fullName = "breakpointEvidenceDir", optional = true)
        public String evidenceDir;

        @Argument(doc = "file for breakpoint intervals output", fullName = "breakpointIntervals", optional = true)
        public String intervalFile;

        @Argument(doc = "file for mapped qname intervals output", fullName = "qnameIntervalsMapped", optional = true)
        public String qNamesMappedFile;

        @Argument(doc = "file for kmer intervals output", fullName = "kmerIntervals", optional = true)
        public String kmerFile;

        @Argument(doc = "file for mapped qname intervals output", fullName = "qnameIntervalsForAssembly", optional = true)
        public String qNamesAssemblyFile;

        @Argument(doc = "output dir for assembled fastqs", fullName = "fastqDir", optional = true)
        public String fastqDir;

        @Argument(doc = "output dir for assemblies", fullName = "gfaDir", optional = true)
        public String gfaDir;

        /**
         * This is a file that calls out the coordinates of intervals in the reference assembly to exclude from
         * consideration when calling putative breakpoints.
         * Each line is a tab-delimited interval with 1-based inclusive coordinates like this:
         *  chr1	124535434	142535434
         */
        @Argument(doc = "file of reference intervals to exclude", fullName = "exclusionIntervals", optional = true)
        public String exclusionIntervalsFile;

        /**
         * This is a path to a file of kmers that appear too frequently in the reference to be usable as probes to localize
         * reads.  We don't calculate it here, because it depends only on the reference.
         * The program FindBadGenomicKmersSpark can produce such a list for you.
         */
        @Argument(doc = "file containing ubiquitous kmer list. see FindBadGenomicKmersSpark to generate it.",
                fullName = "kmersToIgnore")
        public String kmersToIgnoreFile;

        /**
         * This is a path to a text file of contig names (one per line) that will be ignored when looking for inter-contig pairs.
         */
        @Argument(doc = "file containing alt contig names that will be ignored when looking for inter-contig pairs",
                fullName = "crossContigsToIgnore", optional = true)
        public String crossContigsToIgnoreFile;

        @VisibleForTesting
        public static class Locations {
            public final String metadataFile;
            public final String evidenceDir;
            public final String intervalFile;
            public final String qNamesMappedFile;
            public final String kmerFile;
            public final String qNamesAssemblyFile;
            public final String exclusionIntervalsFile;
            public final String crossContigsToIgnoreFile;
            public final String alignerIndexImageFile;

            public Locations( final String metadataFile, final String evidenceDir, final String intervalFile,
                              final String qNamesMappedFile, final String kmerFile, final String qNamesAssemblyFile,
                              final String exclusionIntervalsFile, final String crossContigsToIgnoreFile ,
                              final String alignerIndexImageFile ) {
                this.metadataFile = metadataFile;
                this.evidenceDir = evidenceDir;
                this.intervalFile = intervalFile;
                this.qNamesMappedFile = qNamesMappedFile;
                this.kmerFile = kmerFile;
                this.qNamesAssemblyFile = qNamesAssemblyFile;
                this.exclusionIntervalsFile = exclusionIntervalsFile;
                this.crossContigsToIgnoreFile = crossContigsToIgnoreFile;
                this.alignerIndexImageFile = alignerIndexImageFile;
            }
        }
    }

    public static class DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection implements Serializable {
        private static final long serialVersionUID = 1L;

        // todo: document this better
        // Currently the discovery stage requires a reference parameter in 2bit format (to broadcast) and
        // a reference in FASTA format (to get a good sequence dictionary for sorting variants).
        @Argument(doc = "FASTA formatted reference", shortName = "fastaReference",
                fullName = "fastaReference")
        public String fastaReference;

        @Argument(doc = "Minimum flanking alignment length", shortName = "minAlignLength",
                fullName = "minAlignLength", optional = true)
        public Integer minAlignLength = SVConstants.DiscoveryStepConstants.DEFAULT_MIN_ALIGNMENT_LENGTH;
    }

}
