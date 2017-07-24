package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.argparser.Argument;

import java.io.Serializable;


public class StructuralVariationDiscoveryArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static class FindBreakpointEvidenceSparkArgumentCollection implements Serializable {
        private static final long serialVersionUID = 1L;

        public static final int KMER_SIZE = 51;
        public static final int MAX_DUST_SCORE = KMER_SIZE - 2;

        //--------- parameters ----------

        @Argument(doc = "Kmer size.", fullName = "kSize")
        public int kSize = KMER_SIZE;

        @Argument(doc = "maximum kmer DUST score", fullName = "kmerMaxDUSTScore")
        public int maxDUSTScore = MAX_DUST_SCORE;

        @Argument(doc = "The minimum mapping quality for reads used to gather evidence of breakpoints.",
                fullName = "minEvidenceMapQ", optional = true)
        public int minEvidenceMapQ = 20;

        @Argument(doc = "The minimum length of the matched portion of an interesting alignment.  "+
                "Reads that don't match at least this many reference bases won't be used in gathering evidence.",
                fullName = "minEvidenceMatchLength", optional = true)
        public int minEvidenceMatchLength = 45;

        @Argument(doc = "Proper pairs have the positive strand read upstream of the negative strand read, but "+
                "we allow this much slop for short fragments.",
                fullName = "allowedShortFragmentOverhang", optional = true)
        public int allowedShortFragmentOverhang = 10;

        @Argument(doc = "Largest fragment size that will be explicitly counted in determining " +
                "fragment size statistics.", fullName = "maxTrackedFragmentLength", optional = true)
        public int maxTrackedFragmentLength = 2000;

        @Argument(doc = "Intervals with more than this much coverage are filtered out, because the reads mapped to "+
                "that interval are clearly not exclusively local to the interval.", fullName = "maxIntervalCoverage")
        public int maxIntervalCoverage = 1000;

        @Argument(doc = "Minimum weight of the corroborating read evidence to validate some single piece of evidence.",
                fullName = "minEvidenceCount")
        public int minEvidenceWeight = 15;

        @Argument(doc = "Minimum weight of the evidence that shares a distal target locus to validate the evidence.",
                fullName = "minCoherentEvidenceCount")
        public int minCoherentEvidenceWeight = 7;

        @Argument(doc = "Minimum number of localizing kmers in a valid interval.", fullName="minKmersPerInterval")
        public int minKmersPerInterval = 20;

        @Argument(doc = "KmerCleaner maximum number of intervals for a localizing kmer.", fullName = "cleanerMaxIntervals")
        public int cleanerMaxIntervals = 3;

        @Argument(doc = "KmerCleaner minimum kmer count.", fullName = "cleanerMinKmerCount")
        public int cleanerMinKmerCount = 3;

        @Argument(doc = "KmerCleaner maximum kmer count.", fullName = "cleanerMaxKmerCount")
        public int cleanerMaxKmerCount = 125;

        @Argument(doc = "KmerCleaner unique error-free kmers per partition", fullName = "cleanerKmersPerPartitionGuess")
        public int cleanerKmersPerPartitionGuess = 600000;

        @Argument(doc = "Maximum number of templates containing an assembly kmer.", fullName = "maxQNamesPerKmer")
        public int maxQNamesPerKmer = 500;

        @Argument(doc = "Guess at number of clean kmers per assembly partition.", fullName = "assemblyKmerMapSize")
        public int assemblyKmerMapSize = 250000;

        @Argument(doc = "Guess at the ratio of reads in the final assembly to the number reads mapped to the interval.",
                fullName = "assemblyToMappedSizeRatioGuess")
        public int assemblyToMappedSizeRatioGuess = 7;

        @Argument(doc = "Maximum FASTQ file size.", fullName = "maxFASTQSize")
        public int maxFASTQSize = 3000000;

        @Argument(doc = "Exclusion interval padding.", fullName = "exclusionIntervalPadding")
        public int exclusionIntervalPadding = 0;

        @Argument(doc = "Include read mapping location in FASTQ files.", fullName = "includeMappingLocation")
        public boolean includeMappingLocation = true;

        @Argument(doc = "Don't look for extra reads mapped outside the interval.", fullName = "intervalOnlyAssembly")
        public boolean intervalOnlyAssembly = false;

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

        private static final String OUTPUT_ORDER_SHORT_NAME = "sort";
        private static final String OUTPUT_ORDER_FULL_NAME = "assembliesSortOrder";

        @Argument(doc = "sorting order to be used for the output assembly alignments SAM/BAM file",
                shortName = OUTPUT_ORDER_SHORT_NAME,
                fullName = OUTPUT_ORDER_FULL_NAME,
                optional = true)
        public SAMFileHeader.SortOrder assembliesSortOrder = SAMFileHeader.SortOrder.coordinate;
    }

    public static class DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection implements Serializable {
        private static final long serialVersionUID = 1L;

        public static final int GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY = 50; // alignment with gap of size >= 50 will be broken apart.
        public static final int CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD = 60;
        public static final int MISSING_NM = Integer.MIN_VALUE;
        public static final int ARTIFICIAL_MISMATCH = MISSING_NM;
        public static final int DEFAULT_MIN_ALIGNMENT_LENGTH = 50; // Minimum flanking alignment length filters used when going through contig alignments.

        // todo: document this better
        // Currently the discovery stage requires a reference parameter in 2bit format (to broadcast) and
        // a reference in FASTA format (to get a good sequence dictionary for sorting variants).
        @Argument(doc = "FASTA formatted reference", shortName = "fastaReference",
                fullName = "fastaReference")
        public String fastaReference;

        @Argument(doc = "Minimum flanking alignment length", shortName = "minAlignLength",
                fullName = "minAlignLength", optional = true)
        public Integer minAlignLength = DEFAULT_MIN_ALIGNMENT_LENGTH;
    }

}
