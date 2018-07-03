package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.io.Serializable;


public class StructuralVariationDiscoveryArgumentCollection implements Serializable {
    private static final long serialVersionUID = 1L;

    public static final int STRUCTURAL_VARIANT_SIZE_LOWER_BOUND = 50;

    public static class FindBreakpointEvidenceSparkArgumentCollection implements Serializable {
        private static final long serialVersionUID = 1L;

        public static final int KMER_SIZE = 51;
        public static final int MAX_DUST_SCORE = KMER_SIZE - 2;

        public static final double TRAINING_SET_MEAN_COVERAGE = 42.855164; // mean overlap of data set used to tune evidence filtering
        public static final int TRAINING_SET_OPTIMAL_MIN_OVERLAP = 15; // optimal min overlap of data set used to tune evidence
        public static final int TRAINING_SET_OPTIMAL_MIN_COHERENCE = 7; // optimal min coherence of data set used to tune evidence

        //--------- parameters ----------

        @Argument(doc = "Kmer size.", fullName = "k-size")
        public int kSize = KMER_SIZE;

        @Argument(doc = "Maximum kmer DUST score.", fullName = "kmer-max-dust-score")
        public int maxDUSTScore = MAX_DUST_SCORE;

        @Argument(doc = "The minimum mapping quality for reads used to gather evidence of breakpoints.",
                fullName = "min-evidence-mapq")
        public int minEvidenceMapQ = 20;

        @Argument(doc = "The minimum length of the matched portion of an interesting alignment.  "+
                "Reads that don't match at least this many reference bases won't be used in gathering evidence.",
                fullName = "min-evidence-match-length")
        public int minEvidenceMatchLength = 45;

        @Argument(doc = "Proper pairs have the positive strand read upstream of the negative strand read, but "+
                "we allow this much slop for short fragments.",
                fullName = "allowed-short-fragment-overhang")
        public int allowedShortFragmentOverhang = 10;

        @Argument(doc = "Largest fragment size that will be explicitly counted in determining " +
                "fragment size statistics.", fullName = "max-tracked-fragment-length")
        public int maxTrackedFragmentLength = 2000;

        @Argument(doc = "We filter out contiguous regions of the genome that have coverage of at least high-depth-coverage-factor * avg-coverage and a " +
                "peak coverage of high-depth-coverage-peak-factor * avg-coverage, because the reads mapped to those regions tend to be non-local and high depth prevents accurate assembly.",
                fullName = "high-depth-coverage-peak-factor")
        public int highDepthCoveragePeakFactor = 7;

        @Argument(doc = "We filter out contiguous regions of the genome that have coverage of at least high-depth-coverage-factor * avg-coverage and a " +
                "peak coverage of high-depth-coverage-peak-factor * avg-coverage, because the reads mapped to those regions tend to be non-local and high depth prevents accurate assembly.",
                fullName = "high-depth-coverage-factor")
        public int highDepthCoverageFactor = 3;

        @Argument(doc = "Minimum weight of the corroborating read evidence to validate some single piece of evidence, as a ratio of the mean coverage in the BAM. "
                + "The default value is overlap-count / mean coverage ~ 15 / 42.9 ~ 0.350",
                fullName = "min-evidence-coverage-ratio")
        public double minEvidenceWeightPerCoverage = TRAINING_SET_OPTIMAL_MIN_OVERLAP / TRAINING_SET_MEAN_COVERAGE;

        @Argument(doc = "Minimum weight of the evidence that shares a distal target locus to validate the evidence, as a ratio of the mean coverage in the BAM. "
                + "The default value is coherent-count / mean coverage ~ 7 / 42.9 ~ 0.163",
                fullName = "min-coherent-evidence-coverage-ratio")
        public double minCoherentEvidenceWeightPerCoverage = TRAINING_SET_OPTIMAL_MIN_COHERENCE / TRAINING_SET_MEAN_COVERAGE;

        @Argument(doc = "Minimum number of localizing kmers in a valid interval.", fullName="min-kmers-per-interval")
        public int minKmersPerInterval = 5;

        @Argument(doc = "KmerCleaner maximum number of intervals for a localizing kmer."+
                " If a kmer occurs in too many intervals, it isn't sufficiently local.",
                fullName = "cleaner-max-intervals")
        public int cleanerMaxIntervals = 3;

        @Argument(doc = "KmerCleaner minimum kmer count for a localizing kmer."+
                "  If we see it less often than this many times, we're guessing it's erroneous.",
                fullName = "cleaner-min-kmer-count")
        public int cleanerMinKmerCount = 4;

        @Argument(doc = "KmerCleaner maximum copy number (not count, but copy number) for a kmer."+
                " Kmers observed too frequently are probably mismapped or ubiquitous.",
                fullName = "cleaner-max-copy-number")
        public int cleanerMaxCopyNumber = 4;

        @Argument(doc = "Guess at the ratio of reads in the final assembly to the number reads mapped to the interval.",
                fullName = "assembly-to-mapped-size-ratio-guess")
        public int assemblyToMappedSizeRatioGuess = 7;

        @Argument(doc = "Maximum total bases in FASTQs that can be assembled.", fullName = "max-fastq-size")
        public int maxFASTQSize = 3000000;

        @Argument(doc = "Exclusion interval padding.", fullName = "exclusion-interval-padding")
        public int exclusionIntervalPadding = 0;

        @Argument(doc = "Include read mapping location in FASTQ files.", fullName = "include-mapping-location")
        public boolean includeMappingLocation = true;

        @Argument(doc = "Don't look for extra reads mapped outside the interval.", fullName = "interval-only-assembly")
        public boolean intervalOnlyAssembly = false;

        @Argument(doc = "Weight to give external evidence.", fullName = "external-evidence-weight")
        public int externalEvidenceWeight = 10;

        @Argument(doc = "Uncertainty in location of external evidence.", fullName = "external-evidence-uncertainty")
        public int externalEvidenceUncertainty = 150;

        @Argument(doc = "Adapter sequence.", fullName = "adapter-sequence", optional = true)
        public String adapterSequence;

        @Argument(doc = "Minimum classified probability for a piece of evidence to pass xgboost evidence filter",
                fullName = "sv-evidence-filter-threshold-probability")
        public double svEvidenceFilterThresholdProbability = 0.92;

        @Argument(doc = "Filter method for selecting evidence to group into Assembly Intervals",
                fullName = "sv-evidence-filter-type")
        public SvEvidenceFilterType svEvidenceFilterType = SvEvidenceFilterType.DENSITY;

        // ---------- options -----------

        @Argument(doc = "Write GFA representation of assemblies in fastq-dir.", fullName = "write-gfas")
        public boolean writeGFAs = false;

        @Advanced
        @Argument(doc = "Aggressively simplify local assemblies, ignoring small variants.", fullName = "pop-variant-bubbles")
        public boolean popVariantBubbles = false;

        @Advanced
        @Argument(doc = "Simplify local assemblies by removing contigs shadowed by similar contigs.", fullName = "remove-shadowed-contigs")
        public boolean removeShadowedContigs = true;

        @Advanced
        @Argument(doc = "Traverse assembly graph and produce contigs for all paths.", fullName = "expand-assembly-graph")
        public boolean expandAssemblyGraph = true;

        @Advanced @Argument(doc = "ZDropoff (see Bwa mem manual) for contig alignment.", fullName = "z-dropoff")
        public int zDropoff = 20;

        @Argument(doc = "Allow evidence filter to run without gaps annotation (assume no gaps).", fullName = "run-without-gaps-annotation")
        public boolean runWithoutGapsAnnotation = false;
        @Argument(doc = "Allow evidence filter to run without annotation for single-read mappability of 100-mers (assume all mappable).",
                fullName = "run-without-umap-s100-annotation")
        public boolean runWithoutUmapS100Annotation = false;


        // --------- locations ----------

        @Argument(doc = "bwa-mem index image file", fullName = "aligner-index-image")
        public String alignerIndexImageFile;

        @Argument(doc = "external evidence input file", fullName = "external-evidence", optional = true)
        public String externalEvidenceFile;

        @Argument(doc = "output file for read metadata", fullName = "read-metadata", optional = true)
        public String metadataFile;

        @Argument(doc = "directory for evidence output", fullName = "breakpoint-evidence-dir", optional = true)
        public String evidenceDir;

        @Argument(doc = "directory for evidence output", fullName = "unfiltered-breakpoint-evidence-dir", optional = true)
        public String unfilteredEvidenceDir;

        @Argument(doc = "file for breakpoint intervals output", fullName = "breakpoint-intervals", optional = true)
        public String intervalFile;

        @Argument(doc = "file for high-coverage intervals output", fullName = "high-coverage-intervals", optional = true)
        public String highCoverageIntervalsFile;

        @Argument(doc = "file for mapped qname intervals output", fullName = "qname-intervals-mapped", optional = true)
        public String qNamesMappedFile;

        @Argument(doc = "file for kmer intervals output", fullName = "kmer-intervals", optional = true)
        public String kmerFile;

        @Argument(doc = "file for mapped qname intervals output",
                fullName = "qname-intervals-for-assembly", optional = true)
        public String qNamesAssemblyFile;

        @Argument(doc = "output dir for assembled fastqs", fullName = "fastq-dir", optional = true)
        public String fastqDir;

        @Argument(doc = "output file for non-assembled breakpoints in bedpe format",
                fullName = "target-link-file", optional = true)
        public String targetLinkFile;

        /**
         * This is a file that calls out the coordinates of intervals in the reference assembly to exclude from
         * consideration when calling putative breakpoints.
         * Each line is a tab-delimited interval with 1-based inclusive coordinates like this:
         *  chr1	124535434	142535434
         */
        @Argument(doc = "file of reference intervals to exclude", fullName = "exclusion-intervals", optional = true)
        public String exclusionIntervalsFile;

        /**
         * This is a path to a file of kmers that appear too frequently in the reference to be usable as probes to localize
         * reads.  We don't calculate it here, because it depends only on the reference.
         * The program FindBadGenomicKmersSpark can produce such a list for you.
         */
        @Argument(doc = "file containing ubiquitous kmer list. see FindBadGenomicKmersSpark to generate it.",
                fullName = "kmers-to-ignore")
        public String kmersToIgnoreFile;

        /**
         * This is a path to a text file of contig names (one per line) that will be ignored when looking for inter-contig pairs.
         */
        @Argument(doc = "file containing alt contig names that will be ignored when looking for inter-contig pairs",
                fullName = "cross-contigs-to-ignore", optional = true)
        public String crossContigsToIgnoreFile;

        private static final String OUTPUT_ORDER_SHORT_NAME = "sort";
        private static final String OUTPUT_ORDER_FULL_NAME = "assembled-contigs-output-order";

        @Argument(doc = "sorting order to be used for the output assembly alignments SAM/BAM file (currently only coordinate or query name is supported)",
                shortName = OUTPUT_ORDER_SHORT_NAME,
                fullName = OUTPUT_ORDER_FULL_NAME,
                optional = true)
        public SAMFileHeader.SortOrder assembliesSortOrder = SAMFileHeader.SortOrder.coordinate;

        @Argument(doc = "Path to xgboost classifier model file for evidence filtering",
                fullName = "sv-evidence-filter-model-file", optional=true)
        public String svEvidenceFilterModelFile = null;

        @Argument(doc = "Path to single read 100-mer mappability file in the reference genome, used by classifier to score evidence for filtering. "
                + "To use classifier without specifying mappability file, pass the flag --run-without-umap-s100-annotation",
                fullName = "sv-genome-umap-s100-file", optional = true)
        public String svGenomeUmapS100File = null;

        @Argument(doc = "Path to file enumerating gaps in the reference genome, used by classifier to score evidence for filtering. "
                + "To use classifier without specifying gaps file, pass the flag --run-without-gaps-annotation",
                fullName = "sv-genome-gaps-file", optional = true)
        public String svGenomeGapsFile = null;

        /**
         * Explicit call this method.
         */
        public void validate() {
            ParamUtils.isPositive(kSize, "invalid value provided to kSize: " + kSize);
            ParamUtils.isPositive(maxDUSTScore, "invalid value provided to maxDUSTScore: " + maxDUSTScore);
            ParamUtils.isPositiveOrZero(minEvidenceMapQ, "invalid value provided to minEvidenceMapQ: " + minEvidenceMapQ);
            ParamUtils.isPositive(minEvidenceMatchLength, "invalid value provided to minEvidenceMatchLength: " + minEvidenceMatchLength);
            ParamUtils.isPositiveOrZero(allowedShortFragmentOverhang, "invalid value provided to allowedShortFragmentOverhang: " + allowedShortFragmentOverhang);
            ParamUtils.isPositive(maxTrackedFragmentLength, "invalid value provided to maxTrackedFragmentLength: " + maxTrackedFragmentLength);
            ParamUtils.isPositive(highDepthCoveragePeakFactor, "invalid value provided to highDepthCoveragePeakFactor: " + highDepthCoveragePeakFactor);
            ParamUtils.isPositive(minEvidenceWeightPerCoverage, "invalid value provided to minEvidenceWeightPerCoverage: " + minEvidenceWeightPerCoverage);
            ParamUtils.isPositive(minCoherentEvidenceWeightPerCoverage, "invalid value provided to minCoherentEvidenceWeightPerCoverage: " + minCoherentEvidenceWeightPerCoverage);
            ParamUtils.isPositive(minKmersPerInterval, "invalid value provided to minKmersPerInterval: " + minKmersPerInterval);
            ParamUtils.isPositive(cleanerMaxIntervals, "invalid value provided to cleanerMaxIntervals: " + cleanerMaxIntervals);
            ParamUtils.isPositive(cleanerMinKmerCount, "invalid value provided to cleanerMinKmerCount: " + cleanerMinKmerCount);
            ParamUtils.isPositive(cleanerMaxCopyNumber, "invalid value provided to cleanerMaxCopyNumber: " + cleanerMaxCopyNumber);
            ParamUtils.isPositive(assemblyToMappedSizeRatioGuess, "invalid value provided to assemblyToMappedSizeRatioGuess: " + assemblyToMappedSizeRatioGuess);
            ParamUtils.isPositive(maxFASTQSize, "invalid value provided to maxFASTQSize: " + maxFASTQSize);
            ParamUtils.isPositiveOrZero(exclusionIntervalPadding, "invalid value provided to exclusionIntervalPadding: " + exclusionIntervalPadding);
            ParamUtils.isPositive(externalEvidenceWeight, "invalid value provided to externalEvidenceWeight: " + externalEvidenceWeight);
            ParamUtils.isPositive(externalEvidenceUncertainty, "invalid value provided to externalEvidenceUncertainty: " + externalEvidenceUncertainty);
            if ( !(assembliesSortOrder.equals(SAMFileHeader.SortOrder.coordinate)
                    ||
                    assembliesSortOrder.equals(SAMFileHeader.SortOrder.queryname)) )
                throw new UserException("We currently support only coordinate or query name sort order for assembly alignment SAM output. " +
                        "User provided sort order: " + assembliesSortOrder);
        }
    }

    public enum SvEvidenceFilterType {DENSITY, XGBOOST}

    public static class DiscoverVariantsFromContigAlignmentsSparkArgumentCollection implements Serializable {
        private static final long serialVersionUID = 1L;

        public static final int GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY = STRUCTURAL_VARIANT_SIZE_LOWER_BOUND; // alignment with gap of size >= 50 will be broken apart.

        // TODO: 7/30/18 the following two values essentially perform the same filtering, except one (CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD) is used in a tool (ContigChimericAlignmentIterativeInterpreter) that is about to be phased out, so move it when the kill switch is flipped
        public static final int CHIMERIC_ALIGNMENTS_HIGHMQ_THRESHOLD = 60;
        public static final int ASSEMBLY_ALIGNMENT_MQ_FILTER_DEFAULT = 30;

        public static final int DEFAULT_MIN_ALIGNMENT_LENGTH = 50; // Minimum flanking alignment length filters used when going through contig alignments.
        public static final int DEFAULT_ASSEMBLED_IMPRECISE_EVIDENCE_OVERLAP_UNCERTAINTY = 100;
        public static final int DEFAULT_IMPRECISE_VARIANT_EVIDENCE_THRESHOLD = 7;
        public static final int DEFAULT_TRUTH_INTERVAL_PADDING = 50;
        public static final int DEFAULT_MAX_CALLABLE_IMPRECISE_DELETION_SIZE = 15000;

        @Argument(doc = "Minimum flanking alignment length", fullName = "min-align-length")
        public Integer minAlignLength = DEFAULT_MIN_ALIGNMENT_LENGTH;

        @Argument(doc = "Minimum mapping quality of evidence assembly contig", shortName = "mq", fullName = "min-mq")
        public Integer minMQ = ASSEMBLY_ALIGNMENT_MQ_FILTER_DEFAULT;

        @Hidden
        @Argument(doc = "VCF containing the true breakpoints used only for evaluation (not generation) of calls",
                fullName = "truth-vcf", optional = true)
        public String truthVCF;

        @Argument(doc = "Uncertainty in overlap of assembled breakpoints and evidence target links.",
                fullName = "assembly-imprecise-evidence-overlap-uncertainty")
        public int assemblyImpreciseEvidenceOverlapUncertainty = DEFAULT_ASSEMBLED_IMPRECISE_EVIDENCE_OVERLAP_UNCERTAINTY;

        @Argument(doc = "Number of pieces of imprecise evidence necessary to call a variant in the absence of an assembled breakpoint.",
                fullName = "imprecise-variant-evidence-threshold")
        public int impreciseVariantEvidenceThreshold = DEFAULT_IMPRECISE_VARIANT_EVIDENCE_THRESHOLD;

        @Argument(doc = "External CNV calls file. Should be single sample VCF, and contain only confident autosomal non-reference CNV calls (for now).",
                fullName = "cnv-calls", optional = true)
        public String cnvCallsFile;

        @Argument(doc = "Breakpoint padding for evaluation against truth data.",
                fullName = "truth-interval-padding", optional = true)
        public int truthIntervalPadding = DEFAULT_TRUTH_INTERVAL_PADDING;

        @Argument(doc = "Maximum size deletion to call based on imprecise evidence without corroborating read depth evidence",
                fullName = "max-callable-imprecise-deletion-size", optional=true)
        public int maxCallableImpreciseVariantDeletionSize = DEFAULT_MAX_CALLABLE_IMPRECISE_DELETION_SIZE;

        @Advanced
        @Argument(doc = "Run interpretation tool in debug mode (more information print to screen)", fullName = "debug-mode", optional = true)
        public Boolean runInDebugMode = false;

        /**
         * Explicit call this method.
         */
        public void validate() {
            ParamUtils.isPositive(minAlignLength, "invalid value provided to minAlignLength: " + minAlignLength);
            ParamUtils.isPositive(assemblyImpreciseEvidenceOverlapUncertainty, "invalid value provided to assemblyImpreciseEvidenceOverlapUncertainty: " + assemblyImpreciseEvidenceOverlapUncertainty);
            ParamUtils.isPositive(impreciseVariantEvidenceThreshold, "invalid value provided to impreciseVariantEvidenceThreshold: " + impreciseVariantEvidenceThreshold);
            ParamUtils.isPositive(truthIntervalPadding, "invalid value provided to truthIntervalPadding: " + truthIntervalPadding);
            ParamUtils.isPositive(maxCallableImpreciseVariantDeletionSize, "invalid value provided to maxCallableImpreciseVariantDeletionSize: " + maxCallableImpreciseVariantDeletionSize);
        }
    }
}
