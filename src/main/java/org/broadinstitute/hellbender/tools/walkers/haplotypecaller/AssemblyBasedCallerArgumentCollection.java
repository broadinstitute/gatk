package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAlignmentConstants;

/**
 * Set of arguments for Assembly Based Callers
 */
public abstract class AssemblyBasedCallerArgumentCollection {
    private static final long serialVersionUID = 1L;
    public static final String USE_FILTERED_READS_FOR_ANNOTATIONS_LONG_NAME = "use-filtered-reads-for-annotations";
    public static final String BAM_OUTPUT_LONG_NAME = "bam-output";
    public static final String BAM_OUTPUT_SHORT_NAME = "bamout";
    public static final String BAM_WRITER_TYPE_LONG_NAME = "bam-writer-type";
    public static final String DONT_USE_SOFT_CLIPPED_BASES_LONG_NAME = "dont-use-soft-clipped-bases";
    public static final String DO_NOT_RUN_PHYSICAL_PHASING_LONG_NAME = "do-not-run-physical-phasing";
    public static final String MAX_MNP_DISTANCE_LONG_NAME = "max-mnp-distance";
    public static final String MAX_MNP_DISTANCE_SHORT_NAME = "mnp-dist";

    public static final String MIN_BASE_QUALITY_SCORE_LONG_NAME = "min-base-quality-score";
    public static final String SMITH_WATERMAN_LONG_NAME = "smith-waterman";
    public static final String FORCE_CALL_ALLELES_LONG_NAME = "alleles";
    public static final String FORCE_CALL_FILTERED_ALLELES_LONG_NAME = "force-call-filtered-alleles";
    public static final String FORCE_CALL_FILTERED_ALLELES_SHORT_NAME = "genotype-filtered-alleles";
    public static final String EMIT_REF_CONFIDENCE_LONG_NAME = "emit-ref-confidence";
    public static final String EMIT_REF_CONFIDENCE_SHORT_NAME = "ERC";
    public static final String ALLELE_EXTENSION_LONG_NAME = "allele-informative-reads-overlap-margin";

    public static final String PILEUP_DETECTION_LONG_NAME = "pileup-detection";

    public static final String SMITH_WATERMAN_DANGLING_END_MATCH_VALUE_LONG_NAME = "smith-waterman-dangling-end-match-value";
    public static final String SMITH_WATERMAN_DANGLING_END_MISMATCH_PENALTY_LONG_NAME = "smith-waterman-dangling-end-mismatch-penalty";
    public static final String SMITH_WATERMAN_DANGLING_END_GAP_OPEN_PENALTY_LONG_NAME = "smith-waterman-dangling-end-gap-open-penalty";
    public static final String SMITH_WATERMAN_DANGLING_END_GAP_EXTEND_PENALTY_LONG_NAME = "smith-waterman-dangling-end-gap-extend-penalty";
    public static final String SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_MATCH_VALUE_LONG_NAME = "smith-waterman-haplotype-to-reference-match-value";
    public static final String SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_MISMATCH_PENALTY_LONG_NAME = "smith-waterman-haplotype-to-reference-mismatch-penalty";
    public static final String SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_GAP_OPEN_PENALTY_LONG_NAME = "smith-waterman-haplotype-to-reference-gap-open-penalty";
    public static final String SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_GAP_EXTEND_PENALTY_LONG_NAME = "smith-waterman-haplotype-to-reference-gap-extend-penalty";
    public static final String SMITH_WATERMAN_READ_TO_HAPLOTYPE_MATCH_VALUE_LONG_NAME = "smith-waterman-read-to-haplotype-match-value";
    public static final String SMITH_WATERMAN_READ_TO_HAPLOTYPE_MISMATCH_PENALTY_LONG_NAME = "smith-waterman-read-to-haplotype-mismatch-penalty";
    public static final String SMITH_WATERMAN_READ_TO_HAPLOTYPE_GAP_OPEN_PENALTY_LONG_NAME = "smith-waterman-read-to-haplotype-gap-open-penalty";
    public static final String SMITH_WATERMAN_READ_TO_HAPLOTYPE_GAP_EXTEND_PENALTY_LONG_NAME = "smith-waterman-read-to-haplotype-gap-extend-penalty";

    /**
     * See documentation at {@link SmithWatermanAlignmentConstants#STANDARD_NGS}.
     */
    private static final SWParameters DEFAULT_DANGLING_END_SMITH_WATERMAN_PARAMETERS = SmithWatermanAlignmentConstants.STANDARD_NGS;
    /**
     * See documentation at {@link SmithWatermanAlignmentConstants#NEW_SW_PARAMETERS}.
     */
    private static final SWParameters DEFAULT_HAPLOTYPE_TO_REFERENCE_SMITH_WATERMAN_PARAMETERS = SmithWatermanAlignmentConstants.NEW_SW_PARAMETERS;
    /**
     * See documentation at {@link SmithWatermanAlignmentConstants#ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS}.
     */
    private static final SWParameters DEFAULT_READ_TO_HAPLOTYPE_SMITH_WATERMAN_PARAMETERS = SmithWatermanAlignmentConstants.ALIGNMENT_TO_BEST_HAPLOTYPE_SW_PARAMETERS;

    public ReadThreadingAssembler createReadThreadingAssembler() {
        final ReadThreadingAssembler assemblyEngine = assemblerArgs.makeReadThreadingAssembler();
        assemblyEngine.setDebug(assemblerArgs.debugAssembly);
        assemblyEngine.setMinBaseQualityToUseInAssembly(minBaseQualityScore);

        return assemblyEngine;
    }

    protected abstract ReadThreadingAssemblerArgumentCollection getReadThreadingAssemblerArgumentCollection();

    @ArgumentCollection
    public ReadThreadingAssemblerArgumentCollection assemblerArgs = getReadThreadingAssemblerArgumentCollection();

    @ArgumentCollection
    public LikelihoodEngineArgumentCollection likelihoodArgs = new LikelihoodEngineArgumentCollection();

    @ArgumentCollection
    public PileupDetectionArgumentCollection pileupDetectionArgs = new PileupDetectionArgumentCollection();

    /**
     * The assembled haplotypes and locally realigned reads will be written as BAM to this file if requested.  Really
     * for debugging purposes only. Note that the output here does not include uninformative reads so that not every
     * input read is emitted to the bam.
     * <p>
     * Turning on this mode may result in serious performance cost for the caller.  It's really only appropriate to
     * use in specific areas where you want to better understand why the caller is making specific calls.
     * <p>
     * The reads are written out containing an "HC" tag (integer) that encodes which haplotype each read best matches
     * according to the haplotype caller's likelihood calculation.  The use of this tag is primarily intended
     * to allow good coloring of reads in IGV.  Simply go to "Color Alignments By > Tag" and enter "HC" to more
     * easily see which reads go with these haplotype.
     * <p>
     * Note that the haplotypes (called or all, depending on mode) are emitted as single reads covering the entire
     * active region, coming from sample "HC" and a special read group called "ArtificialHaplotype". This will increase the
     * pileup depth compared to what would be expected from the reads only, especially in complex regions.
     * <p>
     * Note also that only reads that are actually informative about the haplotypes are emitted.  By informative we mean
     * that there's a meaningful difference in the likelihood of the read coming from one haplotype compared to
     * its next best haplotype.
     * <p>
     * If multiple BAMs are passed as input to the tool (as is common for M2), then they will be combined in the bamout
     * output and tagged with the appropriate sample names.
     * <p>
     * The best way to visualize the output of this mode is with IGV.  Tell IGV to color the alignments by tag,
     * and give it the "HC" tag, so you can see which reads support each haplotype.  Finally, you can tell IGV
     * to group by sample, which will separate the potential haplotypes from the reads.  All of this can be seen in
     * <a href="https://www.dropbox.com/s/xvy7sbxpf13x5bp/haplotypecaller%20bamout%20for%20docs.png">this screenshot</a>
     */
    @Advanced
    @Argument(fullName = BAM_OUTPUT_LONG_NAME, shortName = BAM_OUTPUT_SHORT_NAME, doc = "File to which assembled haplotypes should be written", optional = true)
    public String bamOutputPath = null;

    /**
     * The type of BAM output we want to see. This determines whether HC will write out all of the haplotypes it
     * considered (top 128 max) or just the ones that were selected as alleles and assigned to samples.
     */
    @Advanced
    @Argument(fullName = BAM_WRITER_TYPE_LONG_NAME, doc = "Which haplotypes should be written to the BAM", optional = true)
    public HaplotypeBAMWriter.WriterType bamWriterType = HaplotypeBAMWriter.WriterType.CALLED_HAPLOTYPES;

    // -----------------------------------------------------------------------------------------------
    // arguments for debugging / developing
    // -----------------------------------------------------------------------------------------------

    @Advanced
    @Argument(fullName = DONT_USE_SOFT_CLIPPED_BASES_LONG_NAME, doc = "Do not analyze soft clipped bases in the reads", optional = true)
    public boolean dontUseSoftClippedBases = false;

    // Parameters to control read error correction

    /**
     * Bases with a quality below this threshold will not be used for calling.
     */
    @Argument(fullName = MIN_BASE_QUALITY_SCORE_LONG_NAME, shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", optional = true)
    public byte minBaseQualityScore = 10;

    //Annotations

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_LONG_NAME, doc = "Which Smith-Waterman implementation to use, generally FASTEST_AVAILABLE is the right choice", optional = true)
    public SmithWatermanAligner.Implementation smithWatermanImplementation = SmithWatermanAligner.Implementation.JAVA;

    /**
     * The reference confidence mode makes it possible to emit a per-bp or summarized confidence estimate for a site being strictly homozygous-reference.
     * See https://software.broadinstitute.org/gatk/documentation/article.php?id=4017 for information about GVCFs.
     * For Mutect2, this is a BETA feature that functions similarly to the HaplotypeCaller reference confidence/GVCF mode.
     */
    @Advanced
    @Argument(fullName = EMIT_REF_CONFIDENCE_LONG_NAME, shortName = EMIT_REF_CONFIDENCE_SHORT_NAME, doc = "Mode for emitting reference confidence scores (For Mutect2, this is a BETA feature)", optional = true)
    public ReferenceConfidenceMode emitReferenceConfidence = ReferenceConfidenceMode.NONE;

    protected abstract int getDefaultMaxMnpDistance();

    /**
     * Two or more phased substitutions separated by this distance or less are merged into MNPs.
     */
    @Advanced
    @Argument(fullName = MAX_MNP_DISTANCE_LONG_NAME, shortName = MAX_MNP_DISTANCE_SHORT_NAME,
            doc = "Two or more phased substitutions separated by this distance or less are merged into MNPs.", optional = true)
    public int maxMnpDistance = getDefaultMaxMnpDistance();

    @Argument(fullName = FORCE_CALL_ALLELES_LONG_NAME, doc = "The set of alleles to force-call regardless of evidence", optional = true)
    public FeatureInput<VariantContext> alleles;

    @Advanced
    @Argument(fullName = FORCE_CALL_FILTERED_ALLELES_LONG_NAME, shortName = FORCE_CALL_FILTERED_ALLELES_SHORT_NAME, doc = "Force-call filtered alleles included in the resource specified by --alleles", optional = true)
    public boolean forceCallFiltered = false;

    @Advanced
    @Argument(fullName = "soft-clip-low-quality-ends", doc = "If enabled will preserve low-quality read ends as softclips (used for DRAGEN-GATK BQD genotyper model)", optional = true)
    public boolean softClipLowQualityEnds = false;

    @Advanced
    @Argument(fullName = ALLELE_EXTENSION_LONG_NAME,
            doc = "Likelihood and read-based annotations will only take into consideration reads " +
                    "that overlap the variant or any base no further than this distance expressed in base pairs",
            optional = true)
    public int informativeReadOverlapMargin = 2;

    // -----------------------------------------------------------------------------------------------
    // Smith-Waterman parameters for dangling-end recovery
    // -----------------------------------------------------------------------------------------------

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_DANGLING_END_MATCH_VALUE_LONG_NAME,
            doc = "Smith-Waterman match value for dangling-end recovery.",
            minValue = 0,
            optional = true)
    public int smithWatermanDanglingEndMatchValue = DEFAULT_DANGLING_END_SMITH_WATERMAN_PARAMETERS.getMatchValue();

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_DANGLING_END_MISMATCH_PENALTY_LONG_NAME,
            doc = "Smith-Waterman mismatch penalty for dangling-end recovery.",
            maxValue = 0,
            optional = true)
    public int smithWatermanDanglingEndMismatchPenalty = DEFAULT_DANGLING_END_SMITH_WATERMAN_PARAMETERS.getMismatchPenalty();

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_DANGLING_END_GAP_OPEN_PENALTY_LONG_NAME,
            doc = "Smith-Waterman gap-open penalty for dangling-end recovery.",
            maxValue = 0,
            optional = true)
    public int smithWatermanDanglingEndGapOpenPenalty = DEFAULT_DANGLING_END_SMITH_WATERMAN_PARAMETERS.getGapOpenPenalty();

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_DANGLING_END_GAP_EXTEND_PENALTY_LONG_NAME,
            doc = "Smith-Waterman gap-extend penalty for dangling-end recovery.",
            maxValue = 0,
            optional = true)
    public int smithWatermanDanglingEndGapExtendPenalty = DEFAULT_DANGLING_END_SMITH_WATERMAN_PARAMETERS.getGapExtendPenalty();

    // -----------------------------------------------------------------------------------------------
    // Smith-Waterman parameters for haplotype-to-reference alignment
    // -----------------------------------------------------------------------------------------------

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_MATCH_VALUE_LONG_NAME,
            doc = "Smith-Waterman match value for haplotype-to-reference alignment.",
            minValue = 0,
            optional = true)
    public int smithWatermanHaplotypeToReferenceMatchValue = DEFAULT_HAPLOTYPE_TO_REFERENCE_SMITH_WATERMAN_PARAMETERS.getMatchValue();

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_MISMATCH_PENALTY_LONG_NAME,
            doc = "Smith-Waterman mismatch penalty for haplotype-to-reference alignment.",
            maxValue = 0,
            optional = true)
    public int smithWatermanHaplotypeToReferenceMismatchPenalty = DEFAULT_HAPLOTYPE_TO_REFERENCE_SMITH_WATERMAN_PARAMETERS.getMismatchPenalty();

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_GAP_OPEN_PENALTY_LONG_NAME,
            doc = "Smith-Waterman gap-open penalty for haplotype-to-reference alignment.",
            maxValue = 0,
            optional = true)
    public int smithWatermanHaplotypeToReferenceGapOpenPenalty = DEFAULT_HAPLOTYPE_TO_REFERENCE_SMITH_WATERMAN_PARAMETERS.getGapOpenPenalty();

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_HAPLOTYPE_TO_REFERENCE_GAP_EXTEND_PENALTY_LONG_NAME,
            doc = "Smith-Waterman gap-extend penalty for haplotype-to-reference alignment.",
            maxValue = 0,
            optional = true)
    public int smithWatermanHaplotypeToReferenceGapExtendPenalty = DEFAULT_HAPLOTYPE_TO_REFERENCE_SMITH_WATERMAN_PARAMETERS.getGapExtendPenalty();

    // -----------------------------------------------------------------------------------------------
    // Smith-Waterman parameters for read-to-haplotype alignment
    // -----------------------------------------------------------------------------------------------

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_READ_TO_HAPLOTYPE_MATCH_VALUE_LONG_NAME,
            doc = "Smith-Waterman match value for read-to-haplotype alignment.",
            minValue = 0,
            optional = true)
    public int smithWatermanReadToHaplotypeMatchValue = DEFAULT_READ_TO_HAPLOTYPE_SMITH_WATERMAN_PARAMETERS.getMatchValue();

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_READ_TO_HAPLOTYPE_MISMATCH_PENALTY_LONG_NAME,
            doc = "Smith-Waterman mismatch penalty for read-to-haplotype alignment.",
            maxValue = 0,
            optional = true)
    public int smithWatermanReadToHaplotypeMismatchPenalty = DEFAULT_READ_TO_HAPLOTYPE_SMITH_WATERMAN_PARAMETERS.getMismatchPenalty();

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_READ_TO_HAPLOTYPE_GAP_OPEN_PENALTY_LONG_NAME,
            doc = "Smith-Waterman gap-open penalty for read-to-haplotype alignment.",
            maxValue = 0,
            optional = true)
    public int smithWatermanReadToHaplotypeGapOpenPenalty = DEFAULT_READ_TO_HAPLOTYPE_SMITH_WATERMAN_PARAMETERS.getGapOpenPenalty();

    @Advanced
    @Argument(fullName = SMITH_WATERMAN_READ_TO_HAPLOTYPE_GAP_EXTEND_PENALTY_LONG_NAME,
            doc = "Smith-Waterman gap-extend penalty for read-to-haplotype alignment.",
            maxValue = 0,
            optional = true)
    public int smithWatermanReadToHaplotypeGapExtendPenalty = DEFAULT_READ_TO_HAPLOTYPE_SMITH_WATERMAN_PARAMETERS.getGapExtendPenalty();

    public SWParameters getDanglingEndSWParameters() {
        return new SWParameters(
                smithWatermanDanglingEndMatchValue,
                smithWatermanDanglingEndMismatchPenalty,
                smithWatermanDanglingEndGapOpenPenalty,
                smithWatermanDanglingEndGapExtendPenalty);
    }

    public SWParameters getHaplotypeToReferenceSWParameters() {
        return new SWParameters(
                smithWatermanHaplotypeToReferenceMatchValue,
                smithWatermanHaplotypeToReferenceMismatchPenalty,
                smithWatermanHaplotypeToReferenceGapOpenPenalty,
                smithWatermanHaplotypeToReferenceGapExtendPenalty);
    }

    public SWParameters getReadToHaplotypeSWParameters() {
        return new SWParameters(
                smithWatermanReadToHaplotypeMatchValue,
                smithWatermanReadToHaplotypeMismatchPenalty,
                smithWatermanReadToHaplotypeGapOpenPenalty,
                smithWatermanReadToHaplotypeGapExtendPenalty);
    }
}