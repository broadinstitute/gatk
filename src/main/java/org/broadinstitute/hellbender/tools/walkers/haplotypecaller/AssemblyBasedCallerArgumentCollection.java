package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.tools.walkers.genotyper.StandardCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;
import org.broadinstitute.hellbender.utils.smithwaterman.SmithWatermanAligner;

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

    /**
     * The assembled haplotypes and locally realigned reads will be written as BAM to this file if requested.  Really
     * for debugging purposes only. Note that the output here does not include uninformative reads so that not every
     * input read is emitted to the bam.
     *
     * Turning on this mode may result in serious performance cost for the caller.  It's really only appropriate to
     * use in specific areas where you want to better understand why the caller is making specific calls.
     *
     * The reads are written out containing an "HC" tag (integer) that encodes which haplotype each read best matches
     * according to the haplotype caller's likelihood calculation.  The use of this tag is primarily intended
     * to allow good coloring of reads in IGV.  Simply go to "Color Alignments By > Tag" and enter "HC" to more
     * easily see which reads go with these haplotype.
     *
     * Note that the haplotypes (called or all, depending on mode) are emitted as single reads covering the entire
     * active region, coming from sample "HC" and a special read group called "ArtificialHaplotype". This will increase the
     * pileup depth compared to what would be expected from the reads only, especially in complex regions.
     *
     * Note also that only reads that are actually informative about the haplotypes are emitted.  By informative we mean
     * that there's a meaningful difference in the likelihood of the read coming from one haplotype compared to
     * its next best haplotype.
     *
     * If multiple BAMs are passed as input to the tool (as is common for M2), then they will be combined in the bamout
     * output and tagged with the appropriate sample names.
     *
     * The best way to visualize the output of this mode is with IGV.  Tell IGV to color the alignments by tag,
     * and give it the "HC" tag, so you can see which reads support each haplotype.  Finally, you can tell IGV
     * to group by sample, which will separate the potential haplotypes from the reads.  All of this can be seen in
     * <a href="https://www.dropbox.com/s/xvy7sbxpf13x5bp/haplotypecaller%20bamout%20for%20docs.png">this screenshot</a>
     *
     */
    @Advanced
    @Argument(fullName= BAM_OUTPUT_LONG_NAME, shortName= BAM_OUTPUT_SHORT_NAME, doc="File to which assembled haplotypes should be written", optional = true)
    public String bamOutputPath = null;

    /**
     * The type of BAM output we want to see. This determines whether HC will write out all of the haplotypes it
     * considered (top 128 max) or just the ones that were selected as alleles and assigned to samples.
     */
    @Advanced
    @Argument(fullName= BAM_WRITER_TYPE_LONG_NAME, doc="Which haplotypes should be written to the BAM", optional = true)
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
     * (BETA feature) The reference confidence mode makes it possible to emit a per-bp or summarized confidence estimate for a site being strictly homozygous-reference.
     * This is similar to the HaplotypeCaller reference confidence/GVCF mode. See https://software.broadinstitute.org/gatk/documentation/article.php?id=4017 for information about GVCFs.
     */
    @Advanced
    @Argument(fullName="emit-ref-confidence", shortName="ERC", doc="(BETA feature) Mode for emitting reference confidence scores", optional = true)
    public ReferenceConfidenceMode emitReferenceConfidence = ReferenceConfidenceMode.NONE;

    protected abstract int getDefaultMaxMnpDistance();

    /**
     * Two or more phased substitutions separated by this distance or less are merged into MNPs.
     */
    @Advanced
    @Argument(fullName = MAX_MNP_DISTANCE_LONG_NAME, shortName = MAX_MNP_DISTANCE_SHORT_NAME,
            doc = "Two or more phased substitutions separated by this distance or less are merged into MNPs.", optional = true)
    public int maxMnpDistance = getDefaultMaxMnpDistance();
}
