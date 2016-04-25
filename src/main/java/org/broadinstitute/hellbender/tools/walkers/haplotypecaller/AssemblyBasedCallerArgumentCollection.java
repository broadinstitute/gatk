package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.tools.walkers.genotyper.StandardCallerArgumentCollection;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;

import java.nio.file.Path;

/**
 * Set of arguments for Assembly Based Callers
 */
public abstract class AssemblyBasedCallerArgumentCollection extends StandardCallerArgumentCollection {
    private static final long serialVersionUID = 1L;

    @Advanced
    @Argument(fullName="debug", shortName="debug", doc="Print out very verbose debug information about each triggering active region", optional = true)
    public boolean DEBUG;

    @Advanced
    @Argument(fullName="useFilteredReadsForAnnotations", shortName="useFilteredReadsForAnnotations", doc = "Use the contamination-filtered read maps for the purposes of annotating variants", optional=true)
    public boolean USE_FILTERED_READ_MAP_FOR_ANNOTATIONS = false;

    /**
     * The reference confidence mode makes it possible to emit a per-bp or summarized confidence estimate for a site being strictly homozygous-reference.
     * See http://www.broadinstitute.org/gatk/guide/article?id=2940 for more details of how this works.
     * Note that if you set -ERC GVCF, you also need to set -variant_index_type LINEAR and -variant_index_parameter 128000 (with those exact values!).
     * This requirement is a temporary workaround for an issue with index compression.
     */
    @Advanced
    @Argument(fullName="emitRefConfidence", shortName="ERC", doc="Mode for emitting reference confidence scores", optional = true)
    protected ReferenceConfidenceMode emitReferenceConfidence = ReferenceConfidenceMode.NONE;

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
    @Argument(fullName="bamOutput", shortName="bamout", doc="File to which assembled haplotypes should be written", optional = true)
    public String bamWriter = null;

    /**
     * The type of BAM output we want to see. This determines whether HC will write out all of the haplotypes it
     * considered (top 128 max) or just the ones that were selected as alleles and assigned to samples.
     */
    @Advanced
    @Argument(fullName="bamWriterType", shortName="bamWriterType", doc="Which haplotypes should be written to the BAM", optional = true)
    public HaplotypeBAMWriter.WriterType bamWriterType = HaplotypeBAMWriter.WriterType.CALLED_HAPLOTYPES;

    /**
     * If set, certain "early exit" optimizations in HaplotypeCaller, which aim to save compute and time by skipping
     * calculations if an ActiveRegion is determined to contain no variants, will be disabled. This is most likely to be useful if
     * you're using the -bamout argument to examine the placement of reads following reassembly and are interested in seeing the mapping of
     * reads in regions with no variations. Setting the -forceActive and -dontTrimActiveRegions flags may also be necessary.
     */
    @Advanced
    @Argument(fullName = "disableOptimizations", shortName="disableOptimizations", doc="Don't skip calculations in ActiveRegions with no variants",
            optional = true)
    public boolean disableOptimizations = false;

}
