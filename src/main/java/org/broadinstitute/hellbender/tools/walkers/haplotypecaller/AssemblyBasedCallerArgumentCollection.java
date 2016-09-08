package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.Hidden;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.genotyper.StandardCallerArgumentCollection;
import org.broadinstitute.hellbender.utils.haplotype.HaplotypeBAMWriter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * Set of arguments for Assembly Based Callers
 */
public abstract class AssemblyBasedCallerArgumentCollection extends StandardCallerArgumentCollection {
    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    public AssemblyRegionTrimmerArgumentCollection assemblyRegionTrimmerArgs = new AssemblyRegionTrimmerArgumentCollection();

    @ArgumentCollection
    public ReadThreadingAssemblerArgumentCollection assemblerArgs = new ReadThreadingAssemblerArgumentCollection();

    @ArgumentCollection
    public LikelihoodEngineArgumentCollection likelihoodArgs = new LikelihoodEngineArgumentCollection();

    /**
     * rsIDs from this file are used to populate the ID column of the output. Also, the DB INFO flag will be set when appropriate.
     * dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    public DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     * as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field). Records that are
     * filtered in the comp track will be ignored. Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Advanced
    @Argument(fullName = "comp", shortName = "comp", doc = "Comparison VCF file(s)", optional = true)
    public List<FeatureInput<VariantContext>> comps = Collections.emptyList();


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
    public String bamOutputPath = null;

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

    // -----------------------------------------------------------------------------------------------
    // arguments for debugging / developing
    // -----------------------------------------------------------------------------------------------

    @Hidden
    @Argument(fullName = "keepRG", shortName = "keepRG", doc = "Only use reads from this read group when making calls (but use all reads to build the assembly)", optional = true)
    public String keepRG = null;

    /**
     * This argument is intended for benchmarking and scalability testing.
     */
    @Hidden
    @Argument(fullName = "justDetermineActiveRegions", shortName = "justDetermineActiveRegions", doc = "Just determine ActiveRegions, don't perform assembly or calling", optional = true)
    public boolean justDetermineActiveRegions = false;

    /**
     * This argument is intended for benchmarking and scalability testing.
     */
    @Hidden
    @Argument(fullName = "dontGenotype", shortName = "dontGenotype", doc = "Perform assembly but do not genotype variants", optional = true)
    public boolean dontGenotype = false;

    @Advanced
    @Argument(fullName = "dontUseSoftClippedBases", shortName = "dontUseSoftClippedBases", doc = "Do not analyze soft clipped bases in the reads", optional = true)
    public boolean dontUseSoftClippedBases = false;

    @Hidden
    @Argument(fullName = "captureAssemblyFailureBAM", shortName = "captureAssemblyFailureBAM", doc = "Write a BAM called assemblyFailure.bam capturing all of the reads that were in the active region when the assembler failed for any reason", optional = true)
    public boolean captureAssemblyFailureBAM = false;

    // Parameters to control read error correction
    /**
     * Enabling this argument may cause fundamental problems with the assembly graph itself.
     */
    @Hidden
    @Argument(fullName = "errorCorrectReads", shortName = "errorCorrectReads", doc = "Use an exploratory algorithm to error correct the kmers used during assembly", optional = true)
    public boolean errorCorrectReads = false;

    /**
     * As of GATK 3.3, HaplotypeCaller outputs physical (read-based) information (see version 3.3 release notes and documentation for details). This argument disables that behavior.
     */
    @Advanced
    @Argument(fullName = "doNotRunPhysicalPhasing", shortName = "doNotRunPhysicalPhasing", doc = "Disable physical phasing", optional = true)
    public boolean doNotRunPhysicalPhasing = false;

    /**
     * Bases with a quality below this threshold will not be used for calling.
     */
    @Argument(fullName = "min_base_quality_score", shortName = "mbq", doc = "Minimum base quality required to consider a base for calling", optional = true)
    public byte MIN_BASE_QUALTY_SCORE = 10;

    /**
     * Reads with mapping qualities below this threshold will be filtered out
     */
    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for analysis with the HaplotypeCaller", optional = true)
    public int MIN_MAPPING_QUALITY_SCORE = 20;

    //Annotations
    /**
     * Which groups of annotations to add to the output VCF file. The single value 'none' removes the default group. See
     * the VariantAnnotator -list argument to view available groups. Note that this usage is not recommended because
     * it obscures the specific requirements of individual annotations. Any requirements that are not met (e.g. failing
     * to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    @Argument(fullName = "group", shortName = "G", doc = "One or more classes/groups of annotations to apply to variant calls", optional = true)
    public List<String> annotationGroupsToUse = new ArrayList<>(Arrays.asList(new String[]{"StandardAnnotation", "StandardHCAnnotation"}));

    /**
     * Which annotations to exclude from output in the VCF file.  Note that this argument has higher priority than the
     * -A or -G arguments, so these annotations will be excluded even if they are explicitly included with the other
     * options. When HaplotypeCaller is run with -ERC GVCF or -ERC BP_RESOLUTION, some annotations are excluded from the
     * output by default because they will only be meaningful once they have been recalculated by GenotypeGVCFs. As
     * of version 3.3 this concerns ChromosomeCounts, FisherStrand, StrandOddsRatio and QualByDepth.
     *
     */
    @Advanced
    @Argument(fullName = "excludeAnnotation", shortName = "XA", doc = "One or more specific annotations to exclude", optional = true)
    public List<String> annotationsToExclude = new ArrayList<>();

}
