package org.broadinstitute.hellbender.tools.walkers.haplotypecaller;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.DbsnpArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Set of arguments for the {@link HaplotypeCallerEngine}
 */
public class HaplotypeCallerArgumentCollection extends AssemblyBasedCallerArgumentCollection implements Serializable{
    private static final long serialVersionUID = 1L;

    public static final String MAX_MNP_DISTANCE_LONG_NAME = "max-mnp-distance";
    public static final String MAX_MNP_DISTANCE_SHORT_NAME = "mnp-dist";

    /**
     * You can use this argument to specify that HC should process a single sample out of a multisample BAM file. This
     * is especially useful if your samples are all in the same file but you need to run them individually through HC
     * in -ERC GVC mode (which is the recommended usage). Note that the name is case-sensitive.
     */
    @Argument(fullName = StandardArgumentDefinitions.SAMPLE_NAME_LONG_NAME, shortName = StandardArgumentDefinitions.SAMPLE_ALIAS_SHORT_NAME, doc = "Name of single sample to use from a multi-sample bam", optional = true)
    public String sampleNameToUse = null;

    /**
     * rsIDs from this file are used to populate the ID column of the output. Also, the DB INFO flag will be set when appropriate.
     * dbSNP is not used in any way for the calculations themselves.
     */
    @ArgumentCollection
    public DbsnpArgumentCollection dbsnp = new DbsnpArgumentCollection();

    // -----------------------------------------------------------------------------------------------
    // general advanced arguments to control haplotype caller behavior
    // -----------------------------------------------------------------------------------------------

    /**
     * If a call overlaps with a record from the provided comp track, the INFO field will be annotated
     * as such in the output with the track name (e.g. -comp:FOO will have 'FOO' in the INFO field). Records that are
     * filtered in the comp track will be ignored. Note that 'dbSNP' has been special-cased (see the --dbsnp argument).
     */
    @Advanced
    @Argument(fullName = "comp", shortName = "comp", doc = "Comparison VCF file(s)", optional = true)
    public List<FeatureInput<VariantContext>> comps = new ArrayList<>();

    /**
     * The reference confidence mode makes it possible to emit a per-bp or summarized confidence estimate for a site being strictly homozygous-reference.
     * See https://software.broadinstitute.org/gatk/documentation/article.php?id=4017 for more details of how this works.
     */
    @Advanced
    @Argument(fullName="emit-ref-confidence", shortName="ERC", doc="Mode for emitting reference confidence scores", optional = true)
    public ReferenceConfidenceMode emitReferenceConfidence = ReferenceConfidenceMode.NONE;

    /**
     * When HC is run in reference confidence mode with banding compression enabled (-ERC GVCF), homozygous-reference
     * sites are compressed into bands of similar genotype quality (GQ) that are emitted as a single VCF record. See
     * the FAQ documentation for more details about the GVCF format.
     *
     * This argument allows you to set the GQ bands. HC expects a list of strictly increasing GQ values
     * that will act as exclusive upper bounds for the GQ bands. To pass multiple values,
     * you provide them one by one with the argument, as in `-GQB 10 -GQB 20 -GQB 30` and so on
     * (this would set the GQ bands to be `[0, 10), [10, 20), [20, 30)` and so on, for example).
     * Note that GQ values are capped at 99 in the GATK, so values must be integers in [1, 100].
     * If the last value is strictly less than 100, the last GQ band will start at that value (inclusive)
     * and end at 100 (exclusive).
     */
    @Advanced
    @Argument(fullName = "gvcf-gq-bands", shortName = "GQB", doc= "Exclusive upper bounds for reference confidence GQ bands " +
            "(must be in [1, 100] and specified in increasing order)", optional = true)
    public List<Integer> GVCFGQBands = new ArrayList<>(70);
    {
            for (int i=1; i<=60; ++i) {
                GVCFGQBands.add(i);
            }
            GVCFGQBands.add(70); GVCFGQBands.add(80); GVCFGQBands.add(90); GVCFGQBands.add(99);
    };

    /**
     * This parameter determines the maximum size of an indel considered as potentially segregating in the
     * reference model.  It is used to eliminate reads from being indel informative at a site, and determines
     * by that mechanism the certainty in the reference base.  Conceptually, setting this parameter to
     * X means that each informative read is consistent with any indel of size < X being present at a specific
     * position in the genome, given its alignment to the reference.
     */
    @Advanced
    @Argument(fullName = "indel-size-to-eliminate-in-ref-model", doc = "The size of an indel to check for in the reference model", optional = true)
    public int indelSizeToEliminateInRefModel = 10;


    @Advanced
    @Argument(fullName = "use-alleles-trigger", doc = "Use additional trigger on variants found in an external alleles file", optional = true)
    public boolean USE_ALLELES_TRIGGER = false;

    /**
     * If set, certain "early exit" optimizations in HaplotypeCaller, which aim to save compute and time by skipping
     * calculations if an ActiveRegion is determined to contain no variants, will be disabled. This is most likely to be useful if
     * you're using the -bamout argument to examine the placement of reads following reassembly and are interested in seeing the mapping of
     * reads in regions with no variations. Setting the --force-active and --dont-trim-active-regions flags may also be necessary.
     */
    @Advanced
    @Argument(fullName = "disable-optimizations", doc="Don't skip calculations in ActiveRegions with no variants",
            optional = true)
    public boolean disableOptimizations = false;

    @Hidden
    @Argument(fullName = "keep-rg", doc = "Only use reads from this read group when making calls (but use all reads to build the assembly)", optional = true)
    public String keepRG = null;

    /**
     * This argument is intended for benchmarking and scalability testing.
     */
    @Hidden
    @Argument(fullName = "just-determine-active-regions", doc = "Just determine ActiveRegions, don't perform assembly or calling", optional = true)
    public boolean justDetermineActiveRegions = false;

    /**
     * This argument is intended for benchmarking and scalability testing.
     */
    @Hidden
    @Argument(fullName = "dont-genotype", doc = "Perform assembly but do not genotype variants", optional = true)
    public boolean dontGenotype = false;

    /**
     * Two or more phased substitutions separated by this distance or less are merged into MNPs.
     */
    @Advanced
    @Argument(fullName = MAX_MNP_DISTANCE_LONG_NAME, shortName = MAX_MNP_DISTANCE_SHORT_NAME,
            doc = "Two or more phased substitutions separated by this distance or less are merged into MNPs.", optional = true)
    public int maxMnpDistance = 0;
}
