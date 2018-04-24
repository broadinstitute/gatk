package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.DefaultGATKVariantAnnotationArgumentCollection;
import org.broadinstitute.hellbender.cmdline.GATKPlugin.GATKAnnotationArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardMutectAnnotation;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;

import java.util.Arrays;
import java.util.Collections;

public class M2ArgumentCollection extends AssemblyBasedCallerArgumentCollection {
    private static final long serialVersionUID = 9341L;
    public static final String TUMOR_SAMPLE_LONG_NAME = "tumor-sample";
    public static final String TUMOR_SAMPLE_SHORT_NAME = "tumor";
    public static final String NORMAL_SAMPLE_LONG_NAME = "normal-sample";
    public static final String NORMAL_SAMPLE_SHORT_NAME = "normal";
    public static final String PANEL_OF_NORMALS_LONG_NAME = "panel-of-normals";
    public static final String PANEL_OF_NORMALS_SHORT_NAME = "pon";
    public static final String GENOTYPE_PON_SITES_LONG_NAME = "genotype-pon-sites";
    public static final String GENOTYPE_GERMLINE_SITES_LONG_NAME = "genotype-germline-sites";
    public static final String GERMLINE_RESOURCE_LONG_NAME = "germline-resource";
    public static final String DEFAULT_AF_LONG_NAME = "af-of-alleles-not-in-resource";
    public static final String DEFAULT_AF_SHORT_NAME = "default-af";
    public static final String EMISSION_LOG_LONG_NAME = "tumor-lod-to-emit";
    public static final String EMISSION_LOG_SHORT_NAME = "emit-lod";
    public static final String INITIAL_TUMOR_LOD_LONG_NAME = "initial-tumor-lod";
    public static final String INITIAL_TUMOR_LOD_SHORT_NAME = "init-lod";
    public static final String MAX_POPULATION_AF_LONG_NAME = "max-population-af";
    public static final String MAX_POPULATION_AF_SHORT_NAME = "max-af";
    public static final String DOWNSAMPLING_STRIDE_LONG_NAME = "downsampling-stride";
    public static final String DOWNSAMPLING_STRIDE_SHORT_NAME = "stride";
    public static final String MAX_SUSPICIOUS_READS_PER_ALIGNMENT_START_LONG_NAME = "max-suspicious-reads-per-alignment-start";
    public static final String NORMAL_LOD_LONG_NAME = "normal-lod";

    //TODO: HACK ALERT HACK ALERT HACK ALERT
    //TODO: GATK4 does not yet have a way to tag inputs, eg -I:tumor tumor.bam -I:normal normal.bam,
    //TODO: so for now we require the user to specify bams *both* as inputs, with -I tumor.bam -I normal.bam
    //TODO: *and* as sample names e.g. -tumor <tumor sample> -normal <normal sample>

    @Argument(fullName = TUMOR_SAMPLE_LONG_NAME, shortName = TUMOR_SAMPLE_SHORT_NAME, doc = "BAM sample name of tumor.  May be URL-encoded as output by GetSampleName with -encode argument.", optional = false)
    protected String tumorSampleName = null;

    @Argument(fullName = NORMAL_SAMPLE_LONG_NAME, shortName = NORMAL_SAMPLE_SHORT_NAME, doc = "BAM sample name of normal.  May be URL-encoded as output by GetSampleName with -encode argument.", optional = true)
    protected String normalSampleName = null;

    //TODO: END OF HACK ALERT

    /***************************************/
    // Reference Metadata inputs
    /***************************************/

    /**
     * A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants.
     */
    @Argument(fullName= PANEL_OF_NORMALS_LONG_NAME, shortName = PANEL_OF_NORMALS_SHORT_NAME, doc="VCF file of sites observed in normal.", optional = true)
    public FeatureInput<VariantContext> pon;

    /**
     * Usually we exclude sites in the panel of normals from active region determination, which saves time.  Setting this to true
     * causes Mutect to produce a variant call at these sites.  This call will still be filtered, but it shows up in the vcf.
     */
    @Argument(fullName= GENOTYPE_PON_SITES_LONG_NAME, doc="Call sites in the PoN even though they will ultimately be filtered.", optional = true)
    public boolean genotypePonSites = false;

    /**
     * Usually we exclude sites in the panel of normals from active region determination, which saves time.  Setting this to true
     * causes Mutect to produce a variant call at these sites.  This call will still be filtered, but it shows up in the vcf.
     */
    @Argument(fullName= GENOTYPE_GERMLINE_SITES_LONG_NAME, doc="(EXPERIMENTAL) Call all apparent germline site even though they will ultimately be filtered.", optional = true)
    public boolean genotypeGermlineSites = false;

    /**
     * A resource, such as gnomAD, containing population allele frequencies of common and rare variants.
     */
    @Argument(fullName= GERMLINE_RESOURCE_LONG_NAME, doc="Population vcf of germline sequencing containing allele fractions.", optional = true)
    public FeatureInput<VariantContext> germlineResource;

    /**
     * Population allele fraction assigned to alleles not found in germline resource.
     */
    @Argument(fullName= DEFAULT_AF_LONG_NAME, shortName = DEFAULT_AF_SHORT_NAME,
            doc="Population allele fraction assigned to alleles not found in germline resource.  Please see docs/mutect/mutect2.pdf for" +
                    "a derivation of the default value.", optional = true)
    public double afOfAllelesNotInGermlineResource = 5e-8;

    /**
     * Only variants with tumor LODs exceeding this threshold will be written to the VCF, regardless of filter status.
     * Set to less than or equal to tumor_lod. Increase argument value to reduce false positives in the callset.
     * Default setting of 3 is permissive and will emit some amount of negative training data that 
     * {@link FilterMutectCalls} should then filter.
     */
    @Argument(fullName = EMISSION_LOG_LONG_NAME, shortName = EMISSION_LOG_SHORT_NAME, optional = true, doc = "LOD threshold to emit tumor variant to VCF.")
    public double emissionLodThreshold = 3.0;

    /**
     * Only variants with estimated tumor LODs exceeding this threshold will be considered active.
     */
    @Argument(fullName = INITIAL_TUMOR_LOD_LONG_NAME, shortName = INITIAL_TUMOR_LOD_SHORT_NAME, optional = true, doc = "LOD threshold to consider pileup active.")
    public double initialTumorLodThreshold = 2.0;

    /**
     * In tumor-only mode, we discard variants with population allele frequencies greater than this threshold.
     */
    @Argument(fullName = MAX_POPULATION_AF_LONG_NAME, shortName = MAX_POPULATION_AF_SHORT_NAME, optional = true, doc = "Maximum population allele frequency in tumor-only mode.")
    public double maxPopulationAlleleFrequency = 0.01;

    /**
     * Downsample a pool of reads starting within a range of one or more bases.
     */
    @Argument(fullName = DOWNSAMPLING_STRIDE_LONG_NAME, shortName = DOWNSAMPLING_STRIDE_SHORT_NAME, optional = true, doc = "Downsample a pool of reads starting within a range of one or more bases.")
    public int downsamplingStride = 1;

    /**
     * Maximum number of suspicious reads (mediocre mapping quality or too many substitutions) allowed in a downsampling stride.
     */
    @Advanced
    @Argument(fullName = MAX_SUSPICIOUS_READS_PER_ALIGNMENT_START_LONG_NAME, optional = true, doc = "Maximum number of suspicious reads (mediocre mapping quality or too many substitutions) allowed in a downsampling stride.  Set to 0 to disable.")
    public int maxSuspiciousReadsPerAlignmentStart = 0;

    /**
     * This is a measure of the minimum evidence to support that a variant observed in the tumor is not also present in the normal.
     * Applies to normal data in a tumor with matched normal analysis. The default has been tuned for diploid somatic analyses.
     * It is unlikely such analyses will require changing the default value. Increasing the parameter may increase the sensitivity of somatic calling,
     * but may also increase calling false positive, i.e. germline, variants.
     */
    @Argument(fullName = NORMAL_LOD_LONG_NAME, optional = true, doc = "LOD threshold for calling normal variant non-germline.")
    public double normalLodThreshold = 2.2;


    /**
     * Set of annotation arguments to use.
     * Any requirements that are not met, e.g. failing to provide a pedigree file for a pedigree-based annotation, may cause the run to fail.
     */
    @ArgumentCollection
    GATKAnnotationArgumentCollection defaultGATKVariantAnnotationArgumentCollection = new DefaultGATKVariantAnnotationArgumentCollection(
            Arrays.asList(StandardMutectAnnotation.class.getSimpleName()),
            Collections.emptyList(),
            Collections.emptyList());

}
