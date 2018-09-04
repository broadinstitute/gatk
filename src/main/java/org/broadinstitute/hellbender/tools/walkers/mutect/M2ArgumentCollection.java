package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;

import java.io.File;

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
    public static final String EMISSION_LOD_LONG_NAME = "tumor-lod-to-emit";
    public static final String EMISSION_LOG_SHORT_NAME = "emit-lod";
    public static final String INITIAL_TUMOR_LOD_LONG_NAME = "initial-tumor-lod";
    public static final String INITIAL_TUMOR_LOD_SHORT_NAME = "init-lod";
    public static final String INITIAL_PCR_ERROR_QUAL = "initial-pcr-qual";
    public static final String MAX_POPULATION_AF_LONG_NAME = "max-population-af";
    public static final String MAX_POPULATION_AF_SHORT_NAME = "max-af";
    public static final String DOWNSAMPLING_STRIDE_LONG_NAME = "downsampling-stride";
    public static final String DOWNSAMPLING_STRIDE_SHORT_NAME = "stride";
    public static final String MAX_SUSPICIOUS_READS_PER_ALIGNMENT_START_LONG_NAME = "max-suspicious-reads-per-alignment-start";
    public static final String NORMAL_LOD_LONG_NAME = "normal-lod";
    public static final String MAX_MNP_DISTANCE_LONG_NAME = "max-mnp-distance";
    public static final String MAX_MNP_DISTANCE_SHORT_NAME = "mnp-dist";
    public static final String IGNORE_ITR_ARTIFACTS_LONG_NAME = "ignore-itr-artifacts";
    public static final String ARTIFACT_PRIOR_TABLE_NAME = "orientation-bias-artifact-priors";
    public static final String GET_AF_FROM_AD_LONG_NAME = "get-af-from-ad";

    public static final double DEFAULT_AF_FOR_TUMOR_ONLY_CALLING = 5e-8;
    public static final double DEFAULT_AF_FOR_TUMOR_NORMAL_CALLING = 1e-6;

    //TODO: HACK ALERT HACK ALERT HACK ALERT
    //TODO: GATK4 does not yet have a way to tag inputs, eg -I:tumor tumor.bam -I:normal normal.bam,
    //TODO: so for now we require the user to specify bams *both* as inputs, with -I tumor.bam -I normal.bam
    //TODO: *and* as sample names e.g. -tumor <tumor sample> -normal <normal sample>

    @Argument(fullName = TUMOR_SAMPLE_LONG_NAME, shortName = TUMOR_SAMPLE_SHORT_NAME, doc = "BAM sample name of tumor.  May be URL-encoded as output by GetSampleName with -encode argument.", optional = false)
    protected String tumorSample = null;

    @Argument(fullName = NORMAL_SAMPLE_LONG_NAME, shortName = NORMAL_SAMPLE_SHORT_NAME, doc = "BAM sample name of normal.  May be URL-encoded as output by GetSampleName with -encode argument.", optional = true)
    protected String normalSample = null;

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
    public double afOfAllelesNotInGermlineResource = -1;

    public double getDefaultAlleleFrequency() {
        return afOfAllelesNotInGermlineResource >= 0 ? afOfAllelesNotInGermlineResource :
                (normalSample == null ? DEFAULT_AF_FOR_TUMOR_ONLY_CALLING : DEFAULT_AF_FOR_TUMOR_NORMAL_CALLING);
    }

    /**
     * Only variants with tumor LODs exceeding this threshold will be written to the VCF, regardless of filter status.
     * Set to less than or equal to tumor_lod. Increase argument value to reduce false positives in the callset.
     * Default setting of 3 is permissive and will emit some amount of negative training data that 
     * {@link FilterMutectCalls} should then filter.
     */
    @Argument(fullName = EMISSION_LOD_LONG_NAME, shortName = EMISSION_LOG_SHORT_NAME, optional = true, doc = "LOD threshold to emit tumor variant to VCF.")
    public double emissionLod = 3.0;

    /**
     * Only variants with estimated tumor LODs exceeding this threshold will be considered active.
     */
    @Argument(fullName = INITIAL_TUMOR_LOD_LONG_NAME, shortName = INITIAL_TUMOR_LOD_SHORT_NAME, optional = true, doc = "LOD threshold to consider pileup active.")
    public double initialTumorLod = 2.0;

    /**
     * PCR error rate for overlapping fragments in isActive()
     */
    @Argument(fullName = INITIAL_PCR_ERROR_QUAL, optional = true, doc = "PCR error rate for overlapping fragments in isActive()")
    public int initialPCRErrorQual = 40;

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
    public double normalLod = 2.2;

    @Argument(fullName = ARTIFACT_PRIOR_TABLE_NAME, optional = true, doc = "table of prior artifact probabilities for the read orientation filter model")
    public File artifactPriorTable = null;

    /**
     * Two or more phased substitutions separated by this distance or less are merged into MNPs.
     */
    @Advanced
    @Argument(fullName = MAX_MNP_DISTANCE_LONG_NAME, shortName = MAX_MNP_DISTANCE_SHORT_NAME,
            doc = "Two or more phased substitutions separated by this distance or less are merged into MNPs.", optional = true)
    public int maxMnpDistance = 1;

    /**
     * When opposite ends of a fragment are inverted tandem repeats of each other, the sequence past one end may be copied onto the other
     * during library prep.  By default, Mutect2 identifies and clips these artifacts, which are especially prevalent when
     * DNA is damaged as in the case of FFPE samples and ancient DNA.
     */
    @Argument(fullName= IGNORE_ITR_ARTIFACTS_LONG_NAME, doc="Turn off read transformer that clips artifacts associated with end repair insertions near inverted tandem repeats.", optional = true)
    public boolean dontClipITRArtifacts = false;

    @Advanced
    @Argument(fullName = GET_AF_FROM_AD_LONG_NAME, doc="Use allelic depth to calculate tumor allele fraction; recommended for mitochondrial applications", optional = true)
    public boolean calculateAFfromAD = false;

}
