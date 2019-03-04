package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.MutectReadThreadingAssemblerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadThreadingAssemblerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;

import java.io.File;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

public class M2ArgumentCollection extends AssemblyBasedCallerArgumentCollection implements Serializable {
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
    public static final String ANNOTATE_BASED_ON_READS_LONG_NAME = "count-reads";
    public static final String MEDIAN_AUTOSOMAL_COVERAGE_LONG_NAME = "median-autosomal-coverage";
    public static final String MITOCHONDRIA_MODE_LONG_NAME = "mitochondria-mode";

    public static final double DEFAULT_AF_FOR_TUMOR_ONLY_CALLING = 5e-8;
    public static final double DEFAULT_AF_FOR_TUMOR_NORMAL_CALLING = 1e-6;
    public static final double DEFAULT_AF_FOR_MITO_CALLING = 4e-3;
    public static final double DEFAULT_EMISSION_LOD = 3.0;
    public static final double DEFAULT_MITO_EMISSION_LOD = 0;
    public static final double DEFAULT_INITIAL_LOD = 2.0;
    public static final double DEFAULT_MITO_INITIAL_LOD = 0;
    public static final double DEFAULT_GVCF_LOD = Double.NEGATIVE_INFINITY;

    public static final double DEFAULT_MITO_PRUNING_LOG_ODDS_THRESHOLD = -4;

    protected ReadThreadingAssemblerArgumentCollection getReadThreadingAssemblerArgumentCollection(){
        return new MutectReadThreadingAssemblerArgumentCollection();
    }

    @Override
    public ReadThreadingAssembler createReadThreadingAssembler(){
        if(mitochondria ) {
            assemblerArgs.recoverAllDanglingBranches = true;
            if (assemblerArgs.pruningLog10OddsThreshold == ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD) {
                assemblerArgs.pruningLog10OddsThreshold = DEFAULT_MITO_PRUNING_LOG_ODDS_THRESHOLD;
            }
        }

        return super.createReadThreadingAssembler();
    }

    //TODO: HACK ALERT HACK ALERT HACK ALERT
    //TODO: GATK4 does not yet have a way to tag inputs, eg -I:tumor tumor.bam -I:normal normal.bam,
    //TODO: so for now we require the user to specify bams *both* as inputs, with -I tumor.bam -I normal.bam
    //TODO: *and* as sample names e.g. -tumor <tumor sample> -normal <normal sample>

    // As of GATK 4.1, any sample not specified as the normal is considered a tumor sample
    @Deprecated
    @Argument(fullName = TUMOR_SAMPLE_LONG_NAME, shortName = TUMOR_SAMPLE_SHORT_NAME, doc = "BAM sample name of tumor.  May be URL-encoded as output by GetSampleName with -encode argument.", optional = true)
    protected String tumorSample = null;

    @Argument(fullName = NORMAL_SAMPLE_LONG_NAME, shortName = NORMAL_SAMPLE_SHORT_NAME, doc = "BAM sample name of normal(s), if any.  May be URL-encoded as output by GetSampleName with -encode argument.", optional = true)
    protected List<String> normalSamples = new ArrayList<>();

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
    private double afOfAllelesNotInGermlineResource = -1;

    public double getDefaultAlleleFrequency() {
        return afOfAllelesNotInGermlineResource >= 0 ? afOfAllelesNotInGermlineResource :
                (mitochondria ? DEFAULT_AF_FOR_MITO_CALLING:
                (normalSamples.isEmpty() ? DEFAULT_AF_FOR_TUMOR_ONLY_CALLING : DEFAULT_AF_FOR_TUMOR_NORMAL_CALLING));
    }

    /**
     * Mitochondria mode changes default values for --tumor-lod-to-emit and --initial-tumor-lod to 0 to increase sensitivity.
     * --tumor-sample is also not explicitly required in mitochondria mode since a single sample bam is expected as input.
     * Mitochondria mode is also required in FilterMutectCalls if used here.
     */
    @Argument(fullName = MITOCHONDRIA_MODE_LONG_NAME, optional = true, doc="Mitochondria mode sets emission and initial LODs to 0.")
    public Boolean mitochondria = false;

    /**
     * Only variants with tumor LODs exceeding this threshold will be written to the VCF, regardless of filter status.
     * Set to less than or equal to tumor_lod. Increase argument value to reduce false positives in the callset.
     * Default setting of 3 is permissive and will emit some amount of negative training data that 
     * {@link FilterMutectCalls} should then filter.
     */
    @Argument(fullName = EMISSION_LOD_LONG_NAME, shortName = EMISSION_LOG_SHORT_NAME, optional = true, doc = "LOD threshold to emit variant to VCF.")
    private double emissionLodArg = DEFAULT_EMISSION_LOD;

    public double getEmissionLod() {
        if (emitReferenceConfidence != ReferenceConfidenceMode.NONE) {
            return DEFAULT_GVCF_LOD;
        }
        return mitochondria && emissionLodArg == DEFAULT_EMISSION_LOD ? DEFAULT_MITO_EMISSION_LOD : emissionLodArg;
    }

    /**
     * Only variants with estimated tumor LODs exceeding this threshold will be considered active.
     */
    @Argument(fullName = INITIAL_TUMOR_LOD_LONG_NAME, shortName = INITIAL_TUMOR_LOD_SHORT_NAME, optional = true, doc = "LOD threshold to consider pileup active.")
    private double initialLod = DEFAULT_INITIAL_LOD;

    public double getInitialLod() {
        if (emitReferenceConfidence != ReferenceConfidenceMode.NONE) {
            return DEFAULT_GVCF_LOD;
        }
        return mitochondria && initialLod == DEFAULT_INITIAL_LOD ? DEFAULT_MITO_INITIAL_LOD : initialLod;
    }

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

    @Argument(fullName = ARTIFACT_PRIOR_TABLE_NAME, optional = true, doc = "tables of prior artifact probabilities for the read orientation filter model, one per tumor sample")
    public List<File> artifactPriorTables = new ArrayList<>();

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

    /**
     * If set to true, count an overlapping read pair as two separate reads instead of one for {@link StrandArtifact} and {@link StrandBiasBySample} annotations,
     * which is the correct behavior for these annotations. Note that doing so would break the independence assumption of reads and over-count the alt depth in these annotations.
     * On the other hand it could prevent spurious false negatives that could arise if by chance one strand in overlapping pairs is dropped disproportionately
     */
    @Argument(fullName = ANNOTATE_BASED_ON_READS_LONG_NAME, doc = "If true, strand-based annotations use the number of reads, not fragments")
    public boolean annotateBasedOnReads = false;

    /**
     * Used to model autosomal coverage when calling mitochondria. The median tends to be a more robust center statistic.
     */
    @Advanced
    @Argument(fullName = MEDIAN_AUTOSOMAL_COVERAGE_LONG_NAME, doc="For mitochondrial calling only; Annotate possible polymorphic NuMT based on Poisson distribution given median autosomal coverage", optional = true)
    public double autosomalCoverage;

    /**
     * When Mutect2 is run in reference confidence mode with banding compression enabled (-ERC GVCF), homozygous-reference
     * sites are compressed into bands of similar tumor LOD (TLOD) that are emitted as a single VCF record. See
     * the FAQ documentation for more details about the GVCF format.
     * <p>
     * This argument allows you to set the TLOD bands. Mutect2 expects a list of strictly increasing TLOD values
     * that will act as exclusive upper bounds for the TLOD bands. To pass multiple values,
     * you provide them one by one with the argument, as in `-LODB -3.0 -LODB -2.0 -LODB -1.0`
     * (this would set the TLOD bands to be `[-Inf, -3.0), [-3.0, -2.0), [-2.0, -1.0), [-1.0, Inf]`, for example).
     * <p>
     * Note that, unlike the GQ used by HaplotypeCaller GVCFs, here the reference calls with the highest confidence are the most negative.
     */
    @Advanced
    @Argument(fullName = "gvcf-lod-band", shortName = "LODB", doc = "Exclusive upper bounds for reference confidence LOD bands " +
            "(must be specified in increasing order)", optional = true)
    public List<Double> GVCFGQBands = new ArrayList<>(70);
    {
        for (double i = -2.5; i <= 1; i = i + 0.5) {
            GVCFGQBands.add(i);
        }
    }

    @Advanced
    @Argument(fullName = "minimum-allele-fraction", shortName = "min-AF", doc = "Lower bound of variant allele fractions to consider when calculating variant LOD", optional = true)
    public double minAF = 0.00;

    @Argument(fullName="alleles", doc="The set of alleles for which to force genotyping regardless of evidence", optional=true)
    public FeatureInput<VariantContext> alleles;

    @Advanced
    @Argument(fullName = "genotype-filtered-alleles", doc = "Whether to force genotype even filtered alleles", optional = true)
    public boolean genotypeFilteredAlleles = false;
}
