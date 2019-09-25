package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.MutectReadThreadingAssemblerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReadThreadingAssemblerArgumentCollection;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.ReferenceConfidenceMode;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.readthreading.ReadThreadingAssembler;
import org.broadinstitute.hellbender.tools.walkers.mutect.filtering.FilterMutectCalls;
import org.broadinstitute.hellbender.tools.walkers.readorientation.CollectF1R2CountsArgumentCollection;
import org.broadinstitute.hellbender.utils.MathUtils;

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
    public static final String INITIAL_TUMOR_LOG_10_ODDS_LONG_NAME = "initial-tumor-lod";
    public static final String INITIAL_TUMOR_LOG_10_ODDS_SHORT_NAME = "init-lod";
    public static final String MAX_POPULATION_AF_LONG_NAME = "max-population-af";
    public static final String MAX_POPULATION_AF_SHORT_NAME = "max-af";
    public static final String DOWNSAMPLING_STRIDE_LONG_NAME = "downsampling-stride";
    public static final String DOWNSAMPLING_STRIDE_SHORT_NAME = "stride";
    public static final String MAX_SUSPICIOUS_READS_PER_ALIGNMENT_START_LONG_NAME = "max-suspicious-reads-per-alignment-start";
    public static final String NORMAL_LOG_10_ODDS_LONG_NAME = "normal-lod";
    public static final String IGNORE_ITR_ARTIFACTS_LONG_NAME = "ignore-itr-artifacts";
    public static final String MITOCHONDRIA_MODE_LONG_NAME = "mitochondria-mode";
    public static final String CALLABLE_DEPTH_LONG_NAME = "callable-depth";
    public static final String PCR_SNV_QUAL_LONG_NAME = "pcr-snv-qual";
    public static final String PCR_INDEL_QUAL_LONG_NAME = "pcr-indel-qual";
    public static final String F1R2_TAR_GZ_NAME = "f1r2-tar-gz";

    public static final double DEFAULT_AF_FOR_TUMOR_ONLY_CALLING = 5e-8;
    public static final double DEFAULT_AF_FOR_TUMOR_NORMAL_CALLING = 1e-6;
    public static final double DEFAULT_AF_FOR_MITO_CALLING = 4e-3;
    public static final double DEFAULT_EMISSION_LOG_10_ODDS = 3.0;
    public static final double DEFAULT_MITO_EMISSION_LOD = 0;
    public static final double DEFAULT_INITIAL_LOG_10_ODDS = 2.0;
    public static final double DEFAULT_NORMAL_LOG_10_ODDS = 2.2;
    public static final double DEFAULT_MITO_INITIAL_LOG_10_ODDS = 0;
    public static final double DEFAULT_GVCF_LOG_10_ODDS = Double.NEGATIVE_INFINITY;
    public static final int DEFAULT_CALLABLE_DEPTH = 10;

    public static final double DEFAULT_MITO_PRUNING_LOG_ODDS_THRESHOLD = MathUtils.log10ToLog(-4);
    public static final String INDEPENDENT_MATES_LONG_NAME = "independent-mates";
    public static final String MINIMUM_ALLELE_FRACTION_LONG_NAME = "minimum-allele-fraction";
    public static final String MINIMUM_ALLELE_FRACTION_SHORT_NAME = "min-AF";
    public static final String LOD_BAND_LONG_NAME = "gvcf-lod-band";
    public static final String LOD_BAND_SHORT_NAME = "LODB";

    @Override
    protected int getDefaultMaxMnpDistance() { return 1; }

    @Override
    protected ReadThreadingAssemblerArgumentCollection getReadThreadingAssemblerArgumentCollection(){
        return new MutectReadThreadingAssemblerArgumentCollection();
    }

    @Override
    public ReadThreadingAssembler createReadThreadingAssembler(){
        if(mitochondria ) {
            assemblerArgs.recoverAllDanglingBranches = true;
            if (assemblerArgs.pruningLogOddsThreshold == ReadThreadingAssemblerArgumentCollection.DEFAULT_PRUNING_LOG_ODDS_THRESHOLD) {
                assemblerArgs.pruningLogOddsThreshold = DEFAULT_MITO_PRUNING_LOG_ODDS_THRESHOLD;
            }
        }

        return super.createReadThreadingAssembler();
    }

    @ArgumentCollection
    public CollectF1R2CountsArgumentCollection f1r2Args = new CollectF1R2CountsArgumentCollection();

    @Argument(fullName = F1R2_TAR_GZ_NAME, doc = "If specified, collect F1R2 counts and output files into this tar.gz file", optional = true)
    public File f1r2TarGz;

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
    @Argument(fullName = EMISSION_LOD_LONG_NAME, shortName = EMISSION_LOG_SHORT_NAME, optional = true, doc = "Log 10 odds threshold to emit variant to VCF.")
    private double emissionLog10Odds = DEFAULT_EMISSION_LOG_10_ODDS;

    public double getEmissionLogOdds() {
        if (emitReferenceConfidence != ReferenceConfidenceMode.NONE) {
            return MathUtils.log10ToLog(DEFAULT_GVCF_LOG_10_ODDS);
        }
        return MathUtils.log10ToLog(mitochondria && emissionLog10Odds == DEFAULT_EMISSION_LOG_10_ODDS ? DEFAULT_MITO_EMISSION_LOD : emissionLog10Odds);
    }

    /**
     * Only variants with estimated tumor LODs exceeding this threshold will be considered active.
     */
    @Argument(fullName = INITIAL_TUMOR_LOG_10_ODDS_LONG_NAME, shortName = INITIAL_TUMOR_LOG_10_ODDS_SHORT_NAME, optional = true, doc = "Log 10 odds threshold to consider pileup active.")
    private double initialLog10Odds = DEFAULT_INITIAL_LOG_10_ODDS;

    public double getInitialLogOdds() {
        if (emitReferenceConfidence != ReferenceConfidenceMode.NONE) {
            return MathUtils.log10ToLog(DEFAULT_GVCF_LOG_10_ODDS);
        }
        return MathUtils.log10ToLog(mitochondria && initialLog10Odds == DEFAULT_INITIAL_LOG_10_ODDS ? DEFAULT_MITO_INITIAL_LOG_10_ODDS : initialLog10Odds);
    }


    @Argument(fullName = PCR_SNV_QUAL_LONG_NAME, optional = true, doc = "Phred-scaled PCR SNV qual for overlapping fragments")
    public int pcrSnvQual = 40;

    @Argument(fullName = PCR_INDEL_QUAL_LONG_NAME, optional = true, doc = "Phred-scaled PCR SNV qual for overlapping fragments")
    public int pcrIndelQual = 40;

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

    @Argument(fullName = CALLABLE_DEPTH_LONG_NAME, optional = true, doc = "Minimum depth to be considered callable for Mutect stats.  Does not affect genotyping.")
    public int callableDepth = DEFAULT_CALLABLE_DEPTH;

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
    @Argument(fullName = NORMAL_LOG_10_ODDS_LONG_NAME, optional = true, doc = "Log 10 odds threshold for calling normal variant non-germline.")
    public double normalLog10Odds = DEFAULT_NORMAL_LOG_10_ODDS;

    /**
     * When opposite ends of a fragment are inverted tandem repeats of each other, the sequence past one end may be copied onto the other
     * during library prep.  By default, Mutect2 identifies and clips these artifacts, which are especially prevalent when
     * DNA is damaged as in the case of FFPE samples and ancient DNA.
     */
    @Argument(fullName= IGNORE_ITR_ARTIFACTS_LONG_NAME, doc="Turn off read transformer that clips artifacts associated with end repair insertions near inverted tandem repeats.", optional = true)
    public boolean dontClipITRArtifacts = false;

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
    @Argument(fullName = LOD_BAND_LONG_NAME, shortName = LOD_BAND_SHORT_NAME, doc = "Exclusive upper bounds for reference confidence LOD bands " +
            "(must be specified in increasing order)", optional = true)
    public List<Double> GVCFGQBands = new ArrayList<>(70);
    {
        for (double i = -2.5; i <= 1; i = i + 0.5) {
            GVCFGQBands.add(i);
        }
    }

    @Advanced
    @Argument(fullName = MINIMUM_ALLELE_FRACTION_LONG_NAME, shortName = MINIMUM_ALLELE_FRACTION_SHORT_NAME, doc = "Lower bound of variant allele fractions to consider when calculating variant LOD", optional = true)
    public double minAF = 0.00;
    
    @Advanced
    @Argument(fullName = INDEPENDENT_MATES_LONG_NAME, doc = "Allow paired reads to independently support different haplotypes.  Useful for validations with ill-designed synthetic data.", optional = true)
    public boolean independentMates = false;
}
