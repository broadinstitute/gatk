package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.Hidden;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class M2ArgumentCollection extends AssemblyBasedCallerArgumentCollection {
    private static final long serialVersionUID = 9341L;

    //TODO: HACK ALERT HACK ALERT HACK ALERT
    //TODO: GATK4 does not yet have a way to tag inputs, eg -I:tumor tumor.bam -I:normal normal.bam,
    //TODO: so for now we require the user to specify bams *both* as inputs, with -I tumor.bam -I normal.bam
    //TODO: *and* as sample names e.g. -tumor tumorSampleName -normal normalSampleName

    @Argument(fullName = "tumorSampleName", shortName = "tumor", doc = "BAM sample name of tumor", optional = false)
    protected String tumorSampleName = null;

    @Argument(fullName = "normalSampleName", shortName = "normal", doc = "BAM sample name of tumor", optional = true)
    protected String normalSampleName = null;

    //TODO: END OF HACK ALERT

    /***************************************/
    // Reference Metadata inputs
    /***************************************/
    /**
     * Mutect2 has the ability to use COSMIC data in conjunction with dbSNP to adjust the threshold for evidence of a variant
     * in the normal.  If a variant is present in dbSNP, but not in COSMIC, then more evidence is required from the normal
     * sample to prove the variant is not present in germline.
     */
    @Argument(fullName="cosmic", shortName = "cosmic", doc="VCF file of COSMIC sites", optional = true)
    public List<FeatureInput<VariantContext>> cosmicFeatureInput = Collections.emptyList();

    /**
     * A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants.
     */
    @Argument(fullName="normal_panel", shortName = "PON", doc="VCF file of sites observed in normal", optional = true)
    public FeatureInput<VariantContext> normalPanelFeatureInput;

    /**
     * Artifact detection mode is used to prepare a panel of normals. This maintains the specified tumor LOD threshold,
     * but disables the remaining pragmatic filters. See usage examples above for more information.
     */
    @Advanced
    @Argument(fullName = "artifact_detection_mode", optional = true, doc="Enable artifact detection for creating panels of normals")
    public boolean ARTIFACT_DETECTION_MODE = false;

    /**
     * This is the LOD threshold that a variant must pass in the tumor to be emitted to the VCF. Note that the variant may pass this threshold yet still be annotated as FILTERed based on other criteria.
     */
    @Argument(fullName = "initial_tumor_lod", optional = true, doc = "Initial LOD threshold for calling tumor variant")
    public double INITIAL_TUMOR_LOD_THRESHOLD = 4.0;

    /**
     * Only variants with tumor LODs exceeding this threshold can pass filtering.
     */
    @Argument(fullName = "tumor_lod", optional = true, doc = "LOD threshold for calling tumor variant")
    public double TUMOR_LOD_THRESHOLD = 6.3;

    /**
     * This is a measure of the minimum evidence to support that a variant observed in the tumor is not also present in the normal.
     */
    @Argument(fullName = "normal_lod", optional = true, doc = "LOD threshold for calling normal non-germline")
    public double NORMAL_LOD_THRESHOLD = 2.2;

    /**
     * The LOD threshold for the normal is typically made more strict if the variant has been seen in dbSNP (i.e. another
     * normal sample). We thus require MORE evidence that a variant is NOT seen in this tumor's normal if it has been observed as a germline variant before.
     */
    @Argument(fullName = "dbsnp_normal_lod", optional = true, doc = "LOD threshold for calling normal non-variant at dbsnp sites")
    public double NORMAL_DBSNP_LOD_THRESHOLD = 5.5;

    /**
     * This argument is used for the M1-style strand bias filter
     */
    @Argument(fullName="power_constant_qscore", doc="Phred scale quality score constant to use in power calculations", optional = true)
    public int POWER_CONSTANT_QSCORE = 30;

    @Hidden
    @Argument(fullName = "strand_artifact_lod", optional = true, doc = "LOD threshold for calling strand bias")
    public float STRAND_ARTIFACT_LOD_THRESHOLD = 2.0f;

    @Hidden
    @Argument(fullName = "strand_artifact_power_threshold", optional = true, doc = "power threshold for calling strand bias")
    public float STRAND_ARTIFACT_POWER_THRESHOLD = 0.9f;

    @Argument(fullName = "enable_strand_artifact_filter", optional = true, doc = "turn on strand artifact filter")
    public boolean ENABLE_STRAND_ARTIFACT_FILTER = false;

    @Argument(fullName = "enable_clustered_read_position_filter", optional = true, doc = "turn on clustered read position filter")
    public boolean ENABLE_CLUSTERED_READ_POSITION_FILTER = false;

    /**
     * This argument is used for the M1-style read position filter
     */
    @Argument(fullName = "pir_median_threshold", optional = true, doc="threshold for clustered read position artifact median")
    public double PIR_MEDIAN_THRESHOLD = 10;

    /**
     * This argument is used for the M1-style read position filter
     */
    @Argument(fullName = "pir_mad_threshold", optional = true, doc="threshold for clustered read position artifact MAD")
    public double PIR_MAD_THRESHOLD = 3;

    /**
     * Reads with mapping qualities below this threshold will be filtered out
     */
    @Argument(fullName = "min_mapping_quality_score", shortName = "mmq", doc = "Minimum read mapping quality required to consider a read for analysis", optional = true)
    public int MIN_MAPPING_QUALITY_SCORE = 20;

    /**
     * Which annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available annotations.
     * //TODO: port TandemRepeatAnnotator and put it here
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", optional = true)
    protected List<String> annotationsToUse = new ArrayList<>(Arrays.asList(new String[]{"DepthPerAlleleBySample", "BaseQualitySumPerAlleleBySample", "OxoGReadCounts"}));
}
