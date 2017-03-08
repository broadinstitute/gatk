package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.Hidden;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;

import java.util.ArrayList;
import java.util.Arrays;
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
    public FeatureInput<VariantContext> cosmicFeatureInput;

    /**
     * A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants.
     */
    @Argument(fullName="normal_panel", shortName = "PON", doc="VCF file of sites observed in normal", optional = true)
    public FeatureInput<VariantContext> normalPanelFeatureInput;

    /**
     * Minimum number of variant reads in pileup to be considered an active region.
     */
    @Argument(fullName = "min_variants_in_pileup", optional = true, doc = "Minimum number of reads in pileup to be considered active region.")
    public int minVariantsInPileup = 2;

    /**
     * How many standard deviations above the expected number of variant reads due to error we require for a tumor pielup to be considered active
     */
    @Argument(fullName = "tumorStandardDeviationsThreshold", optional = true, doc = "How many standard deviations above the expected number of variant reads due to error we require for a tumor pielup to be considered active.")
    public int tumorStandardDeviationsThreshold = 2;

    /**
     * Minimum allele fraction of variant reads in normal for a pileup to be considered inactive
     */
    @Argument(fullName = "minNormalVariantFraction", optional = true, doc = "Minimum number of reads in pileup to be considered active region.")
    public double minNormalVariantFraction = 0.1;

    /**
     * Only variants with tumor LODs exceeding this threshold can pass filtering.
     */
    @Argument(fullName = "tumor_lod", optional = true, doc = "LOD threshold for calling tumor variant")
    public double TUMOR_LOD_THRESHOLD = 5.3;

    /**
     * This is a measure of the minimum evidence to support that a variant observed in the tumor is not also present in the normal.
     */
    @Argument(fullName = "normal_lod", optional = true, doc = "LOD threshold for calling normal non-germline")
    public double NORMAL_LOD_THRESHOLD = 2.2;

    /**
     * This argument is used for the M1-style strand bias filter
     */
    @Argument(fullName="power_constant_qscore", doc="Phred scale quality score constant to use in power calculations", optional = true)
    public int POWER_CONSTANT_QSCORE = 30;

    /**
     * Which annotations to add to the output VCF file. See the VariantAnnotator -list argument to view available annotations.
     */
    @Advanced
    @Argument(fullName="annotation", shortName="A", doc="One or more specific annotations to apply to variant calls", optional = true)
    protected List<String> annotationsToUse = new ArrayList<>(Arrays.asList(new String[]{"Coverage", "DepthPerAlleleBySample",
            "TandemRepeat", "OxoGReadCounts", "ClippedBases", "ReadPosition", "BaseQuality", "MappingQuality",
            "FragmentLength", "StrandArtifact"}));

    /**
     * Which groups of annotations to add to the output VCF file. The single value 'none' removes the default group. See
     * the VariantAnnotator -list argument to view available groups. Note that this usage is not recommended because
     * it obscures the specific requirements of individual annotations. Any requirements that are not met (e.g. failing
     * to provide a pedigree file for a pedigree-based annotation) may cause the run to fail.
     */
    @Argument(fullName = "group", shortName = "G", doc = "One or more classes/groups of annotations to apply to variant calls", optional = true)
    public List<String> annotationGroupsToUse = new ArrayList<>(Arrays.asList(new String[]{}));

}
