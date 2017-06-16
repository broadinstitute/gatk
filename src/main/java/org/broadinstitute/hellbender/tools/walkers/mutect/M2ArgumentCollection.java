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
     * A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants.
     */
    @Argument(fullName="normal_panel", shortName = "PON", doc="VCF file of sites observed in normal", optional = true)
    public FeatureInput<VariantContext> pon;

    /**
     * A resource, such as gnomAD, containing population allele frequencies of common and rare variants.
     */
    @Argument(fullName="germline_resource", doc="Population vcf of germline sequencing containing allele fractions", optional = true)
    public FeatureInput<VariantContext> germlineResource;

    /**
     * Population allele fraction assigned to alleles not found in germline resource.
     */
    @Argument(fullName="af_of_alleles_not_in_resource", shortName = "default_af",
            doc="Population allele fraction assigned to alleles not found in germline resource.  A reasonable value is" +
                    "1/(2* number of samples in resource) if a germline resource is available; otherwise an average " +
                    "heterozygosity rate such as 0.001 is reasonable.", optional = true)
    public double afOfAllelesNotInGermlineResource = 0.001;

    /**
     * Prior probability that any given site has a somatic allele. Impacts germline probability calculation.
     * For example, -6 translates to one in a million or ~3000 somatic mutations per human genome.
     * Depending on tumor type, mutation rate ranges vary (Lawrence et al. Nature 2013), and so adjust parameter accordingly.
     * For higher expected rate of mutation, adjust number up, e.g. -5. For lower expected rate of mutation, adjust number down, e.g. -7.
     */
    @Argument(fullName="log_somatic_prior",
            doc="Prior probability that a given site has a somatic allele.", optional = true)
    public double log10PriorProbOfSomaticEvent = -6.0;


    /**
     * Minimum number of variant reads in pileup to be considered an active region.
     * In the hypothetical case of extreme high quality alignments, consider adjusting parameter down to 1.
     * In the rare case of extreme deep coverage and where low fraction alleles are not of interest, adjust up to 3.
     * Providing genomic intervals of interest with -L, setting the tumor standard deviation to zero and
     * setting minimum variants in pileup to zero forces the tool to consider all provided regions as active.
     */
    @Argument(fullName = "min_variants_in_pileup", optional = true, doc = "Minimum number of variant reads in pileup to be considered active region.")
    public int minVariantsInPileup = 2;

    /**
     * How many standard deviations above the expected number of variant reads due to error we require for tool to consider a tumor pileup active.
     * Argument sets the z-score. Here, base qualities inform expected error rate.
     * Providing genomic intervals of interest with -L, setting the tumor standard deviation to zero and
     * setting minimum variants in pileup to zero forces the tool to consider all provided regions as active.
     */
    @Argument(fullName = "tumorStandardDeviationsThreshold", optional = true, doc = "How many standard deviations above the expected number of variant reads due to error we require for a tumor pileup to be considered active.")
    public int tumorStandardDeviationsThreshold = 2;

    /**
     * Minimum fraction of variant reads in normal for a pileup to be considered inactive. Applies to normal data in a tumor with matched normal analysis.
     * For value of 0.1, at least one tenth of pileup must be variant for tool to consider it a germline variant site and therefore not a region of interest
     * in the somatic analysis.
     */
    @Argument(fullName = "minNormalVariantFraction", optional = true, doc = "Minimum fraction of variant reads in normal pileup to be considered germline.")
    public double minNormalVariantFraction = 0.1;

    /**
     * Only variants with tumor LODs exceeding this threshold will be written to the VCF, irregardless of filter status.
     * Set to less than or equal to tumor_lod. Increase argument value to reduce false positive calls in the callset.
     * Default setting of 3 is permissive and will emit some amount of false-positive calls that
     * filtering should filter and that allows for downstream training outside of the current workflow, e.g. akin to germline VQSR.
     */
    @Argument(fullName = "tumor_lod_to_emit", optional = true, doc = "LOD threshold to emit tumor variant to results")
    public double emissionLodThreshold = 3.0;

    /**
     * This is a measure of the minimum evidence to support that a variant observed in the tumor is not also present in the normal.
     * Applies to normal data in a tumor with matched normal analysis.
     */
    @Argument(fullName = "normal_lod", optional = true, doc = "LOD threshold for calling normal variant non-germline")
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
