package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.VariantAnnotationArgumentCollection;
import org.broadinstitute.hellbender.engine.FeatureInput;
import org.broadinstitute.hellbender.tools.walkers.annotator.StandardMutectAnnotation;
import org.broadinstitute.hellbender.tools.walkers.haplotypecaller.AssemblyBasedCallerArgumentCollection;

import java.util.Arrays;
import java.util.Collections;

public class M2ArgumentCollection extends AssemblyBasedCallerArgumentCollection {
    private static final long serialVersionUID = 9341L;

    //TODO: HACK ALERT HACK ALERT HACK ALERT
    //TODO: GATK4 does not yet have a way to tag inputs, eg -I:tumor tumor.bam -I:normal normal.bam,
    //TODO: so for now we require the user to specify bams *both* as inputs, with -I tumor.bam -I normal.bam
    //TODO: *and* as sample names e.g. -tumor <tumor sample> -normal <normal sample>

    @Argument(fullName = "tumor-sample", shortName = "tumor", doc = "BAM sample name of tumor", optional = false)
    protected String tumorSampleName = null;

    @Argument(fullName = "normal-sample", shortName = "normal", doc = "BAM sample name of tumor", optional = true)
    protected String normalSampleName = null;

    //TODO: END OF HACK ALERT

    /***************************************/
    // Reference Metadata inputs
    /***************************************/

    /**
     * A panel of normals can be a useful (optional) input to help filter out commonly seen sequencing noise that may appear as low allele-fraction somatic variants.
     */
    @Argument(fullName="panel-of-normals", shortName = "pon", doc="VCF file of sites observed in normal.", optional = true)
    public FeatureInput<VariantContext> pon;

    /**
     * Usually we exclude sites in the panel of normals from active region determination, which saves time.  Setting this to true
     * causes Mutect to produce a variant call at these sites.  This call will still be filtered, but it shows up in the vcf.
     */
    @Argument(fullName="genotype-pon-sites", doc="Whether to call sites in the PoN even though they will ultimately be filtered.", optional = true)
    public boolean genotypePonSites = false;

    /**
     * A resource, such as gnomAD, containing population allele frequencies of common and rare variants.
     */
    @Argument(fullName="germline-resource", doc="Population vcf of germline sequencing containing allele fractions.", optional = true)
    public FeatureInput<VariantContext> germlineResource;

    /**
     * Population allele fraction assigned to alleles not found in germline resource.
     */
    @Argument(fullName="af-of-alleles-not-in-resource", shortName = "default-af",
            doc="Population allele fraction assigned to alleles not found in germline resource.  A reasonable value is" +
                    "1/(2* number of samples in resource) if a germline resource is available; otherwise an average " +
                    "heterozygosity rate such as 0.001 is reasonable.", optional = true)
    public double afOfAllelesNotInGermlineResource = 0.001;

    /**
     * Prior log-10 probability that any given site has a somatic allele. Impacts germline probability calculation.
     * The workflow uses this parameter only towards the germline event filter. It does NOT relate to the LOD threshold.
     * For example, -6 translates to one in a million or ~3000 somatic mutations per human genome.
     * Depending on tumor type, mutation rate ranges vary (Lawrence et al. Nature 2013), and so adjust parameter accordingly.
     * For higher expected rate of mutation, adjust number up, e.g. -5. For lower expected rate of mutation, adjust number down, e.g. -7.
     */
    @Argument(fullName="log-somatic-prior",
            doc="Prior probability that a given site has a somatic allele.", optional = true)
    public double log10PriorProbOfSomaticEvent = -6.0;

    /**
     * Only variants with tumor LODs exceeding this threshold will be written to the VCF, regardless of filter status.
     * Set to less than or equal to tumor_lod. Increase argument value to reduce false positives in the callset.
     * Default setting of 3 is permissive and will emit some amount of negative training data that 
     * {@link FilterMutectCalls} should then filter.
     */
    @Argument(fullName = "tumor-lod-to-emit", shortName = "emit-lod", optional = true, doc = "LOD threshold to emit tumor variant to VCF.")
    public double emissionLodThreshold = 3.0;

    /**
     * Only variants with estimated tumor LODs exceeding this threshold will be considered active.
     */
    @Argument(fullName = "initial-tumor-lod", shortName = "init-lod", optional = true, doc = "LOD threshold to consider pileup active.")
    public double initialTumorLodThreshold = 2.0;

    /**
     * In tumor-only mode, we discard variants with population allele frequencies greater than this threshold.
     */
    @Argument(fullName = "max-population-af", shortName = "max-af", optional = true, doc = "Maximum population allele frequency in tumor-only mode.")
    public double maxPopulationAlleleFrequency = 0.01;

    /**
     * Downsample a pool of reads starting within a range of one or more bases.
     */
    @Argument(fullName = "downsampling-stride", shortName = "stride", optional = true, doc = "Downsample a pool of reads starting within a range of one or more bases.")
    public int downsamplingStride = 1;

    /**
     * Maximum number of suspicious reads (mediocre mapping quality or too many substitutions) allowed in a downsampling stride.
     */
    @Advanced
    @Argument(fullName = "max-suspicious-reads-per-alignment-start", optional = true, doc = "Maximum number of suspicious reads (mediocre mapping quality or too many substitutions) allowed in a downsampling stride.  Set to 0 to disable.")
    public int maxSuspiciousReadsPerAlignmentStart = 0;

    /**
     * This is a measure of the minimum evidence to support that a variant observed in the tumor is not also present in the normal.
     * Applies to normal data in a tumor with matched normal analysis. The default has been tuned for diploid somatic analyses.
     * It is unlikely such analyses will require changing the default value. Increasing the parameter may increase the sensitivity of somatic calling,
     * but may also increase calling false positive, i.e. germline, variants.
     */
    @Argument(fullName = "normal-lod", optional = true, doc = "LOD threshold for calling normal variant non-germline.")
    public double NORMAL_LOD_THRESHOLD = 2.2;


    /**
     * Set of annotation arguments to use.
     * Any requirements that are not met, e.g. failing to provide a pedigree file for a pedigree-based annotation, may cause the run to fail.
     */
    @ArgumentCollection
    VariantAnnotationArgumentCollection variantAnnotationArgumentCollection = new VariantAnnotationArgumentCollection(
            Arrays.asList(StandardMutectAnnotation.class.getSimpleName()),
            Collections.emptyList(),
            Collections.emptyList());

}
