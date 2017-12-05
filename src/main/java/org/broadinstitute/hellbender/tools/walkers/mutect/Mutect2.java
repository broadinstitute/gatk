package org.broadinstitute.hellbender.tools.walkers.mutect;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.tools.walkers.contamination.GetPileupSummaries;
import org.broadinstitute.hellbender.tools.exome.FilterByOrientationBias;
import org.broadinstitute.hellbender.tools.walkers.contamination.CalculateContamination;
import org.broadinstitute.hellbender.utils.downsampling.MutectDownsampler;
import org.broadinstitute.hellbender.utils.downsampling.ReadsDownsampler;

import java.io.File;
import java.util.List;

/**
 * Call somatic short variants, both SNVs and indels, via local assembly of haplotypes
 *
 * <p>
 *     Mutect2 calls somatic single nucleotide (SNV) and insertion and deletion (indel) variants.
 *     The caller combines the DREAM challenge-winning somatic genotyping engine of the original MuTect
 *     (<a href='http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html'>Cibulskis et al., 2013</a>) with the
 *     assembly-based machinery of <a href="https://www.broadinstitute.org/gatk/documentation/tooldocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php">HaplotypeCaller</a>.
 *     Although we present the tool for somatic analyses, it may also apply to other contexts.
 * </p>
 *
 * <h3>How GATK4 Mutect2 differs from GATK3 MuTect2</h3>
 *
 * <dl>
 *     <dd>(i) The filtering functionality is now a separate tool called {@link FilterMutectCalls}.
 *     To filter further based on sequence context artifacts, additionally use {@link FilterByOrientationBias}.</dd>
 *     <dd>(ii) If using a known germline variants resource, then it must contain population allele frequencies, e.g.
 *     from gnomAD or the 1000 Genomes Project. The VCF INFO field contains the allele frequency (AF) tag.
 *     See below or the GATK Resource Bundle for an example.</dd>
 *     <dd>(iii) To create the panel of normals (PoN), call on each normal sample using Mutect2's tumor-only mode and then use GATK4's {@link CreateSomaticPanelOfNormals}.
 *     This contrasts with the GATK3 workflow, which uses an artifact mode in MuTect2 and CombineVariants for PoN creation.
 *     In GATK4, omitting filtering with FilterMutectCalls achieves the same artifact mode.</dd>
 *     <dd>(iv) Instead of using a maximum likelihood estimate, GATK4 Mutect2 marginalizes over allele fractions. 
 *     GATK3 MuTect2 directly uses allele depths (AD) to estimate allele fractions and calculate likelihoods. In contrast, GATK4 Mutect2
 *     factors for the statistical error inherent in allele depths by marginalizing over allele fractions when calculating likelihoods.</dd>
 *     <dd>(v) GATK4 Mutect2 recommends including contamination estimates with the -contaminationFile option from {@link CalculateContamination},
 *     which in turn relies on the results of {@link GetPileupSummaries}.</dd>
 * </dl>
 *
 * <p>
 *     What remains unchanged is that neither tool versions call on seeming loss of heterozygosity (LoH) events.
 *     To detect LoH, see the Copy Number Variant (CNV) and AllelicCNV workflows.
 * </p>
 *
 * <p>Here is an example of a known variants resource with population allele frequencies:</p>
 *
 * <pre>
 *     #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO
 *      1       10067   .       T       TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCC      30.35   PASS    AC=3;AF=7.384E-5
 *      1       10108   .       CAACCCT C       46514.32        PASS    AC=6;AF=1.525E-4
 *      1       10109   .       AACCCTAACCCT    AAACCCT,*       89837.27        PASS    AC=48,5;AF=0.001223,1.273E-4
 *      1       10114   .       TAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTA  *,CAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCCTAACCCTAACCCTAAACCCTA,T      36728.97        PASS    AC=55,9,1;AF=0.001373,2.246E-4,2.496E-5
 *      1       10119   .       CT      C,*     251.23  PASS    AC=5,1;AF=1.249E-4,2.498E-5
 *      1       10120   .       TA      CA,*    14928.74        PASS    AC=10,6;AF=2.5E-4,1.5E-4
 *      1       10128   .       ACCCTAACCCTAACCCTAAC    A,*     285.71  PASS    AC=3,1;AF=7.58E-5,2.527E-5
 *      1       10131   .       CT      C,*     378.93  PASS    AC=7,5;AF=1.765E-4,1.261E-4
 *      1       10132   .       TAACCC  *,T     18025.11        PASS    AC=12,2;AF=3.03E-4,5.049E-5
 * </pre>
 *
 * <h3>How Mutect2 works compared to HaplotypeCaller</h3>
 * <p>Overall, Mutect2 works similarly to HaplotypeCaller, but with a few key differences. </p>
 * <p>
 *     (i) GVCF calling is not a feature of Mutect2.
 *     (ii) While HaplotypeCaller relies on a fixed ploidy assumption to inform its genotype likelihoods that are the basis for genotype probabilities (PL),
 *     Mutect2 allows for varying ploidy in the form of allele fractions for each variant.
 *     Varying allele fractions is often seen within a tumor sample due to fractional purity, multiple subclones and/or copy number variation.
 *     (iii) Mutect2 also differs from the HaplotypeCaller in that it can apply various prefilters to sites and variants depending on the use of
 *     a matched normal (--normalSampleName), a panel of normals (PoN; --normal_panel) and/or a common population variant resource containing allele-specific frequencies (--germline_resource).
 *     If provided, Mutect2 uses the PoN to filter sites and the germline resource and matched normal to filter alleles.
 *     (iv) Mutect2's default variant site annotations differ from those of HaplotypeCaller. See the --annotation parameter description for a list.
 *     (v) Finally, Mutect2 has additional parameters not available to HaplotypeCaller that factor in the decision to reassemble a genomic region,
 *     factor in likelihood calculations that then determine whether to emit a variant, or factor towards filtering.
 *     These parameters include the following and are each described further in the arguments section.
 * </p>
 *
 * <dl>
 *     <dd>--min_variants_in_pileup ==> active region determination</dd>
 *     <dd>--minNormalVariantFraction ==> active region determination</dd>
 *     <dd>--tumorStandardDeviationsThreshold ==> active region determination</dd>
 *     <dd>--af_of_alleles_not_in_resource ==> germline variant prior</dd>
 *     <dd>--log_somatic_prior ==> somatic variant prior</dd>
 *     <dd>--normal_lod ==> filter threshold for variants in tumor not being in the normal, i.e. germline-risk filter</dd>
 *     <dd>--tumor_lod_to_emit ==> cutoff for tumor variants to appear in callset</dd>
 * </dl>
 *
 * <h3>Further points of interest</h3>
 * <p>
 *     Additional parameters that factor towards filtering, including normal_artifact_lod (default threshold 0.0) and
 *     tumor_lod (default threshold 5.3), are available in {@link FilterMutectCalls}. While the tool calculates
 *     normal_lod with a fixed ploidy assumption given by the --sample_ploidy option (default is 2), it calculates
 *     normal_artifact_lod with the same approach it uses for tumor_lod, i.e. with a variable ploidy assumption.
 * </p>
 *
 * <dl>
 *     <dd>If the normal artifact log odds becomes large, then FilterMutectCalls applies the artifact-in-normal filter.
 *     For matched normal samples with tumor contamination, consider increasing the normal_artifact_lod threshold.</dd>
 *     <dd>The tumor log odds, which is calculated independently of any matched normal, determines whether to filter a tumor
 *     variant. Variants with tumor LODs exceeding the threshold pass filtering.</dd>
 * </dl>
 *
 * <p>
 *     If a variant is absent from a given germline resource, then the value for --af_of_alleles_not_in_resource applies. 
 *     For example, gnomAD's 16,000 samples (~32,000 homologs per locus) becomes a probability of one in 32,000 or less.
 *     Thus, an allele's absence from the germline resource becomes evidence that it is not a germline variant.
 * </p>
 *
 * <h3>Examples</h3>
 *
 * <p>Example commands show how to run Mutect2 for typical scenerios.</p>
 *
 * <h4>Tumor with matched normal</h4>
 * <p>
 *     Given a matched normal, Mutect2 is designed to call somatic variants only. The tool includes logic to skip
 *     emitting variants that are clearly present in the germline based on the evidence present in the matched normal. This is done at
 *     an early stage to avoid spending computational resources on germline events. If the variant's germline status is
 *     borderline, then Mutect2 will emit the variant to the callset with a germline-risk filter. Such filtered
 *     emissions enable manual review.
 * </p>
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" Mutect2 \
 *   -R ref_fasta.fa \
 *   -I tumor.bam \
 *   -tumor tumor_sample_name \
 *   -I normal.bam \
 *   -normal normal_sample_name \
 *   --germline_resource af-only-gnomad.vcf.gz \
 *   --normal_panel pon.vcf.gz \
 *   -L intervals.list \
 *   -O tumor_matched_m2_snvs_indels.vcf.gz
 * </pre>
 *
 * <h4>Single tumor sample</h4>
 * <pre>
 *  gatk-launch --javaOptions "-Xmx4g" Mutect2 \
 *   -R ref_fasta.fa \
 *   -I tumor.bam \
 *   -tumor tumor_sample_name \
 *   --germline_resource af-only-gnomad.vcf.gz \
 *   --normal_panel pon.vcf.gz \
 *   -L intervals.list \
 *   -O tumor_unmatched_m2_snvs_indels.vcf.gz
 * </pre>
 *
 * <h4>Single normal sample for panel of normals (PoN) creation</h4>
 * <p>
 *    To create a panel of normals (PoN), call on each normal sample as if a tumor sample. Then use
 *    {@link CreateSomaticPanelOfNormals} to output a PoN of germline and artifactual sites. This contrasts with the
 *    GATK3 workflow, which uses CombineVariants to retain variant sites called in at least two samples and then uses
 *    Picard MakeSitesOnlyVcf to simplify the callset for use as a PoN.
 * </p>
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" Mutect2 \
 *   -R ref_fasta.fa \
 *   -I normal1.bam \
 *   -tumor normal1_sample_name \
 *   --germline_resource af-only-gnomad.vcf.gz \
 *   -L intervals.list \
 *   -O normal1_for_pon.vcf.gz
 * </pre>
 *
 * <h3>Caveats</h3>
 * <p>
 *     Although GATK4 Mutect2 is optimized to accomodate varying coverage depths, further optimization of parameters
 *     is necessary for extreme high depths, e.g. 1000X.
 * </p>
 */
 @CommandLineProgramProperties(
         summary = "Call somatic SNVs and indels via local assembly of haplotypes",
         oneLineSummary = "Call somatic SNVs and indels via local assembly of haplotypes",
         programGroup = VariantProgramGroup.class
 )
@DocumentedFeature
@BetaFeature
public final class Mutect2 extends AssemblyRegionWalker {

    @ArgumentCollection
    protected M2ArgumentCollection MTAC = new M2ArgumentCollection();

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "File to which variants should be written")
    public File outputVCF;

    private VariantContextWriter vcfWriter;

    private Mutect2Engine m2Engine;

    @Override
    protected int defaultReadShardSize() { return NO_INTERVAL_SHARDING; }

    @Override
    protected int defaultReadShardPadding() { return 100; }

    @Override
    protected int defaultMinAssemblyRegionSize() { return 50; }

    @Override
    protected int defaultMaxAssemblyRegionSize() { return 300; }

    @Override
    protected int defaultAssemblyRegionPadding() { return 100; }

    @Override
    protected int defaultMaxReadsPerAlignmentStart() { return 50; }

    @Override
    protected double defaultActiveProbThreshold() { return 0.002; }

    @Override
    protected int defaultMaxProbPropagationDistance() { return 50; }

    @Override
    protected boolean includeReadsWithDeletionsInIsActivePileups() { return true; }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return Mutect2Engine.makeStandardMutect2ReadFilters();
    }

    @Override
    protected ReadsDownsampler createDownsampler() {
        return new MutectDownsampler(maxReadsPerAlignmentStart, MTAC.maxSuspiciousReadsPerAlignmentStart, MTAC.downsamplingStride);
    }

    @Override
    public AssemblyRegionEvaluator assemblyRegionEvaluator() { return m2Engine; }

    @Override
    public void onTraversalStart() {
        m2Engine = new Mutect2Engine(MTAC, createOutputBamIndex, createOutputBamMD5, getHeaderForReads(), referenceArguments.getReferenceFileName());
        final SAMSequenceDictionary sequenceDictionary = getHeaderForReads().getSequenceDictionary();
        vcfWriter = createVCFWriter(outputVCF);
        m2Engine.writeHeader(vcfWriter, sequenceDictionary, getDefaultToolVCFHeaderLines());
    }

    @Override
    public Object onTraversalSuccess() {
        return "SUCCESS";
    }

    @Override
    public void apply(final AssemblyRegion region, final ReferenceContext referenceContext, final FeatureContext featureContext ) {
        m2Engine.callRegion(region, referenceContext, featureContext).forEach(vcfWriter::add);
    }

    @Override
    public void closeTool() {
        if ( vcfWriter != null ) {
            vcfWriter.close();
        }

        if ( m2Engine != null ) {
            m2Engine.shutdown();
        }
    }
}
